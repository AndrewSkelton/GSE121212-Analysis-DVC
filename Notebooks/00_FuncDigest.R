##
## Functional Digest
## Written by Andrew J Skelton
##
## Take a set of enriched functions & DE Genes
## Make a UMAP embedding of Gene to Function Membership
## Assign Clusters to embedding
## Create a Cluster Membership Table
## Create a gene to cluster membership table
##
functional_digest <- function(funcs,genes,genesets,universe,mode_go = F) {
  
  data_out <- list()
  
  data_g2s <- sapply(funcs, simplify = F, function(y) {
                
                genes_tmp <- intersect(genes, genesets[[y]])
                genes_out <- universe %in% genes_tmp
                return(genes_out)
                
              }) %>% 
              purrr::reduce(rbind) %>% 
              `rownames<-`(funcs) %>% 
              `colnames<-`(universe) %>% 
              .[rowSums(.) > 0, colSums(.) > 0, drop=F] %>% 
              {. * 1}
  
  if(mode_go == T) {
    rownames(data_g2s) <- Skin.GeneSets$GO.terms$Term.Info %>% 
                          filter(Type == "biological_process", 
                                 `Term ID` %in% rownames(data_g2s)) %>% 
                          .[match(rownames(data_g2s), .$`Term ID`),] %>% 
                          .$Name
  }
  
  if(nrow(data_g2s) > 1 & ncol(data_g2s) > 1) {
    
    data_g2s_sym   <- data.frame(ensgene = colnames(data_g2s)) %>% 
                      left_join({grch38[,c(1,3)] %>% distinct}) %>% 
                      distinct(ensgene,.keep_all = T)
    data_g2s2      <- data_g2s %>% 
                      .[,colnames(.) %in% data_g2s_sym$ensgene, drop=F] %>% 
                      `colnames<-`(data_g2s_sym$symbol)
    data_g2s2      <- data_g2s2 %>% 
                      .[{data_g2s2 %>% dist %>% hclust %>% .$order},, drop=F] %>% 
                      .[,{data_g2s2 %>% t %>% dist %>% hclust %>% .$order}, drop=F]
    
    data_plotlyg2s <- plot_ly(y = rownames(data_g2s2), x = colnames(data_g2s2),
                              z = data_g2s2, type = "heatmap",
                              colorscale = data.frame(z = c(0,1), col = c("#FFFFFF","#1f77b4")),
                              colorbar = list(ypad = 30, tickvals=c(0.25,0.75), ticktext=c(0,1))) %>% 
                      layout(margin = list(l=360))
    
    data_out[["plotlyg2s"]]    <- data_plotlyg2s

  } else {
    data_out[["plotlyg2s"]]    <- NULL
  }
  
  
  if(!is.null({data_g2s %>% dim}) & (nrow(data_g2s) > 14 & ncol(data_g2s) > 10)) {
      
      ## UMAP of Binary
      data_umap        <- data_g2s %>% umap %>% .$layout %>% 
                          as.data.frame %>% `colnames<-`(c("UMAP1","UMAP2"))
      ##
      
      
      ## Louvain Clustering using KNN of UMAP Embedding
      k          <- 10 
      data_knn   <- data_g2s %>% get.knn(., k = k)
      data_knn2  <- data.frame(from   = rep(1:nrow(data_knn$nn.index),k),
                               to     = {data_knn$nn.index %>% as.vector}, 
                               weight = 1 / (1 + as.vector(data_knn$nn.dist))) %>% 
                    graph_from_data_frame(., directed = F) %>% 
                    simplify %>% cluster_louvain %>% membership
      
      data_umap$Cluster <- data_knn2 %>% factor
      ##
      
      ## ggplot of UMAP Clusters
      data_ggumap    <- data_umap %>% 
                        ggplot(., aes(UMAP1,UMAP2,colour = Cluster)) + 
                          geom_point() + 
                          ggsci::scale_colour_d3("category20", name = "") + 
                          theme_bw(base_size = 15) +
                          theme(legend.position = "bottom") +
                          labs(title = "UMAP Embedding of Enriched Terms based on Gene Membership",
                               x = "UMAP 1", y = "UMAP 2")
      ##
      
      ## plotly of UMAP Embedding
      # data_umap <- data_breakdown$umap
      data_plotly    <- plot_ly(data = {data_umap %>% rownames_to_column("Term")},
                                x = ~UMAP1, y = ~UMAP2, text = ~Term, color = ~Cluster,
                                colors = ggsci::pal_d3("category20")(length(unique(data_umap$Cluster)))) %>% 
                        layout(title = "Gene:Pathway Membership",
                               xaxis = list(zeroline = F),
                               yaxis = list(zeroline = F))
      ##
      
      ## Cluster Membership 
      ##
      data_clusmem   <- data_umap %>% 
                        rownames_to_column("Term") %>% 
                        arrange(Cluster)
      ##
      
      ## Gene Membership 
      ##
      data_genemem <- list()
      for(z in {data_umap$Cluster %>% unique %>% sort}) {
        
        if(mode_go == F) {
          
          data_genemem[[z]]  <- genesets %>% 
                                .[{
                                  data_clusmem %>% 
                                    filter(Cluster == z) %>% 
                                    .$Term
                                }] %>% unlist %>% unique %>% 
                                .[. %in% genes] %>% 
                                data.frame(ensgene = ., Cluster = z)
          
        } else {
          
          data_genemem[[z]]  <- genesets %>% 
                                .[{
                                  Skin.GeneSets$GO.terms$Term.Info %>% 
                                    filter(Type == "biological_process",
                                           Name %in% {
                                             data_clusmem %>% 
                                               filter(Cluster == z) %>% 
                                               .$Term
                                           }) %>% .$`Term ID`
                                  
                                }] %>% unlist %>% unique %>% 
                                .[. %in% genes] %>% 
                                data.frame(ensgene = ., Cluster = z)
          
        }
        
      }
      ##
      
      ## Return
      data_out[["data_g2s"]]     <- data_g2s
      data_out[["umap"]]         <- data_umap
      data_out[["ggumap"]]       <- data_ggumap
      data_out[["plotlyumap"]]   <- data_plotly
      data_out[["plotlyg2s"]]    <- data_plotlyg2s
      data_out[["clusmem"]]      <- data_clusmem
      data_out[["genemem"]]      <- data_genemem
      return(data_out)
      ##
  } else {
    
    ## Return
    data_out[["data_g2s"]]     <- data_g2s
    data_out[["umap"]]         <- NULL
    data_out[["ggumap"]]       <- NULL
    data_out[["plotlyumap"]]   <- NULL
    data_out[["clusmem"]]      <- NULL
    data_out[["genemem"]]      <- NULL
    return(data_out)
    
  }
}