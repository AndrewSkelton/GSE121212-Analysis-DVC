##
## Helper Functions for Analysis 
## With RNA Abundances
##

##
## Get Results from Limma fit
##
## Parameters: 
##    x = Contrast Matrix
##    y = Limma Fit
##    Thres = T (logFC + FDR Filter), or F (Return Everything)
##
getRes <- function(x = NULL, y = NULL, Thres = T) {
  
  results.out.all  <- list()
  results.out.filt <- list()
  
  for(i in colnames(x)) {
    
    results.out.all[[i]]   <- topTable(y, coef = i, number = Inf, confint = T) %>% 
                              as.data.frame %>% rownames_to_column("ensgene") %>% 
                              mutate(ensgene = gsub("_at|\\..*","",ensgene)) %>% 
                              left_join({grch38[,c(1,3,9)] %>% distinct}) %>% 
                              mutate(description = gsub("\\[[^\\]]*\\]", "", description, perl = T)) 
    results.out.filt[[i]]  <- results.out.all[[i]] %>% 
                              filter(abs(logFC) >= log2(1.5), adj.P.Val <= 0.05)
    
    write_csv(results.out.filt[[i]],   file = paste0(var_out,"/Differential_Expression/FDR0_05_FC1_5__",i,".csv"))
    write_csv({results.out.all[[i]]},  file = paste0(var_out,"/Differential_Expression/All_",i,".csv.gz"))
    
  }
  
  if(Thres == T) {
    return(results.out.filt)
  } else {
    return(results.out.all)
  }
  
}
##


## Make Volcano Plot
# X = List results from ucb.getRes
# cont = The contrast string to render a volano plot of
##
volcano_plot <- function(x = NULL, cont = NULL) {
  plot.data    <- x[[cont]] %>% 
                  as.data.frame %>% 
                  dplyr::select(ensgene,logFC,P.Value,FDR=adj.P.Val,symbol) %>% 
                  mutate(`-log10(FDR)` = -log10(FDR),
                         DE           = ifelse((FDR <= 0.05 & abs(logFC) >= log2(1.5)),"DE","Not DE"),
                         Rank         = abs(logFC) * `-log10(FDR)`,
                         Ranku        = logFC * `-log10(FDR)`,
                         Rankd        = -logFC * `-log10(FDR)`) %>% 
                  arrange(desc(Rank)) %>% 
                  mutate(IDX = 1:nrow(.),Label = ifelse(IDX < 35, symbol, "")) %>% 
                  arrange(desc(Ranku)) %>% mutate(IDX1 = 1:nrow(.)) %>% 
                  arrange(desc(Rankd)) %>% mutate(IDX2 = 1:nrow(.)) %>% 
                  mutate(Label = ifelse(IDX1 < 15 | IDX2 < 15, symbol, "")) %>% 
                  mutate(Label = ifelse(FDR < 0.05 & abs(logFC) > 2.5, symbol, Label))
                  
  
  plot.nums    <- c({plot.data %>%  .$ensgene %>% unique %>% length},
                    {plot.data %>% filter(DE == "DE") %>% .$ensgene %>% unique %>% length})
  plot.volcano <- plot.data %>% 
                  ggplot(., aes(logFC,`-log10(FDR)`, colour = DE)) +
                  geom_point(data = . %>% filter(DE != "DE"), alpha = 0.4, stroke = 0, size = 1) +
                  geom_point(data = . %>% filter(DE == "DE"), alpha = 0.6, stroke = 0, size = 2) +
                  geom_hline(yintercept = -log10(0.05), colour = "red", linetype = "dashed", alpha = .4) +
                  geom_vline(xintercept = c(log2(1.5),-log2(1.5)), colour = "red", linetype = "dashed", alpha = .4) +
                  geom_text_repel(data = plot.data %>% filter(IDX < 35), aes(label = Label), show.legend = F, colour = "black", nudge_x = 0.1, nudge_y = 0.1, size = 4) +
                  theme_bw(base_size = 15) +
                  theme(legend.position = "bottom",
                        axis.text  = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        plot.title = element_text(size = 14, face = "bold")) + 
                  scale_colour_manual(name = "", values = c("Not DE"="grey", "DE"="#E41A1C")) +
                  labs(title    = "Volcano Plot", 
                       subtitle = paste0(gsub("_", " ",cont),"\n",
                                         plot.nums[2], " Genes DE (",
                                         round((plot.nums[2]/plot.nums[1]*100),2),"%)"),
                       y = "-log10(FDR)", x = "logFC")
  
  png(paste0(var_out,"/Volcano/VP_",cont,".png"),
      width  = 8,height = 8,
      units  = "in",res    = 300)
  print(plot.volcano); dev.off()
  return(plot.volcano)
}
##

## Enrichment Handler - Enrichment
clean_enrich <- function(x, class = "") {
  x %>% mutate(ID = "No ID", FDR = NA,Class = class) %>% dplyr::select(ID,Term,Pval = pvalue,FDR = qvalue,Class)
}
##



##
## Function to compute MSE from regression for parameter tuning
## 
## sample = Vector of sample gene expression values.
## sigMat = The referent profile signature matrix.
## nu_ks  = Parameter for inferring proportion of support vectors; define several.
##
getMSE <- function(sample, sigMat, nu_ks = seq(0.1, 0.9, 0.05)){
  mses        <- sapply(nu_ks, function(k){
    svr_model    <- svm(sigMat, sample, scale = F, nu = k,
                        kernel = "linear", type = "nu-regression")
    predictedm   <- predict(svr_model, sigMat)
    error        <- sample - predictedm
    mean(error^2)
  })
  names(mses) <- as.character(nu_ks)
  mses
}
##


##
## Function to do the modeling for a given sample
## idx    = sample index
## Y      = sample data frame
## sigMat = signature Matrix
##
modelSample <- function(idx, Y, sigMat){
  
  ## Get best k
  mses                      <- getMSE(Y[,idx,drop = F], sigMat)
  best_k                    <- as.numeric(names(which(min(mses) == mses)))
  ##
  
  ## Run model for best k; get coefficients
  tunedModel                <- svm(sigMat, Y[,idx,drop = F], scale = F, nu = best_k, 
                                   kernel = "linear", type = "nu-regression")
  ##
  
  ## Estimate coefficients; remove negative values; normalize contributions
  coefTuned                 <- t(tunedModel$coefs) %*% tunedModel$SV
  coefTuned[coefTuned < 0 ] <- 0
  coefTunedNormed           <- round(coefTuned /sum(coefTuned), 2)
  ##
  
  return(coefTunedNormed)
}
##