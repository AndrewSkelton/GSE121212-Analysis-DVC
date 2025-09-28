git init
dvc init
git add .dvc .gitignore
git commit -m "Initialize DVC"

dvc add *.RDS
find ./Data ./Workstreams ./Pheno ./Reference -type f -name "*.RDS" | parallel -j1 dvc add {}
dvc add *.gz
find ./Data ./Workstreams ./Pheno ./Reference -type f -name "*.gz" | parallel -j1 dvc add {}

git add .
git remote add origin https://github.com/AndrewSkelton/GSE121212-Analysis-DVC.git
git commit -m "Initialise"
git branch -M main
git push -u origin main


# ## R environment management
# Notebooks/renv/
# !Notebooks/renv.lock

# # R session files
# .Rhistory
# .RData
# .Ruserdata
# .Rproj.user/

# # RStudio project files
# *.Rproj

# # System-specific junk
# .DS_Store
# Thumbs.db

# # DVC Handlers
# *.gz
# *.RDS 
# *.rds