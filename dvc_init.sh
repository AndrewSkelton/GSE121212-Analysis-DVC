## Config for a pull
pip install dvc
pip install dvc-s3
aws configure --profile dvc
# Note, keys and secrets needed for s3 here
git clone https://github.com/AndrewSkelton/GSE121212-Analysis-DVC 
dvc pull
##

## Setting up the repo at first
# Init
git init
dvc init
git add .dvc .gitignore
git commit -m "Initialize DVC"

# Capture binary files for DVC
dvc add *.RDS
find ./Data ./Workstreams ./Pheno ./Reference -type f -name "*.RDS" | parallel -j1 dvc add {}
dvc add *.gz
find ./Data ./Workstreams ./Pheno ./Reference -type f -name "*.gz" | parallel -j1 dvc add {}

# Commit to Git
git add .
git remote add origin https://github.com/AndrewSkelton/GSE121212-Analysis-DVC.git
git commit -m "Initialise"
git branch -M main
git push -u origin main

# Commit to dvcs
dvc push 
##


## .gitignore
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