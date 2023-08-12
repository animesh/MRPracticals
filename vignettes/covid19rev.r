#https://www.r-bloggers.com/2022/08/installation-of-r-4-2-on-ubuntu-22-04-1-lts-and-tips-for-spatial-packages/
#sudo apt-get purge r-base* r-recommended r-cran-*
#sudo apt autoremove
#sudo apt install -y --no-install-recommends software-properties-common dirmngr
#wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
#sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
#sudo apt install -y r-base r-base-core r-recommended r-base-dev

#https://mrcieu.github.io/TwoSampleMR/
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

#http://app.mrbase.org/
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- extract_instruments(c('ebi-a-GCST011075'))
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ebi-a-GCST011075'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

