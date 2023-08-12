#https://www.r-bloggers.com/2022/08/installation-of-r-4-2-on-ubuntu-22-04-1-lts-and-tips-for-spatial-packages/
#sudo apt-get purge r-base* r-recommended r-cran-*
#sudo apt autoremove
#sudo apt install -y --no-install-recommends software-properties-common dirmngr
#wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
#sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

#https://launchpad.net/~c2d4u.team/+archive/ubuntu/c2d4u4.0+
#sudo apt install -y r-base r-base-core r-recommended r-base-dev
#sudo add-apt-repository ppa:c2d4u.team/c2d4u4.0+
#sudo apt update
#sudo apt install -y r-cran-lme4 r-cran-glmnet r-cran-meta 

#https://mrcieu.github.io/TwoSampleMR/
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")

#https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html
library(TwoSampleMR)
library(ggplot2)
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
#> API: public: http://gwas-api.mrcieu.ac.uk/
chd_out_dat <- extract_outcome_data(snps = bmi_exp_dat$SNP, outcomes = 'ieu-a-7')
#> Extracting data for 79 SNP(s) from 1 GWAS(s)
dat <- harmonise_data(bmi_exp_dat, chd_out_dat)
#> Harmonising Body mass index || id:ieu-a-2 (ieu-a-2) and Coronary heart disease || id:ieu-a-7 (ieu-a-7)
res <- mr(dat)
p1 <- mr_scatter_plot(res, dat)
length(p1)
p1[[1]]
ggsave(p1[[1]], file = "filename.png", width = 7, height = 7)
res <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))
#> Analysing 'ieu-a-2' on 'ieu-a-7'
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
ggsave(p1[[1]], file = "filename.png", width = 7, height = 7)

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_single <- mr_singlesnp(dat, all_method = c("mr_ivw", "mr_two_sample_ml"))
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]

mr_method_list()
mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))
mr_heterogeneity(dat)
mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw"))
mr_pleiotropy_test(dat)
res_single <- mr_singlesnp(dat)
res_single <- mr_singlesnp(dat, single_method = "mr_meta_fixed")
res_single <- mr_singlesnp(dat, all_method = "mr_two_sample_ml")
res_loo <- mr_leaveoneout(dat)

#correlations bw SNPs
snplist <- c("rs234", "rs1205")
ld_matrix(snplist)
dat <- harmonise_data(
  exposure_dat = bmi_exp_dat, 
  outcome_dat = chd_out_dat
)
dat2 <- dat_to_MRInput(dat)
MendelianRandomization::mr_ivw(dat2[[1]])
#Alternatively, convert to the MRInput format but also obtaining the LD matrix for the instruments
dat2 <- dat_to_MRInput(dat, get_correlation = TRUE)
MendelianRandomization::mr_ivw(dat2[[1]], correl = TRUE)

#http://app.mrbase.org/
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- extract_instruments(c('ebi-a-GCST011075'))
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ebi-a-GCST011075'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

plot(1:100, 1:100, main = "My Plot")
