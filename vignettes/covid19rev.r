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
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('finn-b-FG_CVD','finn-b-COPD_INSUFFICIENCY','finn-b-J10_CHEMGASRESP','finn-b-R18_SYMPTOMS_SIGNS_INVOLVI_CIRCULATO_RESPI_SYSTEMS','finn-b-ILD_INSUFFICIENCY','finn-b-J10_RESPIRATORY','finn-b-ASTHMA_ACUTE_RESPIRATORY_INFECTIONS','finn-b-J10_ARDS','finn-b-ST19_FOREI_BODY_RESPI_TRACT','finn-b-J10_ACUTELOWERNAS'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 0, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 3)
mr_results <- mr(dat)
p1 <- mr_scatter_plot(mr_results, dat)
length(p1)
p1[[1]]
ggsave(p1[[1]], file = "scatter.png", width = 7, height = 7)
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]], file = "forest.png", width = 7, height = 7)
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
ggsave(p3[[1]], file = "loo.png", width = 7, height = 7)
p4 <- mr_funnel_plot(res_single)
p4[[1]]
ggsave(p4[[1]], file = "funnel.png", width = 7, height = 7)



#https://github.com/rondolab/MR-PRESSO
library(MRPRESSO)
# Load a simulated toy dataset
data(SummaryStats)
# Run MR-PRESSO global method
mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)
# Run MR-PRESSO on a multi-variable MR (MMR) model specifying several exposures
mr_presso(BetaOutcome = "Y_effect", BetaExposure = c("E1_effect", "E2_effect"), SdOutcome = "Y_se", SdExposure = c("E1_se", "E2_se"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)


#replication https://github.com/animesh/mr-base-methods-paper 

#power analysis https://shiny.cnsgenomics.com/mRnd/

#other potential issues from https://www.nature.com/articles/s41598-023-39551-2#Sec2
#Selection for IVs
#The SNPs selected as IVs for BMR were obtained by running the function “extract_instruments” of the R package “TwoSampleMR”. According to previous studies, the IVs selected for exposure used in MR analysis should meet the criteria of p < 5 × E-08 and linkage disequilibrium (LD) r2 < 0.001 within 10 mb23. SNPs were extracted from the summary statistics of the five cardiovascular outcomes, whereas those associated with outcomes with genome-wide significance (p < 5 × E-08) were discarded. The palindromic SNPs with intermediate allele frequencies were then removed using the "harmonise data" function. Subsequently, SNPs were uploaded to the PhenoScanner database24,25 to find and remove those significantly associated with common cardiovascular risk factors (including BMI or obesity, diabetes, smoking, and blood pressure) to reduce the effect of confounding factors (Supplementary Tables S1). Then, the “run_mr_presso” function (set NbDistribution = 10,000, SignifThreshold = 0.05) was used for detecting horizontal pleiotropy26, and the outlier variants were removed when a significant horizontal pleiotropy was detected. Finally, the F statistics of the remaining SNPs used in UVMR were calculated using the following formula: F = R2 × (N–k–1)/[(1–R2) × k], where R2 = 2 × β2 × (1–EAF) × EAF, N is the sample size of BMR, k is the number of SNPs, β is the estimate of genetic effect on BMR, and EAF is the frequency of the effect allele. The SNPs with an F statistic greater than 10 were considered strong IVs of BMR27.
#Statistical analysis
#The inverse variance weighted (IVW) method was used as the primary analysis for the UVMR analysis, with MR-Egger, MR pleiotropy residual sum and outlier (MR-PRESSO), and weighted median (WM) model as supplemental ones. The IVW method can effectively combine the Wald ratio from multiple valid IVs into a single causal estimate for a more precise estimate of causal effect28. If heterogeneity is present, we employed the inverse-variance weighted (IVW) method with a random effects model for analysis, aiming to obtain more conservative outcomes29. When pleiotropy of IVs exists, MR-Egger is a valuable tool for estimating the causal effect30. The MR-PRESSO test can identify and remove any outlier variants, thus obtaining outlier-corrected MR analysis results, which are applicable when < 50% of the IVs are detected with horizontal pleiotropy26. The WM allows up to 50% of the IVs to be invalid and still produce consistent results31. MR-Egger and WM can provide more robust results but with less efficiency, while their consistency with the direction of IVW can increase confidence in estimating causal effects32. Given that five exposures and four analysis methods were included in this study, we consider a significant causal effect at the threshold of p < 0.0025 (Bonferroni correction p = 0.05/20) and nominally significant results at p < 0.05. Following UVMR estimates, we conducted a sensitivity analysis to detect potential heterogeneity and pleiotropy and identify any violations of the basic assumptions of MR. The presence of heterogeneity and pleiotropy can weaken and even invalidate the causal inference process and results. Therefore, we additionally conducted Cochran’s Q-test and MR-Egger intercept analysis. For Cochran’s Q-test, a p-value of < 0.05 indicates the presence of heterogeneity. Further, for MR-Egger intercept analysis, a p-value of < 0.05 signifies the presence of pleiotropy33. Besides, the leave-one-out (LOO) analysis can be used to test the reliability and robustness of the MR analysis. The consistency of LOO analysis and IVW results suggests that the MR analysis process and causal inference results are reliable and robust. Considering the potentially confounding or mediating effects of common cardiovascular risks on the results of MR analysis, we performed MVMR analysis, adjusting for BMI, T2DM, SI, SBP, and DBP, to obtain direct effects of BMR on CVDs as well as to assess possible mediating effects. Similarly, we used the criteria of p < 5 × E-08 and linkage disequilibrium (LD) r2 < 0.01 within 10 mb23 to filter IVs and subsequently removed those significantly associated with outcomes. For the MVMR analysis, we used IVW, WM, and MR-Egger methods as the main analysis approaches. The WM model was superior to the IVW model in detecting heterogeneity because it naturally explained heterogeneity by the bootstrapped variance. In addition, pleiotropy was evaluated by the Egger regression intercept, and a p < 0.05 indicated the presence of pleiotropy.


