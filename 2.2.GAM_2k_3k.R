rm(list=ls())
library(mgcv)
library(binom)
library(ggplot2)
require(mgcv)
library(ggpubr)
library( broom )
library(dplyr)
library(tidyr)

# Download Data and data management
Peru <- read.csv(file = "0.Peru_MFI_RAU_nMFI.csv")
Peru <- Peru[,-1]
PeruRAU <- Peru[, -c(25, 27)] # Work with RAU
PeruRAU_wide <- PeruRAU %>% pivot_wider(names_from = name, values_from = RAU_value)    
PeruRAU_wide<- PeruRAU_wide %>% 
  unnest (cols = c(RES_Adenovirus.T3, RES_Adenovirus.T5, OTH_CMV, VPD_Diptheria.Toxin, OTH_EBV.gp125, OTH_Echovirus, ENT_ENT_erovirus.CoxB3.VP1, OTH_HBeAg, OTH_HCcAg, OTH_HEV.ORF2, VPD_Measles.lysate,
                   VPD_Mumps.lysate, ENT_Norovirus.GII.4.VLP, ENT_Norovirus.GII.4.VP1, ENT_Norovirus.GII.6.VLP, RES_RES_piratory.Syncytial.A, RES_RES_piratory.Syncytial.B, RES_RSV.gG, RES_Rhinovirus.T1A.lysate, OTH_HBsAg,
                   VPD_Tetanus.Toxin, OTH_Varicella.Zoster.Virus, OTH_HAV, OTH_Toxoplasma.gondii.SAG2, LEP_LA0591, HCoV_HKU1_NP, OTH_HPV.16, OTH_HPV.18, LEP_LA3490, HCoV_OC43_NP, HCoV_OC43_S, LEP_RBLs, HCoV_HKU1_S,
                   HCoV_229E_NP, HCoV_229E_S, HCoV_NL63_NP, HCoV_NL63_S, RES_FluA, RES_ADE5, RES_ADE40, VPD_mumps_NP, VPD_measles_NP, VPD_Rub..VLP., HCoV_SARS.CoV.2.Spike.Wuhan, ENT_Rotavirus.VP7,
                   HCoV_SARS.CoV.2.RBD.Wuhan, HCoV_SARS.CoV.2.Spike.Omicron, HCoV_SARS.CoV.2.RBD.Omicron, HCoV_SARS.CoV.2.NP, HCoV_SARS.CoV.2.Spike.S2, HCoV_ME.fusion, VPD_Tet.Tox, ARB_DENV1_DIII, ARB_DENV2_DIII,
                   ARB_DENV3_DIII, ARB_DENV4_DIII, ARB_USUV_NS1, ARB_CHIKV_VLP, ARB_WNV_NS1, ARB_MAYV_E2, ARB_YFV_E, ARB_YFV_NS1, ARB_ZIKV_NS1, ARB_ONNV_VLP, ARB_OROV_GN, ARB_ZIKV_VLP, ARB_DENV1.NS1, ARB_DENV2.NS1,
                   LEP_Ag.LEP_to.CNR, ARB_DENV3.NS1, ARB_DENV4.NS1, ARB_DENV1.VLP, ARB_ZIKVSU.NS1, ARB_CHIKV_E2, ARB_DENV2.VLP, ARB_DENV3.VLP, ARB_DENV4.VLP, ARB_CHIKV_NSP123, ARB_MAYV_NSP1, ARB_WNV_DIII, ARB_ZIKVAS_DIII,
                   ARB_OROV_GC, ARB_JEV_NS1, ARB_JEV_E, MAL_Etramp5.Ag1, MAL_GLURP.R2, MAL_PfAMA1, MAL_PfMSP119, MAL_GEXP18, MAL_Etramp4.Ag2, MAL_MSP2.CH150, MAL_Hyp2, MAL_CSP, MAL_PfMSP1, MAL_PvMSP1, MAL_RAMA,
                   MAL_L02.Pv.fam.a, MAL_DBPSal1, MAL_L34..MSP8, MAL_L01.MSP1.19, MAL_EBPII, MAL_RBP2b.L31, MAL_L54...MSP3.a, MAL_GST, OTH_E.Coli.lysat, MAL_PvMSP8))


vec_entero <- c("ENT_Norovirus.GII.4.VLP", "ENT_Norovirus.GII.4.VP1", "ENT_Norovirus.GII.6.VLP", "ENT_ENT_erovirus.CoxB3.VP1", "ENT_Rotavirus.VP7")
vec_resp <- c("RES_Adenovirus.T3", "RES_Adenovirus.T5","RES_RES_piratory.Syncytial.A","RES_RES_piratory.Syncytial.B", "RES_RSV.gG", "RES_ADE5", "RES_ADE40", "RES_Rhinovirus.T1A.lysate")
vec_HCov <- c("HCoV_HKU1_S","HCoV_HKU1_NP", "HCoV_229E_NP", "HCoV_229E_S", "HCoV_NL63_NP", "HCoV_NL63_S","HCoV_OC43_NP", "HCoV_OC43_S", "HCoV_SARS.CoV.2.Spike.Wuhan", 
              "HCoV_SARS.CoV.2.RBD.Wuhan","HCoV_SARS.CoV.2.NP")
vec_other <- c("OTH_CMV","OTH_EBV.gp125", "OTH_Echovirus", "OTH_HEV.ORF2", "OTH_HAV")

vec_tot <- c(vec_entero, vec_resp, vec_HCov, vec_other)
vec_epidemio <- c(colnames(PeruRAU_wide)[1:23])

# Select cohort from 2021
df_Peru <- PeruRAU_wide[PeruRAU_wide$todays_date> "2021-01-01",] 
df_Peru <- df_Peru[!is.na(df_Peru$Age_num),]
df_Peru <- df_Peru[order(df_Peru$Age_num),]

df_Peru <- as.data.frame(df_Peru[,c(vec_epidemio, vec_tot)])


# Identify antibody responses analysed with GMMxSCM 3k, 2k or for which no model were fitted - estimation dependant on the age stratified seroprevalence estimated by SCMxGMM
vec_3K <- c("ENT_Norovirus.GII.4.VLP", "ENT_Norovirus.GII.6.VLP", "ENT_ENT_erovirus.CoxB3.VP1", "OTH_EBV.gp125", "OTH_Echovirus", "RES_ADE5", "RES_ADE40", "RES_Adenovirus.T3", "RES_Adenovirus.T5")
vec_2K <- c("RES_RSV.gG","ENT_Rotavirus.VP7", "HCoV_229E_S", "HCoV_HKU1_S","HCoV_HKU1_NP", "HCoV_NL63_NP", "HCoV_NL63_S","HCoV_OC43_NP", "HCoV_OC43_S", "OTH_HAV")
vec_noModel_BA <- c("ENT_Norovirus.GII.4.VP1", "HCoV_229E_NP", "HCoV_SARS.CoV.2.Spike.Wuhan", "HCoV_SARS.CoV.2.RBD.Wuhan", "HCoV_SARS.CoV.2.NP", "OTH_CMV")


# Loop on antibody response 
for(Ag in 1: length(vec_3K)){
  # download dfgam
  dfGAM_ag = list.files(pattern = paste0("dfgam_", vec_3K[Ag], "_.csv"))
  list_dfGAM_ag <- lapply(dfGAM_ag, read.csv)
  
  dfGAM_tempo <- as.data.frame(list_dfGAM_ag[[1]])
  dfGAM_tempo <- dfGAM_tempo[-1, -1] # delete empty first line and row

  # data management
  dfGAM_tempo$basalAB <- log(dfGAM_tempo$basalAB)
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$age),]
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$basalAB),]
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$INF),]
  colnames(dfGAM_tempo)[4] <- c('period_time')
  
  
  # Delete L-LIB 0391: data too noisy
  dfGAM_tempo <- dfGAM_tempo[dfGAM_tempo$ID != "L-LIB0391",]
  
  
  
  ###########
  ## MODEL ##
  ###########
  b <- gam(INF~ offset(log(period_time)) + s(basalAB) + s(age), family = binomial(), data=dfGAM_tempo, method = 'REML')
  print(vec_3K[Ag])
  print(summary(b))

  plot(b, pages = 1)

  # Prediction for Figure 4___________________________________________________________________________________________________
  paramtempo = list.files(pattern = paste0("PostDist_", vec_3K[Ag], "_3k_"))
  param_tempo <- read.csv(paramtempo)
  param_tempo$mu <- param_tempo$rho + param_tempo$lambda_bis

  age_seq <- seq(from=0, to=95, by=0.2)

  # Download probabilities of being positive and recovered at each age estimated by the SCMxGMM (df build and save in 1.1.3.Results_Analysis)
  dfribbon2 = list.files(pattern = paste0(vec_3K[Ag], "_3k_pos_est.csv"))
  df_ribbon2 <- read.csv(dfribbon2)

  dfribbon3 = list.files(pattern = paste0(vec_3K[Ag], "_3k_recov_est.csv"))
  df_ribbon3 <- read.csv(dfribbon3)


  # RAU simulation for each age
  rau_est <- as.data.frame(matrix(NA, nrow = length(age_seq), ncol = 1))
  rau_est$age_seq <- rep(age_seq, each = 1)
  for (i in 1: length(age_seq)){
    rau_tempo_med<- c()
    rau_tempo_med <- c(median(param_tempo$mu_N)*max((1 - df_ribbon2$V2[i] - df_ribbon3$V2[i]), 0.0000001) +
                        (median(param_tempo$mu_N) + median(param_tempo$mu_R) + median(param_tempo$mu_P))*df_ribbon2$V2[i] +
                       (median(param_tempo$mu_N) + median(param_tempo$mu_R))*df_ribbon3$V2[i])

    rau_est$V1[rau_est$age_seq == age_seq[i]] <- rau_tempo_med
  }



  # Creation of df with simulated RAU for each age
  colnames(rau_est) <- c('basalAB', 'age')
  rau_est$period_time <- median(dfGAM_tempo$period_time)


  # incidence estimation for each age
  pred_tempo <- predict(b, rau_est, type = "link", se.fit = TRUE)


  rau_est$fit_fig4 <- exp(pred_tempo$fit)
  rau_est$lower_fig4 <- exp(pred_tempo$fit - 1.96 * pred_tempo$se.fit)
  rau_est$upper_fig4 <- exp(pred_tempo$fit + 1.96 * pred_tempo$se.fit)

  # save estimations
  write.csv(rau_est, paste0("Fig4_GAMpred_", vec_3K[Ag], "_.csv"))


  # _________________________________________________________________________________________________________________________________________________________________

  
  dfanalyses <- 0
  dfGAM_ag <- 0
  
  }




# Loop on antibody responses
for(Ag in 1: length(vec_noModel_BA)){
  # download dfgam
  dfGAM_ag = list.files(pattern = paste0("dfgam_", vec_noModel_BA[Ag], "_.csv"))
  list_dfGAM_ag <- lapply(dfGAM_ag, read.csv)
  
  dfGAM_tempo <- as.data.frame(list_dfGAM_ag[[1]])
  dfGAM_tempo <- dfGAM_tempo[-1, -1]                       # delete empty first line and row
  
  
  # data management
  dfGAM_tempo$basalAB <- log(dfGAM_tempo$basalAB)
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$age),]
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$basalAB),]
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$INF),]
  colnames(dfGAM_tempo)[4] <- c('period_time')
  
  
  # Delete L-LIB 0081: Participant with too much noise
  dfGAM_tempo <- dfGAM_tempo[dfGAM_tempo$ID != "L-LIB0391",]
  
  
  
  ###########
  ## MODEL ##
  ###########
  b <- gam(INF~ offset(log(period_time)) + s(basalAB) + s(age), family = binomial(), data=dfGAM_tempo, method = 'REML')
  print(vec_noModel_BA[Ag])
  print(summary(b))
  
  
  
  # Prediction for Figure 4___________________________________________________________________________________________________
  age_seq <- seq(from=0, to=95, by=0.2)
  df_tempo <- as.data.frame(df_Peru[,c(vec_epidemio, vec_noModel_BA[Ag])])
  
  
  # Creation of df with median RAU for each age
  rau_est <- as.data.frame(matrix(NA, nrow = length(age_seq), ncol = 2))
  colnames(rau_est) <- c('basalAB', 'age')
  rau_est$age <- age_seq
  rau_est$basalAB <- log(median(df_tempo[,c(vec_noModel_BA[Ag])], na.rm = TRUE))   # Model with no age pattern: the median RAU titre for each age is the overall median RAU titre
  rau_est$period_time <- median(dfGAM_tempo$period_time)
  
  # incidence estimation for each age
  pred_tempo <- predict(b, rau_est, type = "link", se.fit = TRUE)
  
  
  rau_est$fit_fig4 <- exp(pred_tempo$fit)
  rau_est$lower_fig4 <- exp(pred_tempo$fit - 1.96 * pred_tempo$se.fit)
  rau_est$upper_fig4 <- exp(pred_tempo$fit + 1.96 * pred_tempo$se.fit)
  
  
  # save estimations
  write.csv(rau_est, paste0("Fig4_GAMpred_", vec_noModel_BA[Ag], "_.csv"))
  
  
  # _________________________________________________________________________________________________________________________________________________________________
  
  
  
  dfanalyses <- 0
  dfGAM_ag <- 0
  
}









for(Ag in 1: length(vec_2K)){
  # download dfgam
  dfGAM_ag = list.files(pattern = paste0("dfgam_", vec_2K[Ag], "_.csv"))
  list_dfGAM_ag <- lapply(dfGAM_ag, read.csv)
  
  dfGAM_tempo <- as.data.frame(list_dfGAM_ag[[1]])
  dfGAM_tempo <- dfGAM_tempo[-1, -1]                      # delete empty first line and row
  
  
  # data management
  dfGAM_tempo$basalAB <- log(dfGAM_tempo$basalAB)
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$age),]
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$basalAB),]
  dfGAM_tempo <- dfGAM_tempo[!is.na(dfGAM_tempo$INF),]
  colnames(dfGAM_tempo)[4] <- c('period_time')
  
  
  # Delete L-LIB 0081: Participant with too much noise
  dfGAM_tempo <- dfGAM_tempo[dfGAM_tempo$ID != "L-LIB0391",]
  
  
  
  ###########
  ## MODEL ##
  ###########
  b <- gam(INF~ offset(log(period_time)) + s(basalAB)+ s(age), family = binomial(), data=dfGAM_tempo, method = 'REML')
  print(vec_2K[Ag])
  print(summary(b))



  # Prediction for Figure 4___________________________________________________________________________________________________
  paramtempo = list.files(pattern = paste0("PostDist_", vec_2K[Ag], "_2k_"))
  param_tempo <- read.csv(paramtempo)

  age_seq <- seq(from=0, to=95, by=0.2)

  # Download probabilities of being positive at each age estimated by the SCMxGMM (df build and save in 1.2.3.Results_Analysis)
  dfribbon2 = list.files(pattern = paste0(vec_2K[Ag], "_2k_pos_est.csv"))
  df_ribbon2 <- read.csv(dfribbon2)

  
  # RAU simulation for each age
  rau_est <- as.data.frame(matrix(NA, nrow = length(age_seq), ncol = 1))
  rau_est$age_seq <- rep(age_seq, each = 1)
  for (i in 1: length(age_seq)){
    rau_tempo_med<- c()
    rau_tempo_med <- c(median(param_tempo$mu_N)*max((1 - df_ribbon2$V2[i]), 0.0000001) +
                         (median(param_tempo$mu_N) + median(param_tempo$mu_P))*df_ribbon2$V2[i])

    rau_est$V1[rau_est$age_seq == age_seq[i]] <- rau_tempo_med
  }


  # Creation of df with median simulated RAU for each age
  colnames(rau_est) <- c('basalAB', 'age')
  rau_est$period_time <- median(dfGAM_tempo$period_time)

  
  # incidence estimation for each age
  pred_tempo <- predict(b, rau_est, type = "link", se.fit = TRUE)


  rau_est$fit_fig4 <- exp(pred_tempo$fit)
  rau_est$lower_fig4 <- exp(pred_tempo$fit - 1.96 * pred_tempo$se.fit)
  rau_est$upper_fig4 <- exp(pred_tempo$fit + 1.96 * pred_tempo$se.fit)

  # save estimations
   write.csv(rau_est, paste0("Fig4_GAMpred_", vec_2K[Ag], "_.csv"))
  # _________________________________________________________________________________________________________________________________________________________________




  dfanalyses <- 0
  dfGAM_ag <- 0
  
}





