rm(list = ls())
library(ggplot2)
library(ggpubr)
library(base)
library("deSolve")
library(binom)
library(stringr)
library('tidyr')
library('dplyr')
library(ggridges)
library(scales)



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

vec_3K <- c("ENT_ENT_erovirus.CoxB3.VP1","ENT_Norovirus.GII.4.VLP", "ENT_Norovirus.GII.6.VLP", "OTH_EBV.gp125", "OTH_Echovirus", "RES_ADE5", "RES_ADE40", "RES_Adenovirus.T3", "RES_Adenovirus.T5")
vec_3K_str <- c("Enterovirus CoxB3 VP1", "Norovirus GII4 VLP", "Norovirus GII6 VLP", "EBV gp125", "Echovirus", "Adenovirus 5", "Adenovirus 40", "Adenovirus T3", "Adenovirus T5")
vec_2K <- c("RES_RSV.gG","ENT_Rotavirus.VP7", "HCoV_229E_S", "HCoV_HKU1_S","HCoV_HKU1_NP", "HCoV_NL63_NP", "HCoV_NL63_S","HCoV_OC43_NP", "HCoV_OC43_S", "OTH_HAV")
vec_2k_str <- c("RSV gG","Rotavirus VP7", "HCoV 229E S", "HCoV HKU1 S","HCoV HKU1 NP", "HCoV NL63 NP", "HCoV NL63 S","HCoV OC43 NP", "HCoV OC43 S", "HAV")


vec_tot <- c(vec_entero, vec_resp, vec_HCov, vec_other)
vec_noModel <- vec_tot[!(vec_tot %in% vec_3K)]
vec_noModel <- vec_noModel[!(vec_noModel %in% vec_2K)]
vec_noModel_str <- c("Norovirus GII4 VP1", "Respiratory Syncytial A", "Respiratory Syncytial B", "Rhinovirus T1A lysate",
                     "HCoV 229E NP", "SARS-CoV-2 Spike Wuhan", "SARS-CoV-2 RBD Wuhan", "SARS-CoV-2 NP", "CMV", "HEV ORF2")


vec_epidemio <- c(colnames(PeruRAU_wide)[1:23])

# Select cohort from 2021
df_Peru <- PeruRAU_wide[PeruRAU_wide$todays_date> "2021-01-01",] 
df_Peru <- df_Peru[!is.na(df_Peru$Age_num),]
df_Peru <- df_Peru[order(df_Peru$Age_num),]

df_Peru <- as.data.frame(df_Peru[,c(vec_epidemio, vec_tot)])



# Data management for plotting results
df_Peru_u5 <- df_Peru[df_Peru$Age_num <= 5,]
df_Peru_AgeMid <- df_Peru[df_Peru$Age_num > 5 & df_Peru$Age_num <= 20,]
df_Peru_old <- df_Peru[df_Peru$Age_num > 20,]

df_Peru$Age_group <- "All"
df_Peru_u5$Age_group <- "Children under 5"
df_Peru_AgeMid$Age_group <- "Young adults (5-20)"
df_Peru_old$Age_group <- "Adults (>20)"

cohort <- rbind(df_Peru, df_Peru_u5, df_Peru_AgeMid, df_Peru_old)
cohort$Age_group <- factor(cohort$Age_group, levels=c("All", "Adults (>20)", "Young adults (5-20)", "Children under 5"))





###############
## FIGURE 3K ##
###############
plot_3k <- list()
# Loop for each antibondy response
for (Ag in 1:length(vec_3K)){
  Ag_str <- vec_3K[Ag]
  
  # download posterior distribution
  paramtempo = list.files(pattern = paste0("PostDist_", Ag_str))
  param_tempo <- read.csv(paramtempo)
  param_tempo$mu <- param_tempo$rho + param_tempo$lambda_bis
  
  
  # To plot data ________________________________________________________________________________
  age_bins     <- seq(from=0, to=95, by=2.5)
  age_bins_mid <- seq(from=1.25, to=93.75, by=2.5)
  N_bins <- length(age_bins) - 1
  
  cohort$Age_class <- NA
  for(k in 1:N_bins){
    index <- which( cohort$Age_num>age_bins[k] & cohort$Age_num<=age_bins[k+1] )
    cohort$Age_class[index]  <- age_bins_mid[k]}
  
  cohort$Age_class <- as.numeric(cohort$Age_class)
  cohort$Age_class[cohort$Age_num == 95] <- 93.75
  
  data_boxplot <- as.data.frame(cohort[, c(Ag_str, "Age_class")])
  data_boxplot[,1]<- log(as.numeric(data_boxplot[,1]))
  data_boxplot$data_obs_est <- "Observed"
  colnames(data_boxplot)[1] <- "V1"
  

  # To plot estimations
  age_seq <- seq(from=0, to=95, by=0.2)
   
  # Download age stratified probabilities to be sero positives or recovered (computed in 1.1.3.Results_Analysis)
  dfribbon2 = list.files(pattern = paste0(Ag_str, "_3k_pos_est.csv"))
  df_ribbon2 <- read.csv(dfribbon2)

  dfribbon3 = list.files(pattern = paste0(Ag_str, "_3k_recov_est.csv"))
  df_ribbon3 <- read.csv(dfribbon3)

  # GMM: simulate 100 participants of each ages
  rau_est <- as.data.frame(matrix(NA, nrow = 100*length(age_seq), ncol = 1))
  rau_est$age_seq <- rep(age_seq, each = 100)

  # GMM: simulate RAU value for all simulated participants
  for (i in 1: length(age_seq)){
    tempo <- rmultinom(n=1, size=100, prob=c(max((1 - df_ribbon2$V2[i] - df_ribbon3$V2[i]), 0.0000001),df_ribbon2$V2[i],df_ribbon3$V2[i]))

    rau_tempo_med<- c()
    rau_tempo_med <- c(rnorm(tempo[1], median(param_tempo$mu_N), median(param_tempo$sigma_N)),
                       rnorm(tempo[2], median(param_tempo$mu_N) + median(param_tempo$mu_A) + median(param_tempo$mu_P), median(param_tempo$sigma_P)),
                       rnorm(tempo[3], median(param_tempo$mu_N) + median(param_tempo$mu_A), median(param_tempo$sigma_A)))

    rau_est$V1[rau_est$age_seq == age_seq[i]] <- rau_tempo_med
  }

  # Bin by age class for plotting 
  rau_est$Age_class <- c(NA)
  for(k in 1:N_bins){
    index <- which( rau_est$age_seq>=age_bins[k] & rau_est$age_seq<age_bins[k+1] )
    rau_est$Age_class[index]  <- age_bins_mid[k]}

  rau_est$Age_class <- as.numeric(rau_est$Age_class)
  rau_est$Age_class[rau_est$age_seq == 95] <- 93.75


  rau_est$Age_class[rau_est$age_seq <1] <- 0

  # Build df with median and IC95% of simulated RAU
  rau_est_continuous <- as.data.frame(matrix(NA, nrow = length(unique(rau_est$Age_class))))
  rau_est_continuous$Age_num <- unique(rau_est$Age_class)
  for (i in 1: length(unique(rau_est$Age_class))){
    rau_est_continuous$V1[i] <- median(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]])
    rau_est_continuous$V1_low[i] <- quantile(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]], c(0.025))
    rau_est_continuous$V1_high[i] <- quantile(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]], c(0.975))
  }
  
  
  
  # Plot ________________________________________________________________________________________________________________________
  rau_est_continuous$Antigen <- 1
  rau_est_continuous$age_class_bis <- 1

  rau_est_continuous$V1_low[rau_est_continuous$V1_low<= log(9)] <- log(9)
  rau_est_continuous$V1_high[rau_est_continuous$V1_high>= log(250000)] <- log(250000)

  cohort$age_class_bis <- ifelse(cohort$Age_num <=1, "<= 1yo", "> 1yo")
  cohort_tempo <- cohort[, c("Age_num", Ag_str, "age_class_bis")]
  colnames(cohort_tempo)[2] <- "Antigen"

  if(str_sub(Ag_str, 1, 3) == "RES"){color_est <- "#009933"
                                     color_obs <- "#00CC66"}
  if(str_sub(Ag_str, 1, 3) == "OTH"){color_est <- "#FF9966"
                                     color_obs <- "#FFCC99"}
  if(str_sub(Ag_str, 1, 3) == "ENT"| 
            Ag_str=="OTH_Echovirus"| 
            Ag_str=="OTH_HAV"| 
            Ag_str=="RES_ADE40"){color_est <- "#990000"
                                     color_obs <- "#FF6666"}
  if(str_sub(Ag_str, 1, 3) == "HCo"){color_est <- "#6666FF"
                                     color_obs <- "#9999FF"}

  
  
  plot_3k[[Ag]] <- ggplot(data = cohort_tempo, aes(x = Age_num, y = Antigen)) + geom_point(aes( color = "Observed", fill = "Observed"), alpha = 0.4, size = 1)+   # change Ab
    geom_line(data = rau_est_continuous, mapping = aes(x = Age_num, y =exp(V1), color = "Estimated", fill = "Estimated"), linewidth = 2)+
    geom_ribbon(data = rau_est_continuous, mapping = aes(ymin = exp(V1_low), ymax = exp(V1_high), x = Age_num, fill = "Estimated"), alpha = 0.2)+
    scale_color_manual(name = "Data", values = c("Observed" = color_obs, "Estimated" = color_est)) +
    scale_fill_manual(name = "Data", values = c("Observed" = color_obs, "Estimated" = color_est)) +
    scale_y_log10(limits = c(9,250000), labels=trans_format('log10',math_format(10^.x)))+
    ggtitle(paste0(vec_3K_str[Ag]))+
    xlab("Age")+ylab("RAU")+theme_classic()+
    theme(plot.title = element_text(size=22, face = "bold"),
          axis.text=element_text(size=20),
          axis.title = element_text(size=22, face = "bold"),
          legend.text = element_text(size=24),
          legend.ticks = element_line(linewidth = 16),
          legend.title=element_blank())+
    guides(color = guide_legend(override.aes = list(size = 14)))
}





###############
## FIGURE 2K ##
###############
plot_2k <- list()
# Loop for each antibody response
for (Ag in 1:length(vec_2K)){
  Ag_str <- vec_2K[Ag]
  
  # download posterior distribution
  paramtempo = list.files(pattern = paste0("PostDist_", Ag_str))
  param_tempo <- read.csv(paramtempo)
  
  
  # To plot data ___________________________________________________________________________________
  age_bins     <- seq(from=0, to=95, by=2.5)
  age_bins_mid <- seq(from=1.25, to=93.75, by=2.5)
  N_bins <- length(age_bins) - 1
  
  cohort$Age_class <- NA
  for(k in 1:N_bins){
    index <- which( cohort$Age_num>age_bins[k] & cohort$Age_num<=age_bins[k+1] )
    cohort$Age_class[index]  <- age_bins_mid[k]}
  
  cohort$Age_class <- as.numeric(cohort$Age_class)
  cohort$Age_class[cohort$Age_num == 95] <- 93.75
  
  data_boxplot <- as.data.frame(cohort[, c(Ag_str, "Age_class")])
  data_boxplot[,1]<- log(as.numeric(data_boxplot[,1]))
  data_boxplot$data_obs_est <- "Observed"
  colnames(data_boxplot)[1] <- "V1"
  

  # To plot estimations ________________________________________________________________________________
   age_seq <- seq(from=0, to=95, by=0.2)
  
  # Download age stratified probabilities to be sero positives or recovered (computed in 1.1.3.Results_Analysis)
  dfribbon1 = list.files(pattern = paste0(Ag_str, "_2k_pos_est.csv"))
  df_ribbon1 <- read.csv(dfribbon1)

  
  # GMM: simulate 100 participants of each ages
  rau_est <- as.data.frame(matrix(NA, nrow = 100*length(age_seq), ncol = 1))
  rau_est$age_seq <- rep(age_seq, each = 100)
  
  # GMM: simulate RAU value for all simulated participants
  for (i in 1: length(age_seq)){
    tempo <- runif(100, 0, 1)
    rau_tempo_med<- c()
    for (j in 1: length(tempo)){
      # RAU sampled in neg or positive distribution depending on tpatients probability of being pos 
      rau_tempo_med[j] <- ifelse(tempo[j] < df_ribbon1$V2[i], rnorm(1, median(param_tempo$muPos), median(param_tempo$epsilon_P)), rnorm(1, median(param_tempo$mu_N), median(param_tempo$epsilon_N))) 
    }
    rau_est$V1[rau_est$age_seq == age_seq[i]] <- rau_tempo_med
  }
  
  
  # Bin by age class for plotting 
  rau_est$Age_class <- c(NA)
  for(k in 1:N_bins){
    index <- which( rau_est$age_seq>=age_bins[k] & rau_est$age_seq<age_bins[k+1] )
    rau_est$Age_class[index]  <- age_bins_mid[k]}
  
  rau_est$Age_class <- as.numeric(rau_est$Age_class)
  rau_est$Age_class[rau_est$age_seq == 95] <- 93.75
  
  
  # Build df with median and IC95% of simulated RAU
  rau_est$Age_class[rau_est$age_seq <1] <- 0
  
  rau_est_continuous <- as.data.frame(matrix(NA, nrow = length(unique(rau_est$Age_class))))
  rau_est_continuous$Age_num <- unique(rau_est$Age_class)
  for (i in 1: length(unique(rau_est$Age_class))){
    rau_est_continuous$V1[i] <- median(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]])
    rau_est_continuous$V1_low[i] <- quantile(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]], c(0.025))
    rau_est_continuous$V1_high[i] <- quantile(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]], c(0.975))
  }
  
  
  
  # Plot ___________________________________________________________________________________________________________________
  rau_est_continuous$Antigen <- 1
  rau_est_continuous$age_class_bis <- 1
  
  rau_est_continuous$V1_low[rau_est_continuous$V1_low<= 0] <- 0
  
  cohort$age_class_bis <- ifelse(cohort$Age_num <=1, "<=1 yo", "> 1yo")
  cohort_tempo <- cohort[, c("Age_num", Ag_str, "age_class_bis")]
  colnames(cohort_tempo)[2] <- "Antigen"
  
  if(str_sub(Ag_str, 1, 3) == "RES"){color_est <- "#009933"
                                     color_obs <- "#00CC66"}
  if(str_sub(Ag_str, 1, 3) == "OTH"){color_est <- "#FF9966"
                                     color_obs <- "#FFCC99"}
  if(str_sub(Ag_str, 1, 3) == "ENT"| 
            Ag_str=="OTH_Echovirus"| 
            Ag_str=="OTH_HAV"| 
            Ag_str=="RES_ADE40"){color_est <- "#990000"
                                     color_obs <- "#FF6666"}
  if(str_sub(Ag_str, 1, 3) == "HCo"){color_est <- "#6666FF"
                                     color_obs <- "#9999FF"}
  
  
  
  plot_2k[[Ag]] <- ggplot(data = cohort_tempo, aes(x = Age_num, y = Antigen)) + geom_point(aes( color = "Observed", fill = "Observed"), alpha = 0.4, size = 1)+   # change Ab
    geom_line(data = rau_est_continuous, mapping = aes(x = Age_num, y =exp(V1), color = "Estimated", fill = "Estimated"), linewidth = 2)+
    geom_ribbon(data = rau_est_continuous, mapping = aes(ymin = exp(V1_low), ymax = exp(V1_high), x = Age_num, fill = "Estimated"), alpha = 0.2)+
    scale_color_manual(name = "Data", values = c("Observed" = color_obs, "Estimated" = color_est)) +
    scale_fill_manual(name = "Data", values = c("Observed" = color_obs, "Estimated" = color_est)) +
    scale_y_log10(limits = c(9,250000), labels=trans_format('log10',math_format(10^.x)))+
    # scale_x_log10()+
    ggtitle(paste0(vec_2k_str[Ag]))+
    xlab("Age")+ylab("RAU")+theme_classic()+
    theme(plot.title = element_text(size=22, face = "bold"),
          axis.text=element_text(size=20),
          axis.title = element_text(size=22, face = "bold"),
          legend.text = element_text(size=24),
          legend.ticks = element_line(linewidth = 16),
          legend.title = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size = 14)))
}



# Get legend in black to represent all antigen categories
plot_leg1 <- ggplot(data = cohort_tempo, aes(x = Age_num, y = Antigen)) + geom_point(aes( color = "Observed", fill = "Observed"), alpha = 0.4, size = 1)+   # change Ab
  geom_line(data = rau_est_continuous, mapping = aes(x = Age_num, y =exp(V1), color = "Estimated", fill = "Estimated"), linewidth = 2)+
  geom_ribbon(data = rau_est_continuous, mapping = aes(ymin = exp(V1_low), ymax = exp(V1_high), x = Age_num, fill = "Estimated"), alpha = 0.2)+
  scale_color_manual(name = "Data", values = c("Observed" = "black", "Estimated" = "black")) +
  scale_fill_manual(name = "Data", values = c("Observed" = "black", "Estimated" = "black")) +
  scale_y_log10(limits = c(9,100000))+
  # scale_x_log10()+
  ggtitle(paste0(vec_2k_str[Ag]))+
  xlab("Age")+ylab("RAU")+theme_classic()+
  theme(plot.title = element_text(size=20, face = "bold"),
        axis.text=element_text(size=18),
        axis.title = element_text(size=20, face = "bold"),
        legend.text = element_text(size=24),
        legend.ticks = element_line(linewidth = 20),
        legend.title = element_text(size=22), 
        legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 16)))

leg1 <- get_legend(plot_leg1)





plot_noModel <- list()
# Loop for each antibody response
for (Ag in 1:length(vec_noModel)){
  Ag_str <- vec_noModel[Ag]
  
  
  # To plot data ___________________________________________________________________________________
  age_bins     <- seq(from=0, to=95, by=2.5)
  age_bins_mid <- seq(from=1.25, to=93.75, by=2.5)
  N_bins <- length(age_bins) - 1
  
  cohort$Age_class <- NA
  for(k in 1:N_bins){
    index <- which( cohort$Age_num>age_bins[k] & cohort$Age_num<=age_bins[k+1] )
    cohort$Age_class[index]  <- age_bins_mid[k]}
  
  cohort$Age_class <- as.numeric(cohort$Age_class)
  cohort$Age_class[cohort$Age_num == 95] <- 93.75
  
  data_boxplot <- as.data.frame(cohort[, c(Ag_str, "Age_class")])
  data_boxplot[,1]<- log(as.numeric(data_boxplot[,1]))
  data_boxplot$data_obs_est <- "Observed"
  colnames(data_boxplot)[1] <- "V1"
  
  
  
  # Plot ___________________________________________________________________________________________________________________
  cohort$age_class_bis <- ifelse(cohort$Age_num <=1, "<=1 yo", "> 1yo")
  cohort_tempo <- cohort[, c("Age_num", Ag_str, "age_class_bis")]
  colnames(cohort_tempo)[2] <- "Antigen"
  
  
  if(str_sub(Ag_str, 1, 3) == "RES"){color_obs <- "#00CC66"}
  if(str_sub(Ag_str, 1, 3) == "OTH"){color_obs <- "#FFCC99"}
  if(str_sub(Ag_str, 1, 3) == "ENT"){color_obs <- "#FF6666"}
  if(str_sub(Ag_str, 1, 3) == "HCo"){color_obs <- "#9999FF"}
  
  
  
  
  plot_noModel[[Ag]] <- ggplot(data = cohort_tempo, aes(x = Age_num, y = Antigen)) + geom_point(aes( color = "Observed", fill = "Observed"), alpha = 0.4, size = 1)+   # change Ab
    scale_color_manual(name = "Data", values = c("Observed" = color_obs)) +
    scale_fill_manual(name = "Data", values = c("Observed" = color_obs)) +
    scale_y_log10(limits = c(9,250000), labels=trans_format('log10',math_format(10^.x)))+
    # scale_x_log10()+
    ggtitle(paste0(vec_noModel_str[Ag]))+
    xlab("Age")+ylab("RAU")+theme_classic()+
    theme(plot.title = element_text(size=22, face = "bold"),
          axis.text=element_text(size=20),
          axis.title = element_text(size=22, face = "bold"),
          legend.text = element_text(size=24),
          legend.ticks = element_line(linewidth = 16),
          legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 14)))
}







# Plot leg 2 (To have all colors in the same legend)
GAM_obs = list.files(pattern = paste0("dfgam_sample_"))
list_dfGAM_obs<- lapply(GAM_obs, read.csv)

vec_tot <- c(vec_entero, vec_resp, vec_HCov, vec_other)

dfGAM_obs <- as.data.frame(matrix(ncol = 6, nrow = 0))
for(Ag in 1: length(GAM_obs)){
  dfGAM_tempo = as.data.frame(list_dfGAM_obs[[Ag]])
  dfGAM_tempo <- dfGAM_tempo[, -1]
  
  dfGAM_obs <- rbind(dfGAM_obs, dfGAM_tempo)
}
colnames(dfGAM_obs) <- c("Infection", "BaselineAB", "age", "period_time", 'ID', "Antigen")
dfGAM_obs <- dfGAM_obs[!is.na(dfGAM_obs$Infection),]


dfGAM_obs$Antigen <- factor(dfGAM_obs$Antigen, levels=rev(c(unique(dfGAM_obs$Antigen))))



library(stringr)
dfGAM_obs$antigen_class <- str_sub(dfGAM_obs$Antigen, 1, 3)
dfGAM_obs$AGE_cat <- NA
dfGAM_obs$AGE_cat <- ifelse(dfGAM_obs$age>=0 & dfGAM_obs$age<5, "0-4 yo", dfGAM_obs$AGE_cat)
dfGAM_obs$AGE_cat <- ifelse(dfGAM_obs$age>=5 & dfGAM_obs$age<18, "5-17 yo", dfGAM_obs$AGE_cat)
dfGAM_obs$AGE_cat <- ifelse(dfGAM_obs$age>=18 & dfGAM_obs$age<45, "18-44 yo", dfGAM_obs$AGE_cat)
dfGAM_obs$AGE_cat <- ifelse(dfGAM_obs$age>=45, "45+ yo", dfGAM_obs$AGE_cat)


# Pie 0 to 5 yo
df_pie_age1 <- dfGAM_obs[dfGAM_obs$AGE_cat =="0-4 yo",]
vec_model_single <-  c("ENT_ENT_erovirus.CoxB3.VP1", "ENT_Rotavirus.VP7","ENT_Norovirus.GII.4.VLP","ENT_Norovirus.GII.6.VLP", "HCoV_229E_S", "HCoV_HKU1_S",  "HCoV_NL63_S", 
                       "HCoV_OC43_S", "OTH_EBV.gp125", "OTH_Echovirus", "OTH_HAV", "RES_ADE40", "RES_ADE5", "RES_Adenovirus.T3", "RES_Adenovirus.T5", "RES_RSV.gG")  


df_pie1 <- as.data.frame(matrix(NA, nrow = length(vec_model_single), ncol = 3))
for(Ag in 1: length(vec_model_single)){
  df_tempo <- df_pie_age1[df_pie_age1$Antigen == vec_model_single[Ag],]
  prev_ag <- sum(df_tempo$Infection)/ sum(df_pie_age1$Infection)
  
  df_pie1[Ag,] <- c(prev_ag, vec_model_single[Ag], str_sub(vec_model_single[Ag], 1, 3))
}
colnames(df_pie1) <- c("Proportion",  "Antigen", "antigen cat")

df_pie1 <- df_pie1[order(df_pie1$Antigen), ]   # order rows by character variable
df_pie1$Antigen <- factor(df_pie1$Antigen, levels = rev(df_pie1$Antigen))

df_pie1 <- df_pie1[order(df_pie1$`antigen cat`), ]   # order rows by character variable
df_pie1$`antigen cat 2` <- ifelse(df_pie1$`antigen cat` == "RES", "Respiratory", NA )
df_pie1$`antigen cat 2` <- ifelse(df_pie1$`antigen cat` == "HCo", "Respiratory (Human Corona)", df_pie1$`antigen cat 2` )
df_pie1$`antigen cat 2` <- ifelse(df_pie1$`antigen cat` == "ENT", "Enteric", df_pie1$`antigen cat 2` )
df_pie1$`antigen cat 2` <- ifelse(df_pie1$`antigen cat` == "OTH", "Herpes", df_pie1$`antigen cat 2` )

df_pie1$`antigen cat 2` <- factor(df_pie1$`antigen cat 2`, levels = c('Respiratory', 'Respiratory (Human Corona)', 'Enteric', 'Herpes'))

plot_leg2 <- ggplot(df_pie1, aes(x = "", y = Proportion,  fill = `antigen cat 2`, group = `antigen cat 2`)) +
  geom_col(width = 2, alpha = 0.5, size = 1) +
  scale_fill_manual(name = "Pathogen", values = c("#00CC66",
                                                  "#6666FF",
                                                  "#FF6666", 
                                                  "#FFCC99"))+
  coord_polar("y") +
  theme_void() +
  ggtitle("0 to 4 years old")+ 
  theme(
    plot.title = element_text(size=18, face="bold", hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size=24),
    legend.ticks = element_line(linewidth = 20),
    legend.title = element_text(size=24))+
  guides(fill = guide_legend(override.aes = list(size = 16)))
plot_leg2
leg2 <- get_legend(plot_leg2)














# Build final plot
plot_group <- ggarrange(plot_3k[[6]]+ theme(legend.position = "none",axis.text.x =element_blank(), axis.title.x  = element_blank()), plot_3k[[8]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_3k[[9]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_2k[[1]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_2k[[8]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), 
                        plot_2k[[9]]+ theme(legend.position = "none",axis.text.x =element_blank(), axis.title.x  = element_blank()), plot_2k[[4]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_2k[[5]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_2k[[6]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_2k[[7]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), 
                        plot_2k[[3]]+ theme(legend.position = "none",axis.text.x =element_blank(), axis.title.x  = element_blank()),  plot_noModel[[5]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_noModel[[6]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_noModel[[7]]+ theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_noModel[[8]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), 
                        plot_3k[[4]]+theme(legend.position = "none",axis.text.x =element_blank(), axis.title.x  = element_blank()),   plot_noModel[[9]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_2k[[10]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_3k[[7]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), plot_3k[[5]]+theme(legend.position = "none",axis.text =element_blank(), axis.title  = element_blank()), 
                        plot_3k[[1]]+ theme(legend.position = "none"), plot_3k[[2]]+ theme(legend.position = "none",axis.text.y =element_blank(), axis.title.y  = element_blank()), plot_3k[[3]]+ theme(legend.position = "none",axis.text.y =element_blank(), axis.title.y  = element_blank()),  plot_2k[[2]]+ theme(legend.position = "none",axis.text.y =element_blank(), axis.title.y  = element_blank()), plot_noModel[[1]]+ theme(legend.position = "none",axis.text.y =element_blank(), axis.title.y  = element_blank()),
                        common.legend = F, ncol = 5, nrow = 5, widths = c(1, 0.9, 0.9, 0.9, 0.9))



plotlegend2 <- ggarrange("", leg1, "", leg2, "", ncol = 5, nrow = 1, widths = c(0.07, 0.2, 0.1,  0.33, 0.1))

plot_final <- ggarrange(plot_group, plotlegend2, ncol = 1, nrow = 2, heights = c(0.9, 0.1))


# Save plot
ggsave(filename = "Fig1.jpeg", plot_final, width=22, height=21)
























