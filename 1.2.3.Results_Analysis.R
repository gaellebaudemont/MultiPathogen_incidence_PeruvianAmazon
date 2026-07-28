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



#----------#
#-- Data --#
#----------#
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







# download results
results_2k = list.files(pattern = "_2k.csv")
list_results_2k <- lapply(results_2k, read.csv)



for (Ag in 1:length(results_2k)){
  # Antigens names and place in vec_tot
  Ag_str <- str_split(results_2k[[Ag]], "_")[[1]][3]
  if(Ag_str == "ENT"){ Ag_str <- str_split(results_2k[[Ag]], "_")[[1]][4] }
  if(Ag_str == "RES"){ Ag_str <- str_split(results_2k[[Ag]], "_")[[1]][4] }
  if(Ag_str == "HKU1" |
     Ag_str == "OC43" |
     Ag_str == "229E" |
     Ag_str == "NL63"){ Ag_str <- paste0(str_split(results_2k[[Ag]], "_")[[1]][3], "_", str_split(results_2k[[Ag]], "_")[[1]][4]) }
  
  Ag1 <- which(grepl(Ag_str, vec_tot))
  
  
  # Select results of interest and distinguish chains
  dfanalyses <- list_results_2k[[Ag]]
  dfanalyses$chain <- as.character(c(rep(c(1,2,3,4), each = 7000)))
  dfanalyses$iteration <- c(rep(seq(1,7000, 1), 4))


############################
## POSTERIOR DISTRIBUTION ##
############################
# Posterior distribution
plot_chains1 <- list()
for (param in 1: (dim(dfanalyses)[2] - 2)){
  param_tempo <- colnames(dfanalyses)[param]
  df_tempo <- dfanalyses[,c(param_tempo, "chain")]
  colnames(df_tempo) <- c("parametre", "chain")
  
  plot_chains1[[param]]<- ggplot(df_tempo) + geom_density(aes(x = parametre, group = chain, color = chain), linewidth = 1.5, alpha = 0.2) +
  scale_color_manual(values = c("#0033FF", "#0066FF", "#0099FF", "#0099CC"))+
  geom_vline(aes(xintercept = median(parametre)), color = "black", linetype = "dashed", size = 1.2)+
  geom_text(aes(x =  median(parametre)+0.1*median(parametre), label =  round(median(parametre), 3)), 
            y = 0.1)+
  theme_classic()+ ggtitle(param_tempo)+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16, face = "bold"),
        axis.title.x=element_blank())
  
}

# Chain mixing
plot_chains2 <- list()
for (param in 1: (dim(dfanalyses)[2] - 2)){
  param_tempo <- colnames(dfanalyses)[param]
  df_tempo <- dfanalyses[,c(param_tempo, "chain", "iteration")]
  colnames(df_tempo) <- c("parametre", "chain", "iteration")
  
  plot_chains2[[param]]<- ggplot(df_tempo) + geom_line(aes(x = iteration, y = parametre, group = chain, color = chain), size = 1.5, alpha = 0.2) +
    scale_color_manual(values = c("#0033FF", "#0066FF", "#0099FF", "#0099CC"))+
    theme_classic()+ 
    theme(plot.title = element_text(size=22, face = "bold"),
          axis.text=element_text(size=12),
          axis.title=element_text(size=16, face = "bold"),
          axis.title.x=element_blank())
  
}

plot_chains <- ggarrange(plot_chains1[[1]], plot_chains1[[2]], plot_chains1[[3]], plot_chains1[[4]], plot_chains1[[5]], plot_chains1[[6]], 
                  plot_chains2[[1]], plot_chains2[[2]], plot_chains2[[3]], plot_chains2[[4]], plot_chains2[[5]], plot_chains2[[6]],  
                  ncol = 6, nrow = 2, common.legend = TRUE)

# Save Plot
 ggsave(filename = paste0("Chains_", vec_tot[Ag1], "_.jpeg"), plot_chains, width=17, height=8)






##########################
### PLOT ESTIMATED RAU ###
##########################
# To plot data 
age_bins     <- seq(from=0, to=95, by=2.5)
age_bins_mid <- seq(from=1.25, to=93.75, by=2.5)
N_bins <- length(age_bins) - 1

cohort$Age_class <- NA
for(k in 1:N_bins){
  index <- which( cohort$Age_num>age_bins[k] & cohort$Age_num<=age_bins[k+1] )
  cohort$Age_class[index]  <- age_bins_mid[k]}

cohort$Age_class <- as.numeric(cohort$Age_class)
cohort$Age_class[cohort$Age_num == 95] <- 93.75

data_boxplot <- as.data.frame(cohort[, c(vec_tot[Ag1], "Age_class")])
data_boxplot[,1]<- log(as.numeric(data_boxplot[,1]))
data_boxplot$data_obs_est <- "Observed"
colnames(data_boxplot)[1] <- "V1"


# To plot estimation 
model_M1 = function(a, par_M1){
  lambda <- par_M1[1]
  rho    <- par_M1[2]
  
  Pos <- lambda/(lambda+rho)*(1-exp(-(lambda+rho)*a))
  
  c(Pos)}



age_seq <- seq(from=0, to=95, by=0.2)
par_median <- apply(X=dfanalyses[,1:2], MARGIN=2, FUN=median)
M1_predict_median <- as.data.frame(sapply(age_seq, model_M1, par=par_median))  



# Select 250 random values from the posterior distribution
N_sam = 250
sam_seq = round(seq(from=1, to=nrow(dfanalyses), length=N_sam))


# SCM: Compute model estimations on all ages for 250 parameters estimates from the posterior distribution
M1_predict1 = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
for(k in 1:N_sam){
  tempo <- as.data.frame(sapply(age_seq, model_M1, par = dfanalyses[sam_seq[k],1:2]))
  M1_predict1[k,] = tempo[1,]
}

# SCM: Get the median and CI95%
M1_quant1 = as.data.frame(matrix(NA, nrow=3, ncol=length(age_seq)))
for(j in 1:length(age_seq)){
  M1_quant1[,j] = quantile(as.data.frame(M1_predict1[,j]), prob=c(0.025, 0.5, 0.975), na.rm = T )}

# SCM: Build df with median and CI95%
df_ribbon1 <- as.data.frame(t(M1_quant1))
df_ribbon1$age <- age_seq



# GMM: simulate 100 participants of each ages
rau_est <- as.data.frame(matrix(NA, nrow = 100*length(age_seq), ncol = 1))
rau_est$age_seq <- rep(age_seq, each = 100)

# GMM: simulate RAU value for all simulated participants
for (i in 1: length(age_seq)){
  tempo <- runif(100, 0, 1)
  rau_tempo_med<- c()
  for (j in 1: length(tempo)){
    # RAU sampled in neg or positive distribution depending the probability of being positive at a given age 
    rau_tempo_med[j] <- ifelse(tempo[j] < df_ribbon1$V2[i], rnorm(1, median(dfanalyses$muPos), median(dfanalyses$epsilon_P)), rnorm(1, median(dfanalyses$mu_N), median(dfanalyses$epsilon_N))) 
  }
  rau_est$V1[rau_est$age_seq == age_seq[i]] <- rau_tempo_med
}


# Bin by age class for plotting (boxplot)
rau_est$Age_class <- c(NA)
for(k in 1:N_bins){
  index <- which( rau_est$age_seq>=age_bins[k] & rau_est$age_seq<age_bins[k+1] )
  rau_est$Age_class[index]  <- age_bins_mid[k]}

rau_est$Age_class <- as.numeric(rau_est$Age_class)
rau_est$Age_class[rau_est$age_seq == 95] <- 93.75




#################### PLOT BOXPLOT ###########
rau_est$data_obs_est <- "Estimated"

data_boxplot <- rbind(data_boxplot, rau_est[, c("V1", 'Age_class', 'data_obs_est')])
data_boxplot <- data_boxplot[!is.na(data_boxplot$Age_class),]
data_boxplot$V1[data_boxplot$V1 <= 0] <- 0

# Plot
plot_boxplo <- ggplot(data = data_boxplot, aes(x = factor(Age_class, levels = age_bins_mid), y = V1, colour = data_obs_est, fill = data_obs_est)) + 
  geom_boxplot(aes(colour = data_obs_est, fill = data_obs_est), alpha = 0.6)+
  scale_fill_manual("Data", values = c("black", "#990000")) +
  scale_color_manual("Data", values = c("black", "#990000")) +
  ggtitle("B.")+
  xlab("Age")+ylab("log(RAU)")+theme_classic()+
  theme(plot.title = element_text(size=20, face = "bold"),
        axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title = element_text(size=18, face = "bold"))
plot_boxplo







#################### PLOT CONTINUOUS ###########
rau_est$Age_class[rau_est$age_seq <1] <- 0

rau_est_continuous <- as.data.frame(matrix(NA, nrow = length(unique(rau_est$Age_class))))
rau_est_continuous$Age_num <- unique(rau_est$Age_class)
for (i in 1: length(unique(rau_est$Age_class))){
  rau_est_continuous$V1[i] <- median(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]])
  rau_est_continuous$V1_low[i] <- quantile(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]], c(0.025))
  rau_est_continuous$V1_high[i] <- quantile(rau_est$V1[rau_est$Age_class == rau_est_continuous$Age_num[i]], c(0.975))
}

rau_est_continuous$Antigen <- 1
rau_est_continuous$age_class_bis <- 1

rau_est_continuous$V1_low[rau_est_continuous$V1_low<= 0] <- 0

cohort$age_class_bis <- ifelse(cohort$Age_num <=1, "<=1 yo", "> 1yo")
cohort_tempo <- cohort[, c("Age_num", vec_tot[Ag1], "age_class_bis")]
colnames(cohort_tempo)[2] <- "Antigen"

plot_continuous <- ggplot(data = cohort_tempo, aes(x = Age_num, y = log(Antigen), group = age_class_bis)) + geom_point(aes(alpha = 0.6, color = age_class_bis))+   # change Ab
  scale_fill_manual(values = c("orange", "#CC3333")) +
  scale_color_manual("Data", values = c("orange", "#CC3333")) +
  geom_line(data = rau_est_continuous, mapping = aes(x = Age_num, y =V1 ), linewidth = 1)+
  geom_ribbon(data = rau_est_continuous, mapping = aes(ymin = V1_low, ymax = V1_high, x = Age_num), alpha = 0.3)+
  ggtitle("C.")+
  xlab("Age")+ylab("log(RAU)")+theme_classic()+
  theme(plot.title = element_text(size=20, face = "bold"),
        axis.text=element_text(size=16),
        axis.title = element_text(size=18, face = "bold"),
        legend.position = "none")
plot_continuous




#####################################################
### PLOT ESTIMATED PROBABILITY TO BE SEROPOSITIVE ###
#####################################################
plot_seroprev <- ggplot(data = df_ribbon1, aes(x = age, y = V2)) + geom_line(color = "red", linewidth = 1)+
  geom_ribbon(aes(x = age, ymin = V1, ymax = V3, alpha = 0.5), color = "red", fill = "red")+
  ggtitle("D.")+
  xlab("Age")+ylab("Seroprevalence")+theme_classic()+
  theme(plot.title = element_text(size=20, face = "bold"),
        axis.text=element_text(size=16),
        axis.title = element_text(size=18, face = "bold"),
        legend.position = "none")
plot_seroprev





###################################
### PLOT DISTRIBUTION / DENSITY ###
###################################
# Compute estimated distributions with the median of the posterior distribution
logRAU <- c( seq(min(log(cohort[, c(vec_tot[Ag1])]), na.rm = TRUE), max(log(cohort[, c(vec_tot[Ag1])]), na.rm = TRUE), 0.1))
line_1_estimated <- c((1/median(dfanalyses$epsilon_N*sqrt(2*pi)))*exp(-(1/2)*((logRAU - median(dfanalyses$mu_N))^2)/(median(dfanalyses$epsilon_N))^2))
line_2_estimated <- c((1/median(dfanalyses$epsilon_P*sqrt(2*pi)))*exp(-(1/2)*((logRAU - (median(dfanalyses$muPos)))^2)/(median(dfanalyses$epsilon_P))^2))


# Build df with estimated RAU distribution
estimated_density <- as.data.frame(cbind(logRAU, line_1_estimated, line_2_estimated))
estimated_density$Age_group <- "Estimated distribution"


# Plot estimated RAU distribution compared to overall RAU distribution and RAU distribution from different age categories 
library(ggridges)
dens_plot <- ggplot(data = cohort, aes(x = log(cohort[,vec_tot[Ag1]]), y = Age_group, group = Age_group)) + geom_density_ridges(aes(color = Age_group, fill = Age_group), scale = 2, alpha = 0.5)+
  ggtitle(paste0("A.", vec_tot[Ag1]))+
  scale_y_discrete(limits = c("Children under 5", "Young adults (5-20)", "Adults (>20)", "All"))+
  scale_fill_manual(values = c("#990000", "#CC3333", "#FF6666", "#FF9999")) +
  scale_color_manual(values = c("#990000", "#CC3333", "#FF6666", "#FF9999")) +
  geom_vline(mapping = aes(xintercept = median(dfanalyses$mu_N)), color = "darkgreen", linewidth = 1)+
  geom_vline(mapping = aes(xintercept = median(dfanalyses$mu_N) + median(dfanalyses$mu_P)), color = "red", linewidth = 1)+
  geom_line(data = estimated_density, mapping = aes(y = line_1_estimated, x = logRAU), color = "darkgreen")+
  geom_line(data = estimated_density, mapping = aes(y = line_2_estimated, x = logRAU), color = "red")+
  xlab("log(RAU)")+ylab("density")+theme_classic()+
  theme(plot.title = element_text(size=20, face = "bold"),
        axis.text=element_text(size=16),
        axis.title = element_text(size=18, face = "bold"),
        legend.position = "none")


dens_plot





# Build final plot
plot_tempo <- ggarrange(plot_boxplo, plot_continuous, plot_seroprev, ncol = 1, heights = c(1, 1, 0.6))
plot_final <- ggarrange(dens_plot, plot_tempo, ncol = 2, widths = c(0.75, 1))
plot_final


# Save plot
ggsave(filename = paste0("VPC_", vec_tot[Ag1], "_.jpeg"), plot_final, width=14, height=9)



}









