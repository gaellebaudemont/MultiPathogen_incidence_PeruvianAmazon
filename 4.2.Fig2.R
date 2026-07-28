rm(list=ls())
library('tidyr')
library('dplyr')
library('plotrix')
require(mgcv)
library(ggplot2)
library(ggpubr)
Sys.setlocale("LC_TIME", "English")

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
vec_tot_16 <- c("ENT_ENT_erovirus.CoxB3.VP1", "OTH_EBV.gp125", "OTH_Echovirus", "RES_ADE5", "RES_ADE40", "RES_Adenovirus.T3", "RES_Adenovirus.T5", "RES_RSV.gG",
                "ENT_Rotavirus.VP7", "HCoV_229E_S", "HCoV_HKU1_S","HCoV_HKU1_NP", "HCoV_NL63_NP", "HCoV_NL63_S","HCoV_OC43_NP", "HCoV_OC43_S")
vec_tot_mini <- c("OTH_EBV.gp125", "RES_ADE5", "RES_RSV.gG",
                "HCoV_229E_S")
vec_epidemio <- c(colnames(PeruRAU_wide)[1:23])


# Select cohort from 2021
df_Peru <- PeruRAU_wide[PeruRAU_wide$todays_date> "2021-01-01",] 
df_Peru <- df_Peru[!is.na(df_Peru$Age_num),]
df_Peru <- df_Peru[order(df_Peru$Age_num),]

df_Peru <- as.data.frame(df_Peru[,c(vec_epidemio, vec_tot)])

df_Peru <- df_Peru[, c(vec_epidemio, vec_tot_16)]






########################################
### PLOT INDIVIDUAL ANTIBODY KINETIC ###
########################################

## DATA MANAGEMENT:DELETE BLANKS ________________________________________________________
# (same data management than the one for building the GAM df)


#ID of participants 
vec_patient_ID_tempo <- unique(df_Peru$unique_id_village)
vec_patient_ID_tempo <- vec_patient_ID_tempo[!is.na(vec_patient_ID_tempo)]
  
# ID of participants with several samples in this age category
vec_patient_tempo_several_samples <- vector()
for (patient in 1 : length(vec_patient_ID_tempo)){
    tempo <- df_Peru[df_Peru$unique_id_village == vec_patient_ID_tempo[patient], c("todays_date", "house_id", "Age_num", "sex", "unique_id_village", "Comunidad", "overall_pcr_result")]
    vec_patient_tempo_several_samples[patient] <- ifelse(dim(tempo)[1] > 1, vec_patient_ID_tempo[patient], NA)
}
vec_patient_tempo_several_samples <- vec_patient_tempo_several_samples[!is.na(vec_patient_tempo_several_samples)]
cohort_several_samples_Notempo <- df_Peru[df_Peru$unique_id_village %in% vec_patient_tempo_several_samples,]


# Loop for each antibody response
for (Ab in 1:length(vec_tot_16)){
  cohort_several_samples <- cohort_several_samples_Notempo
  cohort_several_samples$Delete <- 0    # To identify failed samples
  cohort_several_samples$INF <- 0       # To identify potential infection (4-fold increase)
  

  # Infection definition
  for (patient in 1 : length(vec_patient_tempo_several_samples)){
    tempo <- df_Peru[df_Peru$unique_id_village == vec_patient_tempo_several_samples[patient], c("todays_date", "house_id", "Age_num", "sex", vec_tot_16[Ab], "unique_id_village", "Comunidad", "overall_pcr_result")]
    tempo <- tempo[order(tempo$todays_date),]
    tempo <- tempo[!is.na(tempo[,vec_tot_16[Ab]]),]
    
    
    if (dim(tempo)[1] >1){
      tempo$delete <- NA         # Samples identified as a failed sample, to be deleted 
      tempo$drop_date <- NA      # Samples showing a 4 fold drop 
      tempo$medianRAU <- NA      # Median of all previous samples
      tempo$inf <- NA            # Samples Identified as an infection (4-fold boost)
      
      if(dim(tempo)[1] >2){
        for (i in 1:(dim(tempo)[1] - 1)){
          tempo$drop_date[i+1] <- ifelse(tempo[i+1, vec_tot_16[Ab]]< tempo[i, vec_tot_16[Ab]]/4, 1, 0) # Marked as 1 if 4-fold drop is observed compared to previous sample
        }
      }
      
      
      for (i in 1:dim(tempo)[1]){
        tempo$medianRAU[i] <- median(tempo[c(1:i), vec_tot_16[Ab]]) # Median of all previous samples
      }
      
      tempo$inf[1] <- ifelse(tempo[c(1), vec_tot_16[Ab]]> 4*tempo[dim(tempo)[1], "medianRAU"] &                    # 4fold increase compared to the median of all samples
                               tempo[c(1), vec_tot_16[Ab]]> 4*tempo[c(2), vec_tot_16[Ab]] &                        # 4fold increase compared to previous sample (only if the previous is not a failed sample)
                               tempo[c(2), "drop_date"] !=1 , 1, 0)  
      for (i in 2:dim(tempo)[1]){
        tempo$inf[i] <- ifelse(tempo[i, vec_tot_16[Ab]]> 4*tempo[dim(tempo)[1], "medianRAU"] &                     # 4fold increase compared to the median of all samples
                               tempo[i, vec_tot_16[Ab]]> 4*tempo[i-1, vec_tot_16[Ab]] &                            # 4fold increase compared to previous sample
                               tempo[i, "drop_date"] !=1 &                                                         
                               tempo[i-1, "inf"]!=1, 1, 0)                                                         # previous sample is not a blank / failed sample 
        tempo$drop_date[i] <- ifelse(tempo$inf[i-1] == 1, 0, tempo$drop_date[i])                                   # the sample after an infection is not a blank / failed sample
      }
      
      
      
      tempo$delete[1] <- ifelse(tempo[c(1), vec_tot_16[Ab]] < median(tempo[tempo$inf!=1, c(vec_tot_16[Ab])])/4 & 
                                tempo[c(1), vec_tot_16[Ab]] < tempo[c(2), vec_tot_16[Ab]]/4 &
                                tempo$inf[c(2)] != 1, 1,0)
      for (j in 2:dim(tempo)[1]){
        tempo$delete[j] <- ifelse(tempo[j, vec_tot_16[Ab]]< median(tempo[tempo$inf!=1, c(vec_tot_16[Ab])])/4 &      # sample below 4fold decrease of the median of all samples (except samples marked as potential infection)
                                  tempo[j, vec_tot_16[Ab]]< tempo[j-1, vec_tot_16[Ab]]/4 &                          # sample below 4fold decrease of previous sample
                                  tempo$inf[j-1] != 1 ||                                                            # unless previous one is an infection  
                                  tempo$delete[j - 1] ==1 &                                                         # OR previous sample is blank / failed and this one is 4fold decrease compared to median
                                  tempo[j, vec_tot_16[Ab]]< tempo[j-1, "medianRAU"]/4, 1, 0)
      }
  
      
      # Mark infection to each antibody response in the df with all participants
      for (i in 1:dim(tempo)[1]){
        if(tempo$delete[i] == 1){cohort_several_samples_Notempo[cohort_several_samples_Notempo$unique_id_village == vec_patient_tempo_several_samples[patient] & 
                                                                  cohort_several_samples_Notempo$todays_date == tempo$todays_date[i], vec_tot_16[Ab]] <- NA}
        cohort_several_samples_Notempo[cohort_several_samples_Notempo$unique_id_village == vec_patient_tempo_several_samples[patient] &
                                                               cohort_several_samples_Notempo$todays_date == tempo$todays_date[i], c(39+Ab)]<- ifelse(tempo$inf[i] == 1, 1, 0)
      }
 
    }
    
  }
  
}
  

colnames(cohort_several_samples_Notempo)[40:55] <- c("ENT_ENT_erovirus.CoxB3.VP1", "OTH_EBV.gp125", "OTH_Echovirus", "RES_ADE5", "RES_ADE40", "RES_Adenovirus.T3", "RES_Adenovirus.T5", "RES_RSV.gG",
                                                        "ENT_Rotavirus.VP7", "HCoV_229E_S", "HCoV_HKU1_S","HCoV_HKU1_NP", "HCoV_NL63_NP", "HCoV_NL63_S","HCoV_OC43_NP", "HCoV_OC43_S") 
  
  
  
# Create df for plotting
library(data.table)
long1 <- melt(setDT(cohort_several_samples_Notempo[1:39]), id.vars = c(vec_epidemio), variable.name = "Antigen")
long2 <- melt(setDT(cohort_several_samples_Notempo[c(1:23, 40:55)]), id.vars = c(vec_epidemio), variable.name = "Infection")
  
df_long <- cbind(long1, long2[, c(24,25)])
colnames(df_long)[27] <- 'INFECTION'
df_long$INFECTION <- as.factor(df_long$INFECTION)

vec_participants <- vec_patient_tempo_several_samples

  

plot_indiv_mini <- list()
vec_tot_mini <- c("RES_ADE5", "RES_RSV.gG")

# 12 participants selected
vec_participants <- c("L-LIB0081", "L-LIB0170", "U-UMI0696", "L-LIB0294",
                      "L-LIB0334", "L-LIB0103", "U-UMI0674", "U-UMI0588", 
                      "L-LIB0255", "L-LIB0003", "U-UMI0546", "L-LIB0316")


df_long <- df_long[df_long$Antigen %in% vec_tot_mini]
df_long$Antigen <- ifelse(df_long$Antigen== "RES_ADE5","Adenovirus 5","RSV gG")

# Plot ________________________________________________________________________________________________________________________
for (participant in 1: length(vec_participants)){
  tempo2 <- df_long[df_long$unique_id_village == vec_participants[participant],]
  age_tempo <- round(tempo2$Age_num[1], 0)
  
  if(participant == 7 | participant == 11| participant == 12){
    
    plot_indiv_mini[[participant]] <- ggplot(tempo2, aes(x = as.Date(todays_date), y = value, group = Antigen, color = Antigen, shape = INFECTION), na.rm = TRUE)+ 
      geom_point(stroke = 2)+
      geom_line(aes(x = as.Date(todays_date), y = value, group = Antigen, color = Antigen), na.rm = TRUE, size = 0.5)+
      scale_shape_manual("Infection", values = c(16, 4))+
      scale_size_manual("Infection", values = c(3, 5))+
      scale_color_manual(name = "Antigen", values = c("Adenovirus 5" = "#CC6666", "RSV gG" = "#99CCFF")) +
      scale_fill_manual(name = "Antigen", values = c("Adenovirus 5" = "#CC6666", "RSV gG" = "#99CCFF")) +
      scale_y_log10(limits = c(9,100000))+
      theme_classic()+
      ggtitle(paste0("P", participant, ", aged ", age_tempo, "yo"))+ xlab("Date")+ ylab("RAU")+
      theme(plot.title = element_text(size=18, face = "bold", color = "black"),
            axis.text=element_text(size=16),
            axis.title = element_text(size=18, face = "bold"),
            legend.text = element_text(size=20),
            legend.ticks = element_line(linewidth = 22),
            legend.title= element_text(size=20))+
      guides(color = guide_legend(override.aes = list(size = 12)), shape = guide_legend(override.aes = list(size = c(9))))
  }else{
    
      plot_indiv_mini[[participant]] <- ggplot(tempo2, aes(x = as.Date(todays_date), y = value, group = Antigen, color = Antigen, shape = INFECTION), na.rm = TRUE)+ 
    geom_point(stroke = 2)+
    geom_line(aes(x = as.Date(todays_date), y = value, group = Antigen, color = Antigen), na.rm = TRUE, size = 0.5)+
    scale_shape_manual("Infection", values = c(16, 4))+
    scale_size_manual("Infection", values = c(3, 5))+
    scale_color_manual(name = "Antigen", values = c("Adenovirus 5" = "#CC6666", "RSV gG" = "#99CCFF")) +
    scale_fill_manual(name = "Antigen", values = c("Adenovirus 5" = "#CC6666", "RSV gG" = "#99CCFF")) +
    scale_y_log10(limits = c(9,100000))+
    theme_classic()+
    ggtitle(paste0("P", participant, ", aged ", age_tempo, "yo"))+ xlab("Date")+ ylab("RAU")+
    theme(plot.title = element_text(size=18, face = "bold", color = "black"),
          axis.text=element_text(size=16),
          axis.title = element_text(size=18, face = "bold"),
          legend.text = element_text(size=20),
          legend.ticks = element_line(linewidth = 22),
          legend.title= element_text(size=20))+
    guides(color = guide_legend(override.aes = list(size = 12)), shape = guide_legend(override.aes = list(size = c(9,11))))
    
  }

  
}











##############
## PLOT GAM ##
##############
# Download df GAM
df_gam <- read.csv(file = "dfgam_RES_ADE5_.csv")
df_gam <- df_gam[-1, -1]

df_gam$basalAB <- log(df_gam$basalAB)
df_gam <- df_gam[!is.na(df_gam$age),]
colnames(df_gam)[4] <- c('period_time')


# Delete L-LIB 0391: Participant with too much noise
df_gam <- df_gam[df_gam$ID != "L-LIB0391",]


# download results longitudinal
dfanalyses <- read.csv("PostDist_RES_ADE5_3k_.csv")



### PLOT DISTRIBUTION / DENSITY _______________________________________________________________________________________________________________
# 3k estimated serostatus distribution based on baseline ab level 
logRAU <- rep(c( seq(min(df_gam$basalAB, na.rm = TRUE), max(df_gam$basalAB, na.rm = TRUE), 0.1)), 3)
line_estimated <- c(c((1/median(dfanalyses$sigma_N*sqrt(2*pi)))*exp(-(1/2)*((logRAU - median(dfanalyses$mu_N))^2)/(median(dfanalyses$sigma_N))^2)),
                    c((1/median(dfanalyses$sigma_P*sqrt(2*pi)))*exp(-(1/2)*((logRAU - ( median(dfanalyses$mu_N) + median(dfanalyses$mu_A) + median(dfanalyses$mu_P)))^2)/(median(dfanalyses$sigma_P))^2)),
                    c((1/median(dfanalyses$sigma_A*sqrt(2*pi)))*exp(-(1/2)*((logRAU - (median(dfanalyses$mu_N) + median(dfanalyses$mu_A)))^2)/(median(dfanalyses$sigma_A))^2)))
estimated_density <- as.data.frame(cbind(logRAU, line_estimated))
estimated_density$component <- c(rep("Negative", length(logRAU)), 
                                 rep("Positive", length(logRAU)),
                                 rep("Recovered", length(logRAU)))





## GAM __________________________________________________________________________________________________________________________________________
b <- gam(INF~ offset(log(period_time)) + s(basalAB) + s(age), family = binomial(), data=df_gam, method = 'REML')
summary(b)

# Get information from smooth function:
smooth_ADE5 <- plot(b, select = 2, pages = 1, seWithMean = TRUE)
df_smooth_baseline <- data.frame(x = smooth_ADE5[[1]]$x, fit = smooth_ADE5[[1]]$fit, se = smooth_ADE5[[1]]$se)
df_smooth_baseline$ICmin <- df_smooth_baseline$fit - 2*df_smooth_baseline$se
df_smooth_baseline$ICmax <- df_smooth_baseline$fit + 2*df_smooth_baseline$se
df_smooth_age <- data.frame(x = smooth_ADE5[[2]]$x, fit = smooth_ADE5[[2]]$fit, se = smooth_ADE5[[2]]$se)
df_smooth_age$ICmin <- df_smooth_age$fit - 2*df_smooth_age$se
df_smooth_age$ICmax <- df_smooth_age$fit + 2*df_smooth_age$se



# Incidence in person year__________________________________________________________
pred <- predict(b, type = "link", se.fit = TRUE)
df_gam$fit <- exp(pred$fit)
df_gam$lower <- exp(pred$fit - 1.96 * pred$se.fit)
df_gam$upper <- exp(pred$fit + 1.96 * pred$se.fit)


plot(b)
# Obs for plotting
basalAb_cat <- c(min(df_gam$basalAB), 5, 7, 8, 9, max(df_gam$basalAB))
df_gam_BasalAb_obs <- as.data.frame(matrix(ncol = 5, nrow = 5))
colnames(df_gam_BasalAb_obs) <- c("INF", "basalAB", "Ntot_ageclass", "ICmin", "ICmax")

for(BasalAb_tempo in (1:(length(basalAb_cat)-1))){
  df_gam_BasalAb_obs[BasalAb_tempo, ] <- c(sum(df_gam$INF[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/(sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365),
                                           mean(c(basalAb_cat[BasalAb_tempo], basalAb_cat[BasalAb_tempo+1])),
                                           length(df_gam$INF[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]]),
                                           as.numeric(as.vector(poisson.test(sum(df_gam$INF[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]]),
                                                                             (sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365))$conf.int[[1]][1])),
                                           as.numeric(as.vector(poisson.test(sum(df_gam$INF[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]]),
                                                                             (sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365))$conf.int[[2]][1])))
}
df_gam_BasalAb_obs$Data <- "Observed"

# Predictions for plotting
basalAb_cat <- c(min(df_gam$basalAB), 5, 7, 8, 9, max(df_gam$basalAB))
df_gam_BasalAb_pred <- as.data.frame(matrix(ncol = 5, nrow = 5))
colnames(df_gam_BasalAb_pred) <- c("INF", "basalAB", "Ntot_ageclass", "ICmin", "ICmax")

for(BasalAb_tempo in (1:(length(basalAb_cat)-1))){
  df_gam_BasalAb_pred[BasalAb_tempo, ] <- c(sum(df_gam$fit[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/(sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365),
                                            mean(c(basalAb_cat[BasalAb_tempo], basalAb_cat[BasalAb_tempo+1])),
                                            length(df_gam$fit[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]]),
                                            sum(df_gam$lower[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/(sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365),
                                            sum(df_gam$upper[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/(sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365))
}
df_gam_BasalAb_pred$Data <- "Predicted"

df_gam_BasalAb <- rbind(df_gam_BasalAb_obs, df_gam_BasalAb_pred)


# Plot obs
plot_basal_ADE5_tempo <- ggplot(df_gam_BasalAb, aes(x = exp(basalAB), y = INF, shape = Data))+ geom_point(size = 3, alpha = 1, position=position_dodge(0.5), color = "#CC6666")+ 
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), alpha = 0.9, width=0, linewidth = 1,
                position=position_dodge(0.5), color = "#CC6666")+
  scale_x_log10(limits = c(9,100000))+
  ggtitle("C. Adenovirus 5 ")+
  xlab("Baseline Ab level")+ylab("Incidence (person year)")+
    theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=16),
        axis.title = element_text(size=20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x =element_blank())+ ylim(0,6)
plot_basal_ADE5_tempo


plot_basal_ADE5_smooth <- ggplot(data = df_smooth_baseline, mapping = aes(x = exp(x), y = fit))+geom_line()+
  geom_ribbon(mapping =aes(x = exp(x), ymin = ICmin, ymax = ICmax), alpha = 0.1, inherit.aes = FALSE)+
  scale_x_log10(limits = c(9,100000), breaks = c(10, 1000, 100000))+
  xlab("Baseline antibody level")+ylab("smooth component")+
  theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=14),
        axis.title = element_text(size=16, face = "bold"))



library(grid)
plot_basal_ADE5 <- plot_basal_ADE5_tempo + annotation_custom(ggplotGrob(plot_basal_ADE5_smooth), xmin = 3, xmax = 5, 
                                   ymin = 2, ymax = 6)
plot_basal_ADE5

# Legend ________________________________________________________________________________________________________________________________________________________________________
plot_basal_leg <- ggplot(df_gam_BasalAb, aes(x = exp(basalAB), y = INF, shape = Data))+ geom_point(size = 3, alpha = 1, position=position_dodge(0.5), color = "black")+ 
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), alpha = 0.9, width=0, linewidth = 1,
                position=position_dodge(0.5), color = "black")+
  scale_x_log10(limits = c(9,100000))+
  ggtitle("Adenovirus 5 ")+
  xlab("Baseline Ab level")+ylab("Incidence (person year)")+
  theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=16),
        axis.title = element_text(size=20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x =element_blank(),        
        legend.position = "bottom",
        legend.text = element_text(size=20),
        legend.ticks = element_line(linewidth = 16),
        legend.title= element_text(size=20))+
guides(shape = guide_legend(override.aes = list(size = 10)))
plot_basal_leg




plot_component_ADE <- ggplot(estimated_density, aes(y = line_estimated, x = exp(logRAU), colour = (component)), linewidth = 2, alpha = 0.5)+ geom_line()+ 
  geom_area(data = estimated_density, mapping = aes(y = line_estimated, x = exp(logRAU), fill = (component)), alpha = 0.5, position = "identity", stat = "identity")+
  scale_color_manual(name = "Estimated components", values = c( "Smooth component" = "black", "Negative" = "darkgreen", "Positive" = "red", "Recovered" = "darkred"))+
  scale_fill_manual(name = "Estimated components", values = c( "Smooth component" = "black", "Negative" = "darkgreen", "Positive" = "red", "Recovered" = "darkred")) +
  scale_x_log10(limits = c(9,100000))+
  xlab("Baseline Antibody level")+theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=16),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_text(size=18, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size=20),
        legend.ticks = element_line(linewidth = 16),
        legend.title= element_text(size=20))+
  guides(color = guide_legend(override.aes = list(size = 12)))
plot_component_ADE


leg1 <- get_legend(plot_basal_leg)
leg2 <- get_legend(plot_component_ADE)
leg_combined <- ggarrange(leg1, leg2)

plot_component_ADE_tempo <- ggarrange("", plot_component_ADE + theme(legend.position = "none"), "", ncol = 3, nrow = 1, widths = c(0.05, 1, 0.008))

plot_basal_ADE5_tot <- ggarrange(plot_basal_ADE5+ theme(legend.position = "none"), plot_component_ADE_tempo, ncol = 1, nrow = 2, heights = c(1, 0.3))




# Obs for plotting
age_cat <- c(0, 5, 20, 35, 55, 100)
df_gam_ageclass_obs <- as.data.frame(matrix(ncol = 6, nrow = 5))
colnames(df_gam_ageclass_obs) <- c("nb_inf", "age_class", "Ntot_ageclass", "ICmin", "ICmax", "basalAB")

for(age_tempo in (1:(length(age_cat)-1))){
  df_gam_ageclass_obs[age_tempo, ] <- c((sum(df_gam$INF[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))/(sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365),
                                        mean(c(age_cat[age_tempo], age_cat[age_tempo+1])),
                                        length(df_gam$INF[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]),
                                        as.numeric(as.vector(poisson.test((sum(df_gam$INF[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])), 
                                                                          sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365)$conf.int[[1]][1])),
                                        as.numeric(as.vector(poisson.test((sum(df_gam$INF[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])), 
                                                                          sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365)$conf.int[[2]][1])),
                                        mean(df_gam$basalAB[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))
}
df_gam_ageclass_obs$Data <- "Observed"

# predict for plotting
age_cat <- c(0,5, 20, 35, 55, 100)
df_gam_ageclass_predict <- as.data.frame(matrix(ncol = 6, nrow = 5))
colnames(df_gam_ageclass_predict) <- c("nb_inf", "age_class", "Ntot_ageclass", "ICmin", "ICmax", "basalAB")

for(age_tempo in (1:(length(age_cat)-1))){
  df_gam_ageclass_predict[age_tempo, ] <- c((sum(df_gam$fit[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))/(sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365),
                                            mean(c(age_cat[age_tempo], age_cat[age_tempo+1])),
                                            length(df_gam$fit[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]),
                                            (sum(df_gam$lower[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))/(sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365),
                                            (sum(df_gam$upper[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))/(sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365),
                                            mean(df_gam$basalAB[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))
}


df_gam_ageclass_predict$Data <- "Predicted"


df_gam_ageclass <- rbind(df_gam_ageclass_obs, df_gam_ageclass_predict)


# Plot obs
plot_age_ADE5_tempo <- ggplot(df_gam_ageclass, aes(x = age_class, y = nb_inf, shape = Data))+ geom_point(size = 3, alpha = 1, position=position_dodge(3), color = "#CC6666")+ 
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), alpha = 0.9, width=0, linewidth = 1,
                position=position_dodge(3), color = "#CC6666")+
  ggtitle(" ")+
  xlab("Age")+ylab("Incidence (person year)")+theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=16),
        axis.title = element_text(size=18, face = "bold"))+ ylim(0,6)
plot_age_ADE5_tempo


plot_age_ADE5_smooth <- ggplot(data = df_smooth_age, mapping = aes(x = x, y = fit))+geom_line()+
  geom_ribbon(mapping =aes(x = x, ymin = ICmin, ymax = ICmax), alpha = 0.1, inherit.aes = FALSE)+
  xlab("Age")+ylab("smooth component")+
  theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=14),
        axis.title = element_text(size=16, face = "bold"))



library(grid)
plot_age_ADE5 <- plot_age_ADE5_tempo + annotation_custom(ggplotGrob(plot_age_ADE5_smooth), xmin = 35, xmax = 75, 
                                                             ymin = 2, ymax = 6)
plot_age_ADE5




plot_age_ADE5_tot <- ggarrange(plot_age_ADE5+ theme(legend.position = "none"), "", ncol = 1, nrow = 2, heights = c(1, 0.2))



plot_ADE5 <- ggarrange(plot_basal_ADE5_tot, plot_age_ADE5_tot, ncol = 2, nrow = 1, common.legend = TRUE)

















# RSVgG _________________________________________________________________________________________________________________________

df_gam <- read.csv(file = "dfgam_RES_RSV.gG_.csv")
df_gam <- df_gam[-1, -1]

df_gam$basalAB <- log(df_gam$basalAB)
df_gam <- df_gam[!is.na(df_gam$age),]
colnames(df_gam)[4] <- c('period_time')


# Delete L-LIB 0391: Participant with too much noise
df_gam <- df_gam[df_gam$ID != "L-LIB0391",]


# download results longitudinal
dfanalyses <- read.csv("PostDist_RES_RSV.gG_2k.csv")


### PLOT DISTRIBUTION / DENSITY __________________________________________________________________________________________________________
# 3k estimated serostatus distribution based on baseline ab level 
logRAU <- rep(c( seq(min(df_gam$basalAB, na.rm = TRUE), max(df_gam$basalAB, na.rm = TRUE), 0.1)), 3)
line_estimated <- c(c((1/median(dfanalyses$epsilon_N*sqrt(2*pi)))*exp(-(1/2)*((logRAU - median(dfanalyses$mu_N))^2)/(median(dfanalyses$epsilon_N))^2)),
                    c((1/median(dfanalyses$epsilon_P*sqrt(2*pi)))*exp(-(1/2)*((logRAU - (median(dfanalyses$muPos)))^2)/(median(dfanalyses$epsilon_P))^2)))
estimated_density <- as.data.frame(cbind(logRAU, line_estimated))
estimated_density$component <- c(rep("Negative", length(logRAU)), 
                                 rep("Positive", length(logRAU)))


ggplot(df_gam, aes(x = basalAB))+geom_histogram()+scale_x_log10()




## GAM ______________________________________________________________________________________________________________________________
b <- gam(INF~ offset(log(period_time)) + s(basalAB)+ s(age), family = binomial(), data=df_gam, method = 'REML')
summary(b)
plot(b,pages=1,seWithMean=TRUE)


# Get information from smooth function:
smooth_RSV <- plot(b, select = 1, pages = 1, seWithMean = TRUE)
df_smooth_baseline2 <- data.frame(x = smooth_RSV[[1]]$x, fit = smooth_RSV[[1]]$fit, se = smooth_RSV[[1]]$se)
df_smooth_baseline2$ICmin <- (df_smooth_baseline2$fit) - (2*df_smooth_baseline2$se)
df_smooth_baseline2$ICmax <- (df_smooth_baseline2$fit) + (2*df_smooth_baseline2$se)
df_smooth_age2 <- data.frame(x = smooth_RSV[[2]]$x, fit = smooth_RSV[[2]]$fit, se = smooth_RSV[[2]]$se)
df_smooth_age2$ICmin <- df_smooth_age2$fit - 2*df_smooth_age2$se
df_smooth_age2$ICmax <- df_smooth_age2$fit + 2*df_smooth_age2$se


# Incidence in person year__________________________________________________________

pred <- predict(b, type = "link", se.fit = TRUE)
df_gam$fit <- exp(pred$fit)
df_gam$lower <- exp(pred$fit - 1.96 * pred$se.fit)
df_gam$upper <- exp(pred$fit + 1.96 * pred$se.fit)



# Obs for plotting
basalAb_cat <- c(min(df_gam$basalAB), 5, 7, 8, 9, max(df_gam$basalAB))
vec_tot <- c("RES_RSV.gG")
df_gam_BasalAb_obs <- as.data.frame(matrix(ncol = 5, nrow = 5))
colnames(df_gam_BasalAb_obs) <- c("INF", "basalAB", "Ntot_ageclass", "ICmin", "ICmax")

for(BasalAb_tempo in (1:(length(basalAb_cat)-1))){
  df_gam_BasalAb_obs[BasalAb_tempo, ] <- c(sum(df_gam$INF[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/(sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365),
                                           mean(c(basalAb_cat[BasalAb_tempo], basalAb_cat[BasalAb_tempo+1])),
                                           length(df_gam$INF[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]]),
                                           as.numeric(as.vector(poisson.test(sum(df_gam$INF[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]]),
                                                                             (sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365))$conf.int[[1]][1])),
                                           as.numeric(as.vector(poisson.test(sum(df_gam$INF[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]]),
                                                                             (sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365))$conf.int[[2]][1])))
}
df_gam_BasalAb_obs$Data <- "Observed"

# Predictions for plotting
basalAb_cat <- c(min(df_gam$basalAB), 5, 7, 8, 9, max(df_gam$basalAB))
vec_tot <- c("RES_RSV.gG")
df_gam_BasalAb_pred <- as.data.frame(matrix(ncol = 5, nrow = 5))
colnames(df_gam_BasalAb_pred) <- c("INF", "basalAB", "Ntot_ageclass", "ICmin", "ICmax")

for(BasalAb_tempo in (1:(length(basalAb_cat)-1))){
  df_gam_BasalAb_pred[BasalAb_tempo, ] <- c(sum(df_gam$fit[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/(sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365),
                                            mean(c(basalAb_cat[BasalAb_tempo], basalAb_cat[BasalAb_tempo+1])),
                                            length(df_gam$fit[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]]),
                                            sum(df_gam$lower[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/(sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365),
                                            sum(df_gam$upper[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/(sum(df_gam$period_time[df_gam$basalAB > basalAb_cat[BasalAb_tempo] & df_gam$basalAB <= basalAb_cat[BasalAb_tempo +1]])/365))
}
df_gam_BasalAb_pred$Data <- "Predicted"

df_gam_BasalAb <- rbind(df_gam_BasalAb_obs, df_gam_BasalAb_pred)


# Plot obs
plot_basal_RSVgG_tempo <- ggplot(df_gam_BasalAb, aes(x = exp(basalAB), y = INF, shape = Data))+ geom_point(size = 3, alpha = 1, position=position_dodge(0.5), color = "#99CCFF")+ 
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), alpha = 0.9, width=0, linewidth = 1,
                position=position_dodge(0.5), color = "#99CCFF")+
  scale_x_log10(limits = c(9,100000))+
  ggtitle("B. RSV gG")+
  xlab("Baseline Ab level")+ylab("Incidence (person year)")+theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=16),
        axis.title = element_text(size=18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x =element_blank())+ ylim(0,2.1)
plot_basal_RSVgG_tempo




plot_basal_RSVgG_smooth <- ggplot(data = df_smooth_baseline2, mapping = aes(x = exp(x), y = fit))+geom_line()+
  geom_ribbon(mapping =aes(x = exp(x), ymin = ICmin, ymax = ICmax), alpha = 0.1, inherit.aes = FALSE)+
  scale_x_log10(limits = c(9,100000), breaks = c(10, 1000, 100000))+
  xlab("Baseline antibody level")+ylab("smooth component")+
  theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=14),
        axis.title = element_text(size=16, face = "bold"))



plot_basal_RSVgG <- plot_basal_RSVgG_tempo + annotation_custom(ggplotGrob(plot_basal_RSVgG_smooth), xmin = 3, xmax = 5, 
                                                         ymin = 1, ymax = 2)
plot_basal_RSVgG





plot_component_RSV <- ggplot(estimated_density, aes(y = line_estimated, x = exp(logRAU), colour = (component)), linewidth = 2, alpha = 0.5)+ geom_line()+ 
  geom_area(data = estimated_density, mapping = aes(y = line_estimated, x = exp(logRAU), fill = (component)), alpha = 0.5, position = "identity", stat = "identity")+
  scale_color_manual(name = "Components", values = c("darkgreen", "red"))+
  scale_fill_manual(name = "Components",values = c("darkgreen", "red"))+
  scale_x_log10(limits = c(9,100000))+
  xlab("Baseline Antibody level")+theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=16),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_text(size=18, face = "bold"))
plot_component_RSV


plot_component_RSV_tempo <- ggarrange("", plot_component_RSV+ theme(legend.position = "none"), "", ncol = 3, nrow = 1, widths = c(0.05, 1, 0.008))

plot_basal_RSV_tot <- ggarrange(plot_basal_RSVgG+ theme(legend.position = "none"), plot_component_RSV_tempo, ncol = 1, nrow = 2, heights = c(1, 0.3))










# Obs for plotting
age_cat <- c(0, 5, 20, 35, 55, 100)
df_gam_ageclass_obs <- as.data.frame(matrix(ncol = 6, nrow = 5))
colnames(df_gam_ageclass_obs) <- c("nb_inf", "age_class", "Ntot_ageclass", "ICmin", "ICmax", "basalAB")

for(age_tempo in (1:(length(age_cat)-1))){
  df_gam_ageclass_obs[age_tempo, ] <- c((sum(df_gam$INF[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))/(sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365),
                                        mean(c(age_cat[age_tempo], age_cat[age_tempo+1])),
                                        length(df_gam$INF[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]),
                                        as.numeric(as.vector(poisson.test((sum(df_gam$INF[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])), 
                                                                          sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365)$conf.int[[1]][1])),
                                        as.numeric(as.vector(poisson.test((sum(df_gam$INF[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])), 
                                                                          sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365)$conf.int[[2]][1])),
                                        mean(df_gam$basalAB[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))
}
df_gam_ageclass_obs$Data <- "Observed"

# predict for plotting
age_cat <- c(0,5, 20, 35, 55, 100)
df_gam_ageclass_predict <- as.data.frame(matrix(ncol = 6, nrow = 5))
colnames(df_gam_ageclass_predict) <- c("nb_inf", "age_class", "Ntot_ageclass", "ICmin", "ICmax", "basalAB")

for(age_tempo in (1:(length(age_cat)-1))){
  df_gam_ageclass_predict[age_tempo, ] <- c((sum(df_gam$fit[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))/(sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365),
                                            mean(c(age_cat[age_tempo], age_cat[age_tempo+1])),
                                            length(df_gam$fit[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]),
                                            (sum(df_gam$lower[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))/(sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365),
                                            (sum(df_gam$upper[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))/(sum(df_gam$period_time[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]])/365),
                                            mean(df_gam$basalAB[df_gam$age > age_cat[age_tempo] & df_gam$age <= age_cat[age_tempo +1]]))
}


df_gam_ageclass_predict$Data <- "Predicted"


df_gam_ageclass <- rbind(df_gam_ageclass_obs, df_gam_ageclass_predict)


# Plot obs
plot_age_RSVgG_tempo <- ggplot(df_gam_ageclass, aes(x = age_class, y = nb_inf, shape = Data))+ geom_point(size = 3, alpha = 1, position=position_dodge(3), color = "#99CCFF")+ 
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), alpha = 0.9, width=0, linewidth = 1,
                position=position_dodge(3), color = "#99CCFF")+
  ggtitle(" ")+
  xlab("Age")+ylab("Incidence (person year)")+theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=16),
        axis.title = element_text(size=18, face = "bold"))+
  ylim(0,2.1)
plot_age_RSVgG_tempo






plot_age_RSVgG_smooth <- ggplot(data = df_smooth_age2, mapping = aes(x = x, y = fit))+geom_line()+
  geom_ribbon(mapping =aes(x = x, ymin = ICmin, ymax = ICmax), alpha = 0.1, inherit.aes = FALSE)+
  xlab("Age")+ylab("smooth component")+
  theme_classic()+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=14),
        axis.title = element_text(size=16, face = "bold"))



plot_age_RSVgG <- plot_age_RSVgG_tempo + annotation_custom(ggplotGrob(plot_age_RSVgG_smooth), xmin = 35, xmax = 75, 
                                                               ymin = 1, ymax = 2)
plot_basal_RSVgG




plot_age_RSV_tot <- ggarrange(plot_age_RSVgG+ theme(legend.position = "none"), "", ncol = 1, nrow = 2, heights = c(1, 0.2))



plot_RSV <- ggarrange(plot_basal_RSV_tot, plot_age_RSV_tot, ncol = 2, nrow = 1, common.legend = TRUE)







plot_GAM <- ggarrange(plot_RSV, plot_ADE5, leg_combined, ncol = 1, nrow = 3, heights = c(1, 1, 0.05))







# Create a text grob
text_title <- "A."
tgrob <- text_grob(text_title,size = 20, face = "bold")
# Draw the text
plot_text_title <- as_ggplot(tgrob) + theme(plot.margin = margin(0,30,0,0, "cm"))



plot_indiv_mini_all <- ggarrange(plot_indiv_mini[[1]], plot_indiv_mini[[2]], plot_indiv_mini[[3]],
                                 plot_indiv_mini[[4]], plot_indiv_mini[[5]], plot_indiv_mini[[6]],
                                 plot_indiv_mini[[7]], plot_indiv_mini[[8]], plot_indiv_mini[[9]],
                                 plot_indiv_mini[[10]], plot_indiv_mini[[11]], plot_indiv_mini[[12]],
                                 ncol = 3, nrow = 4, common.legend = TRUE, legend = "bottom") 

plot_indiv_mini_all_final <- ggarrange(plot_text_title,plot_indiv_mini_all, ncol = 1,nrow = 2,heights = c(0.05,0.95))



plot_final <- ggarrange(plot_indiv_mini_all_final, "", plot_GAM, ncol= 3, nrow = 1, widths = c(0.9, 0.05, 1.1))
ggsave(filename = "Fig2_5.jpeg", plot_final, width=32, height = 16)



