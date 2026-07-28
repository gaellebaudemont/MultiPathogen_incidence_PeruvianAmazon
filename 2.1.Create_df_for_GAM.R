rm(list=ls())
library('tidyr')
library('dplyr')
library('plotrix')


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








###################################
## DATA MANAGEMENT:DELETE BLANKS ##
###################################
# !! Ultimately failed samples were defined case by case. This code was used to highlight and isolate them more efficiently. 

#ID of all participants 
vec_patient_ID_tempo <-  unique(df_Peru$unique_id_village)[!is.na(unique(df_Peru$unique_id_village))]
  
# ID of participants with multiple samples
vec_patient_tempo_several_samples <- vector()
for (participant in 1 : length(vec_patient_ID_tempo)){
  tempo <- df_Peru[df_Peru$unique_id_village == vec_patient_ID_tempo[participant], c("todays_date", "house_id", "Age_num", "sex", "unique_id_village", "Comunidad")]
  vec_patient_tempo_several_samples[participant] <- ifelse(dim(tempo)[1] > 1, vec_patient_ID_tempo[participant], NA)
}
vec_patient_tempo_several_samples <- vec_patient_tempo_several_samples[!is.na(vec_patient_tempo_several_samples)]
cohort_several_samples_Notempo <- df_Peru[df_Peru$unique_id_village %in% vec_patient_tempo_several_samples,]


# Loop for each antigen  
for (Ab in 1:length(vec_tot)){ 
  cohort_several_samples <- cohort_several_samples_Notempo
  cohort_several_samples$Delete <- 0    # To identify failed samples
  cohort_several_samples$INF <- 0       # To identify potential infection (4-fold increase)
    
  # Blanks/ failed samples definition, loop for each participant
  for (participant in 1 : length(vec_patient_tempo_several_samples)){
    tempo <- df_Peru[df_Peru$unique_id_village == vec_patient_tempo_several_samples[participant], c("todays_date", "house_id", "Age_num", "sex", vec_tot[Ab], "unique_id_village", "Comunidad")]
    tempo <- tempo[order(tempo$todays_date),]
    tempo <- tempo[!is.na(tempo[,vec_tot[Ab]]),]
    
    
    if (dim(tempo)[1] >1){
      tempo$delete <- NA         # Samples identified as a failed sample, to be deleted 
      tempo$drop_date <- NA      # Samples showing a 4 fold drop 
      tempo$medianRAU <- NA      # Median of all previous samples
      tempo$inf <- NA            # Samples Identified as an infection (4-fold boost)
      
      if(dim(tempo)[1] >2){
        for (i in 1:(dim(tempo)[1] - 1)){
          tempo$drop_date[i+1] <- ifelse(tempo[i+1, vec_tot[Ab]]< tempo[i, vec_tot[Ab]]/4, 1, 0) # Marked as 1 if 4-fold drop is observed compared to previous sample
        }
      }
      
      
      for (i in 1:dim(tempo)[1]){
        tempo$medianRAU[i] <- median(tempo[c(1:i), vec_tot[Ab]]) # Median of all previous samples
      }
      
      # For the first sample: compared to the 2nd one (not previous one)
      tempo$inf[1] <- ifelse(tempo[c(1), vec_tot[Ab]]> 4*tempo[dim(tempo)[1], "medianRAU"] &                    # 4fold increase compared to the median of all samples
                               tempo[c(1), vec_tot[Ab]]> 4*tempo[c(2), vec_tot[Ab]] &                           # 4fold increase compared to previous sample (only if the previous is not a failed sample)
                               tempo[c(2), "drop_date"] !=1 , 1, 0)  
      # For all samples
      for (i in 2:dim(tempo)[1]){
        tempo$inf[i] <- ifelse(tempo[i, vec_tot[Ab]]> 4*tempo[dim(tempo)[1], "medianRAU"] &                     # 4fold increase compared to the median of all samples
                               tempo[i, vec_tot[Ab]]> 4*tempo[i-1, vec_tot[Ab]] &                               # 4fold increase compared to previous sample
                               tempo[i, "drop_date"] !=1 , 1, 0)                                                # previous sample is not a blank / failed sample 
        tempo$drop_date[i] <- ifelse(tempo$inf[i-1] == 1, 0, tempo$drop_date[i])                                # the sample after an infection is not a blank / failed sample
      }
      
      
      # For the first sample: compared to the 2nd one (not previous one)
      tempo$delete[1] <- ifelse(tempo[c(1), vec_tot[Ab]] < median(tempo[tempo$inf!=1, c(vec_tot[Ab])])/4 &  
                                tempo[c(1), vec_tot[Ab]] < tempo[c(2), vec_tot[Ab]]/4 &
                                tempo$inf[c(2)] != 1, 1,0)
      # For all samples
      for (j in 2:dim(tempo)[1]){
        tempo$delete[j] <- ifelse(tempo[j, vec_tot[Ab]]< median(tempo[tempo$inf!=1, c(vec_tot[Ab])])/4 &         # sample below 4fold decrease of the median of all samples (except samples marked as potential infection)
                                  tempo[j, vec_tot[Ab]]< tempo[j-1, vec_tot[Ab]]/4 &                             # sample below 4fold decrease of previous sample
                                  tempo$inf[j-1] != 1 ||                                                         # unless previous one is an infection  
                                  tempo$delete[j - 1] ==1 &                                                      # OR previous sample is blank / failed and this one is 4fold decrease compared to median
                                  tempo[j, vec_tot[Ab]]< tempo[j-1, "medianRAU"]/4, 1, 0)
      }
  
      
      # Mark samples to be deleted in the df with all participants
      for (i in 1:dim(tempo)[1]){
        if(tempo$delete[i] == 1){cohort_several_samples$Delete[cohort_several_samples$unique_id_village == vec_patient_tempo_several_samples[participant] & 
                                                               cohort_several_samples$todays_date == tempo$todays_date[i]] <- 1}                                # Identify samples to be deleted in global data base
      }
 
    }
    
  }
  
  
  cohort_several_samples <- cohort_several_samples[cohort_several_samples$Delete == 0,] # delete failed samples from the global database

  
  


  ############
  ## DF GAM ##
  ############
  #ID of participants after deleting failed samples
  vec_patient_ID_tempo2 <- unique(cohort_several_samples$unique_id_village)[!is.na(unique(cohort_several_samples$unique_id_village))]
  
  # ID of participants with several samples (after deleting failed samples)
  vec_patient_tempo_several_samples2 <- vector()
  for (participant in 1 : length(vec_patient_ID_tempo2)){
    tempo <- cohort_several_samples[cohort_several_samples$unique_id_village == vec_patient_ID_tempo2[participant], c("todays_date", "house_id", "Age_num", "sex", vec_tot, "unique_id_village", "Comunidad")]
    vec_patient_tempo_several_samples2[participant] <- ifelse(dim(tempo)[1] > 1, vec_patient_ID_tempo2[participant], NA)
  }
  vec_patient_tempo_several_samples2 <- vec_patient_tempo_several_samples2[!is.na(vec_patient_tempo_several_samples2)]
  cohort_several_samples <- cohort_several_samples[cohort_several_samples$unique_id_village %in% vec_patient_tempo_several_samples2,]
  
  
  # Create df
  df_gam <- as.data.frame(matrix(NA, ncol = 5, nrow = 1))
  colnames(df_gam) <- c("INF", "basalAB", "age", "time", "ID")
  
  # Fill df
  for (participant in 1 : length(vec_patient_tempo_several_samples2)){
    tempo <- cohort_several_samples[cohort_several_samples$unique_id_village == vec_patient_tempo_several_samples2[participant], c("todays_date", "house_id", "Age_num", "sex", vec_tot[Ab], "unique_id_village", "Comunidad")]
    tempo <- tempo[order(tempo$todays_date),]
    tempo$inf <- NA   # Identify infections

    # Define infection
    tempo$inf[1] <- ifelse(tempo[c(1), vec_tot[Ab]]> 4*tempo[c(2), vec_tot[Ab]], 1, 0)                        # For the first sample: compared to 2nd                 
    for (i in 2:dim(tempo)[1]){
      tempo$inf[i] <- ifelse(tempo[i, vec_tot[Ab]]> 4*tempo[i-1, vec_tot[Ab]] & tempo$inf[i-1] !=1, 1, 0)     # 4fold increase compared to previous sample
    }
    
    
    # Fill df_gam
    for (i in 2:dim(tempo)[1]){
      df_gam <- rbind(df_gam, c(tempo$inf[i],                                                             # Infection status for each sample
                                tempo[i-1, vec_tot[Ab]],                                                  # antibody titre of the previous sample
                                tempo$Age_num[i],                                                         # age 
                                as.Date(tempo$todays_date[i]) - as.Date(tempo$todays_date[i-1]),          # number of days since the previous sample
                                vec_patient_tempo_several_samples2[participant]))                         # ID of the participant
    }
  }

df_gam <- df_gam[df_gam$time !=0,] # Delete 22 duplicated samples 
df_gam$Ag <- vec_tot[Ab]
  
df_gam2 <- df_gam[order(df_gam$ID),]

# Save df
write.csv(df_gam2, paste0("dfGAM_", vec_tot[Ab], "_.csv")) 
  
  
  }



  

  
  



