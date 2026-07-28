rm(list=ls())
library(mgcv)
library(binom)
library(ggplot2)
require(mgcv)
library(ggpubr)
library( broom )
library('tidyr')
library('dplyr')
library(stringr)



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


# Select cohort from 2021
df_Peru <- PeruRAU_wide[PeruRAU_wide$todays_date> "2021-01-01",] # cohort form 2021
df_Peru <- df_Peru[!is.na(df_Peru$Age_num),]
df_Peru <- df_Peru[order(df_Peru$Age_num),]



vec_entero <- c("ENT_Norovirus.GII.4.VLP", "ENT_Norovirus.GII.4.VP1", "ENT_Norovirus.GII.6.VLP", "ENT_ENT_erovirus.CoxB3.VP1", "ENT_Rotavirus.VP7")
vec_resp <- c("RES_Adenovirus.T3", "RES_Adenovirus.T5", "RES_RSV.gG", "RES_ADE5", "RES_ADE40")
vec_HCov <- c("HCoV_HKU1_S","HCoV_HKU1_NP", "HCoV_229E_NP", "HCoV_229E_S", "HCoV_NL63_NP", "HCoV_NL63_S","HCoV_OC43_NP", "HCoV_OC43_S", "HCoV_SARS.CoV.2.Spike.Wuhan", 
              "HCoV_SARS.CoV.2.RBD.Wuhan","HCoV_SARS.CoV.2.NP")
vec_other <- c("OTH_CMV","OTH_EBV.gp125", "OTH_Echovirus", "OTH_HAV")

vec_single <- rev(c("OTH_EBV.gp125", "OTH_CMV", "RES_RSV.gG",  "RES_ADE5", 
                "HCoV_OC43_S", "HCoV_NL63_S","HCoV_229E_S", "HCoV_HKU1_S","HCoV_SARS.CoV.2.NP", 
                 "RES_ADE40", "OTH_HAV", "OTH_Echovirus", "ENT_ENT_erovirus.CoxB3.VP1", "ENT_Rotavirus.VP7", "ENT_Norovirus.GII.4.VLP", "ENT_Norovirus.GII.6.VLP"))
vec_tot <- c(vec_entero, vec_resp, vec_HCov, vec_other)

vec_single_str <- rev(c( "EBV gp125", "CMV", "RSV gG",  "Adenovirus 5",
                        "HCoV OC43 S", "HCoV NL63 S", "HCoV 229E S", "HCoV HKU1 S", "SARS-CoV-2 NP",
                        "Adenovirus 40","HAV","Echovirus", "Enterovirus CoxB3 VP1", "Rotavirus VP7","Norovirus GII 4 VLP", "Norovirus GII 6 VLP"))




################################
### PLOT CUMULATED INCIDENCE ###
################################
# Download predictions from GAM
dfGAM_all = list.files(pattern = paste0("Fig4_GAMpred_"))
list_dfGAM_all <- lapply(dfGAM_all, read.csv)

# Create df with all antibody responses
dfGAM <- as.data.frame(matrix(ncol = 7, nrow = 0))
colnames(dfGAM) <- c("basalAB", "age", "period_time", 'fit_fig4', "lower_fig4", "upper_fig4", "Antigen")
for(Ag in 1: length(vec_single)){
  Ag_recup <- which(grepl(vec_single[Ag], dfGAM_all))
  dfGAM_tempo = as.data.frame(list_dfGAM_all[[Ag_recup]])
  dfGAM_tempo <- dfGAM_tempo[, -1]
  dfGAM_tempo$Antigen <- vec_single_str[[Ag]]

  dfGAM <- rbind(dfGAM, dfGAM_tempo)
}

# Data management
dfGAM$incidence <- dfGAM$fit_fig4/(dfGAM$period_time/365)
dfGAM$incidence_loged <- log(dfGAM$incidence)
dfGAM$Antigen <- as.factor(dfGAM$Antigen)

dfGAM$Antigen <- factor(dfGAM$Antigen, levels=rev(vec_single_str))

df <- dfGAM %>%
  group_by(age) %>%
  mutate(value_cum = cumsum(incidence),
         value_log = log10(value_cum))

# Plot cumulated incidence 
plot_inc_cumul <- df %>%
  ggplot(aes(fill = Antigen, color = Antigen, y = value_cum, x = age)) +
  geom_area(position = "identity", stat ="identity", alpha = 0.8 ) +
  scale_color_manual(values = c("#FF9900", "#FFCC99",
                                "#006633",  "#99FF99",
                                "#000033", "#660099", "#6666FF", "#9999FF","#CC99FF",
                                "#660000", "#990000", "#CC0000", "#FF3333", "#FF6666", "#FF9999", "#FFCCCC"))+
  scale_fill_manual(values = c("#FF9900", "#FFCC99",
                               "#006633",  "#99FF99",
                               "#000033", "#660099", "#6666FF", "#9999FF","#CC99FF",
                               "#660000", "#990000", "#CC0000", "#FF3333", "#FF6666", "#FF9999", "#FFCCCC"))+
  scale_alpha_identity()+
  xlab("Age") + ylab("incidence (per person year)") +
  ggtitle("B. Age stratified estimated incidence")+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.title = element_text(size=20, face="bold"),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20),
        legend.ticks = element_line(linewidth = 16),
        legend.title= element_text(size=20))+
  guides(color = guide_legend(override.aes = list(size = 12)))+
  xlim(0.75, 95.1) + ylim(0,12.5)

plot_inc_cumul



# univariate plot
plot_univ_list <- list()
vec_color_univ <- c("#FFCCCC", "#FF9999", "#FF6666", "#FF3333", "#CC0000", "#990000", "#660000",
                    "#CC99FF", "#9999FF", "#6666FF", "#660099", "#000033",
                    "#99FF99", "#006633",
                    "#FFCC99", "#FF9900")

for ( ag in 1:length(vec_single_str)){
  ag_str <- vec_single_str[ag]
  df_univ <- df[df$Antigen == ag_str,] 

  plot_univ_list[[ag]] <- df_univ %>%
  ggplot(aes(fill = Antigen, color = Antigen, y = incidence, x = age)) +
  geom_area(position = "identity", stat ="identity", alpha = 0.6 ) +
  scale_color_manual(values = c(vec_color_univ[ag]))+
  scale_fill_manual(values = c(vec_color_univ[ag]))+
  scale_alpha_identity()+
  xlab("Age") + ylab("incidence (per person year)") +
  ggtitle(paste0(ag_str))+
  theme_bw()+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.title = element_text(size=20, face="bold"),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20),
        legend.ticks = element_line(linewidth = 16),
        legend.title= element_text(size=20), legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 12)))+
  xlim(0.75, 95.1) + ylim(0,2.26)

  print(plot_univ_list[[ag]])
}

plot_univ <- ggarrange(plot_univ_list[[1]], plot_univ_list[[2]], plot_univ_list[[3]], plot_univ_list[[4]],
                       plot_univ_list[[5]], plot_univ_list[[6]], plot_univ_list[[7]], plot_univ_list[[8]],
                       plot_univ_list[[9]], plot_univ_list[[10]], plot_univ_list[[11]], plot_univ_list[[12]],
                       plot_univ_list[[13]], plot_univ_list[[14]], plot_univ_list[[15]], plot_univ_list[[16]])



ggsave(filename = paste0("Supp_incidence_univ.jpeg"), plot_univ, width=20, height=20)



#__________________________________________________________________________________________________________________________________________________________________












#####################################
### PIE CHART: OBSERVED INCIDENCE ###
#####################################
# Dowload all df
GAM_obs = list.files(pattern = paste0("dfgam_"))
list_dfGAM_obs<- lapply(GAM_obs, read.csv)


# Create one df with all markers
dfGAM_obs <- as.data.frame(matrix(ncol = 6, nrow = 0))
for(Ag in 1: length(GAM_obs)){
  dfGAM_tempo = as.data.frame(list_dfGAM_obs[[Ag]])
  dfGAM_tempo <- dfGAM_tempo[, -1]
  
  dfGAM_obs <- rbind(dfGAM_obs, dfGAM_tempo)
}
colnames(dfGAM_obs) <- c("Infection", "BaselineAB", "age", "period_time", 'ID', "Antigen")
dfGAM_obs <- dfGAM_obs[!is.na(dfGAM_obs$Infection),]
dfGAM_obs <- dfGAM_obs[dfGAM_obs$Antigen %in% vec_single,]


dfGAM_obs$antigen_class <- ifelse(dfGAM_obs$Antigen == "OTH_Echovirus" | dfGAM_obs$Antigen == "OTH_HAV" | dfGAM_obs$Antigen == "RES_ADE40", "ENT", str_sub(dfGAM_obs$Antigen, 1, 3))
dfGAM_obs$AGE_cat <- NA
dfGAM_obs$AGE_cat <- ifelse(dfGAM_obs$age>=0 & dfGAM_obs$age<5, "0-4 yo", dfGAM_obs$AGE_cat)
dfGAM_obs$AGE_cat <- ifelse(dfGAM_obs$age>=5 & dfGAM_obs$age<18, "5-17 yo", dfGAM_obs$AGE_cat)
dfGAM_obs$AGE_cat <- ifelse(dfGAM_obs$age>=18 & dfGAM_obs$age<45, "18-44 yo", dfGAM_obs$AGE_cat)
dfGAM_obs$AGE_cat <- ifelse(dfGAM_obs$age>=45, "45+ yo", dfGAM_obs$AGE_cat)



# Pie 0 to 5 yo _____________________________________________________________________________________________________________________________________________________________________
  df_pie_age1 <- dfGAM_obs[dfGAM_obs$AGE_cat =="0-4 yo",]
  
  # Total Follow up time (to compute incidence in infection / person / year) 
  follow_up_time <- 0
  for (participant in 1:length(unique(df_pie_age1$ID))){
    df_tempo <- df_Peru[df_Peru$unique_id_village == unique(df_pie_age1$ID)[participant], ]
    df_tempo <- df_tempo[order(df_tempo$todays_date),]
    
    # Compute follow up time for each participant
    follow_up_time_tempo <- as.Date(df_tempo$todays_date[length(df_tempo$todays_date)]) - as.Date(df_tempo$todays_date[1])
    
    # Compute total follow up time
    follow_up_time <- follow_up_time + follow_up_time_tempo
  }
  
  
  # Build df with proportion of each antigen among infections
  df_pie1 <- as.data.frame(matrix(NA, nrow = length(vec_single), ncol = 3))
  for(Ag in 1: length(vec_single)){
    df_tempo <- df_pie_age1[df_pie_age1$Antigen == vec_single[Ag],]
    prev_ag <- sum(df_tempo$Infection)/ sum(df_pie_age1$Infection)
    
    df_pie1[Ag,] <- c(prev_ag, vec_single[Ag], ifelse(vec_single[Ag]=="OTH_Echovirus"| vec_single[Ag] == "OTH_HAV" | vec_single[Ag] == "RES_ADE40", 'ENT', str_sub(vec_single[Ag], 1, 3)))
  }
  colnames(df_pie1) <- c("Proportion",  "Antigen", "antigen cat")

  
  # order rows by character variable for plotting
  df_pie1 <- df_pie1[order(df_pie1$Antigen), ]   
  df_pie1$Antigen <- factor(df_pie1$Antigen, levels = rev(vec_single))
  
  df_pie1 <- df_pie1[order(df_pie1$`antigen cat`), ]   
  df_pie1$`antigen cat` <- factor(df_pie1$`antigen cat`, levels = c('RES',  'HCo', 'OTH', 'ENT'))
  
  
  # compute incidence
  NInf_0to5 <- sum(df_pie_age1$Infection)/(as.numeric(follow_up_time)/365)

  # Plot
  pie_0to5 <- ggplot(df_pie1, aes(x = "", y = Proportion, color = `antigen cat`, fill = Antigen, group = Antigen)) +
    geom_col(width = 1, alpha = 0.7, size = 1) +
    scale_color_manual(values = c("#009933", 
                                  "#660099",
                                  "#FF9966",
                                  "#990000"))+
    scale_fill_manual(values = c("#FF9900", "#FFCC99",
                                 "#006633",  "#99FF99",
                                 "#000033", "#660099", "#6666FF", "#9999FF","#CC99FF",
                                 "#660000", "#990000", "#CC0000", "#FF3333", "#FF6666", "#FF9999", "#FFCCCC"))+
    coord_polar("y") +
    theme_void() +
    ggtitle(label = "0 to 4 years old", subtitle = paste0(round(NInf_0to5, 0), " Infections / person / year"))+
    theme(
      plot.title = element_text(size=18, face="bold", hjust = 0.5),
      plot.subtitle = element_text(size=16, hjust = 0.5)
    )
    
   
  
# Pie 5 to 17 yo _____________________________________________________________________________________________________________________________
  df_pie_age2 <- dfGAM_obs[dfGAM_obs$AGE_cat =="5-17 yo",]
  
  # Follow up time 
  follow_up_time2 <- 0
  for (participant in 1:length(unique(df_pie_age2$ID))){
    df_tempo <- df_Peru[df_Peru$unique_id_village == unique(df_pie_age2$ID)[participant], ]
    df_tempo <- df_tempo[order(df_tempo$todays_date),]
    
    follow_up_time_tempo <- as.Date(df_tempo$todays_date[length(df_tempo$todays_date)]) - as.Date(df_tempo$todays_date[1])
    follow_up_time2 <- follow_up_time2 + follow_up_time_tempo
  }
  
  N_5to17 <- length(unique(df_pie_age2$ID))
  NInf_5to17 <- sum(df_pie_age2$Infection)/(as.numeric(follow_up_time2)/365)


  
  df_pie2 <- as.data.frame(matrix(NA, nrow = length(vec_single), ncol = 3))
  for(Ag in 1: length(vec_single)){
    df_tempo <- df_pie_age2[df_pie_age2$Antigen == vec_single[Ag],]
    prev_ag <- sum(df_tempo$Infection)/ sum(df_pie_age2$Infection)
    
    df_pie2[Ag,] <- c(prev_ag, vec_single[Ag], ifelse(vec_single[Ag]=="OTH_Echovirus" | vec_single[Ag] == "OTH_HAV" | vec_single[Ag] == "RES_ADE40", 'ENT', str_sub(vec_single[Ag], 1, 3)))
  }
  colnames(df_pie2) <- c("Proportion",  "Antigen", "antigen cat")
  
  df_pie2 <- df_pie2[order(df_pie2$Antigen), ]   # order rows by character variable
  df_pie2$Antigen <- factor(df_pie2$Antigen, levels = rev(vec_single))
  
  df_pie2 <- df_pie2[order(df_pie2$`antigen cat`), ]   # order rows by character variable
  df_pie2$`antigen cat` <- factor(df_pie2$`antigen cat`, levels = c('RES', 'HCo', 'OTH',  'ENT'))
  
  
  pie_5to17 <- ggplot(df_pie2, aes(x = "", y = Proportion, color = `antigen cat`, fill = Antigen, group = Antigen)) +
    geom_col(width = 1, alpha = 0.7, size = 1) +
    scale_color_manual(values = c("#009933", 
                                  "#660099",
                                  "#FF9966",
                                  "#990000"))+
    scale_fill_manual(values = c("#FF9900", "#FFCC99",
                                 "#006633",  "#99FF99",
                                 "#000033", "#660099", "#6666FF", "#9999FF","#CC99FF",
                                 "#660000", "#990000", "#CC0000", "#FF3333", "#FF6666", "#FF9999", "#FFCCCC"))+
    coord_polar("y") +
    theme_void() +
    ggtitle("5 to 17 years old", subtitle = paste0(round(NInf_5to17, 0), " Infections / person / year"))+
    theme(
      plot.title = element_text(size=18, face="bold", hjust = 0.5),
      plot.subtitle = element_text(size=16, hjust = 0.5)
    )

  
  
  
  
  # Pie 18 to 44 yo ____________________________________________________________________________________________________________________
  df_pie_age3 <- dfGAM_obs[dfGAM_obs$AGE_cat =="18-44 yo",]
  
  
  # Follow up time 
  follow_up_time3 <- 0
  for (participant in 1:length(unique(df_pie_age3$ID))){
    df_tempo <- df_Peru[df_Peru$unique_id_village == unique(df_pie_age3$ID)[participant], ]
    df_tempo <- df_tempo[order(df_tempo$todays_date),]
    
    follow_up_time_tempo <- as.Date(df_tempo$todays_date[length(df_tempo$todays_date)]) - as.Date(df_tempo$todays_date[1])
    follow_up_time3 <- follow_up_time3 + follow_up_time_tempo
  }

  N_18to44 <- length(unique(df_pie_age3$ID))
  NInf_18to44 <- sum(df_pie_age3$Infection)/(as.numeric(follow_up_time3)/365)
  
  
  
  df_pie3 <- as.data.frame(matrix(NA, nrow = length(vec_single), ncol = 3))
  for(Ag in 1: length(vec_single)){
    df_tempo <- df_pie_age3[df_pie_age3$Antigen == vec_single[Ag],]
    prev_ag <- sum(df_tempo$Infection)/ sum(df_pie_age3$Infection)
    
    df_pie3[Ag,] <- c(prev_ag, vec_single[Ag], ifelse(vec_single[Ag]=="OTH_Echovirus" | vec_single[Ag] == "OTH_HAV" | vec_single[Ag] == "RES_ADE40", 'ENT', str_sub(vec_single[Ag], 1, 3)))
  }
  colnames(df_pie3) <- c("Proportion",  "Antigen", "antigen cat")
  
  df_pie3 <- df_pie3[order(df_pie3$Antigen), ]   # order rows by character variable
  df_pie3$Antigen <- factor(df_pie3$Antigen, levels = rev(vec_single))
  
  df_pie3 <- df_pie3[order(df_pie3$`antigen cat`), ]   # order rows by character variable
  df_pie3$`antigen cat` <- factor(df_pie3$`antigen cat`, levels = c('RES', 'HCo', 'OTH', 'ENT'))
  
  
  pie_18to44 <- ggplot(df_pie3, aes(x = "", y = Proportion, color = `antigen cat`, fill = Antigen, group = Antigen)) +
    geom_col(width = 1, alpha = 0.7, size = 1) +
    scale_color_manual(values = c("#009933", 
                                  "#660099",
                                  "#FF9966",
                                  "#990000"))+
    scale_fill_manual(values = c("#FF9900", "#FFCC99",
                                 "#006633",  "#99FF99",
                                 "#000033", "#660099", "#6666FF", "#9999FF","#CC99FF",
                                 "#660000", "#990000", "#CC0000", "#FF3333", "#FF6666", "#FF9999", "#FFCCCC"))+
    coord_polar("y") +
    theme_void() +
    ggtitle("18 to 44 years old", subtitle = paste0(round(NInf_18to44, 0), " Infections / person / year"))+
    theme(
      plot.title = element_text(size=18, face="bold", hjust = 0.5),
      plot.subtitle = element_text(size=16, hjust = 0.5)
    )
  

  
  
  
  
  # Pie 45+ yo _______________________________________________________________________________________________________________________________
  df_pie_age4 <- dfGAM_obs[dfGAM_obs$AGE_cat =="45+ yo",]
  
  
  # Follow up time 
  follow_up_time4 <- 0
  for (participant in 1:length(unique(df_pie_age3$ID))){
    df_tempo <- df_Peru[df_Peru$unique_id_village == unique(df_pie_age3$ID)[participant], ]
    df_tempo <- df_tempo[order(df_tempo$todays_date),]
    
    follow_up_time_tempo <- as.Date(df_tempo$todays_date[length(df_tempo$todays_date)]) - as.Date(df_tempo$todays_date[1])
    follow_up_time4 <- follow_up_time4 + follow_up_time_tempo
  }
  
  N_45 <- length(unique(df_pie_age4$ID))
  NInf_45 <- sum(df_pie_age4$Infection)/(as.numeric(follow_up_time4)/365)
  
  
  df_pie4 <- as.data.frame(matrix(NA, nrow = length(vec_single), ncol = 3))
  for(Ag in 1: length(vec_single)){
    df_tempo <- df_pie_age4[df_pie_age4$Antigen == vec_single[Ag],]
    prev_ag <- sum(df_tempo$Infection)/ sum(df_pie_age4$Infection)
    
    df_pie4[Ag,] <- c(prev_ag, vec_single[Ag], ifelse(vec_single[Ag]=="OTH_Echovirus" | vec_single[Ag] == "OTH_HAV" | vec_single[Ag] == "RES_ADE40", 'ENT', str_sub(vec_single[Ag], 1, 3)))
  }
  colnames(df_pie4) <- c("Proportion",  "Antigen", "antigen cat")
  
  df_pie4 <- df_pie4[order(df_pie4$Antigen), ]   # order rows by character variable
  df_pie4$Antigen <- factor(df_pie4$Antigen, levels = rev(vec_single))
  
  df_pie4 <- df_pie4[order(df_pie4$`antigen cat`), ]   # order rows by character variable
  df_pie4$`antigen cat` <- factor(df_pie4$`antigen cat`, levels = c('RES', 'HCo', 'OTH', 'ENT'))
  
  
  pie_45 <- ggplot(df_pie4, aes(x = "", y = Proportion, color = `antigen cat`, fill = Antigen, group = Antigen)) +
    geom_col(width = 1, alpha = 0.7, size = 1) +
    scale_color_manual(values = c("#009933", 
                                  "#660099",
                                  "#FF9966",
                                  "#990000"))+
    scale_fill_manual(values = c("#FF9900", "#FFCC99",
                                 "#006633",  "#99FF99",
                                 "#000033", "#660099", "#6666FF", "#9999FF","#CC99FF",
                                 "#660000", "#990000", "#CC0000", "#FF3333", "#FF6666", "#FF9999", "#FFCCCC"))+
    coord_polar("y") +
    theme_void() +
    ggtitle("45+ years old", subtitle = paste0(round(NInf_45, 0), " Infections / person / year"))+
    theme(
      plot.title = element_text(size=18, face="bold", hjust = 0.5),
      plot.subtitle = element_text(size=16, hjust = 0.5)
    )
  

  
  
  
  
  


# Common title for pie charts
text <- "A. Age stratified observed proportion of antigens among infections"
tgrob <- text_grob(text,size = 20, face = 'bold', hjust = 1.05)
plot_0 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,0,0,0, "cm"))
  

# plot with pies for each age category. widths is proportional to the incidence in the age category
plot_pies <- ggarrange("", pie_0to5, "", pie_5to17, "", pie_18to44, "", pie_45, "", 
                       ncol = 9, nrow = 1, common.legend = TRUE, legend = 'none', widths = c(0.05, NInf_0to5/100, 0.05, NInf_5to17/100, 0.05, NInf_18to44/100, 0.05, NInf_45/100, 0.05))

plot_pies_title <- ggarrange(plot_0, plot_pies, ncol = 1, nrow = 2,  heights = c(0.1, 0.9))


plot_final <- ggarrange(plot_inc_cumul, plot_pies_title, ncol = 1, nrow = 2, common.legend = TRUE, legend = 'right', heights = c(2/3, 1/3))
plot_final <- ggarrange(plot_pies_title, plot_inc_cumul, ncol = 1, nrow = 2, common.legend = TRUE, legend = 'right', heights = c(1/3, 2/3))


# Save plot  
  ggsave(filename = paste0("Fig3.jpeg"), plot_final, width=24, height=14)

  
  
  
  


  