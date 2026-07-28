rm(list=ls())
library('tidyr')
library('dplyr')
library('rstan')


# Argument given to select antigen
args <- commandArgs(trailingOnly=TRUE)
Ab1 =  as.numeric(as.character(args[1]))

# For parallel computation 
options(mc.cores = parallel::detectCores())

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
df_Peru <- df_Peru[df_Peru$Age_num >0.5,] # Delete children under 6 months old to avoid signal from maternal antibodies
df_Peru <- df_Peru[!is.na(df_Peru$Age_num),]
df_Peru <- df_Peru[order(df_Peru$Age_num),]

df_Peru <- as.data.frame(df_Peru[,c(vec_epidemio, vec_tot)])




# Create list for modelling
df_Peru <- df_Peru[!is.na(df_Peru[,Ab1+23]),]
print("creating data for model")

listmodel <- list(N = length(df_Peru$Age_num),           # Number of samples
                  age = df_Peru$Age_num,                 # Age
                  Y = log(df_Peru[,c(Ab1+23)]))          # Antibody titre







#-------------------#
#-- Model fitting --#
#-------------------#
# Stan program to be used
model <- stan_model(file="1.2.2.GMMxSerocat_2k.stan") 


# Initial parameters to be used
init_fun <-function() {
  list(rho = 0.05, lambda = 0.1, mu_N = 4, mu_P = 2, sigma_N = 0.5,  sigma_P = 0.5)
}



# Run model
print("Starting running model")
fit <- sampling(model,                          # Stan program
                data = listmodel,               # named list of data
                chains = 4,                     # number of Markov chains
                warmup = 3000,                  # number of warmup iterations per chain 
                iter = 10000,                   # total number of iterations per chain 
                refresh = 500,                  # show progress every 'refresh' iterations
                init = init_fun,
                control=list(adapt_delta=0.99,
                             max_treedepth=12),
                verbose = TRUE)


print("Model is done")


# Get the estimated parameters and Log Likelihood for each iteration
results_model <- extract(fit)[c(1:7, 9)] 
tempo <- matrix(unlist(results_model), ncol = 8)
colnames(tempo) <- c("lambda", "rho", "mu_N", "mu_P", "epsilon_N", "epsilon_P", "muPos", "LL")

# Save Results
write.csv(tempo, paste0("PostDist_", vec_tot[Ab1],"_2k.csv"), row.names=FALSE) 





