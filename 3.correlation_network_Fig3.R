rm(list=ls())
library(readr)
library(dplyr)
library(stringr)
library(bootnet)
library(igraph)
library(ggraph)
library(ggrepel)
library(patchwork)
library(colorspace)
library(ggplot2)
library(tidyr)
library(huge)


# This code follows the framework implemented in 
# Yman, V., Rosado, J., Ochida, N. et al. Investigating antibody cross-reactivity and transmission dynamics of alphaviruses and flaviviruses using a multiplex serological assay. 
# Nat Commun 17, 3491 (2026). https://doi.org/10.1038/s41467-026-71451-7




# Download Data and data management
Peru <- read_csv("0.Peru_MFI_RAU_nMFI.csv")
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
vec_resp <- c("RES_Adenovirus.T3", "RES_Adenovirus.T5", "RES_RSV.gG", "RES_ADE5", "RES_ADE40")
vec_HCov <- c("HCoV_HKU1_S","HCoV_HKU1_NP", "HCoV_229E_NP", "HCoV_229E_S", "HCoV_NL63_NP", "HCoV_NL63_S","HCoV_OC43_NP", "HCoV_OC43_S", "HCoV_SARS.CoV.2.Spike.Wuhan", 
              "HCoV_SARS.CoV.2.RBD.Wuhan","HCoV_SARS.CoV.2.NP")
vec_other <- c("OTH_CMV","OTH_EBV.gp125", "OTH_Echovirus", "OTH_HAV")

vec_tot <- c(vec_entero, vec_resp, vec_HCov, vec_other)
vec_epidemio <- c(colnames(PeruRAU_wide)[1:23])

# Select cohort from 2021
df_Peru <- PeruRAU_wide[PeruRAU_wide$todays_date> "2021-01-01",]
df_Peru <- df_Peru[!is.na(df_Peru$Age_num),]
df_Peru <- df_Peru[order(df_Peru$Age_num),]
df_Peru <- as.data.frame(df_Peru[,c(vec_epidemio, vec_tot)])


# Analyse antibody response on the log scale
for (i in 1: length(vec_tot)){
  df_Peru[, 23+i] <- log(df_Peru[, 23+i])
  df_Peru <- df_Peru[!is.na(df_Peru[, 23+i]),]
}

# Echovirus, HAV and adenovirus 40 are enteric viruses
colnames(df_Peru)[47] <- "ENT_Echovirus"
vec_tot[24] <- "ENT_Echovirus"
colnames(df_Peru)[48] <- "ENT_HAV"
vec_tot[25] <- "ENT_HAV"
colnames(df_Peru)[33] <- "ENT_ADE40"
vec_tot[10] <- "ENT_ADE40"






##########################################
### CORRELATION NETWORK (EBIC, gLASSO) ###
##########################################
# define parameters for analysis
gamma_ebic     <- 1
incl_thr       <- 0.5
mag_thr        <- 0.1



set.seed(126)

# define colors and legend for plotting 
virus_palette <- c(
  Respiratory="#00CC66", `Respiratory (Human Corona)`="#6666FF", Enteric="#FF6666", Herpes="#FFCC99"
) |>
  desaturate(amount = 0.18) |>
  lighten(amount = 0.2)

virus_order <- c("Respiratory","Respiratory (Human Corona)","Enteric","Herpes")



# estimate correlation network
set.seed(12)
net_ebic <- estimateNetwork(
  df_Peru[, vec_tot],
  default   = "EBICglasso",
  tuning    = gamma_ebic,
  corMethod = "spearman",
  missing   = "pairwise",
  threshold = TRUE
)

# Extract weight matrix from network
P <- net_ebic$graph
P[!is.finite(P)] <- 0   # any non finite correlation is set to 0
diag(P) <- 0            # correlation to itself = 0
P[P < 0] <- 0           # drop negative edges 

# build full graph
g_all <- graph_from_adjacency_matrix(P, mode = "undirected", weighted = TRUE, diag = FALSE)
E(g_all)$w <- E(g_all)$weight

# For plotting: nodes annotation and legend
vnames <- V(g_all)$name
V(g_all)$virus <- case_when(
  str_detect(vnames, "^RES")  ~ "Respiratory",
  str_detect(vnames, "^OTH")  ~ "Herpes",
  str_detect(vnames, "^HCoV")   ~ "Respiratory (Human Corona)",
  str_detect(vnames, "^ENT")   ~ "Enteric"
) |> factor(levels = virus_order)


# Rescale to improve visualization in the plot
rescale_safe <- function(x, to = c(7, 18)) {    # range 7-18 is default       
  if (!length(x) || all(!is.finite(x))) return(rep(mean(to), length(x)))
  rng <- range(x[is.finite(x)])
  if (!is.finite(rng[1]) || rng[1] == rng[2]) return(rep(mean(to), length(x)))
  scales::rescale(x, to = to)
}

V(g_all)$size <- rescale_safe(degree(g_all), to = c(7, 18))


# Compute layout (Fruchterman-Reingold)
set.seed(12)
w_pos <- if (ecount(g_all)) pmax(1e-6, abs(E(g_all)$w)) else numeric(0)     # positive edge weights used for node placement
xy <- if (vcount(g_all)) layout_with_fr(g_all, weights = w_pos, niter = 4000) else matrix(0, 0, 2)   # Compute layout, xy are the coordinates

V(g_all)$x <- if (length(xy)) xy[, 1] else numeric(0)   # attribute the coordinate to the nodes
V(g_all)$y <- if (length(xy)) xy[, 2] else numeric(0)






#################
### BOOTSTRAP ###
#################
set.seed(126)
# Bootstrap
boot_ebic <- bootnet(
  net_ebic,
  nBoots = 1000,
  type   = "nonparametric"
)
boot_sum <- summary(boot_ebic)

edge_boot <- boot_sum %>%
  filter(type == "edge") %>%    
  mutate(
    inclusion = 1 - prop0,      
    key = ifelse(node1 < node2, 
                 paste(node1, node2, sep = "||"),
                 paste(node2, node1, sep = "||"))
  )


# create a unique edge identifier for each edge in the graph
E(g_all)$key <- apply(ends(g_all, E(g_all)), 1, function(x) {
  if (x[1] < x[2]) paste(x[1], x[2], sep = "||") else paste(x[2], x[1], sep = "||")
})

match_idx <- match(E(g_all)$key, edge_boot$key)

E(g_all)$boot_mean  <- edge_boot$mean[match_idx]
E(g_all)$boot_lower <- edge_boot$CIlower[match_idx]
E(g_all)$boot_upper <- edge_boot$CIupper[match_idx]
E(g_all)$boot_incl  <- edge_boot$inclusion[match_idx]

# Define stable edges
is_stable <- (!is.na(E(g_all)$boot_incl)) &
  (E(g_all)$boot_incl >= incl_thr) &
  !(E(g_all)$boot_lower <= 0 & E(g_all)$boot_upper >= 0) &
  (E(g_all)$boot_mean >= 0) &
  (abs(E(g_all)$boot_mean) >= mag_thr)

E(g_all)$is_stable <- is_stable



# stable graph 
g_stable <- delete_edges(g_all, E(g_all)[!is_stable])

# rescale
V(g_stable)$size <- rescale_safe(degree(g_stable), to = c(7, 22))


# use names define previously
idx <- match(V(g_stable)$name, V(g_all)$name)
V(g_stable)$x     <- V(g_all)$x[idx]
V(g_stable)$y     <- V(g_all)$y[idx]
V(g_stable)$virus <- V(g_all)$virus[idx]

# change nodes' label
V(g_stable)$name <-  c("Norovirus GII4 VLP", "Norovirus GII4 VP1","Norovirus GII6 VLP", "Enterovirus CoxB3 VP1", "Rotavirus VP7", "Adenovirus T3", "Adenovirus T5", 
                    "RSV gG", "Adenovirus 5", "Adenovirus 40", "HCoV HKU1 S", "HCoV HKU1 NP", "HCoV 229E NP", "HCoV 229E S", "HCoV NL63 NP", "HCoV NL63 S", "HCoV OC43 NP",
                    "HCoV OC43 S", "SARS-CoV-2 Spike Wuhan", "SARS-CoV-2 RBD Wuhan", "SARS-CoV-2 NP", "CMV", "EBV gp125", "Echovirus", "HAV")

#plot
p_network_stable <-
  ggraph(g_stable, layout = "manual", x = V(g_stable)$x, y = V(g_stable)$y) +
  geom_edge_link0(aes(edge_width = abs(w)), colour = "black", alpha = 0.42) +
  scale_edge_width(range = c(0.25, 5), limits = c(0, 1),
                   breaks = c(0.1, 0.2, 0.4, 0.6, 0.8, 1),
                   name = "Partial correlation") +
  geom_node_point(aes(colour = virus, size = size), stroke = 0.5,
                  show.legend = c(colour = TRUE, size = FALSE)) +
  scale_colour_manual(values = virus_palette, breaks = virus_order, limits = virus_order,
                      name = "Virus",
                      guide = guide_legend(order = 1, override.aes = list(size = 5, alpha = 1))) +
  scale_size_identity() +
  geom_label_repel(aes(x = x, y = y, label = name),
                   size = 2.5, max.overlaps = 200,
                   point.padding = 0.75, box.padding = 0.5,
                   segment.alpha = 0.25) +
  labs(title = "") +
  theme_graph(base_family = "Helvetica") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.justification = "left",
    plot.title = element_text(size = 20),
    text = element_text(family = "Helvetica")
  )

p_network_stable
