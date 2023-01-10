# This code was originally written by Olival et al. (2017)
# and was adapted (lightly) by Liam Shaw 2019 (liam.philip.shaw at gmail dot com)
# for this project.

# See: https://zenodo.org/record/807517 for the original code repository this code was sourced from

# This file is based on: scripts/04-fit-models.R from that repository

# I am grateful to Olival et al. for making their original code available under an MIT License, which also applies here. 
# https://opensource.org/licenses/MIT
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

library(mgcv)
library(dplyr)
library(stringi)
library(parallel)
library(purrr)
library(ggplot2)
library(viridis)
library(knitr)
library(svglite)
library(MuMIn)
library(corrplot)

rm(list = ls())
if(!is.null(dev.list())) dev.off()
base_dir <- "C:/Users/tomos/OneDrive/PhD/ExtinctHosts_PathogenRichness/"
source(paste0(base_dir, "Olival-functions3.R"))
#source(paste0(base_dir, "Shaw_et_al_2020/scripts"))
host.df <- as.data.frame(read.csv(paste0(base_dir, "intermediates/olival_hosts_CLOVER_Farrell_Plourde_nameEdit4.csv")))
host.df2 <- host.df <- host.df[!is.na(host.df$TotVirusPerHost == T),]

# correlation between Olival viruses per host and Clover viruses per host
max_virus_no <- max(host.df2$TotVirusPerHost)
increments <- as.data.frame(c(1:max_virus_no))
colnames(increments) <- "inc"

olival_clover_correlation_plot <- ggplot() + 
  theme_classic() +
  geom_point(data = host.df2, 
             aes(x = host.df2$TotVirusPerHost, 
                 y = host.df2$n_virus_CLOVER),
             color = "coral",
             position=position_jitter(h=0.0,w=0.4)) +
  xlab("Viruses per host - Olival et al (2017)") +
  ylab("Viruses per host - CLOVER") +
  theme(legend.position = "none") + 
  geom_bar(data = increments, aes(x = inc, y = inc), stat = 'identity', alpha = 0, color = "black") +
  geom_vline(xintercept = seq(from=0.5, to=max_virus_no+0.5, by = 1), color = "lightgrey")
#geom_abline(slope = 1, intercept = 0)
plot(olival_clover_correlation_plot)


host.df2 <- host.df[!is.na(host.df$n_virus_Shaw == T),]
# correlation between Shaw viruses per host and Clover viruses per host
max_virus_no <- max(host.df2$n_virus_Shaw)
increments <- as.data.frame(c(1:max_virus_no))
colnames(increments) <- "inc"

shaw_clover_correlation_plot <- ggplot() + 
  theme_classic() +
  geom_point(data = host.df2, 
             aes(x = host.df2$n_virus_Shaw, 
                 y = host.df2$n_virus_CLOVER),
             color = "skyblue",
             position=position_jitter(h=0.0,w=0.4)) +
  xlab("Viruses per host - Shaw et al (2020)") +
  ylab("Viruses per host - CLOVER") +
  theme(legend.position = "none") + 
  geom_bar(data = increments, aes(x = inc, y = inc), stat = 'identity', alpha = 0, color = "black") +
  geom_vline(xintercept = seq(from=0.5, to=max_virus_no+0.5, by = 1), color = "lightgrey")
#geom_abline(slope = 1, intercept = 0)
plot(shaw_clover_correlation_plot)



############## ADD METRICS OF EXTINCT MAMMALS TO THE CURRNET data_set
### ADD METRICS ASSOCIATED WITH PAST EXTINCTIONS to Albery et al (2021) trait data

#host.df <- as.data.frame(readRDS(paste0(base_dir, "albery_hosts_CLOVER_Farrell_Plourde_nameEdit.rds")))
#host.df <- as.data.frame(read.csv(paste0(base_dir, "olival_hosts_CLOVER_Farrell_Plourde.csv")))

phylacine_extant_species_table <- read.csv(paste0(base_dir, "intermediates/Extinct_sympatric_metrics.csv"))
phylacine_extant_species_df <- as.data.frame(phylacine_extant_species_table)

host.df$extinct.host.no <- phylacine_extant_species_df$sympatric_extinct_no[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.total.phylo.dist <- phylacine_extant_species_df$total_phylo_dist[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
#host.df$extinct.log.total.phylo.dist <- log(phylacine_extant_species_df$total_phylo_dist[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)])
host.df$extinct.mean.phylo.dist <- phylacine_extant_species_df$average_phylo_dist[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
#host.df$extinct.log.mean.phylo.dist <- log(phylacine_extant_species_df$average_phylo_dist[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)])
host.df$extinct.min.phylo.dist <- phylacine_extant_species_df$min_phylo_dist[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
#host.df$extinct.log.min.phylo.dist <- log(phylacine_extant_species_df$min_phylo_dist[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)])
host.df$extinct.median.phylo.dist <- phylacine_extant_species_df$median_phylo_dist[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
#host.df$extinct.log.median.phylo.dist <- log(phylacine_extant_species_df$median_phylo_dist[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)])
host.df$extinct.mean.HR <- phylacine_extant_species_df$mean_HR[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.median.HR <- phylacine_extant_species_df$median_HR[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.mean.FD <- phylacine_extant_species_df$mean_FD[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.median.FD <- phylacine_extant_species_df$median_FD[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.mean.DR <- phylacine_extant_species_df$mean_DR[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.median.DR <- phylacine_extant_species_df$median_DR[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.mean.PT <- phylacine_extant_species_df$mean_PT[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.median.PT <- phylacine_extant_species_df$median_PT[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.mean.biomass <- phylacine_extant_species_df$mean_biomass[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.median.biomass <- phylacine_extant_species_df$median_biomass[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.n.order <- phylacine_extant_species_df$no_overlapping_order[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.n.family <- phylacine_extant_species_df$no_overlapping_family[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df$extinct.n.genera <- phylacine_extant_species_df$no_overlapping_genera[match(host.df$hHostNameFinal_Phylacine, phylacine_extant_species_df$Extant_species...1.)]
host.df <- host.df[!is.na(host.df$extinct.total.phylo.dist) == T, ]
host.df$hOrder <- as.factor(host.df$hOrder)

#host.df$extinct.n.order <- host.df$extinct.n.order/host.df$extinct.host.no
#host.df$extinct.n.family <- host.df$extinct.n.family/host.df$extinct.host.no

##host.df$extinct.mean.biomass <- log(host.df$extinct.mean.biomass)
##host.df$extinct.median.biomass <- log(host.df$extinct.median.biomass)
##host.df$extinct.mean.DR <- log(host.df$extinct.mean.DR)
##host.df$extinct.median.DR <- log(host.df$extinct.median.DR)
##host.df$extinct.mean.FD <- log(host.df$extinct.mean.FD)
##host.df$extinct.median.FD <- log(host.df$extinct.median.FD)
##host.df$extinct.mean.HR <- log(host.df$extinct.mean.HR)
##host.df$extinct.median.HR <- log(host.df$extinct.median.HR)

host.df$extinct.median.phylo.dist <- log(host.df$extinct.median.phylo.dist)
host.df$extinct.mean.phylo.dist <- log(host.df$extinct.mean.phylo.dist)
host.df$extinct.min.phylo.dist <- log(host.df$extinct.min.phylo.dist)

host.df$Sympat_n_Phylacine <-  host.df$extant_sym_decrease
host.df$Sympat_n_Phylacine <- log(host.df$Sympat_n_Phylacine) #######


### Only focusing on clover viruses, bacteria, helminths and arthropods
### As these are the groups with hosts that have a parasite count greater than 1
hosts.virus.clover <- host.df[!is.na(host.df$n_virus_CLOVER),]
hosts.bacteria.clover <- host.df[!is.na(host.df$n_bact_CLOVER),]
#hosts.fungi.clover <- host.df[!is.na(host.df$n_fungi_CLOVER),]
#hosts.protozoa.clover <- host.df[!is.na(host.df$n_protozoa_CLOVER),]
hosts.helminth.clover <- host.df[!is.na(host.df$n_helminth_CLOVER),]
hosts.helminth.clover <- hosts.helminth.clover[!is.na(hosts.helminth.clover$S100),]
hosts.helminth.gmpd2 <-host.df[!is.na(host.df$n_helminth_GMPD2),]
hosts.helminth.eid2 <- host.df[!is.na(host.df$n_helminth_EID2),]
hosts.arthropod.farrell <- host.df[!is.na(host.df$n_arthropods_Farrell),]
hosts.virus.shaw <- host.df[!is.na(host.df$n_virus_Shaw),]
hosts.virus.olival <- host.df[!is.na(host.df$TotVirusPerHost),]
hosts.bacteria.shaw <- host.df[!is.na(host.df$n_bact_Shaw),]

hosts.virus.shaw.zoo <- host.df[!is.na(host.df$n_virus_Shaw_zoo),]
hosts.bacteria.shaw.zoo <- host.df[!is.na(host.df$n_bact_Shaw_zoo),]

hosts.virus.shaw.rna <- host.df[!is.na(host.df$n_virus_Shaw_RNA),]
hosts.virus.shaw.dna <- host.df[!is.na(host.df$n_virus_Shaw_DNA),]

#host.df$helminth_sum <- host.df$n_helminth_EID2 + host.df$n_helminth_GMPD2
#host.df$helminth_diff <- host.df$helminth_sum - host.df$n_helminth_CLOVER
#unique(host.df$helminth_diff)


clover_virus <- readRDS(paste0(base_dir, "intermediates/clover_hostvirus"))
clover_bacteria <- readRDS(paste0(base_dir, "intermediates/clover_hostbacteria"))
clover_helminth4 <- readRDS(paste0(base_dir, "intermediates/clover_hosthelminth"))
shaw_hostvirus <- readRDS(paste0(base_dir, "intermediates/shaw_hostvirus"))
shaw_hostbacteria <- readRDS(paste0(base_dir, "intermediates/shaw_hostbacteria"))
gmpd2_helminth <- readRDS(paste0(base_dir, "intermediates/gmpd2_hosthelminth"))  




####---- CLOVER Virus GAM - All Associations ----
data_set <- hosts.virus.clover %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))

data_set <- data_set %>%
  tidyr::drop_na(c(LnAreaHost,
                   hMassGramsPVR,
                   S100, 
                   S80, 
                   S50,
                   S40,
                   S20,
                   S,
                   hDiseaseZACitesLn,
                   extinct.host.no,
                   Plourde_PC1,
                   Plourde_PC2,
                   #extinct.total.phylo.dist,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   #extinct.mean.biomass,
                   #extinct.median.biomass,
                   #extinct.mean.DR,
                   #extinct.median.DR,
                   #extinct.mean.FD,
                   #extinct.median.FD,
                   #extinct.mean.HR,
                   #extinct.median.HR,
                   #CR, CR_prop,
                   #EN, EN_prop,
                   #VU, VU_prop,
                   #total_threat, #total_threat_prop,
                   Sympat_n_Phylacine,
                   #Sympat_min_PhyloDist,
                   #Sympat_mean_PhyloDist,
                   #Sympat_median_PhyloDist,
                   #Phylacine_extinct_over_extantnatural,
                   #extinct.n.order,
                   #extinct.n.genera
                   extinct.n.family
                   
                   ))


clover_virushost_in_olival <- clover_virus[which(gsub(" ","_", clover_virus$Host) %in% tolower(data_set$hHostNameFinal)),]
print(paste0("Number of unique hosts - virus associations, from CLOVER, used in this study: ", nrow(clover_virushost_in_olival)))
print(paste0("Unique hosts: ", nrow(data_set), "; Unique viruses: ", length(unique(clover_virushost_in_olival$Pathogen))))



#  Create dummy variables for orders to use as random effects
# But note that these become non-conformable
# Could just leave them out?
data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
# The problem (from pairs plot) is hOrderSCANDENTIA and hOrderPERAMELEMORPHIA which both have zero
dummys = dummys[, colSums(dummys) >= 5]

data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)


## Create data.frame of all possible models
############# NEEDS ADJUSTING TO INCLUDE EXTINCT MAMMAL METRICS
terms = list(
  s1 = c("s(LnAreaHost, bs='cs', k = 10)"),
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(extinct.host.no, bs='cs', k = 10)"),
  s5 = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
         "s(extinct.median.phylo.dist, bs='cs', k = 10)",
         "s(extinct.min.phylo.dist, bs='cs', k = 10)"##,
         ##         "s(extinct.host.no, bs='cs', k = 10)",
         ##         "s(extinct.n.family, bs='cs', k = 10)",
         ##         "s(extinct.n.order, bs='cs', k = 10)",
         ##         "s(extinct.n.genera, bs='cs', k = 10)"
  ),
  #s6 = c(#"s(extinct.mean.HR, bs='cs', k = 10)",
  #"s(extinct.median.HR, bs='cs', k = 10)",
  #"s(extinct.mean.FD, bs='cs', k = 10)",
  #"s(extinct.median.FD, bs='cs', k = 10)",
  #"s(extinct.mean.DR, bs='cs', k = 10)",
  #"s(extinct.median.DR, bs='cs', k = 10)"),
  s6 = c("s(extinct.host.no, bs='cs', k = 10)",
         "s(extinct.n.family, bs='cs', k = 10)"#, 
         #"s(extinct.n.order, bs='cs', k = 10)",
         #"s(extinct.n.genera, bs='cs', k = 10)"
  ),
  s7 = c("s(Plourde_PC1, bs='cs', k = 10)"), ###
  s8 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  #  s10 = c("s(Sympat_min_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_mean_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_median_PhyloDist, bs='cs', k = 7)" 
  #        ),
  #s8 = c("s(Phylacine_extinct_over_extantnatural, bs='cs', k = 7)"),
  #s9 = c(
  #"s(CR, bs='cs', k = 10)", ###
  #"s(EN, bs='cs', k = 10)", ###
  #"s(VU, bs='cs', k = 10)", ###
  #"s(total_threat, bs='cs', k = 20)"
  #),
  s10 = c("s(Sympat_n_Phylacine, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

saveRDS(data_set, paste0(base_dir, "intermediates/CLOVER_virus_dataset.rds"))

#all_virus <- fit_all_gams_qpoisson(data_set,
#                                outcome_variable = "n_virus_CLOVER",
#                                  terms)
#saveRDS(all_virus, paste0(base_dir, "intermediates/clover_virus_models_qpoisson.rds"))

all_virus <- fit_all_gams_poisson(data_set,
                                   outcome_variable = "n_virus_CLOVER",
                                   terms)
saveRDS(all_virus, paste0(base_dir, "intermediates/clover_virus_models_poisson.rds"))

all_virus <- fit_all_gams_tw(data_set,
                                  outcome_variable = "n_virus_CLOVER",
                                  terms)
saveRDS(all_virus, paste0(base_dir, "intermediates/clover_virus_models_tw.rds"))

#all_virus <- fit_all_gams_nb(data_set,
#                                  outcome_variable = "n_virus_CLOVER",
#                                  terms)
#saveRDS(all_virus, paste0(base_dir, "intermediates/clover_virus_models_nb.rds"))

#all_virus <- fit_all_gams_scat(data_set,
#                             outcome_variable = "n_virus_CLOVER",
#                             terms)
#saveRDS(all_virus, paste0(base_dir, "intermediates/clover_virus_models_scat.rds"))

#all_virus <- fit_all_gams_gamma(data_set,
#                             outcome_variable = "n_virus_CLOVER",
#                             terms)
#saveRDS(all_virus, paste0(base_dir, "intermediates/clover_virus_models_gamma.rds"))

#all_virus <- fit_all_gams_gau(data_set,
#                                outcome_variable = "n_virus_CLOVER",
#                                terms)
#saveRDS(all_virus, paste0(base_dir, "intermediates/clover_virus_models_gau.rds"))


#### Version of model without extinction variables:

terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 10)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(Plourde_PC1, bs='cs', k = 10)"),
  #s5 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

#all_virus_noExtinction <- fit_all_gams_qpoisson(data_set,
#                                               outcome_variable = "n_virus_CLOVER",
#                                               terms)
#saveRDS(all_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_virus_models_qpoisson.rds"))

all_virus_noExtinction <- fit_all_gams_poisson(data_set,
                                                outcome_variable = "n_virus_CLOVER",
                                                terms)
saveRDS(all_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_virus_models_poisson.rds"))

all_virus_noExtinction <- fit_all_gams_tw(data_set,
                                               outcome_variable = "n_virus_CLOVER",
                                               terms)
saveRDS(all_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_virus_models_tw.rds"))

#all_virus_noExtinction <- fit_all_gams_nb(data_set,
#                                               outcome_variable = "n_virus_CLOVER",
#                                               terms)
#saveRDS(all_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_virus_models_nb.rds"))

#all_virus_noExtinction <- fit_all_gams_scat(data_set,
#                                          outcome_variable = "n_virus_CLOVER",
#                                          terms)
#saveRDS(all_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_virus_models_scat.rds"))

#all_virus_noExtinction <- fit_all_gams_gamma(data_set,
#                                          outcome_variable = "n_virus_CLOVER",
#                                          terms)
#saveRDS(all_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_virus_models_gamma.rds"))

#all_virus_noExtinction <- fit_all_gams_gau(data_set,
#                                             outcome_variable = "n_virus_CLOVER",
#                                             terms)
#saveRDS(all_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_virus_models_gau.rds"))




####---- CLOVER Bacteria GAM - All Associations ----
data_set <- hosts.bacteria.clover %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order),
         !is.na(extinct.host.no))

data_set <- data_set %>%
  tidyr::drop_na(c(LnAreaHost,
                   hMassGramsPVR,
                   S100, 
                   S80, 
                   S50,
                   S40,
                   S20,
                   S,
                   hDiseaseZACitesLn,
                   extinct.host.no,
                   Plourde_PC1,
                   Plourde_PC2,
                   #extinct.total.phylo.dist,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   #extinct.mean.biomass,
                   #extinct.median.biomass,
                   #extinct.mean.DR,
                   #extinct.median.DR,
                   #extinct.mean.FD,
                   #extinct.median.FD,
                   #extinct.mean.HR,
                   #extinct.median.HR,
                   #CR, CR_prop,
                   #EN, EN_prop,
                   #VU, VU_prop,
                   #total_threat, #total_threat_prop,
                   Sympat_n_Phylacine,
                   #Sympat_min_PhyloDist,
                   #Sympat_mean_PhyloDist,
                   #Sympat_median_PhyloDist,
                   #Phylacine_extinct_over_extantnatural,
                   #extinct.n.order,
                   #extinct.n.genera
                   extinct.n.family
                   
  ))



clover_bacthost_in_olival <- clover_bacteria[which(gsub(" ","_", clover_bacteria$Host) %in% tolower(data_set$hHostNameFinal)),]
print(paste0("Number of unique hosts - bacteria associations, from CLOVER, used in this study: ", nrow(clover_bacthost_in_olival)))
print(paste0("Unique hosts: ", nrow(data_set), "; Unique viruses: ", length(unique(clover_bacthost_in_olival$Pathogen))))




#  Create dummy variables for orders to use as random effects
# But note that these become non-conformable
# Could just leave them out?
data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
# The problem (from pairs plot) is hOrderSCANDENTIA and hOrderPERAMELEMORPHIA which both have zero
dummys = dummys[, colSums(dummys) >= 5]

data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

data_set$hOrder <- as.factor(data_set$hOrder) ###

## Create data.frame of all possible models
############# NEEDS ADJUSTING TO INCLUDE EXTINCT MAMMAL METRICS
terms = list(
  s1 = c("s(LnAreaHost, bs='cs', k = 10)"),
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"#,
    #"s(Sympat_n_Phylacine, bs='cs', k = 10)"
  ),
  #s4 = c("s(extinct.host.no, bs='cs', k = 10)"),
  s5 = c(#"s(extinct.mean.phylo.dist, bs='cs', k = 10)",
         #"s(extinct.median.phylo.dist, bs='cs', k = 10)",
         "s(extinct.min.phylo.dist, bs='cs', k = 10)",
         "s(extinct.host.no, bs='cs', k = 10)",
         "s(extinct.n.family, bs='cs', k = 10)"#,
         #"s(extinct.median.HR, bs='cs', k = 10)",
         #"s(extinct.median.FD, bs='cs', k = 10)",
         #"s(extinct.median.DR, bs='cs', k = 10)"
  ),
  #s6 = c(#"s(extinct.mean.HR, bs='cs', k = 10)",
  #"s(extinct.median.HR, bs='cs', k = 10)",
  #"s(extinct.mean.FD, bs='cs', k = 10)",
  #"s(extinct.median.FD, bs='cs', k = 10)",
  #"s(extinct.mean.DR, bs='cs', k = 10)",
  #"s(extinct.median.DR, bs='cs', k = 10)"),
  #s6 = c(#"s(extinct.host.no, bs='cs', k = 10)",
  #"s(extinct.n.family, bs='cs', k = 10)"#, 
  #"s(extinct.n.order, bs='cs', k = 10)",
  #"s(extinct.n.genera, bs='cs', k = 10)"
  #),
  s7 = c("s(Plourde_PC1, bs='cs', k = 10)"), ###
  s8 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  #  s10 = c("s(Sympat_min_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_mean_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_median_PhyloDist, bs='cs', k = 7)" 
  #        ),
  #s8 = c("s(Phylacine_extinct_over_extantnatural, bs='cs', k = 7)"),
##  s9 = c(  ########################################
  #"s(CR, bs='cs', k = 10)", ###
  #"s(EN, bs='cs', k = 10)", ###
  #"s(VU, bs='cs', k = 10)", ###
##      "s(total_threat, bs='cs', k = 10)" 
##    ), 
  s10 = c("s(Sympat_n_Phylacine, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)



# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

saveRDS(data_set, paste0(base_dir, "intermediates/CLOVER_bacteria_dataset.rds"))

#all_bacteria <- fit_all_gams_qpoisson(data_set,
#                                     outcome_variable = "n_bact_CLOVER",
#                                     terms)
#saveRDS(all_bacteria, paste0(base_dir, "intermediates/clover_bacteria_models_qpoisson.rds"))

all_bacteria <- fit_all_gams_poisson(data_set,
                                      outcome_variable = "n_bact_CLOVER",
                                      terms)
saveRDS(all_bacteria, paste0(base_dir, "intermediates/clover_bacteria_models_poisson.rds"))

all_bacteria <- fit_all_gams_tw(data_set,
                                     outcome_variable = "n_bact_CLOVER",
                                     terms)
saveRDS(all_bacteria, paste0(base_dir, "intermediates/clover_bacteria_models_tw.rds"))

#all_bacteria <- fit_all_gams_nb(data_set,
#                                     outcome_variable = "n_bact_CLOVER",
#                                     terms)
#saveRDS(all_bacteria, paste0(base_dir, "intermediates/clover_bacteria_models_nb.rds"))


#all_bacteria <- fit_all_gams_scat(data_set,
#                               outcome_variable = "n_bact_CLOVER",
#                               terms)
#saveRDS(all_bacteria, paste0(base_dir, "intermediates/clover_bacteria_models_scat.rds"))

#all_bacteria <- fit_all_gams_gamma(data_set,
#                               outcome_variable = "n_bact_CLOVER",
#                               terms)
#saveRDS(all_bacteria, paste0(base_dir, "intermediates/clover_bacteria_models_gamma.rds"))

#all_bacteria <- fit_all_gams_gau(data_set,
#                                   outcome_variable = "n_bact_CLOVER",
#                                   terms)
#saveRDS(all_bacteria, paste0(base_dir, "intermediates/clover_bacteria_models_gau.rds"))

#### Version of model without extinction variables:

terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 10)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(Plourde_PC1, bs='cs', k = 10)"),
  #s5 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

#all_bacteria_noExtinction <- fit_all_gams_qpoisson(data_set,
#                                                  outcome_variable = "n_bact_CLOVER",
#                                                  terms)
#saveRDS(all_bacteria_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_bacteria_models_qpoisson.rds"))

all_bacteria_noExtinction <- fit_all_gams_poisson(data_set,
                                                   outcome_variable = "n_bact_CLOVER",
                                                   terms)
saveRDS(all_bacteria_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_bacteria_models_poisson.rds"))

all_bacteria_noExtinction <- fit_all_gams_tw(data_set,
                                                  outcome_variable = "n_bact_CLOVER",
                                                  terms)
saveRDS(all_bacteria_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_bacteria_models_tw.rds"))

#all_bacteria_noExtinction <- fit_all_gams_nb(data_set,
#                                                  outcome_variable = "n_bact_CLOVER",
#                                                  terms)
#saveRDS(all_bacteria_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_bacteria_models_nb.rds"))

#all_bacteria_noExtinction <- fit_all_gams_scat(data_set,
#                                             outcome_variable = "n_bact_CLOVER",
#                                             terms)
#saveRDS(all_bacteria_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_bacteria_models_scat.rds"))

#all_bacteria_noExtinction <- fit_all_gams_gamma(data_set,
#                                             outcome_variable = "n_bact_CLOVER",
#                                             terms)
#saveRDS(all_bacteria_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_bacteria_models_gamma.rds"))

#all_bacteria_noExtinction <- fit_all_gams_gau(data_set,
#                                                outcome_variable = "n_bact_CLOVER",
#                                                terms)
#saveRDS(all_bacteria_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_bacteria_models_gau.rds"))



####---- CLOVER Helminth GAM - All Associations ----
data_set <- hosts.helminth.clover %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order),
         !is.na(extinct.host.no))

data_set <- data_set %>%
  tidyr::drop_na(c(LnAreaHost,
                   hMassGramsPVR,
                   S100, 
                   S80, 
                   S50,
                   S40,
                   S20,
                   S,
                   hDiseaseZACitesLn,
                   extinct.host.no,
                   Plourde_PC1,
                   Plourde_PC2,
                   #extinct.total.phylo.dist,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   #extinct.mean.biomass,
                   #extinct.median.biomass,
                   #extinct.mean.DR,
                   #extinct.median.DR,
                   #extinct.mean.FD,
                   #extinct.median.FD,
                   #extinct.mean.HR,
                   #extinct.median.HR,
                   #CR, CR_prop,
                   #EN, EN_prop,
                   #VU, VU_prop,
                   #total_threat, #total_threat_prop,
                   Sympat_n_Phylacine,
                   #Sympat_min_PhyloDist,
                   #Sympat_mean_PhyloDist,
                   #Sympat_median_PhyloDist,
                   #Phylacine_extinct_over_extantnatural,
                   #extinct.n.order,
                   #extinct.n.genera
                   extinct.n.family,
                   prop_veg,
                   prop_invert,
                   prop_vert
                   
                   
  ))


clover_helhost_in_olival <- clover_helminth4[which(gsub(" ","_", clover_helminth4$Host) %in% tolower(data_set$hHostNameFinal)),]
print(paste0("Number of unique hosts - helminth associations, from CLOVER, used in this study: ", nrow(clover_helhost_in_olival)))
print(paste0("Unique hosts: ", nrow(data_set), "; Unique viruses: ", length(unique(clover_helhost_in_olival$Pathogen))))


#  Create dummy variables for orders to use as random effects
# But note that these become non-conformable
# Could just leave them out?
data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
# The problem (from pairs plot) is hOrderSCANDENTIA and hOrderPERAMELEMORPHIA which both have zero
dummys = dummys[, colSums(dummys) >= 5]

data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
############# NEEDS ADJUSTING TO INCLUDE EXTINCT MAMMAL METRICS
terms = list(
  s1 = c("s(LnAreaHost, bs='cs', k = 10)"),
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"#,
    #"s(Sympat_n_Phylacine, bs='cs', k = 10)"
  ),
  #s4 = c("s(extinct.host.no, bs='cs', k = 10)"),
  s5 = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
         "s(extinct.median.phylo.dist, bs='cs', k = 10)",
         "s(extinct.min.phylo.dist, bs='cs', k = 10)",
         "s(extinct.host.no, bs='cs', k = 10)",
         "s(extinct.n.family, bs='cs', k = 10)"#,
         ##         "s(extinct.n.order, bs='cs', k = 10)",
         ##         "s(extinct.n.genera, bs='cs', k = 10)"
  ),
  #s6 = c(#"s(extinct.mean.HR, bs='cs', k = 10)",
  #"s(extinct.median.HR, bs='cs', k = 10)",
  #"s(extinct.mean.FD, bs='cs', k = 10)",
  #"s(extinct.median.FD, bs='cs', k = 10)",
  #"s(extinct.mean.DR, bs='cs', k = 10)",
  #"s(extinct.median.DR, bs='cs', k = 10)"),
#  s6 = c("s(extinct.host.no, bs='cs', k = 10)",
##         "s(extinct.n.family, bs='cs', k = 10)"#, 
         #"s(extinct.n.order, bs='cs', k = 10)",
         #"s(extinct.n.genera, bs='cs', k = 10)"
##  ),
  s7 = c("s(Plourde_PC1, bs='cs', k = 10)"),
  s8 = c("s(Plourde_PC2, bs='cs', k = 10)"), 
  #  s10 = c("s(Sympat_min_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_mean_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_median_PhyloDist, bs='cs', k = 7)" 
  #        ),
  #s8 = c("s(Phylacine_extinct_over_extantnatural, bs='cs', k = 7)"),
  s9 = c("s(prop_veg, bs='cs', k = 10)",
       "s(prop_vert, bs='cs', k = 10"),
  #s9 = c(
    #"s(CR, bs='cs', k = 10)", 
    #"s(EN, bs='cs', k = 10)", 
    #"s(VU, bs='cs', k = 10)", 
  #  "s(total_threat, bs='cs', k = 10)"
  #),
  s10 = c("s(Sympat_n_Phylacine, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

saveRDS(data_set, paste0(base_dir, "intermediates/CLOVER_helminth_dataset.rds"))

all_helminth <- fit_all_gams_poisson(data_set,
                                     outcome_variable = "n_helminth_CLOVER",
                                     terms)
saveRDS(all_helminth, paste0(base_dir, "intermediates/clover_helminth_models_poisson.rds"))

#all_helminth <- fit_all_gams_qpoisson(data_set,
#                                     outcome_variable = "n_helminth_CLOVER",
#                                     terms)
#saveRDS(all_helminth, paste0(base_dir, "intermediates/clover_helminth_models_qpoisson.rds"))

all_helminth <- fit_all_gams_tw(data_set,
                                     outcome_variable = "n_helminth_CLOVER",
                                     terms)
saveRDS(all_helminth, paste0(base_dir, "intermediates/clover_helminth_models_tw.rds"))

all_helminth <- fit_all_gams_nb(data_set,
                                     outcome_variable = "n_helminth_CLOVER",
                                     terms)
saveRDS(all_helminth, paste0(base_dir, "intermediates/clover_helminth_models_nb.rds"))

#all_helminth <- fit_all_gams_gau(data_set,
#                                outcome_variable = "n_helminth_CLOVER",
#                                terms)
#saveRDS(all_helminth, paste0(base_dir, "intermediates/clover_helminth_models_gau.rds"))
 
#all_helminth <- fit_all_gams_gamma(data_set,
#                                 outcome_variable = "n_helminth_CLOVER",
#                                 terms)
#saveRDS(all_helminth, paste0(base_dir, "intermediates/clover_helminth_models_gamma.rds"))

#### Version of model without extinction variables:

terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 10)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(Plourde_PC1, bs='cs', k = 10)"),
  #s5 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  #s6 = c("s(prop_veg, bs='cs', k = 10)",
  #       "s(prop_vert, bs='cs', k = 10"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)


# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

all_helminth_noExtinction <- fit_all_gams_poisson(data_set,
                                                  outcome_variable = "n_helminth_CLOVER",
                                                  terms)
saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_helminth_models_poisson.rds"))


#all_helminth_noExtinction <- fit_all_gams_qpoisson(data_set,
#                                                  outcome_variable = "n_helminth_CLOVER",
#                                                  terms)
#saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_helminth_models_qpoisson.rds"))


all_helminth_noExtinction <- fit_all_gams_tw(data_set,
                                                  outcome_variable = "n_helminth_CLOVER",
                                                  terms)
saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_helminth_models_tw.rds"))


all_helminth_noExtinction <- fit_all_gams_nb(data_set,
                                                  outcome_variable = "n_helminth_CLOVER",
                                                  terms)
saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_helminth_models_nb.rds"))

#all_helminth_noExtinction <- fit_all_gams_gau(data_set,
#                                               outcome_variable = "n_helminth_CLOVER",
#                                               terms)
#saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_helminth_models_gau.rds"))

#all_helminth_noExtinction <- fit_all_gams_gamma(data_set,
#                                              outcome_variable = "n_helminth_CLOVER",
#                                              terms)
#saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_clover_helminth_models_gamma.rds"))


####---- CLOVER Helminth GAM - Just GMPD2 ----
data_set <- hosts.helminth.gmpd2 %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order),
         !is.na(extinct.host.no))

data_set <- data_set %>%
  tidyr::drop_na(c(LnAreaHost,
                   hMassGramsPVR,
                   S100, 
                   S80, 
                   S50,
                   S40,
                   S20,
                   S,
                   hDiseaseZACitesLn,
                   extinct.host.no,
                   Plourde_PC1,
                   Plourde_PC2,
                   #extinct.total.phylo.dist,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   #extinct.mean.biomass,
                   #extinct.median.biomass,
                   #extinct.mean.DR,
                   #extinct.median.DR,
                   #extinct.mean.FD,
                   #extinct.median.FD,
                   #extinct.mean.HR,
                   #extinct.median.HR,
                   #CR, CR_prop,
                   #EN, EN_prop,
                   #VU, VU_prop,
                   #total_threat, #total_threat_prop,
                   Sympat_n_Phylacine,
                   #Sympat_min_PhyloDist,
                   #Sympat_mean_PhyloDist,
                   #Sympat_median_PhyloDist,
                   #Phylacine_extinct_over_extantnatural,
                   #extinct.n.order,
                   #extinct.n.genera
                   extinct.n.family,
                   prop_veg,
                   prop_invert,
                   prop_vert
                   
                   
  ))


gmpd2_helminth <- gmpd2_helminth[which(gsub(" ","_", gmpd2_helminth$Host) %in% tolower(data_set$hHostNameFinal)),]
print(paste0("Number of unique hosts - helminth associations, from GMPD2, used in this study: ", nrow(gmpd2_helminth)))
print(paste0("Unique hosts: ", nrow(data_set), "; Unique helminths: ", length(unique(gmpd2_helminth$Pathogen))))




#  Create dummy variables for orders to use as random effects
# But note that these become non-conformable
# Could just leave them out?
data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
# The problem (from pairs plot) is hOrderSCANDENTIA and hOrderPERAMELEMORPHIA which both have zero
dummys = dummys[, colSums(dummys) >= 5]

data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
############# NEEDS ADJUSTING TO INCLUDE EXTINCT MAMMAL METRICS
terms = list(
  s1 = c("s(LnAreaHost, bs='cs', k = 10)"),
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"#,
    #"s(Sympat_n_Phylacine, bs='cs', k = 10)"
  ),
  #s4 = c("s(extinct.host.no, bs='cs', k = 10)"),
  s5 = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
         "s(extinct.median.phylo.dist, bs='cs', k = 10)",
         "s(extinct.min.phylo.dist, bs='cs', k = 10)"#,
         #"s(extinct.host.no, bs='cs', k = 10)",
         #"s(extinct.n.family, bs='cs', k = 10)"#,
         ##         "s(extinct.n.order, bs='cs', k = 10)",
         ##         "s(extinct.n.genera, bs='cs', k = 10)"
  ),
  #s6 = c(#"s(extinct.mean.HR, bs='cs', k = 10)",
  #"s(extinct.median.HR, bs='cs', k = 10)",
  #"s(extinct.mean.FD, bs='cs', k = 10)",
  #"s(extinct.median.FD, bs='cs', k = 10)",
  #"s(extinct.mean.DR, bs='cs', k = 10)",
  #"s(extinct.median.DR, bs='cs', k = 10)"),
 s6 = c("s(extinct.host.no, bs='cs', k = 10)",
        "s(extinct.n.family, bs='cs', k = 10)"#, 
  #"s(extinct.n.order, bs='cs', k = 10)",
  #"s(extinct.n.genera, bs='cs', k = 10)"
    ),
  s7 = c("s(Plourde_PC1, bs='cs', k = 10)"),
  s8 = c("s(Plourde_PC2, bs='cs', k = 10)"), 
  #  s10 = c("s(Sympat_min_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_mean_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_median_PhyloDist, bs='cs', k = 7)" 
  #        ),
  #s8 = c("s(Phylacine_extinct_over_extantnatural, bs='cs', k = 7)"),
  s9 = c("s(prop_veg, bs='cs', k = 10)",
         "s(prop_vert, bs='cs', k = 10"),
  #s9 = c(
  #"s(CR, bs='cs', k = 10)", 
  #"s(EN, bs='cs', k = 10)", 
  #"s(VU, bs='cs', k = 10)", 
  #  "s(total_threat, bs='cs', k = 10)"
  #),
  s10 = c("s(Sympat_n_Phylacine, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

saveRDS(data_set, paste0(base_dir, "intermediates/GMPD2_helminth_dataset.rds"))

all_helminth <- fit_all_gams_poisson(data_set,
                                     outcome_variable = "n_helminth_GMPD2",
                                     terms)
saveRDS(all_helminth, paste0(base_dir, "intermediates/gmpd2_helminth_models_poisson.rds"))

#all_helminth <- fit_all_gams_qpoisson(data_set,
#                                      outcome_variable = "n_helminth_GMPD2",
#                                      terms)
#saveRDS(all_helminth, paste0(base_dir, "intermediates/gmpd2_helminth_models_qpoisson.rds"))

all_helminth <- fit_all_gams_tw(data_set,
                                outcome_variable = "n_helminth_GMPD2",
                                terms)
saveRDS(all_helminth, paste0(base_dir, "intermediates/gmpd2_helminth_models_tw.rds"))

all_helminth <- fit_all_gams_nb(data_set,
                                outcome_variable = "n_helminth_GMPD2",
                                terms)
saveRDS(all_helminth, paste0(base_dir, "intermediates/gmpd2_helminth_models_nb.rds"))

#all_helminth <- fit_all_gams_gau(data_set,
#                                 outcome_variable = "n_helminth_GMPD2",
#                                 terms)
#saveRDS(all_helminth, paste0(base_dir, "intermediates/gmpd2_helminth_models_gau.rds"))

#all_helminth <- fit_all_gams_scat(data_set,
#                                  outcome_variable = "n_helminth_GMPD2",
#                                  terms)
#saveRDS(all_helminth, paste0(base_dir, "intermediates/gmpd2_helminth_models_scat.rds"))

#all_helminth <- fit_all_gams_gamma(data_set,
#                                  outcome_variable = "n_helminth_GMPD2",
#                                  terms)
#saveRDS(all_helminth, paste0(base_dir, "intermediates/gmpd2_helminth_models_gamma.rds"))


#### Version of model without extinction variables:

terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 10)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(Plourde_PC1, bs='cs', k = 10)"), ###
  #s5 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  #s6 = c("s(prop_veg, bs='cs', k = 10)",
  #       "s(prop_vert, bs='cs', k = 10"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)


# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

all_helminth_noExtinction <- fit_all_gams_poisson(data_set,
                                                  outcome_variable = "n_helminth_GMPD2",
                                                  terms)
saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_gmpd2_helminth_models_poisson.rds"))

#all_helminth_noExtinction <- fit_all_gams_qpoisson(data_set,
#                                                   outcome_variable = "n_helminth_GMPD2",
#                                                   terms)
#saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_gmpd2_helminth_models_qpoisson.rds"))


all_helminth_noExtinction <- fit_all_gams_tw(data_set,
                                             outcome_variable = "n_helminth_GMPD2",
                                             terms)
saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_gmpd2_helminth_models_tw.rds"))


all_helminth_noExtinction <- fit_all_gams_nb(data_set,
                                             outcome_variable = "n_helminth_GMPD2",
                                             terms)
saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_gmpd2_helminth_models_nb.rds"))

#all_helminth_noExtinction <- fit_all_gams_scat(data_set,
#                                               outcome_variable = "n_helminth_GMPD2",
#                                               terms)
#saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_gmpd2_helminth_models_scat.rds"))

#all_helminth_noExtinction <- fit_all_gams_gau(data_set,
#                                               outcome_variable = "n_helminth_GMPD2",
#                                               terms)
#saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_gmpd2_helminth_models_gau.rds"))

#all_helminth_noExtinction <- fit_all_gams_gamma(data_set,
#                                              outcome_variable = "n_helminth_GMPD2",
#                                              terms)
#saveRDS(all_helminth_noExtinction, paste0(base_dir, "intermediates/noExtinction_gmpd2_helminth_models_gamma.rds"))




####---- Shaw Virus GAM - All Associations ----
data_set <- hosts.virus.shaw %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))

data_set <- data_set %>%
  tidyr::drop_na(c(LnAreaHost,
                   hMassGramsPVR,
                   S100, 
                   S80, 
                   S50,
                   S40,
                   S20,
                   S,
                   hDiseaseZACitesLn,
                   extinct.host.no,
                   Plourde_PC1,
                   Plourde_PC2,
                   #extinct.total.phylo.dist,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   #extinct.mean.biomass,
                   #extinct.median.biomass,
                   #extinct.mean.DR,
                   #extinct.median.DR,
                   #extinct.mean.FD,
                   #extinct.median.FD,
                   #extinct.mean.HR,
                   #extinct.median.HR,
                   #CR, CR_prop,
                   #EN, EN_prop,
                   #VU, VU_prop,
                   #total_threat, #total_threat_prop,
                   Sympat_n_Phylacine,
                   #Sympat_min_PhyloDist,
                   #Sympat_mean_PhyloDist,
                   #Sympat_median_PhyloDist,
                   #Phylacine_extinct_over_extantnatural,
                   #extinct.n.order,
                   #extinct.n.genera,
                   extinct.n.family
                   
  ))

shaw_hostvirus <- shaw_hostvirus[which(shaw_hostvirus$HostSpecies %in% data_set$hHostNameFinal),]
print(paste0("Number of unique hosts - virus associations, from Shaw et al. (2020), used in this study: ", nrow(shaw_hostvirus)))
print(paste0("Unique hosts: ", nrow(data_set), "; Unique viruses: ", length(unique(shaw_hostvirus$Species))))


data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
dummys = dummys[, colSums(dummys) >= 5]
data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
terms = list(
  s1 = c("s(LnAreaHost, bs='cs', k = 10)"),
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(extinct.host.no, bs='cs', k = 10)"),
  s5 = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
         "s(extinct.median.phylo.dist, bs='cs', k = 10)",
         "s(extinct.min.phylo.dist, bs='cs', k = 10)"##,
##         "s(extinct.host.no, bs='cs', k = 10)",
##         "s(extinct.n.family, bs='cs', k = 10)",
##         "s(extinct.n.order, bs='cs', k = 10)",
##         "s(extinct.n.genera, bs='cs', k = 10)"
  ),
  #s6 = c(#"s(extinct.mean.HR, bs='cs', k = 10)",
  #"s(extinct.median.HR, bs='cs', k = 10)",
  #"s(extinct.mean.FD, bs='cs', k = 10)",
  #"s(extinct.median.FD, bs='cs', k = 10)",
  #"s(extinct.mean.DR, bs='cs', k = 10)",
  #"s(extinct.median.DR, bs='cs', k = 10)"),
  s6 = c("s(extinct.host.no, bs='cs', k = 10)",
         "s(extinct.n.family, bs='cs', k = 10)"#, 
         #"s(extinct.n.order, bs='cs', k = 10)",
         #"s(extinct.n.genera, bs='cs', k = 10)"
  ),
  s7 = c("s(Plourde_PC1, bs='cs', k = 10)"), ###
  s8 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  #  s10 = c("s(Sympat_min_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_mean_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_median_PhyloDist, bs='cs', k = 7)" 
  #        ),
  #s8 = c("s(Phylacine_extinct_over_extantnatural, bs='cs', k = 7)"),
  #s9 = c(
    #"s(CR, bs='cs', k = 10)", ###
    #"s(EN, bs='cs', k = 10)", ###
    #"s(VU, bs='cs', k = 10)", ###
    #"s(total_threat, bs='cs', k = 20)"
  #),
  s10 = c("s(Sympat_n_Phylacine, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

terms = c(dummy_terms, terms)

saveRDS(data_set, paste0(base_dir, "intermediates/Shaw_virus_dataset.rds"))

#shaw_virus <- fit_all_gams_qpoisson(data_set,
#                                   outcome_variable = "n_virus_Shaw",
#                                   terms)
#saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_models_qpoisson.rds"))

shaw_virus <- fit_all_gams_poisson(data_set,
                                   outcome_variable = "n_virus_Shaw",
                                   terms)
saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_models_poisson.rds"))

shaw_virus <- fit_all_gams_tw(data_set,
                                   outcome_variable = "n_virus_Shaw",
                                   terms)
saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_models_tw.rds"))

#shaw_virus <- fit_all_gams_nb(data_set,
#                                   outcome_variable = "n_virus_Shaw",
#                                   terms)
#saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_models_nb.rds"))

#shaw_virus <- fit_all_gams_scat(data_set,
#                              outcome_variable = "n_virus_Shaw",
#                              terms)
#saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_models_scat.rds"))

#shaw_virus <- fit_all_gams_gamma(data_set,
#                                outcome_variable = "n_virus_Shaw",
#                                terms)
#saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_models_gamma.rds"))

terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 10)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(Plourde_PC1, bs='cs', k = 10)"), 
  #s5 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

terms = c(dummy_terms, terms)

shaw_virus_noExtinction <- fit_all_gams_poisson(data_set,
                                                outcome_variable = "n_virus_Shaw",
                                                terms)
saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_models_poisson.rds"))

#shaw_virus_noExtinction <- fit_all_gams_qpoisson(data_set,
#                                                outcome_variable = "n_virus_Shaw",
#                                                terms)
#saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_models_qpoisson.rds"))

shaw_virus_noExtinction <- fit_all_gams_tw(data_set,
                                              outcome_variable = "n_virus_Shaw",
                                              terms)
saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_models_tw.rds"))

#shaw_virus_noExtinction <- fit_all_gams_nb(data_set,
#                                           outcome_variable = "n_virus_Shaw",
#                                           terms)
#saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_models_nb.rds"))

#shaw_virus_noExtinction <- fit_all_gams_scat(data_set,
#                                                outcome_variable = "n_virus_Shaw",
#                                                terms)
#saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_models_scat.rds"))

#shaw_virus_noExtinction <- fit_all_gams_gamma(data_set,
#                                             outcome_variable = "n_virus_Shaw",
#                                             terms)
#saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_models_gamma.rds"))






####---- Shaw Virus RNA GAM - All Associations ----
data_set <- hosts.virus.shaw.rna %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))

data_set <- data_set %>%
  tidyr::drop_na(c(LnAreaHost,
                   hMassGramsPVR,
                   S100, 
                   S80, 
                   S50,
                   S40,
                   S20,
                   S,
                   hDiseaseZACitesLn,
                   extinct.host.no,
                   Plourde_PC1,
                   Plourde_PC2,
                   #extinct.total.phylo.dist,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   #extinct.mean.biomass,
                   #extinct.median.biomass,
                   #extinct.mean.DR,
                   #extinct.median.DR,
                   #extinct.mean.FD,
                   #extinct.median.FD,
                   #extinct.mean.HR,
                   #extinct.median.HR,
                   #CR, CR_prop,
                   #EN, EN_prop,
                   #VU, VU_prop,
                   #total_threat, #total_threat_prop,
                   Sympat_n_Phylacine,
                   #Sympat_min_PhyloDist,
                   #Sympat_mean_PhyloDist,
                   #Sympat_median_PhyloDist,
                   #Phylacine_extinct_over_extantnatural,
                   #extinct.n.order,
                   #extinct.n.genera,
                   extinct.n.family
                   
  ))


data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
dummys = dummys[, colSums(dummys) >= 5]
data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
terms = list(
  s1 = c("s(LnAreaHost, bs='cs', k = 10)"),
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(extinct.host.no, bs='cs', k = 10)"),
  s5 = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
         "s(extinct.median.phylo.dist, bs='cs', k = 10)",
         "s(extinct.min.phylo.dist, bs='cs', k = 10)"##,
         ##         "s(extinct.host.no, bs='cs', k = 10)",
         ##         "s(extinct.n.family, bs='cs', k = 10)",
         ##         "s(extinct.n.order, bs='cs', k = 10)",
         ##         "s(extinct.n.genera, bs='cs', k = 10)"
  ),
  #s6 = c(#"s(extinct.mean.HR, bs='cs', k = 10)",
  #"s(extinct.median.HR, bs='cs', k = 10)",
  #"s(extinct.mean.FD, bs='cs', k = 10)",
  #"s(extinct.median.FD, bs='cs', k = 10)",
  #"s(extinct.mean.DR, bs='cs', k = 10)",
  #"s(extinct.median.DR, bs='cs', k = 10)"),
  s6 = c("s(extinct.host.no, bs='cs', k = 10)",
         "s(extinct.n.family, bs='cs', k = 10)"#, 
         #"s(extinct.n.order, bs='cs', k = 10)",
         #"s(extinct.n.genera, bs='cs', k = 10)"
  ),
  s7 = c("s(Plourde_PC1, bs='cs', k = 10)"), ###
  s8 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  #  s10 = c("s(Sympat_min_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_mean_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_median_PhyloDist, bs='cs', k = 7)" 
  #        ),
  #s8 = c("s(Phylacine_extinct_over_extantnatural, bs='cs', k = 7)"),
  #s9 = c(
  #"s(CR, bs='cs', k = 10)", ###
  #"s(EN, bs='cs', k = 10)", ###
  #"s(VU, bs='cs', k = 10)", ###
  #"s(total_threat, bs='cs', k = 20)"
  #),
  s10 = c("s(Sympat_n_Phylacine, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

terms = c(dummy_terms, terms)


#shaw_virus <- fit_all_gams_qpoisson(data_set,
#                                   outcome_variable = "n_virus_Shaw",
#                                   terms)
#saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_models_qpoisson.rds"))

shaw_virus <- fit_all_gams_poisson(data_set,
                                   outcome_variable = "n_virus_Shaw_RNA",
                                   terms)
saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_RNA_models_poisson.rds"))

shaw_virus <- fit_all_gams_tw(data_set,
                              outcome_variable = "n_virus_Shaw_RNA",
                              terms)
saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_RNA_models_tw.rds"))


terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 10)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  s4 = c("s(Plourde_PC1, bs='cs', k = 10)"), 
  s5 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

terms = c(dummy_terms, terms)

shaw_virus_noExtinction <- fit_all_gams_poisson(data_set,
                                                outcome_variable = "n_virus_Shaw_RNA",
                                                terms)
saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_RNA_models_poisson.rds"))

shaw_virus_noExtinction <- fit_all_gams_tw(data_set,
                                           outcome_variable = "n_virus_Shaw_RNA",
                                           terms)
saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_RNA_models_tw.rds"))





####---- Shaw Virus DNA GAM - All Associations ----
data_set <- hosts.virus.shaw.dna %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))

data_set <- data_set %>%
  tidyr::drop_na(c(LnAreaHost,
                   hMassGramsPVR,
                   S100, 
                   S80, 
                   S50,
                   S40,
                   S20,
                   S,
                   hDiseaseZACitesLn,
                   extinct.host.no,
                   Plourde_PC1,
                   Plourde_PC2,
                   #extinct.total.phylo.dist,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   #extinct.mean.biomass,
                   #extinct.median.biomass,
                   #extinct.mean.DR,
                   #extinct.median.DR,
                   #extinct.mean.FD,
                   #extinct.median.FD,
                   #extinct.mean.HR,
                   #extinct.median.HR,
                   #CR, CR_prop,
                   #EN, EN_prop,
                   #VU, VU_prop,
                   #total_threat, #total_threat_prop,
                   Sympat_n_Phylacine,
                   #Sympat_min_PhyloDist,
                   #Sympat_mean_PhyloDist,
                   #Sympat_median_PhyloDist,
                   #Phylacine_extinct_over_extantnatural,
                   #extinct.n.order,
                   #extinct.n.genera,
                   extinct.n.family
                   
  ))


data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
dummys = dummys[, colSums(dummys) >= 5]
data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
terms = list(
  s1 = c("s(LnAreaHost, bs='cs', k = 10)"),
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(extinct.host.no, bs='cs', k = 10)"),
  s5 = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
         "s(extinct.median.phylo.dist, bs='cs', k = 10)",
         "s(extinct.min.phylo.dist, bs='cs', k = 10)"##,
         ##         "s(extinct.host.no, bs='cs', k = 10)",
         ##         "s(extinct.n.family, bs='cs', k = 10)",
         ##         "s(extinct.n.order, bs='cs', k = 10)",
         ##         "s(extinct.n.genera, bs='cs', k = 10)"
  ),
  #s6 = c(#"s(extinct.mean.HR, bs='cs', k = 10)",
  #"s(extinct.median.HR, bs='cs', k = 10)",
  #"s(extinct.mean.FD, bs='cs', k = 10)",
  #"s(extinct.median.FD, bs='cs', k = 10)",
  #"s(extinct.mean.DR, bs='cs', k = 10)",
  #"s(extinct.median.DR, bs='cs', k = 10)"),
  s6 = c("s(extinct.host.no, bs='cs', k = 10)",
         "s(extinct.n.family, bs='cs', k = 10)"#, 
         #"s(extinct.n.order, bs='cs', k = 10)",
         #"s(extinct.n.genera, bs='cs', k = 10)"
  ),
  s7 = c("s(Plourde_PC1, bs='cs', k = 10)"), ###
  s8 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  #  s10 = c("s(Sympat_min_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_mean_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_median_PhyloDist, bs='cs', k = 7)" 
  #        ),
  #s8 = c("s(Phylacine_extinct_over_extantnatural, bs='cs', k = 7)"),
  #s9 = c(
  #"s(CR, bs='cs', k = 10)", ###
  #"s(EN, bs='cs', k = 10)", ###
  #"s(VU, bs='cs', k = 10)", ###
  #"s(total_threat, bs='cs', k = 20)"
  #),
  s10 = c("s(Sympat_n_Phylacine, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

terms = c(dummy_terms, terms)


shaw_virus <- fit_all_gams_nb(data_set,
                                   outcome_variable = "n_virus_Shaw_DNA",
                                   terms)
saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_DNA_models_nb.rds"))

shaw_virus <- fit_all_gams_poisson(data_set,
                                   outcome_variable = "n_virus_Shaw_DNA",
                                   terms)
saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_DNA_models_poisson.rds"))

shaw_virus <- fit_all_gams_tw(data_set,
                              outcome_variable = "n_virus_Shaw_DNA",
                              terms)
saveRDS(shaw_virus, paste0(base_dir, "intermediates/shaw_virus_DNA_models_tw.rds"))


terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 10)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  s4 = c("s(Plourde_PC1, bs='cs', k = 10)"), 
  s5 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

terms = c(dummy_terms, terms)

shaw_virus_noExtinction <- fit_all_gams_nb(data_set,
                                                outcome_variable = "n_virus_Shaw_DNA",
                                                terms)
saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_DNA_models_nb.rds"))

shaw_virus_noExtinction <- fit_all_gams_poisson(data_set,
                                                outcome_variable = "n_virus_Shaw_DNA",
                                                terms)
saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_DNA_models_poisson.rds"))

shaw_virus_noExtinction <- fit_all_gams_tw(data_set,
                                           outcome_variable = "n_virus_Shaw_DNA",
                                           terms)
saveRDS(shaw_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_virus_DNA_models_tw.rds"))





####---- Shaw Bacteria GAM - All Associations ----
data_set <- hosts.bacteria.shaw %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))

data_set <- data_set %>%
  tidyr::drop_na(c(LnAreaHost,
                   hMassGramsPVR,
                   S100, 
                   S80, 
                   S50,
                   S40,
                   S20,
                   S,
                   hDiseaseZACitesLn,
                   extinct.host.no,
                   #Plourde_PC1,
                   #Plourde_PC2,
                   #extinct.total.phylo.dist,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   #extinct.mean.biomass,
                   #extinct.median.biomass,
                   #extinct.mean.DR,
                   #extinct.median.DR,
                   #extinct.mean.FD,
                   #extinct.median.FD,
                   #extinct.mean.HR,
                   #extinct.median.HR,
                   #CR, CR_prop,
                   #EN, EN_prop,
                   #VU, VU_prop,
                   #total_threat, #total_threat_prop,
                   Sympat_n_Phylacine,
                   #Sympat_min_PhyloDist,
                   #Sympat_mean_PhyloDist,
                   #Sympat_median_PhyloDist,
                   #Phylacine_extinct_over_extantnatural,
                   #extinct.n.order,
                   #extinct.n.genera
                   extinct.n.family
                   
  ))


shaw_hostbacteria <- shaw_hostbacteria[which(shaw_hostbacteria$HostSpecies %in% data_set$hHostNameFinal),]
print(paste0("Number of unique hosts - bacteria associations, from Shaw et al. (2020), used in this study: ", nrow(shaw_hostbacteria)))
print(paste0("Unique hosts: ", nrow(data_set), "; Unique bacteria: ", length(unique(shaw_hostbacteria$Species))))



data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
# The problem (from pairs plot) is hOrderSCANDENTIA and hOrderPERAMELEMORPHIA which both have zero
dummys = dummys[, colSums(dummys) >= 5]

data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
############# NEEDS ADJUSTING TO INCLUDE EXTINCT MAMMAL METRICS
terms = list(
  s1 = c("s(LnAreaHost, bs='cs', k = 20)"),
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 20)",
    "s(S80, bs='cs', k = 20)",
    "s(S50, bs='cs', k = 20)",
    "s(S40, bs='cs', k = 20)",
    "s(S20, bs='cs', k = 20)",
    "s(S, bs='cs', k = 20)"
  ),
  #s4 = c("s(extinct.host.no, bs='cs', k = 10)"),
  s5 = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
         "s(extinct.median.phylo.dist, bs='cs', k = 10)",
         "s(extinct.min.phylo.dist, bs='cs', k = 10)"##,
         ##         "s(extinct.host.no, bs='cs', k = 10)",
         ##         "s(extinct.n.family, bs='cs', k = 10)",
         ##         "s(extinct.n.order, bs='cs', k = 10)",
         ##         "s(extinct.n.genera, bs='cs', k = 10)"
  ),
  #s6 = c(#"s(extinct.mean.HR, bs='cs', k = 10)",
  #"s(extinct.median.HR, bs='cs', k = 10)",
  #"s(extinct.mean.FD, bs='cs', k = 10)",
  #"s(extinct.median.FD, bs='cs', k = 10)",
  #"s(extinct.mean.DR, bs='cs', k = 10)",
  #"s(extinct.median.DR, bs='cs', k = 10)"),
  s6 = c("s(extinct.host.no, bs='cs', k = 10)",
         "s(extinct.n.family, bs='cs', k = 10)"#, 
         #"s(extinct.n.order, bs='cs', k = 10)",
         #"s(extinct.n.genera, bs='cs', k = 10)"
  ),
  s7 = c("s(Plourde_PC1, bs='cs', k = 10)"), ###
  s8 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  #  s10 = c("s(Sympat_min_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_mean_PhyloDist, bs='cs', k = 7)", 
  #        "s(Sympat_median_PhyloDist, bs='cs', k = 7)" 
  #        ),
  #s8 = c("s(Phylacine_extinct_over_extantnatural, bs='cs', k = 7)"),
  #s9 = c(
  #"s(CR, bs='cs', k = 10)", ###
  #"s(EN, bs='cs', k = 10)", ###
  #"s(VU, bs='cs', k = 10)", ###
      #"s(total_threat, bs='cs', k = 20)"
  #),
  s10 = c("s(Sympat_n_Phylacine, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

saveRDS(data_set, paste0(base_dir, "intermediates/Shaw_bacteria_dataset.rds"))

shaw_bact <- fit_all_gams_poisson(data_set,
                                  outcome_variable = "n_bact_Shaw",
                                  terms)
saveRDS(shaw_bact, paste0(base_dir, "intermediates/shaw_bact_models_poisson.rds"))

#shaw_bact <- fit_all_gams_qpoisson(data_set,
#                                  outcome_variable = "n_bact_Shaw",
#                                  terms)
#saveRDS(shaw_bact, paste0(base_dir, "intermediates/shaw_bact_models_qpoisson.rds"))

shaw_bact <- fit_all_gams_tw(data_set,
                                  outcome_variable = "n_bact_Shaw",
                                  terms)
saveRDS(shaw_bact, paste0(base_dir, "intermediates/shaw_bact_models_tw.rds"))

#shaw_bact <- fit_all_gams_nb(data_set,
#                                  outcome_variable = "n_bact_Shaw",
#                                  terms)
#saveRDS(shaw_bact, paste0(base_dir, "intermediates/shaw_bact_models_nb.rds"))

#shaw_bact <- fit_all_gams_scat(data_set,
#                             outcome_variable = "n_bact_Shaw",
#                             terms)
#saveRDS(shaw_bact, paste0(base_dir, "intermediates/shaw_bact_models_scat.rds"))

#shaw_bact <- fit_all_gams_gamma(data_set,
#                               outcome_variable = "n_bact_Shaw",
#                               terms)
#saveRDS(shaw_bact, paste0(base_dir, "intermediates/shaw_bact_models_gamma.rds"))


terms = list(
  s1 = "s(LnAreaHost, bs='cs', k = 20)",
  s2 = "s(hMassGramsPVR, bs='cs', k = 10)",
  s3 = c(
    "s(S100, bs='cs', k = 10)",
    "s(S80, bs='cs', k = 10)",
    "s(S50, bs='cs', k = 10)",
    "s(S40, bs='cs', k = 10)",
    "s(S20, bs='cs', k = 10)",
    "s(S, bs='cs', k = 10)"
  ),
  #s4 = c("s(Plourde_PC1, bs='cs', k = 10)"), 
  #s5 = c("s(Plourde_PC2, bs='cs', k = 10)"),
  bias = c("s(hDiseaseZACitesLn, bs='cs', k = 20)"),
  stringsAsFactors = FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)


shaw_bact_noExtinction <- fit_all_gams_tw(data_set,
                                          outcome_variable = "n_bact_Shaw",
                                          terms)
saveRDS(shaw_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_bact_models_tw.rds"))

shaw_bact_noExtinction <- fit_all_gams_poisson(data_set,
                                          outcome_variable = "n_bact_Shaw",
                                          terms)
saveRDS(shaw_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_bact_models_poisson.rds"))

#shaw_bact_noExtinction <- fit_all_gams_qpoisson(data_set,
#                                               outcome_variable = "n_bact_Shaw",
#                                               terms)
#saveRDS(shaw_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_bact_models_qpoisson.rds"))

#shaw_bact_noExtinction <- fit_all_gams_nb(data_set,
#                                               outcome_variable = "n_bact_Shaw",
#                                               terms)
#saveRDS(shaw_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_bact_models_nb.rds"))

#shaw_bact_noExtinction <- fit_all_gams_scat(data_set,
#                                          outcome_variable = "n_bact_Shaw",
#                                          terms)
#saveRDS(shaw_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_bact_models_scat.rds"))

#shaw_bact_noExtinction <- fit_all_gams_scat(data_set,
#                                          outcome_variable = "n_bact_Shaw",
#                                          terms)
#saveRDS(shaw_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_bact_models_scat.rds"))




####---- Shaw ZOONOTIC Virus GAM - variables used by Olival et al. (2017) ----
data_set <- hosts.virus.shaw.zoo %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))

data_set$ln_virus_Shaw <- log(data_set$n_virus_Shaw)

data_set <- data_set %>%
  tidyr::drop_na(c(hMassGramsPVR,
                   HabAreaCropLn,
                   HabAreaCropChgLn,
                   HabAreaUrbanLn,
                   HabAreaUrbanChgLn,
                   HabAreaGrassLn,
                   HabAreaGrassChgLn,
                   HabInhabitedLn,
                   HabInhabitedChgLn,
                   TotHumPopLn, 
                   TotHumPopChgLn,
                   UrbRurPopRatioLn, 
                   UrbRurPopRatioChg, 
                   HumPopDensLn,
                   HumPopDensLnChg,
                   UrbRurPopRatioLn,
                   UrbRurPopRatioChg,
                   hHuntedIUCN,
                   hArtfclHbttUsrIUCN,
                   PdHoSa.cbCst,
                   PdHoSaSTPD,
                   hDiseaseZACitesLn,
                   ln_virus_Shaw,
                   LnTotNumVirus,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   extinct.host.no,
                   extinct.n.family
                   
  ))



#  Create dummy variables for orders to use as random effects
# But note that these become non-conformable
# Could just leave them out?
data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
# The problem (from pairs plot) is hOrderSCANDENTIA and hOrderPERAMELEMORPHIA which both have zero
dummys = dummys[, colSums(dummys) >= 5]


data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
############# NEEDS ADJUSTING TO INCLUDE EXTINCT MAMMAL METRICS

terms = list(
  mass = "s(hMassGramsPVR, bs = 'tp', k=7)",
  interaction = c(
    "s(HabAreaCropLn, bs = 'tp', k=7) + s(HabAreaCropChgLn, bs = 'tp', k=7)",
    "s(HabAreaGrassLn, bs = 'tp', k=7)  + s(HabAreaGrassChgLn, bs = 'tp', k=7)",
    "s(HabAreaUrbanLn, bs = 'tp', k=7)  + s(HabAreaUrbanChgLn, bs = 'tp', k=7)",
    "s(HabInhabitedLn, bs = 'tp', k=7)  + s(HabInhabitedChgLn, bs = 'tp', k=7)",
    "s(TotHumPopLn, bs = 'tp', k=7) + s(TotHumPopChgLn, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)",
    "s(HumPopDensLn, bs = 'tp', k=7) + s(HumPopDensLnChg, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)"),
  interaction2 = "s(hHuntedIUCN, bs='re')",
  interaction3 = "s(hArtfclHbttUsrIUCN, bs='re')",
  phylo_distance = c("s(PdHoSa.cbCst, bs = 'tp', k=7)", "s(PdHoSaSTPD, bs = 'tp', k=7)"),
  extinction = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
                 "s(extinct.median.phylo.dist, bs='cs', k = 10)",
                 "s(extinct.min.phylo.dist, bs='cs', k = 10)",
                 "s(extinct.host.no, bs='cs', k = 10)",
                 "s(extinct.n.family, bs='cs', k = 10)"
  ),
  bias = c("s(hDiseaseZACitesLn, bs = 'tp', k=7)"),
  offset = "offset(LnTotNumVirus)",
  stringsAsFactors=FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

shaw_virus_zoo <- fit_all_gams_poisson(data_set,
                                       outcome_variable = "n_virus_Shaw_zoo",
                                       terms)
saveRDS(shaw_virus_zoo, paste0(base_dir, "intermediates/shaw_zoo_virus_models2_poisson.rds"))

shaw_virus_zoo <- fit_all_gams_tw(data_set,
                                  outcome_variable = "n_virus_Shaw_zoo",
                                  terms)
saveRDS(shaw_virus_zoo, paste0(base_dir, "intermediates/shaw_zoo_virus_models2_tw.rds"))

shaw_virus_zoo <- fit_all_gams_qpoisson(data_set,
                                        outcome_variable = "n_virus_Shaw_zoo",
                                        terms)
saveRDS(shaw_virus_zoo, paste0(base_dir, "intermediates/shaw_zoo_virus_models2_qpoisson.rds"))

shaw_virus_zoo <- fit_all_gams_nb(data_set,
                                  outcome_variable = "n_virus_Shaw_zoo",
                                  terms)
saveRDS(shaw_virus_zoo, paste0(base_dir, "intermediates/shaw_zoo_virus_models2_nb.rds"))

shaw_virus_zoo <- fit_all_gams_gau(data_set,
                                   outcome_variable = "n_virus_Shaw_zoo",
                                   terms)
saveRDS(shaw_virus_zoo, paste0(base_dir, "intermediates/shaw_zoo_virus_models2_gau.rds"))

shaw_virus_zoo <- fit_all_gams_scat(data_set,
                                    outcome_variable = "n_virus_Shaw_zoo",
                                    terms)
saveRDS(shaw_virus_zoo, paste0(base_dir, "intermediates/shaw_zoo_virus_models2_scat.rds"))


terms = list(
  mass = "s(hMassGramsPVR, bs = 'tp', k=7)",
  interaction = c(
    "s(HabAreaCropLn, bs = 'tp', k=7)   + s(HabAreaCropChgLn, bs = 'tp', k=7)",
    "s(HabAreaGrassLn, bs = 'tp', k=7)  + s(HabAreaGrassChgLn, bs = 'tp', k=7)",
    "s(HabAreaUrbanLn, bs = 'tp', k=7)  + s(HabAreaUrbanChgLn, bs = 'tp', k=7)",
    "s(HabInhabitedLn, bs = 'tp', k=7)  + s(HabInhabitedChgLn, bs = 'tp', k=7)",
    "s(TotHumPopLn, bs = 'tp', k=7) + s(TotHumPopChgLn, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)",
    "s(HumPopDensLn, bs = 'tp', k=7) + s(HumPopDensLnChg, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)"),
  interaction2 = "s(hHuntedIUCN, bs='re')",
  interaction3 = "s(hArtfclHbttUsrIUCN, bs='re')",
  phylo_distance = c("s(PdHoSa.cbCst, bs = 'tp', k=7)", "s(PdHoSaSTPD, bs = 'tp', k=7)"),
  bias = c("s(hDiseaseZACitesLn, bs = 'tp', k=7)"),
  offset = "offset(LnTotNumVirus)",
  stringsAsFactors=FALSE)
# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

shaw_zoo_virus_noExtinction <- fit_all_gams_tw(data_set,
                                               outcome_variable = "n_virus_Shaw_zoo",
                                               terms)
saveRDS(shaw_zoo_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_virus_models2_tw.rds"))

shaw_zoo_virus_noExtinction <- fit_all_gams_poisson(data_set,
                                                    outcome_variable = "n_virus_Shaw_zoo",
                                                    terms)
saveRDS(shaw_zoo_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_virus_models2_poisson.rds"))

shaw_zoo_virus_noExtinction <- fit_all_gams_qpoisson(data_set,
                                                     outcome_variable = "n_virus_Shaw_zoo",
                                                     terms)
saveRDS(shaw_zoo_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_virus_models2_qpoisson.rds"))

shaw_zoo_virus_noExtinction <- fit_all_gams_nb(data_set,
                                               outcome_variable = "n_virus_Shaw_zoo",
                                               terms)
saveRDS(shaw_zoo_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_virus_models2_nb.rds"))

shaw_zoo_virus_noExtinction <- fit_all_gams_gau(data_set,
                                                outcome_variable = "n_virus_Shaw_zoo",
                                                terms)
saveRDS(shaw_zoo_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_virus_models2_gau.rds"))

shaw_zoo_virus_noExtinction <- fit_all_gams_scat(data_set,
                                                 outcome_variable = "n_virus_Shaw_zoo",
                                                 terms)
saveRDS(shaw_zoo_virus_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_virus_models2_scat.rds"))

#ggplot(data = data_set) +
#  geom_point(aes(x = LnTotNumVirus, y = ln_virus_Shaw, color = hHostNameFinal),
#             position=position_jitter(h=0.2,w=0)) +
#  geom_smooth(aes(x = LnTotNumVirus, y = ln_virus_Shaw)) +
#  theme(legend.position = "none")







####---- Shaw ZOONOTIC Virus GAM - variables used by Olival et al. (2017) ----
data_set <- hosts.bacteria.shaw.zoo %>%
  filter(hMarOTerr == "Terrestrial",
         hWildDomFAO == "wild",
         !is.na(PdHoSa.cbCst_order))

data_set$ln_bact_Shaw <- log(data_set$n_bact_Shaw)

data_set <- data_set %>%
  tidyr::drop_na(c(hMassGramsPVR,
                   HabAreaCropLn,
                   HabAreaCropChgLn,
                   HabAreaUrbanLn,
                   HabAreaUrbanChgLn,
                   HabAreaGrassLn,
                   HabAreaGrassChgLn,
                   HabInhabitedLn,
                   HabInhabitedChgLn,
                   TotHumPopLn, 
                   TotHumPopChgLn,
                   UrbRurPopRatioLn, 
                   UrbRurPopRatioChg, 
                   HumPopDensLn,
                   HumPopDensLnChg,
                   UrbRurPopRatioLn,
                   UrbRurPopRatioChg,
                   hHuntedIUCN,
                   hArtfclHbttUsrIUCN,
                   PdHoSa.cbCst,
                   PdHoSaSTPD,
                   hDiseaseZACitesLn,
                   ln_bact_Shaw,
                   extinct.mean.phylo.dist,
                   extinct.median.phylo.dist,
                   extinct.min.phylo.dist,
                   extinct.host.no,
                   extinct.n.family
                   
  ))



#  Create dummy variables for orders to use as random effects
# But note that these become non-conformable
# Could just leave them out?
data_set$hOrder <- droplevels(data_set$hOrder)
dummys = as.data.frame(with(data_set, model.matrix(~hOrder))[,-1])
# The problem (from pairs plot) is hOrderSCANDENTIA and hOrderPERAMELEMORPHIA which both have zero
dummys = dummys[, colSums(dummys) >= 5]


data_set = cbind(data_set, dummys)
dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
names(dummy_terms) <- names(dummys)

## Create data.frame of all possible models
############# NEEDS ADJUSTING TO INCLUDE EXTINCT MAMMAL METRICS

terms = list(
  mass = "s(hMassGramsPVR, bs = 'tp', k=7)",
  interaction = c(
    "s(HabAreaCropLn, bs = 'tp', k=7) + s(HabAreaCropChgLn, bs = 'tp', k=7)",
    "s(HabAreaGrassLn, bs = 'tp', k=7)  + s(HabAreaGrassChgLn, bs = 'tp', k=7)",
    "s(HabAreaUrbanLn, bs = 'tp', k=7)  + s(HabAreaUrbanChgLn, bs = 'tp', k=7)",
    "s(HabInhabitedLn, bs = 'tp', k=7)  + s(HabInhabitedChgLn, bs = 'tp', k=7)",
    "s(TotHumPopLn, bs = 'tp', k=7) + s(TotHumPopChgLn, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)",
    "s(HumPopDensLn, bs = 'tp', k=7) + s(HumPopDensLnChg, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)"),
  interaction2 = "s(hHuntedIUCN, bs='re')",
  interaction3 = "s(hArtfclHbttUsrIUCN, bs='re')",
  phylo_distance = c("s(PdHoSa.cbCst, bs = 'tp', k=7)", "s(PdHoSaSTPD, bs = 'tp', k=7)"),
  extinction = c("s(extinct.mean.phylo.dist, bs='cs', k = 10)",
                 "s(extinct.median.phylo.dist, bs='cs', k = 10)",
                 "s(extinct.min.phylo.dist, bs='cs', k = 10)",
                 "s(extinct.host.no, bs='cs', k = 10)",
                 "s(extinct.n.family, bs='cs', k = 10)"
  ),
  bias = c("s(hDiseaseZACitesLn, bs = 'tp', k=7)"),
  offset = "offset(ln_bact_Shaw)",
  stringsAsFactors=FALSE)

# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

shaw_bact_zoo <- fit_all_gams_poisson(data_set,
                                      outcome_variable = "n_bact_Shaw_zoo",
                                      terms)
saveRDS(shaw_bact_zoo, paste0(base_dir, "intermediates/shaw_zoo_bact_models2_poisson.rds"))

shaw_bact_zoo <- fit_all_gams_tw(data_set,
                                 outcome_variable = "n_bact_Shaw_zoo",
                                 terms)
saveRDS(shaw_bact_zoo, paste0(base_dir, "intermediates/shaw_zoo_bact_models2_tw.rds"))

shaw_bact_zoo <- fit_all_gams_qpoisson(data_set,
                                       outcome_variable = "n_bact_Shaw_zoo",
                                       terms)
saveRDS(shaw_bact_zoo, paste0(base_dir, "intermediates/shaw_zoo_bact_models2_qpoisson.rds"))

shaw_bact_zoo <- fit_all_gams_nb(data_set,
                                 outcome_variable = "n_bact_Shaw_zoo",
                                 terms)
saveRDS(shaw_bact_zoo, paste0(base_dir, "intermediates/shaw_zoo_bact_models2_nb.rds"))

shaw_bact_zoo <- fit_all_gams_gau(data_set,
                                  outcome_variable = "n_bact_Shaw_zoo",
                                  terms)
saveRDS(shaw_bact_zoo, paste0(base_dir, "intermediates/shaw_zoo_bact_models2_gau.rds"))

shaw_bact_zoo <- fit_all_gams_scat(data_set,
                                   outcome_variable = "n_bact_Shaw_zoo",
                                   terms)
saveRDS(shaw_bact_zoo, paste0(base_dir, "intermediates/shaw_zoo_bact_models2_scat.rds"))


terms = list(
  mass = "s(hMassGramsPVR, bs = 'tp', k=7)",
  interaction = c(
    "s(HabAreaCropLn, bs = 'tp', k=7)   + s(HabAreaCropChgLn, bs = 'tp', k=7)",
    "s(HabAreaGrassLn, bs = 'tp', k=7)  + s(HabAreaGrassChgLn, bs = 'tp', k=7)",
    "s(HabAreaUrbanLn, bs = 'tp', k=7)  + s(HabAreaUrbanChgLn, bs = 'tp', k=7)",
    "s(HabInhabitedLn, bs = 'tp', k=7)  + s(HabInhabitedChgLn, bs = 'tp', k=7)",
    "s(TotHumPopLn, bs = 'tp', k=7) + s(TotHumPopChgLn, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)",
    "s(HumPopDensLn, bs = 'tp', k=7) + s(HumPopDensLnChg, bs = 'tp', k=7) + s(UrbRurPopRatioLn, bs = 'tp', k=7) + s(UrbRurPopRatioChg, bs = 'tp', k=7)"),
  interaction2 = "s(hHuntedIUCN, bs='re')",
  interaction3 = "s(hArtfclHbttUsrIUCN, bs='re')",
  phylo_distance = c("s(PdHoSa.cbCst, bs = 'tp', k=7)", "s(PdHoSaSTPD, bs = 'tp', k=7)"),
  bias = c("s(hDiseaseZACitesLn, bs = 'tp', k=7)"),
  offset = "offset(ln_bact_Shaw)",
  stringsAsFactors=FALSE)
# The dummy terms might be causing problems - so leave it out
terms = c(dummy_terms, terms)

shaw_zoo_bact_noExtinction <- fit_all_gams_tw(data_set,
                                              outcome_variable = "n_bact_Shaw_zoo",
                                              terms)
saveRDS(shaw_zoo_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_bact_models2_tw.rds"))

shaw_zoo_bact_noExtinction <- fit_all_gams_poisson(data_set,
                                                   outcome_variable = "n_bact_Shaw_zoo",
                                                   terms)
saveRDS(shaw_zoo_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_bact_models2_poisson.rds"))

shaw_zoo_bact_noExtinction <- fit_all_gams_qpoisson(data_set,
                                                    outcome_variable = "n_bact_Shaw_zoo",
                                                    terms)
saveRDS(shaw_zoo_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_bact_models2_qpoisson.rds"))

shaw_zoo_bact_noExtinction <- fit_all_gams_nb(data_set,
                                              outcome_variable = "n_bact_Shaw_zoo",
                                              terms)
saveRDS(shaw_zoo_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_bact_models2_nb.rds"))

shaw_zoo_bact_noExtinction <- fit_all_gams_gau(data_set,
                                               outcome_variable = "n_bact_Shaw_zoo",
                                               terms)
saveRDS(shaw_zoo_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_bact_models2_gau.rds"))

shaw_zoo_bact_noExtinction <- fit_all_gams_scat(data_set,
                                                outcome_variable = "n_bact_Shaw_zoo",
                                                terms)
saveRDS(shaw_zoo_bact_noExtinction, paste0(base_dir, "intermediates/noExtinction_shaw_zoo_bact_models2_scat.rds"))







