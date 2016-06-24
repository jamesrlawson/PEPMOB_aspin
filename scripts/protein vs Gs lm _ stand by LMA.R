require(readr)
require(ggplot2)

source('scripts/functions.R')

phys <- read_csv('data/Aspinwall_leafphys_160623.csv')
prot <- read_csv('data/Aspinwall_proteomics_160623.csv')
prot <- prot[-97,]

phys_prot <- merge(phys, prot, by = c('id', 'temp_treatment', 'time', 'species'))

phys_prot$temp_treatment <- as.factor(phys_prot$temp_treatment)
phys_prot$time <- factor(phys_prot$time, levels = c('pre_heatwave', 'heatwave', 'post_1day', 'post_6day'))
phys_prot$species <- factor(phys_prot$species, levels = c('eucbot', 'eucsmi', 'eucter', 'euccam'))
phys_prot$log10gs <- log10(phys_prot$gs)
#phys_prot$log10HSPs <- log10(phys_prot$HSPs)

phys_prot_pre <- phys_prot[phys_prot$time == 'pre_heatwave',]
phys_prot_heatwave <- phys_prot[phys_prot$time == 'heatwave',]
phys_prot_post_1day <- phys_prot[phys_prot$time == 'post_1day',]
phys_prot_post_6day <- phys_prot[phys_prot$time == 'post_6day',]

phys_prot[,24:52] <- lapply(phys_prot[,24:52], function(x) x/phys_prot$LMA_g_m2)


##### HSPs #####

# single regression per timepoint

regression <- regression_output(phys_prot, c('time'), log10gs, HSPs)
write.csv(regression, 'output/hsps_vs_log10gs_temp_spp_lumped_170623.csv')
#plot_panel(phys_prot, c('time'), log10gs, HSPs)
write_panel(phys_prot, c('time'), log10gs, HSPs)


# one regression per species, temp treatments lumped
regression <- regression_output(phys_prot, c('time', 'species'), log10gs, HSPs)
write.csv(regression, 'output/hsps_vs_log10gs_temp_lumped_170623.csv')
#plot_panel(phys_prot, c('time', 'species'), log10gs, HSPs)
write_panel(phys_prot, c('time', 'species'), log10gs, HSPs)

# one regression per temp treatment, species lumped

regression <- regression_output(phys_prot, c('time', 'temp_treatment'), log10gs, HSPs)
write.csv(regression, 'output/hsps_vs_log10gs_spp_lumped_170623.csv')
#plot_panel(phys_prot, c('time', 'temp_treatment'), log10gs, HSPs)
write_panel(phys_prot, c('time', 'temp_treatment'), log10gs, HSPs)

# regressions for every combination of species and treatment

regression <- regression_output(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, HSPs)
write.csv(regression, 'output/hsps_vs_log10gs_all_models_170623.csv')
#plot_panel(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, HSPs)
write_panel(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, HSPs)

##### HSPs_protein #####

phys_prot$HSPs_protein <- phys_prot$HSPs + phys_prot$protein

# single regression per timepoint

regression <- regression_output(phys_prot, c('time'), log10gs, HSPs_protein)
write.csv(regression, 'output/HSPs_protein_vs_log10gs_temp_spp_lumped_170623.csv')
#plot_panel(phys_prot, c('time'), log10gs, HSPs_protein)
write_panel(phys_prot, c('time'), log10gs, HSPs_protein)

# one regression per species, temp treatments lumped
regression <- regression_output(phys_prot, c('time', 'species'), log10gs, HSPs_protein)
write.csv(regression, 'output/HSPs_protein_vs_log10gs_temp_lumped_170623.csv')
#plot_panel(phys_prot, c('time', 'species'), log10gs, HSPs_protein)
write_panel(phys_prot, c('time', 'species'), log10gs, HSPs_protein)

# one regression per temp treatment, species lumped

regression <- regression_output(phys_prot, c('time', 'temp_treatment'), log10gs, HSPs_protein)
write.csv(regression, 'output/HSPs_protein_vs_log10gs_spp_lumped_170623.csv')
#plot_panel(phys_prot, c('time', 'temp_treatment'), log10gs, HSPs_protein)
write_panel(phys_prot, c('time', 'temp_treatment'), log10gs, HSPs_protein)

# regressions for every combination of species and treatment

regression <- regression_output(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, HSPs_protein)
write.csv(regression, 'output/HSPs_protein_vs_log10gs_all_models_170623.csv')
#plot_panel(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, HSPs_protein)
write_panel(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, HSPs_protein)

#### abiotic stress + protein folding ####

phys_prot$abioticStress_protein <- phys_prot$abiotic_stress + phys_prot$protein

# single regression per timepoint

regression <- regression_output(phys_prot, c('time'), log10gs, abioticStress_protein)
write.csv(regression, 'output/abioticStress_protein_vs_log10gs_temp_spp_lumped_170623.csv')
#plot_panel(phys_prot, c('time'), log10gs, abioticStress_protein)
write_panel(phys_prot, c('time'), log10gs, abioticStress_protein)

# one regression per species, temp treatments lumped
regression <- regression_output(phys_prot, c('time', 'species'), log10gs, abioticStress_protein)
write.csv(regression, 'output/abioticStress_protein_vs_log10gs_temp_lumped_170623.csv')
#plot_panel(phys_prot, c('time', 'species'), log10gs, abioticStress_protein)
write_panel(phys_prot, c('time', 'species'), log10gs, abioticStress_protein)

# one regression per temp treatment, species lumped

regression <- regression_output(phys_prot, c('time', 'temp_treatment'), log10gs, abioticStress_protein)
write.csv(regression, 'output/abioticStress_protein_vs_log10gs_spp_lumped_170623.csv')
#plot_panel(phys_prot, c('time', 'temp_treatment'), log10gs, abioticStress_protein)
write_panel(phys_prot, c('time', 'temp_treatment'), log10gs, abioticStress_protein)

# regressions for every combination of species and treatment

regression <- regression_output(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, abioticStress_protein)
write.csv(regression, 'output/abioticStress_protein_vs_log10gs_all_models_170623.csv')
#plot_panel(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, abioticStress_protein)
write_panel(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, abioticStress_protein)

#### stress + protein folding ####

phys_prot$stress_protein <- phys_prot$stress + phys_prot$protein

# single regression per timepoint

regression <- regression_output(phys_prot, c('time'), log10gs, abioticStress_protein)
write.csv(regression, 'output/abioticStress_protein_vs_log10gs_temp_spp_lumped_170623.csv')
plot_panel(phys_prot, c('time'), log10gs, abioticStress_protein)
write_panel(phys_prot, c('time'), log10gs, abioticStress_protein)

# one regression per species, temp treatments lumped
regression <- regression_output(phys_prot, c('time', 'species'), log10gs, abioticStress_protein)
write.csv(regression, 'output/abioticStress_protein_vs_log10gs_temp_lumped_170623.csv')
#plot_panel(phys_prot, c('time', 'species'), log10gs, abioticStress_protein)
write_panel(phys_prot, c('time', 'species'), log10gs, abioticStress_protein)

# one regression per temp treatment, species lumped

regression <- regression_output(phys_prot, c('time', 'temp_treatment'), log10gs, abioticStress_protein)
write.csv(regression, 'output/abioticStress_protein_vs_log10gs_spp_lumped_170623.csv')
#plot_panel(phys_prot, c('time', 'temp_treatment'), log10gs, abioticStress_protein)
write_panel(phys_prot, c('time', 'temp_treatment'), log10gs, abioticStress_protein)

# regressions for every combination of species and treatment

regression <- regression_output(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, abioticStress_protein)
write.csv(regression_output, 'output/abioticStress_protein_vs_log10gs_all_models_170623.csv')
#plot_panel(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, abioticStress_protein)
write_panel(phys_prot, c('time', 'temp_treatment', 'species'), log10gs, abioticStress_protein)
