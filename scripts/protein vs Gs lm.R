require(readr)
require(ggplot2)

source('scripts/regression_outputs.R')

phys <- read_csv('data/Aspinwall_leafphys_170623.csv')
prot <- read_csv('data/Aspinwall_proteomics_170623.csv')
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

##### HSPs #####

# single regression per timepoint

ggplot(phys_prot, aes(y=HSPs, x=log10gs)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Heat Shock Protein (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

temp_spp_lumped_HSPs <- regression_output(phys_prot, c('time'), HSPs, log10gs)
write.csv(temp_spp_lumped_HSPs, 'output/hsps_vs_log10gs_temp_spp_lumped_170623.csv')

# one regression per species, temp treatments lumped

ggplot(phys_prot, aes(y=HSPs, x=log10gs, color = species)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Heat Shock Protein (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

temp_lumped_HSPs <- regression_output(phys_prot, c('time', 'species'), HSPs, log10gs)
write.csv(temp_lumped_HSPs, 'output/hsps_vs_log10gs_temp_lumped_170623.csv')

# one regression per temp treatment, species lumped

ggplot(phys_prot, aes(y=HSPs, x=log10gs, shape = temp_treatment, linetype = temp_treatment)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Heat Shock Protein (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

spp_lumped_HSPs <- regression_output(phys_prot, c('time', 'temp_treatment'), HSPs, log10gs)
write.csv(spp_lumped, 'output/hsps_vs_log10gs_spp_lumped_170623.csv')

# regressions for every combination of species and treatment

ggplot(phys_prot, aes(y=HSPs, x=gs, color=species, shape = temp_treatment, linetype = temp_treatment)) +
  geom_point(size = 2) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Heat Shock Protein (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)')

all_models_HSPs <- regression_output(phys_prot, c('time', 'temp_treatment', 'species'), HSPs, log10gs)
write.csv(spp_lumped, 'output/hsps_vs_log10gs_all_models_170623.csv')

# HSPs vs log10gs, data = heatwave minus preheatwave

phys_prot1 <- phys_prot
HSPs_stand <- subset(phys_prot1, time == 'heatwave')$HSPs - subset(phys_prot1, time == 'pre_heatwave')$HSPs
log10gs_stand <- subset(phys_prot1, time == 'heatwave')$log10gs - subset(phys_prot1, time == 'pre_heatwave')$log10gs
plot(HSPs_stand ~ log10gs_stand, main = "HSPs vs log10gs, heatwave vals standardised by pre_heatwave vals")
model <- lm(HSPs_stand ~ log10gs_stand)
abline(model)
summary(model)
plot(model)

####### HSPs + protein folding #######

phys_prot$HSPs_protein <- phys_prot$HSPs + phys_prot$protein

# single regression per timepoint

ggplot(phys_prot, aes(y=HSPs_protein, x=log10gs)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Heat Shock Protein + protein folding (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

temp_spp_lumped_HSPs_protein <- regression_output(phys_prot, c('time'), HSPs_protein, log10gs)
write.csv(temp_spp_lumped_HSPs_protein, 'output/HSPs_protein_vs_log10gs_temp_spp_lumped_170623.csv')

# one regression per species, temp treatments lumped

ggplot(phys_prot, aes(y=HSPs_protein, x=log10gs, color = species)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Heat Shock Protein + protein folding (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

temp_lumped_HSPs_protein <- regression_output(phys_prot, c('time', 'species'), HSPs_protein, log10gs)
write.csv(temp_lumped_HSPs_protein, 'output/HSPs_proteinvs_log10gs_temp_lumped_170623.csv')

# one regression per temp treatment, species lumped

ggplot(phys_prot, aes(y=HSPs_protein, x=log10gs, shape = temp_treatment, linetype = temp_treatment)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Heat Shock Protein + protein folding (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

spp_lumped_HSPs_protein <- regression_output(phys_prot, c('time', 'temp_treatment'), HSPs_protein, log10gs)
write.csv(spp_lumped_HSPs_protein, 'output/HSPs_protein_vs_log10gs_spp_lumped_170623.csv')

# regressions for every combination of species and treatment

ggplot(phys_prot, aes(y=HSPs_protein, x=gs, color=species, shape = temp_treatment, linetype = temp_treatment)) +
  geom_point(size = 2) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Heat Shock Protein + protein folding (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)')

all_models_HSPs_protein <- regression_output(phys_prot, c('time', 'temp_treatment', 'species'), HSPs, log10gs)


#### abiotic stress + protein folding ####

phys_prot$abioticStress_protein <- phys_prot$abiotic_stress + phys_prot$protein

# single regression per timepoint

ggplot(phys_prot, aes(y=abioticStress_protein, x=log10gs)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Protein folding + abiotic stress response (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

temp_spp_lumped_protein_abiotic <- regression_output(phys_prot, c('time'), abioticStress_protein, log10gs)
write.csv(temp_spp_lumped_abioticStress_protein, 'output/abioticStress_protein_vs_log10gs_temp_spp_lumped_170623.csv')

# one regression per species, temp treatments lumped

ggplot(phys_prot, aes(y=abioticStress_protein, x=log10gs, color = species)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Protein folding + abiotic stress response (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

temp_lumped_protein_abiotic <- regression_output(phys_prot, c('time', 'species'), abioticStress_protein, log10gs)
write.csv(temp_lumped_HSPs_protein, 'output/abioticStress_protein_vs_log10gs_temp_lumped_170623.csv')

# one regression per temp treatment, species lumped

ggplot(phys_prot, aes(y=abioticStress_protein, x=log10gs, shape = temp_treatment, linetype = temp_treatment)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Protein folding + abiotic stress response (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

spp_lumped_protein_abiotic <- regression_output(phys_prot, c('time', 'temp_treatment'), abioticStress_protein, log10gs)
write.csv(spp_lumped_HSPs_protein, 'output/abioticStress_protein_vs_log10gs_spp_lumped_170623.csv')

# regressions for every combination of species and treatment

ggplot(phys_prot, aes(y=abioticStress_protein, x=gs, color=species, shape = temp_treatment, linetype = temp_treatment)) +
  geom_point(size = 2) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Protein folding + abiotic stress response (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)')

all_models_protein_abiotic <- regression_output(phys_prot, c('time', 'temp_treatment', 'species'), HSPs, log10gs)


#### stress + protein folding ####

phys_prot$stress_protein <- phys_prot$stress + phys_prot$protein

# single regression per timepoint

ggplot(phys_prot, aes(y=stress_protein, x=log10gs)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Protein folding + stress response (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

temp_spp_lumped_protein_stress <- regression_output(phys_prot, c('time'), stress_protein, log10gs)
write.csv(temp_spp_lumped_stress_protein, 'output/stress_protein_vs_log10gs_temp_spp_lumped_170623.csv')

# one regression per species, temp treatments lumped

ggplot(phys_prot, aes(y=stress_protein, x=log10gs, color = species)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Protein folding + stress response (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

temp_lumped_protein_stress <- regression_output(phys_prot, c('time', 'species'), stress_protein, log10gs)
write.csv(temp_lumped_HSPs_protein, 'output/stress_protein_vs_log10gs_temp_lumped_170623.csv')

# one regression per temp treatment, species lumped

ggplot(phys_prot, aes(y=stress_protein, x=log10gs, shape = temp_treatment, linetype = temp_treatment)) +
  geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Protein folding + stress response (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)') 

spp_lumped_protein_stress <- regression_output(phys_prot, c('time', 'temp_treatment'), stress_protein, log10gs)
write.csv(spp_lumped_HSPs_protein, 'output/stress_protein_vs_log10gs_spp_lumped_170623.csv')

# regressions for every combination of species and treatment

ggplot(phys_prot, aes(y=stress_protein, x=gs, color=species, shape = temp_treatment, linetype = temp_treatment)) +
  geom_point(size = 2) +
  geom_smooth(method=lm,se=FALSE) +
  scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
  facet_grid(. ~ time, scales = "free_x") +
  ylab('Protein folding + stress response (mg / m2)') +
  xlab('log10(Gs) (mmol H20 / s / m2)')

all_models_protein_stress <- regression_output(phys_prot, c('time', 'temp_treatment', 'species'), HSPs, log10gs)

