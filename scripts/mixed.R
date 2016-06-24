require(lme4)

?lme


mod <- lmer(log10gs ~ HSPs * temp_treatment + (1|time) + (1|species), data = phys_prot, REML=FALSE)
summary(mod)

mod1 <- lmer(log10gs ~ HSPs + temp_treatment + (1|time) + (1|species), data = phys_prot, REML=FALSE)
summary(mod1)

mod2 <- lmer(log10gs ~ HSPs + temp_treatment + (1|species), data = phys_prot, REML=FALSE)
summary(mod2)

mod3 <- lmer(log10gs ~ HSPs + temp_treatment + (1|time), data = phys_prot, REML=FALSE)
summary(mod3)

mod4 <- lmer(log10gs ~ HSPs + temp_treatment + (temp_treatment|species) + (1|time), data = phys_prot, REML=FALSE)
summary(mod4)

mod5 <- lmer(log10gs ~ HSPs + temp_treatment + (1|species) + (temp_treatment|time), data = phys_prot, REML=FALSE)

mod6 <-  lmer(log10gs ~ HSPs + temp_treatment + (temp_treatment|species), data = phys_prot, REML=FALSE)
  
  
AIC(mod,mod1,mod2,mod3,mod4,mod5,mod6)

anova(mod,mod1) # protein/temp treatment interaction is not significant
anova(mod1,mod2) # time has significant effect
anova(mod1,mod3) # species has significant effect
anova(mod1,mod4) # sig
anova(mod3,mod4) # highly sig
anova(mod1,mod5) # sig
anova(mod4,mod5) # no diff
anova(mod4,mod6) # time itself not important
