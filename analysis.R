library(dplyr)
library(lme4)
library(glmmTMB)
library(ggplot2)
library(ggeffects)
library(performance)

options(scipen=999)


#######################################################################
# Create datasets
#######################################################################

# all conventional flocks
df <- read.csv("data_conv.csv") %>%
    filter(
        FirstThinningAge > 20 | is.na(FirstThinningAge)
    )
nrow(df)

df_thinned <- df %>%
    filter(
        TimesThinned > 0
    )

df_goodFPL <- df %>%
    filter(
        PrevYearlyFPLGroup == "Good"
    )

df_mort <- df %>%
    mutate(
        TotalDeadAtFirstThinning = round(
            MortalityAtFirstThinning/100 * StartFlockSize
        )
    ) %>%
    filter(
        MortalityAtFirstThinning < Mortality,  #also removes NA
        MortalityAtFirstThinning > 0,
    ) %>%
    group_by(Farm)
nrow(df_mort)

df_mort_goodFPL <- df_mort %>%
    filter(PrevYearlyFPLGroup == "Good")

#number of rows in FPL models
nrow(df %>% filter(!is.na(FootpadLesions)))


#######################################################################
# Create models
#######################################################################

# ANTIBIOTICS WEEK 1

abwk1 <- glmer(AntibioticsWeek1 ~ 
                PlacementStockDens +
                factor(HatchYear) +
                factor(HatchQuarter) +
                (1|Farm) + (1|VetId),
            data = df,
            family = binomial,
            control= lme4::glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=100000)))
summary(abwk1)
save(abwk1, file = "FinalModel_abwk1.rda")

check_collinearity(abwk1)

# Without previous FPL scores
abwk1_goodFPL <- glmer(AntibioticsWeek1 ~ 
                   PlacementStockDens +
                   factor(HatchYear) +
                   factor(HatchQuarter) +
                   (1|Farm) + (1|VetId),
               data = df_goodFPL,
               family = binomial,
               control= lme4::glmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=100000)))
summary(abwk1_goodFPL)
save(abwk1_goodFPL, file = "AltModel_abwk1.rda")


###############################################################

# ANTIBIOTICS AFTER WEEK 1

ab <- glmer(AntibioticsAfterWeek1 ~ 
                PlacementStockDens + 
                factor(TimesThinned) +
                scale(FractionThinned) + 
                scale(AgeAtSlaughter) + 
                factor(HatchYear) +
                factor(HatchQuarter) +
                (1|Farm) + (1|VetId),
            data = df,
            family = binomial,
            control= lme4::glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=100000)))
summary(ab)
save(ab, file = "FinalModel_ab_all.rda")
check_collinearity(ab)

ab_goodFPL <- glmer(AntibioticsAfterWeek1 ~ 
                        PlacementStockDens + 
                        factor(TimesThinned) +
                        scale(FractionThinned) + 
                        scale(AgeAtSlaughter) + 
                        factor(HatchYear) +
                        factor(HatchQuarter) +
                        (1|Farm) + (1|VetId),
                    data = df_goodFPL,
                    family = binomial,
                    control= lme4::glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun=100000)))
summary(ab_goodFPL)
save(ab_goodFPL, file = "AltModel_ab.rda")


#########################################################

# MORTALITY (UNTIL THINNING)

mort <- glmmTMB::glmmTMB(TotalDeadAtFirstThinning ~ 
                                 PlacementStockDens +
                                 FirstThinningAge +
                                 factor(HatchYear) +
                                 factor(HatchQuarter) +
                                 offset(log(StartFlockSize)) +
                                 offset(log(FirstThinningAge)) +
                                 (1|Farm) + (1|VetId),
                             data = df_mort,
                             family = glmmTMB::nbinom2,
                             control = glmmTMB::glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
)
summary(mort)
check_collinearity(mort)
save(mort, file = "FinalModel_mort.rda")

# Only good prev_year
mort_alt <- glmmTMB::glmmTMB(TotalDeadAtFirstThinning ~ 
                                 PlacementStockDens +
                                 FirstThinningAge +
                                 factor(HatchYear) +
                                 factor(HatchQuarter) +
                                 offset(log(StartFlockSize)) +
                                 offset(log(FirstThinningAge)) +
                                 (1|Farm) + (1|VetId),
                         data = df_mort_goodFPL,
                         family = glmmTMB::nbinom2,
                         control = glmmTMB::glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
)
summary(mort_alt)
save(mort_alt, file = "AltModel_mort.rda")


######################################################################    

# FOOTPAD LESIONS

fpl = glmmTMB::glmmTMB(FootpadLesions ~ 
                            PlacementStockDens + 
                            factor(TimesThinned) +
                            scale(FractionThinned) +
                            scale(AgeAtSlaughter) +
                            factor(HatchYear) +
                            factor(HatchQuarter) +
                            (1|Farm) + (1|VetId),
                            data = df,
                            family = glmmTMB::nbinom2,
                            control = glmmTMB::glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
)
summary(fpl)
save(fpl, file = 'FinalModelFPL.rda')
check_collinearity(fpl)


# only flocks with good FPL scores

fpl_goodFPL = glmmTMB::glmmTMB(FootpadLesions ~ 
                           PlacementStockDens +
                           factor(TimesThinned) +
                           scale(FractionThinned) + 
                           scale(AgeAtSlaughter) + 
                           factor(HatchYear) +
                           factor(HatchQuarter) +
                           (1|Farm) + (1|VetId),
                       data = df_goodFPL,
                       family = glmmTMB::nbinom2,
                       control = glmmTMB::glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
)
summary(fpl_goodFPL)
save(fpl_goodFPL, file = 'AltModelFPL.rda')



#####################################################
# Results
#####################################################

###### AB week 1

load(file = "FinalModel_abwk1.rda") #is named 'abwk1'

# Get odds ratios and Wald 95% CI
coefs <- summary(abwk1)$coefficients
OR <- exp(coefs[, "Estimate"])
lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
p <- coefs[, "Pr(>|z|)"]

results_abwk1 <- data.frame(
    Estimate = coefs[, "Estimate"],
    OR = OR,
    CI_low = lower,
    CI_high = upper,
    p_value = p
) %>%
    mutate(
        across(where(is.numeric), ~ round(.x, 2))
    )
results_abwk1



load(file = "AltModel_abwk1.rda") #is named 'abwk1_goodFPL'

# Get odds ratios and Wald 95% CI
coefs <- summary(abwk1_goodFPL)$coefficients
OR <- exp(coefs[, "Estimate"])
lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
p <- coefs[, "Pr(>|z|)"]

results_abwk1_goodFPL <- data.frame(
    Estimate = coefs[, "Estimate"],
    OR = OR,
    CI_low = lower,
    CI_high = upper,
    p_value = p
) %>%
    mutate(
        across(where(is.numeric), ~ round(.x, 2))
    )
results_abwk1_goodFPL



##### AB After Week 1

load(file = "FinalModel_ab_all.rda") #is named 'ab'

pred_ab <- ggpredict(
    ab,
    terms = "PlacementStockDens [10:25, by = 1]",
    bias_correction = FALSE
) %>%
    as.data.frame() %>%
    mutate(Model = "Antibiotics after week 1")

save(pred_ab, file = "pred_ab.csv")

ggplot(pred_ab,
       aes(x = x, y = predicted)) +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.15,
                linewidth = 0,
                color = "#de8500", fill = "#de8500"
                ) +
    geom_line(
        linewidth = 2,
        color = "#de8500"
        ) +
    labs(
        x = "Placement stocking density (birds/m2)",
        y = "Predicted treatment probability",
        title = "A"
    ) +
    theme_classic() +
    theme(
        text=element_text(size = 18),
        axis.title.x = element_text(size = 16, margin = margin(t = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    ) +
    scale_y_continuous(
        labels = scales::percent,
        limits = c(0, NA))

# Get odds ratios and Wald 95% CI
coefs <- summary(ab)$coefficients
OR <- exp(coefs[, "Estimate"])
lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
p <- coefs[, "Pr(>|z|)"]

results_ab <- data.frame(
    Estimate = coefs[, "Estimate"],
    OR = OR,
    CI_low = lower,
    CI_high = upper,
    p_value = p
) %>%
    mutate(
        across(where(is.numeric), ~ round(.x, 2))
    )
results_ab

# random effects
r2_nakagawa(ab)
summary(ab)


load(file = "AltModel_ab.rda") #is named 'ab_goodFPL'

# Get odds ratios and Wald 95% CI
coefs <- summary(ab_goodFPL)$coefficients
OR <- exp(coefs[, "Estimate"])
lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
p <- coefs[, "Pr(>|z|)"]

results_ab_goodFPL <- data.frame(
    Estimate = coefs[, "Estimate"],
    OR = OR,
    CI_low = lower,
    CI_high = upper,
    p_value = p
) %>%
    mutate(
        across(where(is.numeric), ~ round(.x, 2))
    )
results_ab_goodFPL

#scale
sd(df$FractionThinned)
sd(df$AgeAtSlaughter)
sd(df_goodFPL$FractionThinned)
sd(df_goodFPL$AgeAtSlaughter)


##### Mortality until thinning

load(file = "FinalModel_mort.rda") #is named 'mort'

# Get RATE ratios and Wald 95% CI
coefs <- summary(mort)$coefficients$cond
RR <- exp(coefs[, "Estimate"])
lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
p <- coefs[, "Pr(>|z|)"]
results_mort <- data.frame(
    Estimate = coefs[, "Estimate"],
    RR = RR,
    CI_low = lower,
    CI_high = upper,
    p_value = p
) %>%
    mutate(
        across(where(is.numeric), ~ round(.x, 2))
    )

results_mort
#confint(mort, method="profile")  #doesn't work!!
confint(mort, method="boot")


##### Footpad lesions

min(df$PlacementStockDens)
max(df$PlacementStockDens)

load(file = "FinalModelFPL.rda") #is named 'fpl'

summary(fpl)
sd(df$PlacementStockDens)


# Get RATE ratios and Wald 95% CI
coefs <- summary(fpl)$coefficients$cond
IRR <- exp(coefs[, "Estimate"])
lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
p <- coefs[, "Pr(>|z|)"]
results_fpl <- data.frame(
    Estimate = coefs[, "Estimate"],
    IRR = IRR,
    CI_low = lower,
    CI_high = upper,
    p_value = p
) %>%
    mutate(
        across(where(is.numeric), ~ round(.x, 2))
    )

results_fpl

# effect of farm vs fixed effects:
icc(fpl)
r2_nakagawa(fpl)
# ! fpl model is not meant to explain all variance. so possible that much of 
# unexplained variance goes into 'farm variation'



load(file = "AltModelFPL.rda") #is named 'fpl_goodFPL'

# Get RATE ratios and Wald 95% CI
coefs <- summary(fpl_goodFPL)$coefficients$cond
IRR <- exp(coefs[, "Estimate"])
lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
p <- coefs[, "Pr(>|z|)"]
results_fpl_goodFPL <- data.frame(
    Estimate = coefs[, "Estimate"],
    IRR = IRR,
    CI_low = lower,
    CI_high = upper,
    p_value = p
) %>%
    mutate(
        across(where(is.numeric), ~ round(.x, 2))
    )

results_fpl_goodFPL



pred_fpl <- ggpredict(
    fpl,
    terms = "PlacementStockDens [10:25, by = 1]",
    bias_correction = FALSE
) %>%
    as.data.frame()

ggplot(pred_fpl,
       aes(x = x, y = predicted)) +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.15,
                linewidth = 0,
                color = "#de8500", fill = "#de8500"
    ) +
    geom_line(
        linewidth = 2,
        color = "#de8500"
    ) +
    labs(
        x = "Placement stocking density (birds/m2)",
        y = "Predicted FPL score",
        title = "B"
    ) +
    theme_classic() +
    theme(
        text=element_text(size = 18),
        axis.title.x = element_text(size = 16, margin = margin(t = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    ) +
    scale_y_continuous(
        limits = c(0, NA)
    )



##################################
# Descriptives
#################################

nrow(df)

# n houses
nrow(unique(df[, c("Farm", "Pen")]))

# n farms
nrow(df[unique(df$Farm), ])

# n flocks with previous score
nrow(df_goodFPL)
nrow(df %>% filter(!is.na(FootpadLesions)))
nrow(df) - nrow(df_goodFPL)
(nrow(df) - nrow(df_goodFPL)) / nrow(df)

# thinning
mean(df_thinned$FractionThinned)
sd(df_thinned$FractionThinned)

table(df$TimesThinned)
prop.table(table(df$TimesThinned))

# age
mean(df$AgeAtSlaughter)
summary(df$AgeAtSlaughter)

df_all <- df %>%
    summarise(
        median_place_sd = median(PlacementStockDens),
        place_min = min(PlacementStockDens),
        place_max= max(PlacementStockDens),
        place_sd_25 = quantile(PlacementStockDens, 0.25),
        place_sd_75 = quantile(PlacementStockDens, 0.75),
        
        median_final_sd = median(FinalStockDens),
        max_final_sd = max(FinalStockDens),
        final_sd_25 = quantile(FinalStockDens, 0.25),
        final_sd_75 = quantile(FinalStockDens, 0.75),
        
        median_age = median(AgeAtSlaughter),
        age_25 = quantile(AgeAtSlaughter, 0.25),
        age_75 = quantile(AgeAtSlaughter, 0.75),
        
        mean_kg = mean(SlaughterWeightPerChicken),
        sd_kg = sd(SlaughterWeightPerChicken),
        fraction_mean = mean(FractionThinned, na.rm = TRUE),
        fraction_sd = sd(FractionThinned, na.rm = TRUE)
    ) %>%
    t()
df_all

df_per_t <- df %>%
    group_by(TimesThinned) %>%
    summarise(
        median_place_sd = median(PlacementStockDens),
        place_sd_25 = quantile(PlacementStockDens, 0.25),
        place_sd_75 = quantile(PlacementStockDens, 0.75),
        
        median_final_sd = median(FinalStockDens),
        final_sd_25 = quantile(FinalStockDens, 0.25),
        final_sd_75 = quantile(FinalStockDens, 0.75),
        
        median_age = median(AgeAtSlaughter),
        age_25 = quantile(AgeAtSlaughter, 0.25),
        age_75 = quantile(AgeAtSlaughter, 0.75),
        
        first_thinning = median(FirstThinningAge, na.rm=TRUE),
        f_age_25 = quantile(FirstThinningAge, 0.25, na.rm=TRUE),
        f_age_75 = quantile(FirstThinningAge, 0.75, na.rm=TRUE),
        
        second_thinning = median(SecondThinningAge, na.rm=TRUE),
        s_age_25 = quantile(SecondThinningAge, 0.25, na.rm=TRUE),
        s_age_75 = quantile(SecondThinningAge, 0.75, na.rm=TRUE),
        
        mean_kg = mean(SlaughterWeightPerChicken),
        sd_kg = sd(SlaughterWeightPerChicken),
        fraction_mean = mean(FractionThinned, na.rm = TRUE),
        fraction_sd = sd(FractionThinned, na.rm = TRUE)
    ) %>%
    t()
df_per_t

# weight
mean(df$SlaughterWeightPerChicken)
summary(df$SlaughterWeightPerChicken)


# histogram stocking density per times thinned

df %>%
    ggplot( aes(x=PlacementStockDens, fill=factor(TimesThinned))) +
    geom_histogram(
        color="#e9ecef",
        alpha=0.5,
        position = 'identity',
        binwidth = 0.5) +
    scale_fill_manual(
        values= c('#E50545', '#6297BE', '#022D4A')) +
    theme_classic() +
    labs(y="Count",
         x = "Placement stocking density (birds/m2)",
         fill = "Times thinned") +
    theme(
        text=element_text(size = 18),
        axis.title.x = element_text(size = 16, margin = margin(t = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
        legend.title = element_text(size = 14)
    )

# one without color:
df %>%
    ggplot( aes(x=PlacementStockDens)) +
    geom_histogram(
        color="#e9ecef",
        fill = "#4781AE",
        alpha=1,
        position = 'identity',
        binwidth = 0.5) +
    theme_classic() +
    labs(y="Count",
         x = "Placement stocking density (birds/m2)",
         fill = "Times thinned") +
    theme(
        text=element_text(size = 18),
        axis.title.x = element_text(size = 16, margin = margin(t = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
        legend.title = element_text(size = 14)
    )

df %>%
    ggplot(
        aes(
            x = FinalStockDens, 
            fill = forcats::fct_rev(factor(TimesThinned))
        )
    ) +
    geom_histogram(
        color="#e9ecef",
        alpha=0.5,
        position = 'identity',
        binwidth = 0.5) +
    scale_fill_manual(
        values= c('#022D4A', '#6297BE', '#E50545')) +
    theme_classic() +
    labs(y="Count",
         x = "End stocking density (kg/m2)",
         fill = "Times thinned") +
    theme(
        text=element_text(size = 18),
        axis.title.x = element_text(size = 16, margin = margin(t = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
        legend.title = element_text(size = 14)
    )

df %>%
    ggplot( aes(x=FinalStockDens)) +
    geom_histogram(
        color="#e9ecef",
        fill = "#4781AE",
        alpha=1,
        position = 'identity',
        binwidth = 0.5) +
    theme_classic() +
    labs(y="Count",
         x = "End stocking density (kg/m2)",
         fill = "Times thinned") +
    theme(
        text=element_text(size = 18),
        axis.title.x = element_text(size = 16, margin = margin(t = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
        legend.title = element_text(size = 14)
    )