## Libraries
library(haven)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GGally)
library(corrplot)
library(nlme)
library(lmtest)
library(DescTools)
library(effects)
library(ggeffects)
library(MASS)

## Import and fix the data
alz <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")

head(alz)
summary(alz)

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)

## Create baseline values
alz$ab_base <- alz$abpet0
alz$tau_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0

summary(alz)

## Create longitudinal dataset
alz_df <- data.frame(alz)

alz_long <- alz_df %>%
  pivot_longer(

    cols = matches("^(bprs|cdrsb|abpet|taupet)\\d+$"),
    
    names_to = c(".value", "year"),
 
    names_pattern = "(bprs|cdrsb|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                            
    sample = factor(rep(1:nrow(alz_df), each = 7)) # ID per ogni paziente
  )

## Discretize variables

## Maybe any 5 years?
alz_long$age_disc <- (alz_long$age %/% 5) * 5
alz_long$age_disc <- as.factor(alz_long$age_disc)

## bmi any 4
alz_long$bmi_disc <- (alz_long$bmi %/% 4) * 4
alz_long$bmi_disc <- as.factor(alz_long$bmi_disc)

## inkomen any 500
alz_long$inkomen_disc <- (alz_long$inkomen %/% 500) * 500
alz_long$inkomen_disc <- as.factor(alz_long$inkomen_disc)

## adl any 5
alz_long$adl_disc <- (as.numeric(alz_long$adl) %/% 5) * 5
alz_long$adl_disc <- as.factor(alz_long$adl_disc)

## cdrsb any 5
alz_long$cdrsb_disc <- (alz_long$cdrsb %/% 5) * 5
alz_long$cdrsb_disc <- as.factor(alz_long$cdrsb_disc)

## abpet any 0.2
alz_long$abpet_disc <- (alz_long$abpet %/% 0.2) * 0.2
alz_long$abpet_disc <- as.factor(alz_long$abpet_disc)

## taupet any 0.2
alz_long$taupet_disc <- (alz_long$taupet %/% 0.2) * 0.2
alz_long$taupet_disc <- as.factor(alz_long$taupet_disc)

## adl any 5
alz$adl_disc <- (as.numeric(alz$adl) %/% 5) * 5
alz$adl_disc <- as.factor(alz$adl_disc)

## year discrete
alz_long$year_seq <- ave(alz_long$year, alz_long$sample, FUN = function(x) as.integer(factor(x)))



#### MULTIVARIATE MODEL ####

## Start by fitting one of the most general model we can think of

## Maybe try to use heterogeneous AR(1) as a starting point

mult_model_1 <- gls(
  bprs ~ (trial + age + edu + bmi + inkomen + adl_num + wzc + cdrsb_base + ab_base + tau_base +
            sex + job) * year ,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

## This is our benchmark model

## Now we would like to reduce the mean structure


### Computationally long

mult_model_naive <- gls(
  bprs ~ 1,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

n <- nrow(alz)

step_forw <- stepAIC(mult_model_naive, 
                     direction = "forward",
                     scope = list(lower = ~1, 
                                  upper = formula(mult_model_1)),
                     k = log(n),  ## BIC Criterion
                     trace = TRUE)

step_back <- stepAIC(mult_model_1, 
                     direction = "backward",
                     scope = list(lower = ~1, 
                                  upper = formula(mult_model_1)),
                     k = log(n),
                     trace = TRUE)

step_both <- stepAIC(mult_model_1, 
                     direction = "both",
                     scope = list(lower = ~1, 
                                  upper = formula(mult_model_1)),
                     k = log(n),
                     trace = TRUE)

## Retrieve explicitly the models

mult_model_forw <- gls(
  bprs ~ year + age + trial + adl_num + wzc + cdrsb_base + 
    year:cdrsb_base,
  correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

mult_model_back <- gls(
  model = bprs ~ trial + age + bmi + adl_num + wzc + cdrsb_base + 
    job + year + age:year + adl_num:year + cdrsb_base:year + 
    job:year, 
  data = alz_long,
  correlation = corAR1(form = ~year | sample),
  weights = varIdent(form = ~1 | year), 
  method = "ML", 
  na.action = na.exclude
)

mult_model_both <- gls(
  model = bprs ~ trial + age + bmi + adl_num + wzc + cdrsb_base + 
    job + year + age:year + adl_num:year + cdrsb_base:year + 
    job:year, 
  data = alz_long,
  correlation = corAR1(form = ~year | sample),
  weights = varIdent(form = ~1 | year), 
  method = "ML", 
  na.action = na.exclude
)

## Both and Backward are the same

## Now we can compare the two models: they are nested so we can use lr test

anova(mult_model_back, mult_model_forw)

## Take the forward one

mult_model_final <- mult_model_back


### DIFFERENT COVARIANCE ###

## Simple diagonal cov matrix

mult_model_final_naive <- gls(
  bprs ~ trial + age + bmi + adl_num + wzc + cdrsb_base + 
    job + year + age:year + adl_num:year + cdrsb_base:year + 
    job:year,
  #correlation = corAR1(form = ~ year | sample),
  #weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_final, mult_model_final_naive)

## Definitely no

## Different elements on the diagonal

mult_model_final_naive <- gls(
  bprs ~ trial + age + bmi + adl_num + wzc + cdrsb_base + 
    job + year + age:year + adl_num:year + cdrsb_base:year + 
    job:year,
  #correlation = corAR1(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_final, mult_model_final_naive)

## No

## Exponential?

mult_model_final_un <- gls(
  bprs ~ trial + age + bmi + adl_num + wzc + cdrsb_base + 
    job + year + age:year + adl_num:year + cdrsb_base:year + 
    job:year,
  correlation = corExp(form = ~ year | sample),
  weights = varIdent(form = ~ 1 | year),
  method = "ML",
  data = alz_long,
  na.action = na.exclude
)

anova(mult_model_final, mult_model_final_un)

## Basically the same

