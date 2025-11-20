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
library(Matrix)
library(MASS)
library(metafor)

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
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))

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


#### 2 STAGE MODEL ####

### STAGE 1 MODEL ###

## Start by fitting any time a linear regression

count <- 0
coeff_stage1 <- matrix(nrow = length(alz$patid), ncol = 2)
sigma_stage_1 <- matrix(nrow = length(alz$patid), ncol = 1)
r2_stage1 <- matrix(nrow = length(alz$patid), ncol = 1)
r2_stage1_quad <- matrix(nrow = length(alz$patid), ncol = 1)
sum_squares <- matrix(nrow = length(alz$patid), ncol = 3)

# We assume basically that bprs_i = beta_0i + beta_1i * year_i + eps_i
# We should define a structure for the errors
# Usually it is reasonable to consider esp_i ~ N(0, Sigma)
# and Sigma = sigma^2 * I

for (i in 1:(length(alz$patid))) {
  idx <- (7 * i + 1):(7 * i + 7)
  bprs_values <- alz_long$bprs[idx]
  mean_bprs <- mean(bprs_values, na.rm = TRUE)
  sum_squares[i, 1] <- sum((bprs_values - mean_bprs)^2, na.rm = TRUE)
  mod_prova <- lm(bprs ~ year, 
                  data = alz_long[c((7*count + 1) : (7*count + 7)), ])
  mod_prova_quad <- lm(bprs ~ year + I(year^2), 
                       data = alz_long[c((7*count + 1) : (7*count + 7)), ])
  coeff_stage1[i, ] <- mod_prova$coefficients
  r2_stage1[i] <- summary(mod_prova)$r.squared
  r2_stage1_quad[i] <- summary(mod_prova_quad)$r.squared
  sigma_prov <- sqrt(sum(residuals(mod_prova)^2) / df.residual(mod_prova))
  sigma_stage_1[i] <- sigma_prov
  sum_squares[i, 2] <- sum(residuals(mod_prova)^2)
  sum_squares[i, 3] <- sum(residuals(mod_prova_quad)^2)
  count = count + 1
}

## Quick visualization of the data

r_squared_meta <- 1 - (sum(sum_squares[, 2], na.rm = TRUE) / sum(sum_squares[, 1], na.rm = TRUE))
r_squared_meta_quad <- 1 - (sum(sum_squares[, 3], na.rm = TRUE) / sum(sum_squares[, 1], na.rm = TRUE))

## Value really high, not that bad

# Visualization

r_squared_visual <- data.frame(
  n_obs = alz$n_obs_data,
  r2_stage1 = r2_stage1
)

r_squared_visual_quad <- data.frame(
  n_obs = alz$n_obs_data,
  r2_stage1 = r2_stage1_quad
)

ggplot(r_squared_visual, aes(x = n_obs, y = r2_stage1)) +
  geom_point(alpha = 0.6, color = "blue", size = 2) +
  geom_hline(yintercept = r_squared_meta, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "Scatterplot of R² under Linear Model",
    x = "Number n_i of measurements",
    y = "Coefficient Ri²"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggplot(r_squared_visual_quad, aes(x = n_obs, y = r2_stage1)) +
  geom_point(alpha = 0.6, color = "blue", size = 2) +
  geom_hline(yintercept = r_squared_meta_quad, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "Scatterplot of R² under Quadratic Model",
    x = "Number n_i of measurements",
    y = "Coefficient Ri²"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

### Model Comparison - F TEST 

## Prova veloce

SSE_L <- sum(sum_squares[, 2], na.rm = TRUE)
SSE_Q <- sum(sum_squares[, 3], na.rm = TRUE)

df_L <- sum(alz$n_obs_data - 2)               
df_Q <- sum(alz$n_obs_data - 3)               

# F-test 
F_meta <- ((SSE_L - SSE_Q) / (df_L - df_Q)) / (SSE_Q / df_Q)
p_meta <- 1 - pf(F_meta, df_L - df_Q, df_Q)

cat("F_meta =", F_meta, "  p-value =", p_meta, "\n")

# Here the stage 1 model is over


### STAGE 2 MODEL ###

res.list <- lmList(bprs ~ year | sample, data = alz_long, na.action = na.exclude)
b <- lapply(res.list, coef)
V <- lapply(res.list, vcov)

estm <- rep(c("intercept","slope"), length(b))
subj <- rep(names(b), each=2)


b <- unlist(b)
V <- bldiag(V)

covariate_data <- alz_long[!duplicated(alz_long$sample), 
                           c("sample", "trial", "sex", "age", "edu", 
                             "bmi", "inkomen", "job", "adl_num", "wzc", 
                             "abpet", "taupet", "cdrsb")]

subj_names_ordered <- names(res.list)
covariate_data_aligned <- covariate_data[match(subj_names_ordered, covariate_data$sample), ]

trial_rep   <- rep(covariate_data_aligned$trial, each = 2)
sex_rep     <- rep(covariate_data_aligned$sex, each = 2)
age_rep     <- rep(covariate_data_aligned$age, each = 2)
edu_rep     <- rep(covariate_data_aligned$edu, each = 2)
bmi_rep     <- rep(covariate_data_aligned$bmi, each = 2)
inkomen_rep <- rep(covariate_data_aligned$inkomen, each = 2)
job_rep     <- rep(covariate_data_aligned$job, each = 2)
adl_num_rep     <- rep(covariate_data_aligned$adl_num, each = 2)
wzc_rep     <- rep(covariate_data_aligned$wzc, each = 2)
abpet_rep   <- rep(covariate_data_aligned$abpet, each = 2)
taupet_rep  <- rep(covariate_data_aligned$taupet, each = 2)
cdrsb_rep   <- rep(covariate_data_aligned$cdrsb,  each = 2)
length(coef(res.list[[1]]))

subj <- rep(subj_names_ordered, each = 2)

length(b)           # 2 * n
length(trial_rep)   
length(subj) 


res2_full <- rma.mv(b ~ estm + 
                      estm:trial_rep + estm:edu_rep + estm:inkomen_rep +
                      estm:sex_rep + estm:age_rep + 
                      estm:bmi_rep + estm:job_rep + estm:adl_num_rep + 
                      estm:wzc_rep + estm:abpet_rep + estm:taupet_rep + estm:cdrsb_rep - 1, 
                    V = V, 
                    random = ~ estm | subj, 
                    struct = "UN",
                    method = "ML")

res2_full

## Get rid of the non-significant variables and re-estimate

# Indicatori
is_intercept <- as.numeric(estm == "intercept")
is_slope     <- as.numeric(estm == "slope")

# TRIAL: versione solo per l'intercept
trial_intercept <- trial_rep


if (is.factor(trial_intercept)) {
  levels(trial_intercept) <- c(levels(trial_intercept), "none")
  trial_intercept[is_slope == 1] <- "none"  # ora ok
  trial_intercept <- droplevels(trial_intercept)  # pulisce eventuali livelli non usati
} else {
  trial_intercept[is_slope == 1] <- 0
}

bmi_intercept <- bmi_rep * is_intercept

# Final Model
res2_prova <- rma.mv(
  b ~ estm + trial_intercept + estm:age_rep + estm:bmi_intercept + 
    estm:job_rep + estm:adl_num_rep + estm:wzc_rep+ estm:cdrsb_rep - 1,
  V = V,
  random = ~ estm | subj,
  struct = "UN",
  method = "ML",
  control = list(stepadj = 0.5, optimizer = "optim")
)

summary(res2_prova)