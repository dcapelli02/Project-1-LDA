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
    # 1️⃣ Seleziona tutte le colonne che iniziano con cdrsb, abpet o taupet,
    #     seguite da un numero (0–6)
    cols = matches("^(bprs|cdrsb|abpet|taupet)\\d+$"),
    
    # 2️⃣ Crea due nuove colonne:
    #     - ".value" → il prefisso (cdrsb, abpet, taupet)
    #     - "year"   → il numero alla fine
    names_to = c(".value", "year"),
    
    # 3️⃣ Divide i nomi in due parti: (prefisso, numero)
    names_pattern = "(bprs|cdrsb|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                            # converte in numerico
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

### COMPARAZIONE MODELLI - F TEST 

## Prova veloce

# Totali aggregati
SSE_L <- sum(sum_squares[, 2], na.rm = TRUE)
SSE_Q <- sum(sum_squares[, 3], na.rm = TRUE)

df_L <- sum(alz$n_obs_data - 2)               # due parametri: beta0, beta1
df_Q <- sum(alz$n_obs_data - 3)               # tre parametri: beta0, beta1, beta2

# F-test aggregato
F_meta <- ((SSE_L - SSE_Q) / (df_L - df_Q)) / (SSE_Q / df_Q)
p_meta <- 1 - pf(F_meta, df_L - df_Q, df_Q)

cat("F_meta =", F_meta, "  p-value =", p_meta, "\n")

# Here the stage 1 model is over


### STAGE 2 MODEL ###

## We can group beta0 and beta1

## First of all store the results
beta0 <- coeff_stage1[, 1]
beta1 <- coeff_stage1[, 2]


## Tutto numerico
alz_prova <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")

## Group them into a single vector
y <- as.vector((cbind(beta0, beta1)))

Z <- model.matrix(~ trial + sex + job + age + inkomen + 
                    adl + wzc + abpet0 + taupet0 +
                    edu + bmi + cdrsb0, data = alz_prova)

block_Z <- kronecker(diag(2), Z)

colnames(block_Z) <- c(
  paste0(colnames(Z), "_int"),   # per il primo blocco (intercept)
  paste0(colnames(Z), "_time")   # per il secondo blocco (slope)
)

alz_long2 <- data.frame(
  sample = rep(alz$patid, 2),
  coef_type = factor(c(rep("beta0", times = nrow(alz)), rep("beta1", times = nrow(alz)))),
  y = y,
  block_Z = block_Z  # o puoi ricostruire con model.matrix
)

names(alz_long2) <- sub("^block_Z\\.", "", names(alz_long2))

alz_long2$trial_int[alz_long2$coef_type == "beta1"] <- alz_long2$trial_time[alz_long2$coef_type == "beta1"] + 25
alz_long2$sex_int[alz_long2$coef_type == "beta1"] <- alz_long2$sex_time[alz_long2$coef_type == "beta1"] + 2
alz_long2$edu_int[alz_long2$coef_type == "beta1"] <- alz_long2$edu_time[alz_long2$coef_type == "beta1"] + 4
alz_long2$job_int[alz_long2$coef_type == "beta1"] <- alz_long2$job_time[alz_long2$coef_type == "beta1"] + 2
alz_long2$wzc_int[alz_long2$coef_type == "beta1"] <- alz_long2$wzc_time[alz_long2$coef_type == "beta1"] + 2

alz_long2$trial_int <- as.factor(alz_long2$trial_int)
alz_long2$sex_int <- as.factor(alz_long2$sex_int)
alz_long2$edu_int <- as.factor(alz_long2$edu_int)
alz_long2$job_int <- as.factor(alz_long2$job_int)
alz_long2$wzc_int <- as.factor(alz_long2$wzc_int)

alz_long2$trial_time <- as.factor(alz_long2$trial_time)
alz_long2$sex_time <- as.factor(alz_long2$sex_time)
alz_long2$edu_time <- as.factor(alz_long2$edu_time)
alz_long2$job_time <- as.factor(alz_long2$job_time)
alz_long2$wzc_time <- as.factor(alz_long2$wzc_time)

alz_long2 <- alz_long2[order(alz_long2$sample, alz_long2$coef_type), ]

contrasts(alz_long2$job_time) <- matrix(c(0,1), ncol=1)
colnames(contrasts(alz_long2$job_time)) <- "job_time1"

## Questo è da ridurre perchè non fitta tutto a quanto pare...

contrasts(alz_long2$trial_int) <- contr.sum(26)
contrasts(alz_long2$trial_time) <- contr.sum(26)

View(model.matrix(y ~ sex_time + job_time + edu_time - 1,
             correlation = corCompSymm(form = ~ 1 | sample),
             weights = varIdent(form = ~ 1 | coef_type),
             method = "ML",
             na.action = na.exclude,
             data = alz_long2
  
))

rank(model.matrix(y ~ sex_time + job_time + edu_time - 1,
                correlation = corCompSymm(form = ~ 1 | sample),
                weights = varIdent(form = ~ 1 | coef_type),
                method = "ML",
                na.action = na.exclude,
                data = alz_long2
                
))

mod_2_full <- gls(y ~ job_time + edu_time + sex_time - 1,
                  correlation = corCompSymm(form = ~ 1 | sample),
                  weights = varIdent(form = ~ 1 | coef_type),
                  method = "ML",
                  na.action = na.exclude,
                  data = alz_long2)

mod_2_naive <- gls(y ~ .Intercept._time + .Intercept._int - 1,
                                 correlation = corSymm(form = ~ 1 | sample),
                                 weights = varIdent(form = ~ 1 | coef_type),
                                 method = "ML",
                                 na.action = na.exclude,
                                 data = alz_long2)

## Stepwise procedure

step_stage2 <- stepAIC(mod_2_naive, 
                       direction = "forward",
                       scope = list(lower = ~ .Intercept._time + .Intercept._int - 1, 
                                    upper = ~ . - 1 - sample - coef_type - trial_time - edu_int - edu_time
                                    - trial_int - age_int),
                                    #upper = formula(mod_2_full)),
                       trace = TRUE)





###################### FINO QUA ###################################















#nuovo metodo
library(metafor)

res.list <- lmList(bprs ~ year | sample, data = alz_long, na.action = na.exclude)
b <- lapply(res.list, coef)
V <- lapply(res.list, vcov)

estm <- rep(c("intercept","slope"), length(b))
subj <- rep(names(b), each=2)


b <- unlist(b)
V <- bldiag(V)

# Prendi solo la prima osservazione di ogni soggetto
covariate_data <- alz_long[!duplicated(alz_long$sample), 
                           c("sample", "trial", "sex", "age", "edu", 
                             "bmi", "inkomen", "job", "adl_num", "wzc", 
                             "abpet", "taupet", "cdrsb")]

# Allinea all'ordine dei soggetti di lmList
subj_names_ordered <- names(res.list)
covariate_data_aligned <- covariate_data[match(subj_names_ordered, covariate_data$sample), ]

# Crea vettori ripetuti per ogni stima (intercept/slope)
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
length(trial_rep)   # deve essere uguale a length(b)
length(subj) 


is_intercept <- as.numeric(estm == "intercept")
## Per i fattori
trial_intercept <- trial_rep
trial_intercept[estm == "slope"] <- NA   # rimuove trial per le slope
#trial_intercept[is_intercept == 0] <- "none"  # livello fittizio
#trial_intercept <- factor(trial_intercept) # ricrea il factor correttamente

## Per le variabili numeriche
age_intercept <- age_rep * is_intercept

res2_prova <- rma.mv(b ~ estm + 
                       trial_intercept + estm:age_rep - 1, 
                     V = V, 
                     random = ~ estm | subj, 
                     struct = "UN",
                     method = "ML",
                     na.action = na.exclude)

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

## Dimenticato un po' di roba ma va beh

## Leviamo slope*age, slope*bmi, abpet direttamente, taupet direttamente,,
## anche sex*interc andrebbe tolta

# Indicatori
is_intercept <- as.numeric(estm == "intercept")
is_slope     <- as.numeric(estm == "slope")

# TRIAL: versione solo per l'intercept
trial_intercept <- trial_rep


if (is.factor(trial_intercept)) {
  # aggiungi esplicitamente il nuovo livello
  levels(trial_intercept) <- c(levels(trial_intercept), "none")
  trial_intercept[is_slope == 1] <- "none"  # ora ok
  trial_intercept <- droplevels(trial_intercept)  # pulisce eventuali livelli non usati
} else {
  # se è numerico
  trial_intercept[is_slope == 1] <- 0
}


# BMI: variabile numerica, interazione esplicita con estm
bmi_intercept <- bmi_rep * is_intercept

# Modello finale
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









# Indicatori
is_intercept <- as.numeric(estm == "intercept")
is_slope     <- as.numeric(estm == "slope")

# TRIAL: versione solo per l'intercept
trial_intercept <- trial_rep


if (is.factor(trial_intercept)) {
  # aggiungi esplicitamente il nuovo livello
  levels(trial_intercept) <- c(levels(trial_intercept), "none")
  trial_intercept[is_slope == 1] <- "none"  # ora ok
  trial_intercept <- droplevels(trial_intercept)  # pulisce eventuali livelli non usati
} else {
  # se è numerico
  trial_intercept[is_slope == 1] <- 0
}


# AGE: variabile numerica, interazione esplicita con estm
age_intercept <- age_rep * is_intercept

# Modello finale
res2_prova <- rma.mv(
  b ~ estm + trial_intercept + estm:age_rep - 1,
  V = V,
  random = ~ estm | subj,
  struct = "UN",
  method = "ML",
  control = list(stepadj = 0.5, optimizer = "optim")
)

summary(res2_prova)











#### ALTRO TENTATIVO ####
## Now create a new dataset with the beta0 and beta1 values

alz_long_beta <- cbind(
  cbind(alz, coef_type = "int",  y = beta0),
  cbind(alz, coef_type = "time", y = beta1)
)

alz_long_beta$coef_type <- as.factor(alz_long_beta$coef_type)

## Then fit a gls to this model

mod2_full <- gls(
  y ~ coef_type * (trial + sex + job + age + inkomen + 
                     adl_num + wzc + ab_base + tau_base +
                     edu + bmi + cdrsb_base) - 1,
  correlation = corSymm(form = ~ 1 | patid),
  weights = varIdent(form = ~ 1 | coef_type),
  data = alz_long_beta,
  na.action = na.exclude,
  method = "ML"
)


## Then stepwise procedure

mod2_naive <- gls(
  y ~ 1,
  correlation = corSymm(form = ~ 1 | patid),
  weights = varIdent(form = ~ 1 | coef_type),
  data = alz_long_beta,
  na.action = na.exclude,
  method = "ML"
)


step2 <- stepAIC(mod2_naive, 
                 direction = "forward",
                 scope = list(lower = ~1, 
                              upper = formula(mod2_full)),
                 trace = TRUE)











## Then create a new block matrix
Z <- model.matrix(~ trial + sex + job + age + inkomen + 
                    adl_num + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base, data = alz)

block_Z <- kronecker(diag(2), Z)

colnames(block_Z) <- c(
  paste0(colnames(Z), "_int"),   # per il primo blocco (intercept)
  paste0(colnames(Z), "_time")   # per il secondo blocco (slope)
)

alz_long2 <- data.frame(
  id = rep(alz$patid, each = 2),
  coef_type = factor(rep(c("beta0", "beta1"), times = nrow(alz))),
  y = y,
  block_Z = block_Z  # o puoi ricostruire con model.matrix
)

## Then we can fit a gls to this model,
## we want to impose a generic covariance structure

mod_2_full <- gls(y ~ block_Z - 1,
                  correlation = corSymm(form = ~ 1 | id),
                  weights = varIdent(form = ~ 1 | coef_type),
                  method = "ML",
                  na.action = na.exclude,
                  data = alz_long2)

mod_2_naive <- mod_2_full <- gls(y ~ 1,
                                 correlation = corSymm(form = ~ 1 | id),
                                 weights = varIdent(form = ~ 1 | coef_type),
                                 method = "ML",
                                 na.action = na.exclude,
                                 data = alz_long2)

## Stepwise procedure

step_stage2 <- stepAIC(mod_2_naive, 
                      direction = "forward",
                      scope = list(lower = ~1, 
                                   upper = formula(mod_2_full)),
                      trace = TRUE)






## First of all store the results
beta0 <- coeff_stage1[, 1]
beta1 <- coeff_stage1[, 2]

model_beta0 <- lm(beta0 ~ trial + sex + job + age + inkomen + 
                    adl_num + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base,
                  #weights = weights_lm,
                  data = alz)

model_beta1 <- lm(beta1 ~ trial + sex + job + age + inkomen + 
                    adl_num + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base,
                  #weights = weights_lm,
                  data = alz)

## So these are the starting points
## We want to reduce the models

step_beta_0 <- step(model_beta0)

model_beta0_final <- step_beta_0$call

## Se servisse
# lm(formula = beta0 ~ job + age + inkomen + adl_disc + wzc + ab_base + 
# tau_base + edu + bmi + cdrsb_base, data = alz)

step_beta_1 <- step(model_beta1)

model_beta1_final <- step_beta_1$call

## Se servisse
# lm(formula = beta1 ~ adl_disc + wzc + edu + cdrsb_base, data = alz)

coeff_beta0_final <- model_beta0_final$coefficients
coeff_beta1_final <- model_beta1_final$coefficients






#### PROVA ####
library(Matrix)

# 1️⃣ Definisci le risposte (stacked vector)
y <- as.vector(t(cbind(beta0, beta1)))  # [β0_1, β1_1, β0_2, β1_2, ...]

# 2️⃣ Definisci la matrice delle covariate individuali Z (una per paziente)
Z <- model.matrix(~ trial + sex + job + age + inkomen + 
                    adl_num + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base, data = alz)

# 3️⃣ Crea la grande matrice a blocchi diagonali
#    Ogni blocco corrisponde a Z_i ripetuto per ciascun coefficiente (β0, β1)
bigZ <- kronecker(diag(2), Z)

# 4️⃣ Fitta il modello combinato
fit_big <- gls(y ~ bigZ - 1,
               na.action = na.exclude)  # "-1" per evitare un intercept duplicato
summary(fit_big)






#### POSSIBILE MODIFICA DELLO STAGE 1 ####

## We need to decide the form of the error matrix Sigma_i

## Look at the correlation plot
plot(alz[, c(18:24)])

## There is clearly correlation among the measures
## It is not reasonable to assume that the error matrix is sigma^2 * I
## as conditional independence is not reasonable
## (the measures are correlated over time)

## What is a reasonable structure?

## Measurement error + serial correlation
## Measurement: sigma^2 * I (ragionevole sia sempre lo stesso)
## Serial correlation: exponential (but we need enough data)

## If we don't have enough data we simply fit a gls model

coeff_stage1 <- matrix(nrow = length(alz$patid), ncol = 2)
r2_stage1 <- matrix(nrow = length(alz$patid), ncol = 1)
sum_squares <- matrix(nrow = length(alz$patid), ncol = 2)


for (i in 1:(length(alz$patid))) {
  idx <- (7 * (i-1) + 1):(7 * (i-1) + 7)
  dati_i <- alz_long[idx, ]
  
  dati_i <- dati_i[complete.cases(dati_i[, c("bprs", "year")]), ]
  
  # Compute only if we have data
  
  if (nrow(dati_i) > 0) {
    mean_bprs <- mean(dati_i$bprs, na.rm = TRUE)
    sum_squares[i, 1] <- sum((dati_i$bprs - mean_bprs)^2, na.rm = TRUE)
  } else {
    sum_squares[i, 1] <- NA
    next
  }
  
  # We need at least two values
  if (length(unique(dati_i$year)) < 2) {
    coeff_stage1[i, ] <- NA
    r2_stage1[i] <- NA
    sum_squares[i, 2] <- NA
    next
  }
  
  # Model
  try({
    if (alz$n_obs_data[i] > 3) {
      mod_prova <- gls(
        bprs ~ year,
        correlation = corExp(form = ~ year),
        data = dati_i,
        method = "ML",
        na.action = na.exclude
      )
    } else {
      mod_prova <- lm(bprs ~ year, data = dati_i)
    }
    
    coeff_stage1[i, ] <- mod_prova$coefficients
    sum_squares[i, 2] <- sum(residuals(mod_prova)^2, na.rm = TRUE)
    r2_stage1[i] <- 1 - (sum_squares[i, 2]/sum_squares[i, 1])
    
  }, silent = TRUE)
}

r_squared_meta <- 1 - (sum(sum_squares[, 2], na.rm = TRUE) / sum(sum_squares[, 1], na.rm = TRUE))

# Visualization

r_squared_visual <- data.frame(
  n_obs = alz$n_obs_data,
  r2_stage1 = r2_stage1
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


### STAGE 2 MODEL ###

## First of all store the results
beta0 <- coeff_stage1[, 1]
beta1 <- coeff_stage1[, 2]

## Then we can perform a regression

## Here usually we consider the matrix of errors D to be unstructured
## so we need to use gls and not lm

## In order to use it we need to construct a new dataframe
## (otherwise gls does not work properly)

alz$beta0 <- beta0
alz$beta1 <- beta1

alz_long_beta <- alz %>%
  mutate(id = row_number()) %>%
  select(id, sex, job, age, inkomen, adl_disc, wzc,
         ab_base, tau_base, edu, bmi, cdrsb_base, beta0, beta1) %>%
  pivot_longer(cols = c(beta0, beta1),
               names_to = "param",
               values_to = "beta")

alz_long_beta$param <- as.factor(alz_long_beta$param)

model_2stage <- gls(beta ~ param + sex + job + age + inkomen +
                      adl_disc + wzc + ab_base + tau_base +
                      edu + bmi + cdrsb_base,
                    data = alz_long_beta,
                    correlation = corSymm(form = ~ 1 | id),
                    weights     = varIdent(form = ~ 1 | param),
                    method      = "ML",
                    na.action   = na.exclude)

model_2stage$coefficients

## Then I should reduce the model?
## But if I do like this I do not take into consideration that I would like
## to exclude different covariates for beta0 and beta1...

## Actually, if we proceed in this way we need to take into consideration
## the same regressors both for beta0 and beta1

## So I prefer to run a lm regression for both the variables, maybe by simply
## considering a different weight for the estimations carried on by the 
## fact that we are considering estimates from a model

## Fix the NA values
# sostituisci NA/0 con un piccolo valore, ad esempio
sigma_stage_1_clean <- ifelse(is.na(sigma_stage_1) | sigma_stage_1 == 0,
                              min(sigma_stage_1[sigma_stage_1 > 0], na.rm=TRUE),
                              sigma_stage_1)

weights_lm <- 1 / sigma_stage_1_clean^2

model_beta0 <- lm(beta0 ~ trial + sex + job + age + inkomen + 
                    adl_disc + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base,
                  #weights = weights_lm,
                  data = alz)

model_beta1 <- lm(beta1 ~ trial + sex + job + age + inkomen + 
                    adl_disc + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base,
                  #weights = weights_lm,
                  data = alz)

## So these are the starting points
## We want to reduce the models

step_beta_0 <- step(model_beta0)

model_beta0_final <- step_beta_0$call

## Se servisse
# lm(formula = beta0 ~ job + age + inkomen + adl_disc + wzc + ab_base + 
# tau_base + edu + bmi + cdrsb_base, data = alz)

step_beta_1 <- step(model_beta1)

model_beta1_final <- step_beta_1$call

## Se servisse
# lm(formula = beta1 ~ adl_disc + wzc + edu + cdrsb_base, data = alz)

coeff_beta0_final <- model_beta0_final$coefficients
coeff_beta1_final <- model_beta1_final$coefficients



#### PROVA CON TRIAL ####

## Stage 1 the same

### STAGE 2 MODEL ###

model_beta0_trial <- lm(beta0 ~ sex + job + age + inkomen + 
                    adl_num + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base + trial,
                  #weights = weights_lm,
                  data = alz)

model_beta1_trial <- lm(beta1 ~ sex + job + age + inkomen + 
                    adl_num + wzc + ab_base + tau_base +
                    edu + bmi + cdrsb_base + trial,
                  #weights = weights_lm,
                  data = alz)

step0_trial <- step(model_beta0_trial)
step1_trial <- step(model_beta1_trial)

model_beta0_trial_final <- step0_trial$call
model_beta1_trial_final <- step1_trial$call



### E CON ADL CONTINUA? ###

model_beta0_num <- lm(beta0 ~ sex + job + age + inkomen + 
                          adl_num + wzc + ab_base + tau_base +
                          edu + bmi + cdrsb_base,
                        #weights = weights_lm,
                        data = alz)

model_beta1_num <- lm(beta1 ~ sex + job + age + inkomen + 
                          adl_num + wzc + ab_base + tau_base +
                          edu + bmi + cdrsb_base,
                        #weights = weights_lm,
                        data = alz)

step0_num <- step(model_beta0_num)
step1_num <- step(model_beta1_num)

model_beta0_num <- step0_num$call
model_beta1_num <- step1_num$call

#> step0_num$call
#lm(formula = beta0 ~ sex + job + age + adl_num + wzc + ab_base + 
#     bmi + cdrsb_base, data = alz)
#> step1_num$call
#lm(formula = beta1 ~ wzc + edu + cdrsb_base, data = alz)

