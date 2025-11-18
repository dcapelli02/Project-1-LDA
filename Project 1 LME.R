library(nlme)
library(ggplot2)
library(dplyr)
library(ggrepel)  
library(broom)
library(lme4)
library(emdbook)
library(RLRsim)
library(robustlmm)

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


#### 1. Preliminary Mean Structure - OLS ####
mod_ols <- lm(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  na.action = na.exclude
)


#### 2. Preliminary Random-eﬀects Structure ####
# Residuals OLS
df_ols <- alz_long %>%
  mutate(
    resid = resid(mod_ols),        # OLS
    resid2 = resid(mod_ols)^2,              # squared
    sample = as.factor(sample)     # ID paziente come fattore
  ) %>%
  dplyr::select(sample, year, resid, resid2) %>%
  filter(is.finite(resid), is.finite(year))

set.seed(123)  # reproducibility
top30 <- sample(unique(df_ols$sample), 30)

df_ols_30 <- df_ols %>% filter(sample %in% top30)

# Plot
ggplot(df_ols_30, aes(x = year, y = resid, group = sample, color = sample)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  labs(title = "Profili residui OLS (30 pazienti)",
       x = "Anno", y = "Residui OLS") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",        
    strip.text = element_text(face = "bold")
  )

# ---> RANDOM INTERCEPT

# Smoothing
ggplot(df_ols, aes(x = year, y = resid2)) +
  geom_point(alpha = 0.15, color = "grey40", size = 1) +
  geom_smooth(se = FALSE, method = "loess",
              color = "#1f77b4", linewidth = 1.3) +
  labs(title = "Smoothed variance OLS",
       x = "Year", y = "OLS Residuals^2") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ---> RANDOM SLOPES

mod_lme <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  random = ~ year | sample,
  na.action = na.exclude
)

# Estrazione D e sigma^2 dall'LME

vc <- VarCorr(mod_lme)

# Variance
d00 <- as.numeric(vc["(Intercept)", "Variance"])
d11 <- as.numeric(vc["year", "Variance"])
sigma2 <- as.numeric(vc["Residual", "Variance"])

# Correlation intercept - slope
rho <- as.numeric(vc["year", "Corr"])

# Covariance
d01 <- rho * sqrt(d00 * d11)

# Matrice D
D <- matrix(c(d00, d01, d01, d11), nrow = 2)

# Components of D
d00 <- D[1, 1]                 # Var(random intercept)
d11 <- D[2, 2]                 # Var(random slope year)
d01 <- D[1, 2]                 # Cov(intercept, year)

# Grid
t_seq <- sort(unique(model.frame(mod_ols)$year))

var_fun <- d00 + 2 * d01 * t_seq + d11 * (t_seq^2) + sigma2

df_var <- data.frame(
  year = t_seq,
  resid2 = var_fun,
  modello = "Fitted variance (LME)"
)

ggplot() +
  geom_smooth(data = df_ols, aes(x = year, y = resid2),
              method = "loess", se = FALSE,
              linewidth = 1.3, color = "#1f77b4") +
  geom_line(data = df_var, aes(x = year, y = resid2),
            linewidth = 1.2, linetype = "dashed", color = "#d62728") +
  geom_text_repel(
    data = df_ols %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "Smoothed OLS"),
    color = "#1f77b4", hjust = 0, nudge_x = 1, size = 4
  ) +
  geom_text_repel(
    data = df_var %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "Fitted Variance LME"),
    color = "#d62728", hjust = 0, nudge_x = 1, size = 4
  ) +
  labs(
    x = "Year",
    y = "Squared Residuals"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# Model LME only random intercept

mod_lme1 <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  random = ~ 1 | sample,
  na.action = na.exclude
)


vc_1 <- VarCorr(mod_lme1)

d00_1 <- as.numeric(vc_1["(Intercept)", "Variance"])   # Var(random intercept)
sigma2_1 <- as.numeric(vc_1["Residual", "Variance"])   # Var(residuo)

var_fun_1 <- d00_1 + sigma2_1 

df_var_1 <- data.frame(
  year = sort(unique(df_ols$year)),
  resid2 = var_fun_1,
  modello = "Fitted variance (LME)"
)


ggplot() +
  geom_smooth(data = df_ols, aes(x = year, y = resid2),
              method = "loess", se = FALSE,
              linewidth = 1.3, color = "#1f77b4") +
  geom_hline(data = df_var_1, yintercept = var_fun_1,
             linetype = "dashed", color = "#d62728", linewidth = 1.2) +
  labs(x = "Year", y = "Squared Residuals") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

## Comparison

ggplot() +
  geom_smooth(data = df_ols, aes(x = year, y = resid2),
              method = "loess", se = FALSE,
              linewidth = 1.3, color = "#1f77b4") +
  geom_line(data = df_var, aes(x = year, y = resid2),
            linewidth = 1.2, linetype = "dashed", color = "#d62728") +
  geom_line(data = df_var_1, aes(x = year, y = resid2),
            linewidth = 1.2, linetype = "dashed", color = "orange") +
  geom_text_repel(
    data = df_ols %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "OLS (smoothed)"),
    color = "#1f77b4", hjust = 0, nudge_x = 1, size = 4
  ) +
  geom_text_repel(
    data = df_var %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "Variance LME"),
    color = "#d62728", hjust = 0, nudge_x = 1, size = 4
  ) +
  geom_text_repel(
    data = df_var_1 %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "Random intercept"),
    color = "orange", hjust = 0, nudge_x = 1, size = 4
  ) +
  labs(
    x = "Year",
    y = "Sqrd Residuals"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


#### 3. Residual Covariance Structure ####
mod_lme_hetero <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base + ab_base + tau_base) * year,
  data = alz_long,
  random = ~ 1 + year | sample,                
  #correlation = corGaus(form = ~ year | sample),# AR(1) serial correlation
  weights = varIdent(form = ~ 1 | year_seq),
  na.action = na.exclude,
  method = "REML"
)

mod_lme_gaus <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base + ab_base + tau_base) * year,
  data = alz_long,
  random = ~ 1 + year | sample,                # intercetta + slope casuali
  correlation = corGaus(form = ~ year | sample),# AR(1) serial correlation
  na.action = na.exclude,
  method = "REML"
)

mod_lme_exp <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base + ab_base + tau_base) * year,
  data = alz_long,
  random = ~ 1 + year | sample,                # intercetta + slope casuali
  correlation = corExp(form = ~ year | sample),# AR(1) serial correlation
  na.action = na.exclude,
  method = "REML"
)

anova(mod_lme_hetero, mod_lme_gaus,
      mod_lme_exp, mod_lme) 

## Comparison

vc_gaus <- VarCorr(mod_lme_gaus)
vc_exp <- VarCorr(mod_lme_exp)
vc_het <- VarCorr(mod_lme_hetero)

# Variances
d00_gaus <- as.numeric(vc_gaus["(Intercept)", "Variance"])
d11_gaus <- as.numeric(vc_gaus["year", "Variance"])
sigma2_gaus <- as.numeric(vc_gaus["Residual", "Variance"])

d00_exp <- as.numeric(vc_exp["(Intercept)", "Variance"])
d11_exp <- as.numeric(vc_exp["year", "Variance"])
sigma2_exp <- as.numeric(vc_exp["Residual", "Variance"])

d00_het <- as.numeric(vc_het["(Intercept)", "Variance"])
d11_het <- as.numeric(vc_het["year", "Variance"])
sigma2_het <- as.numeric(vc_het["Residual", "Variance"])

## Hetweroskedasticity

weights_het <- coef(mod_lme_hetero$modelStruct$varStruct, unconstrained = FALSE)
weights_het
levels_year <- names(weights_het)
weights_full <- c("0" = 1, weights_het)

sigma2_het <- sigma2_het * weights_full

# Correlation intercetta-slope
rho_gaus <- as.numeric(vc_gaus["year", "Corr"])
rho_exp <- as.numeric(vc_exp["year", "Corr"])
rho_het <- as.numeric(vc_het["year", "Corr"])

# Covariance = correlation * sqrt(var_intercetta * var_slope)
d01_gaus <- rho_gaus * sqrt(d00_gaus * d11_gaus)
d01_exp <- rho_exp * sqrt(d00_exp * d11_exp)
d01_het <- rho_het * sqrt(d00_het * d11_het)

# Matrice D
D_gaus <- matrix(c(d00_gaus, d01_gaus, d01_gaus, d11_gaus), nrow = 2)
D_exp <- matrix(c(d00_exp, d01_exp, d01_exp, d11_exp), nrow = 2)
D_het <- matrix(c(d00_het, d01_het, d01_het, d11_het), nrow = 2)

# Components of D
d00_gaus <- D_gaus[1, 1]                 # Var(random intercept)
d11_gaus <- D_gaus[2, 2]                 # Var(random slope year)
d01_gaus <- D_gaus[1, 2]                 # Cov(intercept, year)

d00_exp <- D_exp[1, 1]                 # Var(random intercept)
d11_exp <- D_exp[2, 2]                 # Var(random slope year)
d01_exp <- D_exp[1, 2]                 # Cov(intercept, year)

d00_het <- D_het[1, 1]                 # Var(random intercept)
d11_het <- D_het[2, 2]                 # Var(random slope year)
d01_het <- D_het[1, 2]                 # Cov(intercept, year)

t_seq <- sort(unique(model.frame(mod_ols)$year))

var_fun_gaus <- d00_gaus + 2 * d01_gaus * t_seq + d11_gaus * (t_seq^2) + sigma2_gaus
var_fun_exp <- d00_exp + 2 * d01_exp * t_seq + d11_exp * (t_seq^2) + sigma2_exp
var_fun_het <- d00_het + 2 * d01_het * t_seq + d11_het * (t_seq^2) + sigma2_het

df_var_gaus <- data.frame(
  year = t_seq,
  resid2 = var_fun_gaus,
  modello = "Fitted variance (LME)"
)

df_var_exp <- data.frame(
  year = t_seq,
  resid2 = var_fun_exp,
  modello = "Fitted variance (LME)"
)

df_var_het <- data.frame(
  year = t_seq,
  resid2 = var_fun_het,
  modello = "Fitted variance (LME)"
)

## First comparison

ggplot() +
  geom_smooth(data = df_ols, aes(x = year, y = resid2),
              method = "loess", se = FALSE,
              linewidth = 1.3, color = "#1f77b4") +
  geom_line(data = df_var, aes(x = year, y = resid2),
            linewidth = 1.2, linetype = "dashed", color = "#d62728") +
  geom_line(data = df_var_gaus, aes(x = year, y = resid2),
            linewidth = 1.2, linetype = "dashed", color = "orange") +
  geom_line(data = df_var_exp, aes(x = year, y = resid2),
            linewidth = 1.2, linetype = "dashed", color = "purple") +
  geom_line(data = df_var_het, aes(x = year, y = resid2),
            linewidth = 1.2, linetype = "dashed", color = "green") +
  geom_text_repel(
    data = df_ols %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "OLS (smussato)"),
    color = "#1f77b4", hjust = 0, nudge_x = 1, size = 4
  ) +
  geom_text_repel(
    data = df_var %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "LME"),
    color = "#d62728", hjust = 0, nudge_x = 1, size = 4
  ) +
  geom_text_repel(
    data = df_var_gaus %>% slice_tail(n = 1),
    aes(x = year, y = resid2, label = "LME gauss"),
    color = "orange", hjust = 0, nudge_x = 1, size = 4
  ) +
  labs(
    x = "Anno",
    y = "Varianza / residui al quadrato"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )



#### 4. Reduction of Preliminary Random-eﬀects Structure ####

mod_lme1 <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  random = ~ 1 | sample,
  na.action = na.exclude
)
mod_lme2 <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  random = list(
    trial = pdIdent(~ 1),        # intercept per centro
    sample = pdSymm(~ 1 + year)   # intercept+slope per paziente
  ),
  na.action = na.exclude,
  method = "REML"
)


#test slopes+intercept vs intercept
#mixture of chis-squared

lrt <- anova(mod_lme1, mod_lme) 
LRT_value <- lrt$"L.Ratio"[2]
pval <- dchibarsq(LRT_value, df = c(1,2))
cat("p-value mixture of chi-squared:", pval)


#### 5. Reduction of Preliminary Mean Structure ####
summary(mod_lme)


mod_full_lm <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  random = ~ year | sample,
  na.action = na.exclude,
  method = "ML" 
)  ######## qua si usa ML

mod_full_reml <- lme(
  bprs ~ (trial + sex + age + edu + bmi + inkomen + job +
            adl_num + wzc + cdrsb_base +
            ab_base + tau_base) * year,
  data = alz_long,
  random = ~ year | sample,
  na.action = na.exclude,
  method = "REML" 
) 


## Drop parameters one at a time

obs <- sum(alz$n_obs_data)

drop1(mod_full_lm, 
      test = "Chisq",
      k = log(obs))

## We can try to remove the terms with higher p-value?

red1 <- update(mod_full_lm, . ~ . - trial:year - bmi:year - tau_base:year)
red1$call

anova(mod_full_lm, red1)

## Go on

drop1(red1, test = "Chisq", k = log(obs))

## Forse ab_base:year e tau_base

red2 <- update(red1, . ~ . - ab_base:year - tau_base)
red2$call

anova(red1, red2)   ## abbasranza ok

## Avanti
drop1(red2, test = "Chisq", k = log(obs))

## Forse si puo togliere ab_base (però forzata interpretazione biologica)
## hanno p-valu 0.11 circa anche sex:year, inkomen:year, wzc:year

## Proviamo entrambi

red3 <- update(red2, . ~ . - ab_base)
red3$call

anova(red2, red3) ## ok teniamo red3

red4 <- update(red2, . ~ . - sex:year - inkomen:year - wzc:year - ab_base)

anova(red2, red4)   ## ok keep red4
anova(red3, red4)   ## ok keep red4


drop1(red4, test = "Chisq", k = log(obs))

## Togliamo inkomen (sex forse ma io non lo toglierei)

red5 <- update(red4, . ~ . - inkomen - sex)

anova(red4, red5)  ## yes

drop1(red5, test = "Chisq", k = log(obs))   ## ok dovremmo essere arrivati alla fine

red5$call

## Final comparison

anova(mod_full_lm, red5)  ## ok




#### 6. Analysis of random effects ####

final_model_mixed <- lme(
  fixed = bprs ~ trial + age + edu + bmi + job + adl_num + 
    wzc + cdrsb_base + year + age:year + edu:year + job:year + 
    adl_num:year + cdrsb_base:year,
  data = alz_long, 
  random = ~year | sample,
  method = "REML", 
  na.action = na.exclude)


# Slope and Intercept

re <- ranef(final_model_mixed)
head(re)
df_re <- data.frame(id = rownames(re),
                    intercept = re[, "(Intercept)"],
                    slope_time = re[, "year"])

hist(df_re$intercept, main="Intercepts", xlab="Intercepts")

plot(df_re$intercept, df_re$slope_time,
     xlab="Intercepts", ylab="Slopes time")
text(df_re$intercept, df_re$slope_time, labels=df_re$id, pos=4, cex=0.7)
#strong correlation between slopes and intercept

hist(df_re$slope_time, main="Slopes time", xlab="Slopes time")

## Interessante forse vedere qualche esempio di come si fitta

set.seed(1234)
casual_sample <- sample(1:length(alz$patid), 5)
alz_long$pred <- predict(final_model_mixed, level = 1)  ## con random effect
alz_long$pred_no_re <- predict(final_model_mixed, level = 0)  ## senza random effect

rdeff_plot <- alz_long %>% filter(sample %in% casual_sample)

ggplot(rdeff_plot, aes(x = year, y = bprs, group = patid, 
                       color = sample)) + 
  geom_line(alpha = 1, size = 1, aes(linetype = "Observed")) +
  geom_line(data = rdeff_plot, aes(x = year, y = pred, group= patid,
                                   color = sample, linetype = "Pred RE"),
            size = 0.7) +
  geom_line(data = rdeff_plot, aes(x = year, y = pred_no_re, group= patid,
                                   color = sample, linetype = "Pred no RE"),
            size = 0.7) +
  scale_linetype_manual(values = c("Observed" = "solid",
                                   "Pred RE" = "dashed",
                                   "Pred no RE" = "dotdash")) +
  theme_bw() +
  labs(title = "Time Evolution of BPRS", linetype = "") +
  guides(color = "none") +   # 🔹 nasconde la legenda del colore (sample)
  theme(legend.key.width = unit(2, "lines"))  # 🔹 restringe la colonna della legenda



## Prendiamo quello che sembra un outlier
rdeff_plot <- alz_long %>% filter(sample == 1116)

ggplot(rdeff_plot, aes(x = year, y = bprs, group = patid, 
                       color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 1, show.legend = FALSE, size = 1) +
  geom_line(data = rdeff_plot, aes(x = year, y = pred, group= patid,
                                   color = sample, show.legend = FALSE),
            linetype = "dashed", size = 0.7
  ) +
  geom_line(data = rdeff_plot, aes(x = year, y = pred_no_re, group= patid,
                                   color = sample, show.legend = FALSE),
            linetype = "dotdash", size = 0.7
  ) +
  theme_bw() +
  labs(title = "BPRS against year of measuring")

rdeff_plot <- alz_long %>% filter(sample == 930)

ggplot(rdeff_plot, aes(x = year, y = bprs, group = patid, 
                       color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 1, show.legend = FALSE, size = 1) +
  geom_line(data = rdeff_plot, aes(x = year, y = pred, group= patid,
                                   color = sample, show.legend = FALSE),
            linetype = "dashed", size = 0.7
  ) +
  geom_line(data = rdeff_plot, aes(x = year, y = pred_no_re, group= patid,
                                   color = sample, show.legend = FALSE),
            linetype = "dotdash", size = 0.7
  ) +
  theme_bw() +
  labs(title = "BPRS against year of measuring")

## Forse da analizzare un po' i random effect


## Heterogeneity model? Io lo salterei...




