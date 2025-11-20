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



##### EXPLORATORY DATA ANALYSIS ####

## For the spaghetti plot we need a random sample

casual1 <- sample(1:length(alz$patid), 20)

## Now we can start by looking at random values for the mean and see
## if we can work on the mean and so on

alz_rist1 <- alz[casual1, ]
alz_rist1_df <- data.frame(alz_rist1)

alz_rist1_long <- alz_rist1_df %>%
  pivot_longer(

    cols = matches("^(bprs|cdrsb|abpet|taupet)\\d+$"),
    
    names_to = c(".value", "year"),
    
    names_pattern = "(bprs|cdrsb|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                            # numeric
    sample = factor(rep(1:nrow(alz_rist1_df), each = 7)) # ID per ogni paziente
  )


### SPAGHETTI PLOT

ggplot(alz_rist1_long, aes(x = year, y = bprs, group = patid, 
                           color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 1, show.legend = FALSE, size = 1.1) +
  theme_bw() +
  labs(title = "Time Evolution of BPRS")


### BEHAVIOUR WRT DIFFERENT GROUPS
## General behavior
ggplot(alz_long, aes(x = year, y = bprs)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior")

ggplot(alz_long, aes(x = year, y = bprs, group = sex, color = sex)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt sex")

ggplot(alz_long, aes(x = year, y = bprs, group = trial, color = trial)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt trial")

ggplot(alz_long, aes(x = year, y = bprs, group = trial, color = trial)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt trial without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = age, color = age)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age")

ggplot(alz_long, aes(x = year, y = bprs, group = age, color = age)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = edu, color = edu)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt edu")

ggplot(alz_long, aes(x = year, y = bprs, group = bmi, color = bmi)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt bmi")

ggplot(alz_long, aes(x = year, y = bprs, group = bmi, color = bmi)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt bmi without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = inkomen, color = inkomen)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt inkomen")

ggplot(alz_long, aes(x = year, y = bprs, group = inkomen, color = inkomen)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt inkomen without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = job, color = job)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt job")

ggplot(alz_long, aes(x = year, y = bprs, group = adl, color = adl)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt adl")

ggplot(alz_long, aes(x = year, y = bprs, group = wzc, color = wzc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt wzc")


## Discretized categorical variables
ggplot(alz_long, aes(x = year, y = bprs, group = age_disc, color = age_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age_disc without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = bmi_disc, color = bmi_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt bmi_disc without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = inkomen_disc, color = inkomen_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt inkomen_disc without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = adl_disc, color = adl_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt adl_disc without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = cdrsb_disc, color = cdrsb_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt cdrsb_disc without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = abpet_disc, color = abpet_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt abpet_disc without variance")

ggplot(alz_long, aes(x = year, y = bprs, group = taupet_disc, color = taupet_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt taupet_disc without variance")

### EMPIRICAL VARIANCE
var_by_year <- alz_long %>%
  group_by(year) %>%
  summarise(var_bprs = var(bprs, na.rm = TRUE))

ggplot(var_by_year, aes(x = year, y = var_bprs)) +
  geom_line(size = 1.5) + geom_point(size=3) +
  theme_minimal() +
  labs(title = "Variance of BPRS over time")

### COVARIANCE MATRIX
cov_matrix_bprs <- cov(alz[, c(18:24)], use = "pairwise.complete.obs")
round(cov_matrix_bprs, 2)

heatmap(cov_matrix_bprs, main = "Covariance matrix of BPRS")


### CORRELATION MATRIX
cor_matrix_bprs <- cor(alz[, c(18:24)], use = "pairwise.complete.obs")
round(cor_matrix_bprs, 2)

heatmap(cor_matrix_bprs, main = "Correlation matrix of BPRS")


### INFORMATIVE DROPOUT
#I find the latest reading for each patient
ultima_rilevazione <- alz_long %>%
  group_by(patid) %>%
  summarise(last_time = max(year))

#  indidcator 1 = there is untile the last one observation
dati <- alz_long %>%
  left_join(ultima_rilevazione, by = "patid") %>%
  mutate(last_visit = ifelse(year == last_time, 1, 0))

ultimo_anno <- max(alz_long$year, na.rm = TRUE)

bprs_baseline <- alz_long %>%
  group_by(patid) %>%
  summarise(
    # first value of BPRS
    BPRS_start = bprs[!is.na(bprs) & year == min(year[!is.na(bprs)])][1],
    # latest reading
    last_year = max(year[!is.na(bprs)], na.rm = TRUE)
  ) %>%
  mutate(
    completed = ifelse(last_year == ultimo_anno, 1, 0)
  )

bprs_baseline %>%
  group_by(completed) %>%
  summarise(
    mean_BPRS = mean(BPRS_start, na.rm = TRUE),
    sd_BPRS = sd(BPRS_start, na.rm = TRUE),
    n = n()
  )

t_test <- t.test(BPRS_start ~ completed, data = bprs_baseline)
t_test
#conclusion:
#there is an informative drop out, who start the study with high levels of bprs
#drop out earlier than those who start with low levels of bprs

