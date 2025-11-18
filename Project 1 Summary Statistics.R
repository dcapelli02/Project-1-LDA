library(ggplot2)
library(gridExtra)

# Define the outcome variable
outcome <- "bprs"

alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))

alz_long <- alz_long %>%
  left_join(alz %>% select(patid, n_obs_data), by = "patid")


# Compute summary statistics: AUC, endpoint, increment
summary_stats <- alz_long %>%
  group_by(patid) %>%
  summarise(
    trial = first(trial),
    sex = first(sex),
    age = first(age),
    edu = first(edu),
    bmi = first(bmi),
    inkomen = first(inkomen),
    job = first(job),
    adl = first(adl),
    wzc = first(wzc),
    abpet_base = first(abpet_base),
    taupet_base = first(taupet_base),
    cdrsb_base = first(cdrsb_base),
    
    # Baseline Value
    y0 = first(.data[[outcome]]),
    
    # Last Available Measurement
    yini = last(na.omit(.data[[outcome]])),
    
    # Last Time Point Observed
    tmax = last(n_obs_data),
    
    # Increment
    increment = yini - y0,
    
    # Area Under the Curve (AUC) using trapezoidal rule
    AUC = sum(diff(year) * (head(.data[[outcome]], -1) + tail(.data[[outcome]], -1)) / 2, na.rm = TRUE),
    
    # AUC Normalised
    nAUC = AUC / tmax,
    
    # Rate of Increment
    rate = increment / tmax
  )

summary(summary_stats)


#Build a linear model to see how baseline factors explain each person's overall summary statistics:

#AUC: Area Under the Curve normalised
model_nauc <- lm(nAUC ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                data = summary_stats)

summary(model_nauc)

cat("The most significant ones are: sex, edu, inkomen, job, adl, wzc, abpet_base and taupet_base (***)")

#Endpoints
model_endpoint <- lm(yini ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                    data = summary_stats)

summary(model_endpoint)

cat("The most significant ones are: age, HigherEdu, inkomen, job, adl, wzc, abpet_base and cdrsb_base(***)")

#Covariance (Depending also on the bprs baseline) Linear Regression including bprs as a predictor
model_ancova <- lm(yini ~ y0 + trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                   data = summary_stats)

summary(model_ancova)

cat("The most significant ones are: inkomen, job, adl, wzc, abpet_base and cdrsb_base(***) [y0 highly significant (**)]")

#Increment
model_increment <- lm(increment ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                      data = summary_stats)

summary(model_increment)

cat("The most significant ones are: inkomen, job, adl, wzc, abpet_base and cdrsb_base(***)")

# Rate of Increment
model_rate <- lm(rate ~ trial + sex + age + edu + bmi + inkomen + job + adl + wzc + abpet_base + taupet_base + cdrsb_base, 
                      data = summary_stats)

summary(model_rate)

cat("The most significant ones are: sex, edu, inkomen, job, adl, wzc, abpet_base and taupet_base (***)")


AIC(model_nauc, model_endpoint, model_increment, model_ancova, model_rate)


# Create predicted vs observed plots for all models
p1 <- ggplot(summary_stats, aes(x = nAUC, y = fitted(model_nauc))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: nAUC",
       x = "Observed nAUC", 
       y = "Predicted nAUC") +
  theme_minimal()

p2 <- ggplot(summary_stats, aes(x = yini, y = fitted(model_endpoint))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: Endpoint",
       x = "Observed Endpoint", 
       y = "Predicted Endpoint") +
  theme_minimal()

p3 <- ggplot(summary_stats, aes(x = yini, y = fitted(model_ancova))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: Endpoint (ANCOVA)",
       x = "Observed ANCOVA", 
       y = "Predicted ANCOVA") +
  theme_minimal()

p4 <- ggplot(summary_stats, aes(x = increment, y = fitted(model_increment))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: Increment",
       x = "Observed Increment", 
       y = "Predicted Increment") +
  theme_minimal()

p5 <- ggplot(summary_stats, aes(x = rate, y = fitted(model_rate))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed: Rate of Change",
       x = "Observed Rate",
       y = "Predicted Rate") +
  theme_minimal()


grid.arrange(p1, p2, p3, p4, p5, ncol = 2)

# Residuals vs Fitted plots for all models
par(mfrow = c(3, 2))

plot(fitted(model_nauc), residuals(model_nauc),
     main = "Residuals vs Fitted: nAUC",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

plot(fitted(model_endpoint), residuals(model_endpoint),
     main = "Residuals vs Fitted: Endpoint",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

plot(fitted(model_ancova), residuals(model_ancova),
     main = "Residuals vs Fitted: Endpoint (ANCOVA)",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

plot(fitted(model_increment), residuals(model_increment),
     main = "Residuals vs Fitted: Increment",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

plot(fitted(model_rate), residuals(model_rate),
     main = "Residuals vs Fitted: Rate",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")





