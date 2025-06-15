library(lmtest)
library(forecast)
library(stargazer)
library(fredr)
library(ggplot2)
library(urca)
library(vars)
library(knitr)
library(gridExtra)
library(dplyr)
library(reshape2)

fredr_set_key("68930503314ecaad758dc1a52090ca2f")

start_date <- as.Date("1968-01-01")
end_date <- as.Date("2025-05-31")

# DATA COLLECTION AND CLEANING
cpi_data <- fredr(series_id = "CPIAUCSL", observation_start = start_date, observation_end = end_date)
ffr_data <- fredr(series_id = "FEDFUNDS", observation_start = start_date, observation_end = end_date)
unemployment_data <- fredr(series_id = "UNRATE", observation_start = start_date, observation_end = end_date)
gdp_data <- fredr(series_id = "GDPC1", observation_start = start_date, observation_end = end_date)
commodity_data <- fredr(series_id = "WTISPLC", observation_start = start_date, observation_end = end_date)
money_data <- fredr(series_id = "M2REAL", observation_start = start_date, observation_end = end_date)

inflation_ts <- ts(cpi_data$value, start = c(1968, 1), frequency = 12)
interest_ts <- ts(ffr_data$value, start = c(1968, 1), frequency = 12)
unemployment_ts <- ts(unemployment_data$value, start = c(1968, 1), frequency = 12)
commodity_ts <- ts(commodity_data$value, start = c(1968, 1), frequency = 12)
money_ts <- ts(money_data$value, start = c(1968, 1), frequency = 12)

gdp_ts <- ts(gdp_data$value, start = c(1968, 1), frequency = 4)
gdp_monthly <- approx(time(gdp_ts), gdp_ts, xout = time(inflation_ts))$y
gdp_monthly_ts <- ts(gdp_monthly, start = c(1968, 1), frequency = 12)
gdp_monthly_ts <- na.approx(gdp_monthly_ts)

# STATIONARITY TESTING AND TRANSFORMATIONS
diff_inflation <- 100 * diff(log(inflation_ts), lag = 12)
gdp_monthly_diff <- diff(gdp_monthly_ts)
money_diff <- diff(money_ts)

adf_test_inf <- ur.df(diff_inflation, type = "drift", selectlags = "AIC")
adf_test_int <- ur.df(interest_ts, type = "drift", selectlags = "AIC")
adf_test_emp <- ur.df(unemployment_ts, type = "drift", selectlags = "AIC")
adf_test_gdp_diff <- ur.df(gdp_monthly_diff, type = "drift", selectlags = "AIC")
adf_test_com <- ur.df(commodity_ts, type = "trend", selectlags = "AIC")
adf_test_m2_diff <- ur.df(money_diff, type = "drift", selectlags = "AIC")


summary(adf_test_inf)
summary(adf_test_int)
summary(adf_test_emp)
summary(adf_test_gdp_diff)
summary(adf_test_com)
summary(adf_test_m2_diff)


# BIVARIATE VAR MODEL
combined_data_bivar <- na.omit(ts.union(diff_inflation, interest_ts))
colnames(combined_data_bivar) <- c("Inflation", "Interest")

VAR_lag_selection_bivar <- VARselect(combined_data_bivar, lag.max = 30, type = "const")
bic_lag_bivar <- VAR_lag_selection_bivar$selection["SC(n)"]
var_model_bivariate <- VAR(combined_data_bivar, p = bic_lag_bivar, type = "const")
summary(var_model_bivariate)

# MULTIVARIATE VAR MODEL
model_data_multi <- ts.union(
  commodity_ts,
  diff_inflation,
  gdp_monthly_diff,
  unemployment_ts,
  interest_ts,
  money_diff
)
colnames(model_data_multi) <- c(
  "Commodity_Prices",
  "Inflation",
  "GDP_Growth",
  "Unemployment",
  "Interest",
  "Money_Supply_Growth"
)
model_data_multi <- na.omit(model_data_multi)

VAR_lag_selection_multi <- VARselect(model_data_multi, lag.max = 30, type = "const")
bic_lag_multi <- VAR_lag_selection_multi$selection["SC(n)"]
var_model_multi <- VAR(model_data_multi, p = bic_lag_multi, type = "const")
summary(var_model_multi)


# DIAGNOSTIC CHECKS

#Bivariate Model
var_residuals_bivar <- residuals(var_model_bivariate)
ccf(var_residuals_bivar[, "Inflation"], var_residuals_bivar[, "Interest"],
    main = "Bivariate VAR: CCF - Inflation vs. Interest Residuals",
    sub = "Check spike at lag 0 for significance")
serial.test(var_model_bivariate, lags.pt = 24 , type = "PT.asymptotic")


#Multivariate Model
var_residuals_multi <- residuals(var_model_multi)

residual_correlation_matrix_multi <- cor(var_residuals_multi)
print(round(residual_correlation_matrix_multi, 3))

serial_test_lags <- max(bic_lag_multi + 10, 2 * bic_lag_multi, 24)
serial.test(var_model_multi, lags.pt = serial_test_lags, type = "PT.asymptotic")



# GRANGER CAUSALITY ANALYSIS
granger_int_to_inf_bivar <- causality(var_model_bivariate, cause = "Interest")
granger_int_to_inf_multi <- causality(var_model_multi, cause = "Interest")

print(granger_int_to_inf_bivar$Granger)
print(granger_int_to_inf_multi$Granger)




# IMPULSE RESPONSE ANALYSIS
horizon <- 48

irf_int_to_inf_bivar <- irf(var_model_bivariate, impulse = "Interest", response = "Inflation",
                            n.ahead = horizon, boot = TRUE, seed = 12345, runs = 1000, ci = 0.95)
irf_int_to_inf_multi <- irf(var_model_multi, impulse = "Interest", response = "Inflation",
                            n.ahead = horizon, boot = TRUE, seed = 12345, runs = 1000, ci = 0.95)

irf_data_bivar <- data.frame(
  Horizon = 0:horizon,
  Response = c(irf_int_to_inf_bivar$irf$Interest),
  Lower = c(irf_int_to_inf_bivar$Lower$Interest),
  Upper = c(irf_int_to_inf_bivar$Upper$Interest),
  Model = "Bivariate VAR"
)

irf_data_multi <- data.frame(
  Horizon = 0:horizon,
  Response = c(irf_int_to_inf_multi$irf$Interest),
  Lower = c(irf_int_to_inf_multi$Lower$Interest),
  Upper = c(irf_int_to_inf_multi$Upper$Interest),
  Model = "Multivariate VAR"
)

irf_combined <- rbind(irf_data_bivar, irf_data_multi)

ggplot(irf_combined, aes(x = Horizon, y = Response, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("Bivariate VAR" = "#1f77b4", "Multivariate VAR" = "#ff7f0e")) +
  scale_fill_manual(values = c("Bivariate VAR" = "#1f77b4", "Multivariate VAR" = "#ff7f0e")) +
  labs(x = "Horizon (Months)", y = "Response of Inflation") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# CONDITIONAL FORECASTING
train_data_bivar <- window(combined_data_bivar, end = c(2022,2))
train_data_multi <- window(model_data_multi, end = c(2022,2))

var_train_bivar <- VAR(train_data_bivar, p = bic_lag_bivar, type = "const")
var_train_multi <- VAR(train_data_multi, p = bic_lag_multi, type = "const")

forecast_start_date <- c(2022, 3)
forecast_end_date <- c(2024, 2)
num_forecast_steps <- 24

# Bivariate scenario
scenario_ts_dates <- ts(NA, start = forecast_start_date, end = forecast_end_date, frequency = 12)
scenario_matrix_bivar <- matrix(NA, nrow = length(scenario_ts_dates), ncol = ncol(combined_data_bivar))
colnames(scenario_matrix_bivar) <- colnames(combined_data_bivar)
scenario_matrix_bivar[,"Interest"] <- 0.08
scenario_matrix_bivar[,"Inflation"] <- NA

forecast_cond_bivar <- predict(var_train_bivar, n.ahead = num_forecast_steps, dumvar = scenario_matrix_bivar)
pred_inflation_cond_bivar <- forecast_cond_bivar$fcst$Inflation[, "fcst"]

# Multivariate scenario
scenario_matrix_multi <- matrix(NA, nrow = length(scenario_ts_dates), ncol = ncol(model_data_multi))
colnames(scenario_matrix_multi) <- colnames(model_data_multi)
scenario_matrix_multi[,"Interest"] <- 0.08

for (col_name in colnames(scenario_matrix_multi)) {
  if (col_name != "Inflation" && col_name != "Interest") {
    scenario_matrix_multi[, col_name] <- window(model_data_multi[, col_name], start = forecast_start_date, end = forecast_end_date)
  }
}

forecast_cond_multi <- predict(var_train_multi, n.ahead = num_forecast_steps, dumvar = scenario_matrix_multi)
pred_inflation_cond_multi <- forecast_cond_multi$fcst$Inflation[, "fcst"]

actual_inflation_counterfactual <- window(combined_data_bivar[, "Inflation"], start = forecast_start_date, end = forecast_end_date)

forecast_data <- data.frame(
  Date = as.Date(time(actual_inflation_counterfactual)),
  Actual = as.numeric(actual_inflation_counterfactual),
  Bivariate_Predicted = as.numeric(pred_inflation_cond_bivar),
  Multivariate_Predicted = as.numeric(pred_inflation_cond_multi)
)

forecast_long <- melt(forecast_data, id.vars = "Date", variable.name = "Series", value.name = "Inflation")

forecast_long$Series <- factor(forecast_long$Series, 
                               levels = c("Actual", "Bivariate_Predicted", "Multivariate_Predicted"),
                               labels = c("Actual", "Bivariate VAR", "Multivariate VAR"))

ggplot(forecast_long, aes(x = Date, y = Inflation, color = Series, linetype = Series)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Actual" = "black", "Bivariate VAR" = "blue", "Multivariate VAR" = "red")) +
  scale_linetype_manual(values = c("Actual" = "solid", "Bivariate VAR" = "dashed", "Multivariate VAR" = "dotted")) +
  labs(x = "Date", y = "Inflation (YoY %)") +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# SENSITIVITY ANALYSIS - MODEL COMPARISON
lag_minus_1 <- max(1, bic_lag_multi - 1)
lag_plus_1 <- bic_lag_multi + 1
num_vars <- ncol(model_data_multi)
num_obs_total <- nrow(model_data_multi)

var_model_minus_1 <- VAR(model_data_multi, p = lag_minus_1, type = "const")
var_model_plus_1 <- VAR(model_data_multi, p = lag_plus_1, type = "const")

# Extract inflation equation coefficients and statistics
extract_inflation_stats <- function(var_model, model_name) {
  inf_eq <- var_model$varresult$Inflation
  p <- var_model$p
  T <- num_obs_total - p  # Adjusted number of observations
  k <- num_vars
  total_coeffs <- k * (k * p + 1)
  obs_per_coeff <- T / total_coeffs
  
  data.frame(
    Model = model_name,
    Lags = p,
    R_squared = summary(inf_eq)$r.squared,
    Adj_R_squared = summary(inf_eq)$adj.r.squared,
    F_statistic = summary(inf_eq)$fstatistic[1],
    Log_likelihood = logLik(inf_eq)[1],
    Obs_per_Coeff = obs_per_coeff
  )
}

model_comparison <- rbind(
  extract_inflation_stats(var_model_minus_1, paste0("VAR(", lag_minus_1, ")")),
  extract_inflation_stats(var_model_multi, paste0("VAR(", bic_lag_multi, ")")),
  extract_inflation_stats(var_model_plus_1, paste0("VAR(", lag_plus_1, ")"))
)

# Plot Adjusted R²
p1 <- ggplot(model_comparison, aes(x = factor(Lags), y = Adj_R_squared, fill = factor(Lags))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(Adj_R_squared, 3)), vjust = -0.5, size = 6) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  ylim(0, max(model_comparison$Adj_R_squared) * 1.25) +
  labs(x = "Number of Lags", y = "Adjusted R²") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Plot F-Statistic
p2 <- ggplot(model_comparison, aes(x = factor(Lags), y = F_statistic, fill = factor(Lags))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(F_statistic, 1)), vjust = -0.5, size = 6) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  ylim(0, max(model_comparison$F_statistic) * 1.25) +
  labs(x = "Number of Lags", y = "F-Statistic") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Plot Observations per Coefficient
p3 <- ggplot(model_comparison, aes(x = factor(Lags), y = Obs_per_Coeff, fill = factor(Lags))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(Obs_per_Coeff, 2)), vjust = -0.5, size = 6) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  ylim(0, max(model_comparison$Obs_per_Coeff) * 1.25) +
  labs(x = "Number of Lags", y = "Obs per Coefficient") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Arrange plots
grid.arrange(p1, p2, p3, ncol = 3)
