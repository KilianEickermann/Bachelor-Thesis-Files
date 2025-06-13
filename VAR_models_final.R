# Load required packages
library(lmtest)
library(forecast)
library(stargazer)
library(fredr)
library(ggplot2)
library(urca)
library(vars)
library(knitr)

# Set FRED API key
fredr_set_key("68930503314ecaad758dc1a52090ca2f")

# Define common start and end dates
start_date <- as.Date("1968-01-01")
end_date <- as.Date("2025-01-06")

# --- SECTION 1: BIVARIATE VAR MODEL (Inflation & Federal Funds Rate) ---

# --- 1.1 Data Collection for Bivariate Model ---
# Download CPI data (All Urban Consumers)
cpi_data <- fredr(series_id = "CPIAUCSL", observation_start = start_date, observation_end = end_date)
# Download Federal Funds Rate data
ffr_data <- fredr(series_id = "FEDFUNDS", observation_start = start_date, observation_end = end_date)

# Create time series objects
inflation_ts <- ts(cpi_data$value, start = c(1968, 1), frequency = 12)
interest_ts <- ts(ffr_data$value, start = c(1968, 1), frequency = 12)

# Plot original series
par(mfrow = c(1, 2))
plot(inflation_ts, main = "Inflation (CPI)", ylab = "Index")
plot(interest_ts, main = "Federal Funds Rate", ylab = "Percent")
par(mfrow = c(1, 1)) # Reset plot layout

# --- 1.2 Data Transformation and Stationarity Testing for Bivariate Model ---
# Calculate annual inflation rate (YoY % change)
diff_inflation <- 100 * diff(log(inflation_ts), lag = 12)

# Test for stationarity
adf_test_inf_bivar <- ur.df(diff_inflation, type = "drift", selectlags = "AIC")
summary(adf_test_inf_bivar)

adf_test_int_bivar <- ur.df(interest_ts, type = "drift", selectlags = "AIC")
summary(adf_test_int_bivar)

2.# Create combined dataset for VAR, handling NAs from differencing
combined_data_bivar <- na.omit(ts.union(diff_inflation, interest_ts))
colnames(combined_data_bivar) <- c("Inflation", "Interest")

# Plot combined stationary series
plot(combined_data_bivar, main = "Bivariate VAR: Inflation and Interest Rate Over Time (Stationary)", xlab = "Year", lwd = 1)

# --- 1.3 VAR Model Estimation for Bivariate Model (using BIC/SC for lag selection) ---
VAR_lag_selection_bivar <- VARselect(combined_data_bivar, lag.max = 30, type = "const")
bic_lag_bivar <- VAR_lag_selection_bivar$selection["SC(n)"]

# Plot BIC/SC criteria for visual inspection of lag selection
plot(VAR_lag_selection_bivar$criteria[3,], main = "Bivariate VAR: Optimal Lag Selection (SC/BIC)", ylab = "SC/BIC Value", xlab = "Lag Order (p)", type = "b", pch = 19, col = "darkgreen")
abline(v = bic_lag_bivar, col = "red", lty = 2)

# Estimate the Bivariate VAR model with BIC-selected lag (p=14)
var_model_bivariate <- VAR(combined_data_bivar, p = bic_lag_bivar, type = "const")

# Model summary for Bivariate VAR
summary(var_model_bivariate)

# Display regression results using stargazer
stargazer(var_model_bivariate$varresult$Inflation, var_model_bivariate$varresult$Interest,
          type = "text",
          column.labels = c("Equation for Inflation", "Equation for Federal Funds Rate"),
          dep.var.labels = c("Inflation (YoY %)", "Federal Funds Rate"),
          title = paste0("Bivariate VAR(", bic_lag_bivar, ") Model Estimation Results"),
          align = TRUE,
          no.space = TRUE,
          notes.align = "l",
          font.size = "small")

# --- 1.4 Diagnostic Checks for Bivariate Model ---
var_residuals_bivar <- residuals(var_model_bivariate)

# ACF and PACF of residuals
par(mfrow = c(2, 2), oma = c(0,0,3,0))
acf(var_residuals_bivar[, "Inflation"], main = "ACF: Bivariate Inflation Residuals")
pacf(var_residuals_bivar[, "Inflation"], main = "PACF: Bivariate Inflation Residuals")
acf(var_residuals_bivar[, "Interest"], main = "ACF: Bivariate Interest Residuals")
pacf(var_residuals_bivar[, "Interest"], main = "PACF: Bivariate Interest Residuals")
mtext("Bivariate VAR: ACF and PACF of Residuals", side = 3, line = 1, outer = TRUE, cex = 1.2)
par(mfrow = c(1, 1), oma = c(0,0,0,0)) # Reset plotting parameters

# Check cross-correlation of residuals
ccf(var_residuals_bivar[, "Inflation"], var_residuals_bivar[, "Interest"],
    main = "Bivariate VAR: CCF - Inflation vs. Interest Residuals",
    sub = "Check spike at lag 0 for significance")

# Portmanteau Test (Ljung-Box) for residual autocorrelation
serial.test(var_model_bivariate, lags.pt = 24 , type = "PT.asymptotic")


# --- 1.5 Granger Causality Tests for Bivariate Model ---
granger_inf_to_int_bivar <- causality(var_model_bivariate, cause = "Inflation")
print(granger_inf_to_int_bivar$Granger)

granger_int_to_inf_bivar <- causality(var_model_bivariate, cause = "Interest")
print(granger_int_to_inf_bivar$Granger)

# --- 1.6 Impulse Response Analysis for Bivariate Model ---
horizon_bivar <- 48 # 48 months = 4 year horizon

# IRF: Response of Inflation to Interest Rate Shock (Bivariate)
irf_int_to_inf_bivar <- irf(var_model_bivariate, impulse = "Interest", response = "Inflation",
                            n.ahead = horizon_bivar, boot = TRUE, seed = 12345, runs = 1000, ci = 0.95)
plot(irf_int_to_inf_bivar, main = "Response of Inflation to FFR Shock", ylab = "Response", xlab = "Horizon (Months)")


# --- 1.7 Conditional Forecasting for Bivariate Model ---
# Use data until Feb 2022 for training for bivariate model
train_data_bivar <- window(combined_data_bivar, end = c(2022,2))

# Estimate VAR model on training data for bivariate model (using BIC-selected lag)
var_train_bivar <- VAR(train_data_bivar, p = bic_lag_bivar, type = "const")

# Create scenario matrix (fixed interest rate at 0.08% from March 2022 to Feb 2024)
forecast_start_date_bivar <- c(2022, 3)
forecast_end_date_bivar <- c(2024, 2)
num_forecast_steps_bivar <- 24

# Generate a time series object for the scenario period to get correct dates/freq
scenario_ts_dates_bivar <- ts(NA, start = forecast_start_date_bivar, end = forecast_end_date_bivar, frequency = 12)
scenario_matrix_bivar <- matrix(NA, nrow = length(scenario_ts_dates_bivar), ncol = ncol(combined_data_bivar))
colnames(scenario_matrix_bivar) <- colnames(combined_data_bivar)
rownames(scenario_matrix_bivar) <- time(scenario_ts_dates_bivar)

# Fix Interest Rate at 0.08% for the scenario period
scenario_matrix_bivar[,"Interest"] <- 0.08

# Inflation is the variable to be predicted in this scenario, so it remains NA
scenario_matrix_bivar[,"Inflation"] <- NA

# Predict under the conditional scenario for bivariate model
forecast_cond_bivar <- predict(var_train_bivar, n.ahead = num_forecast_steps_bivar, dumvar = scenario_matrix_bivar)
pred_inflation_cond_bivar <- forecast_cond_bivar$fcst$Inflation[, "fcst"]

# Get actual inflation for comparison
actual_inflation_counterfactual_bivar <- window(combined_data_bivar[, "Inflation"], start = forecast_start_date_bivar, end = forecast_end_date_bivar)

# Plotting the counterfactual using ggplot2
plot_df_counterfactual_bivar <- data.frame(
  Date = as.Date(time(actual_inflation_counterfactual_bivar)),
  Predicted = as.numeric(pred_inflation_cond_bivar),
  Actual = as.numeric(actual_inflation_counterfactual_bivar)
)

ggplot(data = plot_df_counterfactual_bivar, aes(x = Date)) +
  geom_line(aes(y = Predicted, color = "Predicted (Counterfactual)"), size = 1) +
  geom_line(aes(y = Actual, color = "Actual Inflation"), size = 1) +
  labs(title = "Actual vs. Predicted Inflation",
       x = "Date", y = "Inflation (YoY %)") +
  scale_color_manual(name = "Series", values = c("Predicted (Counterfactual)" = "red", "Actual Inflation" = "blue")) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

# --- SECTION 2: MULTIPLE VAR MODEL (Comprehensive Macroeconomic System) ---

# --- 2.1 Data Collection for Multiple VAR Model ---
unemployment_data <- fredr(series_id = "UNRATE", observation_start = start_date, observation_end = end_date)
gdp_data <- fredr(series_id = "GDPC1", observation_start = start_date, observation_end = end_date)
commodity_data <- fredr(series_id = "WTISPLC", observation_start = start_date, observation_end = end_date)
money_data <- fredr(series_id = "M2REAL", observation_start = start_date, observation_end = end_date)

# --- 2.2 Create Time Series Objects and Transformations for Multiple VAR Model ---
unemployment_ts <- ts(unemployment_data$value, start = c(1968, 1), frequency = 12)
commodity_ts <- ts(commodity_data$value, start = c(1968, 1), frequency = 12)
money_ts <- ts(money_data$value, start = c(1968, 1), frequency = 12)

# Convert quarterly GDP to monthly (linear interpolation and na.approx for leading NAs)
gdp_ts <- ts(gdp_data$value, start = c(1968, 1), frequency = 4)
gdp_monthly <- approx(time(gdp_ts), gdp_ts, xout = time(inflation_ts))$y
gdp_monthly_ts <- ts(gdp_monthly, start = c(1968, 1), frequency = 12)
gdp_monthly_ts <- na.approx(gdp_monthly_ts)

# Test for stationarity for all multiple VAR variables (and apply differencing where needed)
adf_test_inf <- ur.df(diff_inflation, type = "drift", selectlags = "AIC")
summary(adf_test_inf)

adf_test_int <- ur.df(interest_ts, type = "drift", selectlags = "AIC")
summary(adf_test_int)

adf_test_emp <- ur.df(unemployment_ts, type = "drift", selectlags = "AIC")
summary(adf_test_emp)

adf_test_gdp <- ur.df(gdp_monthly_ts, type = "trend", selectlags = "AIC")
summary(adf_test_gdp)
gdp_monthly_diff <- diff(gdp_monthly_ts)
adf_test_gdp_diff <- ur.df(gdp_monthly_diff, type = "drift", selectlags = "AIC")
summary(adf_test_gdp_diff)

adf_test_com <- ur.df(commodity_ts, type = "trend", selectlags = "AIC")
summary(adf_test_com)

adf_test_m2 <- ur.df(money_ts, type = "trend", selectlags = "AIC")
summary(adf_test_m2)
money_diff <- diff(money_ts)
adf_test_m2_diff <- ur.df(money_diff, type = "drift", selectlags = "AIC")
summary(adf_test_m2_diff)


# --- 2.3 Combine Variables into Dataset for Multiple VAR Model ---
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

# --- 2.4 VAR Model Estimation for Multiple VAR Model (Lag Selection using BIC/SC only) ---
VAR_lag_selection_multi <- VARselect(model_data_multi, lag.max = 30, type = "const")

bic_lag_multi <- VAR_lag_selection_multi$selection["SC(n)"]

# Plot BIC/SC criteria for visual inspection (single plot)
par(mfrow = c(1, 1))
plot(VAR_lag_selection_multi$criteria[3,], main = "Multiple VAR: Optimal Lag Selection (SC/BIC)",
     ylab = "SC/BIC Value", xlab = "Lag Order (p)", type = "b", pch = 19, col = "darkgreen")
abline(v = bic_lag_multi, col = "red", lty = 2)

# Estimate Multiple VAR model with BIC-selected lag (p=4)
var_model_multi <- VAR(model_data_multi, p = bic_lag_multi, type = "const")

# Model summary
summary(var_model_multi)

# Display regression results
stargazer(var_model_multi$varresult[[1]], var_model_multi$varresult[[2]],
          type = "text",
          column.labels = colnames(model_data_multi)[1:2],
          dep.var.labels = colnames(model_data_multi)[1:2],
          title = paste0("Multiple VAR(", bic_lag_multi, ") Model Estimation Results (Partial)"),
          align = TRUE,
          no.space = TRUE,
          notes.align = "l",
          font.size = "small")


# --- 2.5 Diagnostic Checks for Multiple VAR Model ---
var_residuals_multi <- residuals(var_model_multi)

# Check Cross-Correlation Function (CCF) of VAR Residuals
residual_correlation_matrix_multi <- cor(var_residuals_multi)
print(round(residual_correlation_matrix_multi, 3))

# Plot CCF for all pairs (using your original loop for visualization)
k_multi <- ncol(var_residuals_multi)
par(mfrow = c(3, 5), oma = c(0,0,3,0))
for (i in 1:(k_multi - 1)) {
  for (j in (i + 1):k_multi) {
    ccf(var_residuals_multi[, i], var_residuals_multi[, j],
        main = paste0("CCF: ", colnames(var_residuals_multi)[i], " & ", colnames(var_residuals_multi)[j]),
        sub = "Dashed lines indicate significance bounds")
  }
}
mtext("Multiple VAR: Cross-Correlation Functions of Residuals", side = 3, line = 1, outer = TRUE, cex = 1.2)
par(mfrow = c(1, 1), oma = c(0,0,0,0)) # Reset plotting parameters


# Portmanteau Test (Ljung-Box) for residual autocorrelation (joint test)
serial_test_lags <- max(bic_lag_multi + 10, 2 * bic_lag_multi, 24)
serial.test(var_model_multi, lags.pt = serial_test_lags, type = "PT.asymptotic")


# --- 2.6 Granger Causality Tests for Multiple VAR Model ---
granger_inf_to_int_multi <- causality(var_model_multi, cause = "Inflation")
print(granger_inf_to_int_multi$Granger)

granger_int_to_inf_multi <- causality(var_model_multi, cause = "Interest")
print(granger_int_to_inf_multi$Granger)

granger_com_to_inf_multi <- causality(var_model_multi, cause = "Commodity_Prices")
print(granger_com_to_inf_multi$Granger)


# --- 2.7 Impulse Response Analysis ---
horizon <- 48


# IRF: Response of Inflation to Interest Rate Shock
irf_int_to_inf <- irf(var_model_multi, impulse = "Interest", response = "Inflation",
                      n.ahead = horizon, boot = TRUE, seed = 12345, runs = 1000, ci = 0.95)
plot(irf_int_to_inf, main = "Response of Inflation to FFR Shock", ylab = "Response", xlab = "Horizon (Months)")


# --- 2.8 Conditional Forecasting for Multiple VAR Model ---
# Use data until Feb 2022 for training
train_data <- window(model_data_multi, end = c(2022,2))

# Estimate VAR model on training data
var_train <- VAR(train_data, p = bic_lag_multi, type = "const")

# Create scenario matrix (fixed interest rate)
forecast_start_date <- c(2022, 3)
forecast_end_date <- c(2024, 2)
num_forecast_steps <- 24

# Generate a time series object for the scenario period to get correct dates/freq
scenario_ts_dates <- ts(NA, start = forecast_start_date, end = forecast_end_date, frequency = 12)
scenario_matrix <- matrix(NA, nrow = length(scenario_ts_dates), ncol = ncol(model_data_multi))
colnames(scenario_matrix) <- colnames(model_data_multi)
rownames(scenario_matrix) <- time(scenario_ts_dates)

# Fix Interest Rate at 0.08% for the scenario period
scenario_matrix[,"Interest"] <- 0.08

# Populate actual values for other variables that are NOT being forecasted conditional on
for (col_name in colnames(scenario_matrix)) {
  if (col_name != "Inflation" && col_name != "Interest") {
    scenario_matrix[, col_name] <- window(model_data_multi[, col_name], start = forecast_start_date, end = forecast_end_date)
  }
}

# Predict under the conditional scenario
forecast_cond_multi <- predict(var_train, n.ahead = num_forecast_steps, dumvar = scenario_matrix)
pred_inflation_cond_multi <- forecast_cond_multi$fcst$Inflation[, "fcst"]

# Get actual inflation for comparison
actual_inflation_counterfactual <- window(model_data_multi[, "Inflation"], start = c(2022, 3), end = c(2024, 2))

# Plotting the counterfactual using ggplot2
plot_df_counterfactual <- data.frame(
  Date = as.Date(time(actual_inflation_counterfactual)),
  Predicted = as.numeric(pred_inflation_cond_multi),
  Actual = as.numeric(actual_inflation_counterfactual)
)

ggplot(data = plot_df_counterfactual, aes(x = Date)) +
  geom_line(aes(y = Predicted, color = "Predicted (Counterfactual)"), size = 1) +
  geom_line(aes(y = Actual, color = "Actual Inflation"), size = 1) +
  labs(title = "Actual vs. Predicted Inflation Under Counterfactual Interest Rate",
       x = "Date", y = "Inflation (YoY %)") +
  scale_color_manual(name = "Series", values = c("Predicted (Counterfactual)" = "red", "Actual Inflation" = "blue")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m")


# --- SECTION 3: SENSITIVITY ANALYSIS FOR MULTIPLE VAR MODEL ---

# Define alternative lag lengths
lag_minus_1 <- max(1, bic_lag_multi - 1)
lag_plus_1 <- bic_lag_multi + 1

# --- Scenario 1.1: VAR(p = bic_lag_multi) (Original Model) ---
cat(paste("\n### Scenario 1.1: VAR(p = ", bic_lag_multi, ") - Original Model Results ###\n"))
cat("Portmanteau Test for Residual Autocorrelation:\n")
serial_test_lags_current <- max(bic_lag_multi + 10, 2 * bic_lag_multi, 24)
print(serial.test(var_model_multi, lags.pt = serial_test_lags_current, type = "PT.asymptotic"))
cat("\nGranger Causality (Interest -> Inflation):\n")
print(causality(var_model_multi, cause = "Interest")$Granger)


# --- Scenario 1.2: VAR(p = lag_minus_1) ---

# Estimate the alternative VAR model with p-1 lags
var_model_minus_1 <- VAR(model_data_multi, p = lag_minus_1, type = "const")

cat("Portmanteau Test for Residual Autocorrelation:\n")
serial_test_lags_minus_1 <- max(lag_minus_1 + 10, 2 * lag_minus_1, 24)
print(serial.test(var_model_minus_1, lags.pt = serial_test_lags_minus_1, type = "PT.asymptotic"))

cat("\nGranger Causality (Interest -> Inflation):\n")
print(causality(var_model_minus_1, cause = "Interest")$Granger)


# --- Scenario 1.3: VAR(p = lag_plus_1) ---
cat(paste("\n### Scenario 1.3: VAR(p = ", lag_plus_1, ") Results ###\n"))

# Estimate the alternative VAR model with p+1 lags
var_model_plus_1 <- VAR(model_data_multi, p = lag_plus_1, type = "const")

cat("Portmanteau Test for Residual Autocorrelation:\n")
serial_test_lags_plus_1 <- max(lag_plus_1 + 10, 2 * lag_plus_1, 24)
print(serial.test(var_model_plus_1, lags.pt = serial_test_lags_plus_1, type = "PT.asymptotic"))

cat("\nGranger Causality (Interest -> Inflation):\n")
print(causality(var_model_plus_1, cause = "Interest")$Granger)

# --- 3.2 Sensitivity to Variable Ordering (for Impulse Response Functions) ---

alternative_order_colnames <- c(
  "Interest",
  "Inflation",
  "Commodity_Prices",
  "GDP_Growth",
  "Unemployment",
  "Money_Supply_Growth"
)

# Convert model_data_multi to a data frame for reliable column reordering by name
model_data_multi_df <- as.data.frame(model_data_multi)
# Reorder columns by name
model_data_multi_alt_order_df <- model_data_multi_df[, alternative_order_colnames]
# Convert back to a time-series matrix for the VAR function
model_data_multi_alt_order_ts <- ts(model_data_multi_alt_order_df,
                                    start = start(model_data_multi),
                                    frequency = frequency(model_data_multi))


# Estimate VAR model with alternative ordering and the same BIC-selected lag
var_model_multi_alt_order <- VAR(model_data_multi_alt_order_ts, p = bic_lag_multi, type = "const")

cat("\n### Scenario 2.1: IRFs - Original Ordering (Interest shock) ###\n")
par(mfrow = c(2, 1))
irf_int_to_inf_orig_order <- irf(var_model_multi, impulse = "Interest", response = "Inflation",
                                 n.ahead = horizon, boot = TRUE, seed = 12345, runs = 1000, ci = 0.95)
plot(irf_int_to_inf_orig_order, main = "IRF: Original Order (Interest Shock on Inflation)",
     ylab = "Response", xlab = "Horizon (Months)")


cat("\n### Scenario 2.2: IRFs - Alternative Ordering (Interest shock) ###\n")
irf_int_to_inf_alt_order <- irf(var_model_multi_alt_order, impulse = "Interest", response = "Inflation",
                                n.ahead = horizon, boot = TRUE, seed = 12345, runs = 1000, ci = 0.95)
plot(irf_int_to_inf_alt_order, main = "IRF: Alternative Order (Interest Shock on Inflation)",
     ylab = "Response", xlab = "Horizon (Months)")
