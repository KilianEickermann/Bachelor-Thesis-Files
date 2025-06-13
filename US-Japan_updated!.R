# Load required packages
library(fredr)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(lmtest)    
library(sandwich)   
library(tidyr)    
library(tseries)

# Set FRED API key
fredr_set_key("68930503314ecaad758dc1a52090ca2f")

# --- Data Download and Initial Preparation ---
start_date <- as.Date("2000-01-01")
end_date <- as.Date("2021-06-01")

# Retrieve CPI and Interest Rate Data from FRED
us_cpi <- fredr(series_id = "CPIAUCSL", observation_start = start_date)
japan_cpi <- fredr(series_id = "CPALTT01JPM661S", observation_start = start_date)
fed_rate <- fredr(series_id = "FEDFUNDS", observation_start = start_date)

# Load Japan interest rate data
japan_rate <- read.csv("Interest_japan.csv")
japan_rate <- japan_rate[2:nrow(japan_rate), 1:2]
colnames(japan_rate) <- c("date", "interest_rate")

# Convert 'date' column to Date format
japan_rate$date <- as.Date(paste0(japan_rate$date, "/01"), format = "%Y/%m/%d")

# Remove row names (index)
rownames(japan_rate) <- NULL

# Convert interest_rate to numeric explicitly
japan_rate$interest_rate <- as.numeric(as.character(japan_rate$interest_rate))

# Convert to time series format
us_cpi_ts <- ts(us_cpi$value, start = c(year(min(us_cpi$date)), month(min(us_cpi$date))), frequency = 12)
japan_cpi_ts <- ts(japan_cpi$value, start = c(year(min(japan_cpi$date)), month(min(japan_cpi$date))), frequency = 12)
fed_rate_ts <- ts(fed_rate$value, start = c(year(min(fed_rate$date)), month(min(fed_rate$date))), frequency = 12)
japan_rate_ts <- ts(japan_rate$interest_rate, start = c(year(min(japan_rate$date)), month(min(japan_rate$date))), frequency = 12)

# Compute Year-over-Year Inflation Rate
us_inflation <- 100 * diff(log(us_cpi_ts), lag = 12)
japan_inflation <- 100 * diff(log(japan_cpi_ts), lag = 12)

# Combine all series into a dataframe for consistent date ranges
all_data_ts_aligned <- ts.union(us_inflation, japan_inflation, fed_rate_ts, japan_rate_ts)
all_data_df <- data.frame(
  date = as.Date(time(zoo::as.zoo(all_data_ts_aligned))),
  us_inflation = coredata(us_inflation)[match(time(all_data_ts_aligned), time(us_inflation))],
  japan_inflation = coredata(japan_inflation)[match(time(all_data_ts_aligned), time(japan_inflation))],
  fed_rate = coredata(fed_rate_ts)[match(time(all_data_ts_aligned), time(fed_rate_ts))],
  japan_rate = coredata(japan_rate_ts)[match(time(all_data_ts_aligned), time(japan_rate_ts))]
) %>% na.omit()


# --- ANALYSIS: US (Treated) vs. Japan (Control) ---

# Define analysis period and event date
start_analysis_date <- as.Date("2003-05-01")
end_analysis_date <- as.Date("2006-05-31")
event_date <- as.Date("2004-05-01") # Start of Fed rate hike cycle

# Filter data for the analysis period and create DiD variables
filtered_data <- all_data_df %>%
  filter(date >= start_analysis_date & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= event_date, 1, 0),
    period = case_when(
      date < event_date ~ "pre",
      date < as.Date("2005-05-01") ~ "period1",
      TRUE ~ "period2"
    ),
    period_num = case_when(
      period == "pre" ~ 0,
      period == "period1" ~ 1,
      period == "period2" ~ 2
    )
  )

# Create panel dataset for DiD analysis
panel_data <- bind_rows(
  filtered_data %>% mutate(
    region = "US (Treated)",
    treatment = 1,
    inflation_rate = us_inflation
  ),
  filtered_data %>% mutate(
    region = "Japan (Control)",
    treatment = 0,
    inflation_rate = japan_inflation
  )
) %>% na.omit()


# --- PART 1: Parallel Trends Assumption Check ---

# Calculate optimal lag for Newey-West HAC standard errors
T_obs <- nrow(panel_data)
L_nw <- floor(4 * (T_obs / 100)^(2/9))

# Visual Check for Parallel Trends in pre-period
p1_parallel_trends <- ggplot(filtered_data %>% filter(date < event_date), aes(x = date)) +
  geom_line(aes(y = us_inflation, color = "US Inflation")) +
  geom_line(aes(y = japan_inflation, color = "Japan Inflation")) +
  geom_smooth(aes(y = us_inflation, color = "US Inflation"), method = "lm", se = FALSE, linetype = "dashed") +
  geom_smooth(aes(y = japan_inflation, color = "Japan Inflation"), method = "lm", se = FALSE, linetype = "dashed") +
  labs(title = "Pre-Treatment Trends: US vs. Japan Inflation (2003-2004)",
       subtitle = paste0("Data points in pre-period: ", nrow(filtered_data %>% filter(date < event_date))),
       x = "Date", y = "Inflation Rate (%)",
       color = "Series") +
  theme_minimal()
print(p1_parallel_trends)

# Statistical Check for Parallel Trends: Test for differential trend in pre-period
pre_treatment_data <- filtered_data %>% filter(date < event_date)
pre_diff_model <- lm(us_inflation - japan_inflation ~ date, data = pre_treatment_data)
cat("--- Statistical Test for Differential Pre-Trend ---\n")
print(summary(pre_diff_model))


# --- PART 2: Standard DiD Analysis (Canonical Form) ---
# Estimates one average treatment effect for the entire post-period.

standard_did_canonical <- lm(inflation_rate ~ treatment * post, data = panel_data)

# Applying HAC standard errors using calculated L_nw
robust_se_standard_did <- coeftest(standard_did_canonical, vcov = NeweyWest(standard_did_canonical, lag = L_nw, prewhite = FALSE))
cat("\n==== STANDARD DiD MODEL (US vs. JAPAN) ====\n")
print(robust_se_standard_did)

# Diagnostic Checks for Standard DiD Model
cat("\n--- Diagnostic Checks for Standard DiD Model ---\n")
cat("Durbin-Watson Test for Autocorrelation:\n")
print(dwtest(standard_did_canonical))
cat("Breusch-Pagan Test for Heteroskedasticity:\n")
print(bptest(standard_did_canonical))


# --- PART 3: Sequential DiD Analysis (Event Study) ---
# Provides insights into the *dynamics* of the treatment effect across specific post-treatment periods.

# Model 1: Using factor(period_num) directly for period-specific interaction terms
seq_did_model1 <- lm(inflation_rate ~ treatment * factor(period_num), data = panel_data)
robust_se_seq1 <- coeftest(seq_did_model1, vcov = NeweyWest(seq_did_model1, lag = L_nw, prewhite = FALSE))
cat("\n==== SEQUENTIAL DiD MODEL (US vs. JAPAN) - Model 1 ====\n")
print(robust_se_seq1)

# Model 2: Explicitly adding interaction terms for each period for clarity
panel_data_seq_explicit <- panel_data %>%
  mutate(
    treat_period1 = ifelse(period == "period1" & treatment == 1, 1, 0),
    treat_period2 = ifelse(period == "period2" & treatment == 1, 1, 0)
  )
seq_did_model2 <- lm(inflation_rate ~ treatment + factor(period_num) + treat_period1 + treat_period2, data = panel_data_seq_explicit)
robust_se_seq2 <- coeftest(seq_did_model2, vcov = NeweyWest(seq_did_model2, lag = L_nw, prewhite = FALSE))
cat("\n==== SEQUENTIAL DiD MODEL (US vs. JAPAN) - Model 2 ====\n")
print(robust_se_seq2)

# Diagnostic Checks for Sequential DiD Model (using Model 2 for example)
cat("\n--- Diagnostic Checks for Sequential DiD Model ---\n")
cat("Durbin-Watson Test for Autocorrelation:\n")
print(dwtest(seq_did_model2))
cat("Breusch-Pagan Test for Heteroskedasticity:\n")
print(bptest(seq_did_model2))


# --- PART 4: Visualization ---
# Combined visualization of US and Japan inflation trends with period markers.
p_combined_did <- ggplot() +
  geom_rect(aes(xmin = as.Date("2004-05-01"), xmax = as.Date("2005-04-30"), ymin = -Inf, ymax = Inf), fill = "lightgray", alpha = 0.3) +
  geom_rect(aes(xmin = as.Date("2005-05-01"), xmax = max(filtered_data$date), ymin = -Inf, ymax = Inf), fill = "darkgray", alpha = 0.3) +
  geom_line(data = filtered_data, aes(x = date, y = us_inflation, color = "US Inflation")) +
  geom_line(data = filtered_data, aes(x = date, y = japan_inflation, color = "Japan Inflation")) +
  geom_vline(xintercept = as.Date("2004-05-01"), linetype = "dashed", color = "red") +
  geom_vline(xintercept = as.Date("2005-05-01"), linetype = "dotted", color = "blue") +
  annotate("text", x = as.Date("2004-06-01"), y = max(filtered_data$us_inflation, na.rm = TRUE) - 0.5, label = "Period 1 Begins", hjust = 0, size = 3) +
  annotate("text", x = as.Date("2005-06-01"), y = max(filtered_data$us_inflation, na.rm = TRUE) - 0.5, label = "Period 2 Begins", hjust = 0, size = 3) +
  labs(title = "US-Japan DiD Analysis: Inflation Trends",
       subtitle = "Effect of Fed rate hike on inflation rates",
       x = "Date", y = "Inflation Rate (%)",
       color = "Series") +
  theme_minimal()
print(p_combined_did)


# --- PART 5: Basic Sensitivity Analysis ---

# This section explores the robustness of the Standard DiD model's results
# to minor changes in the definition of the pre-treatment period or event date.

cat("\n\n==== PART 5: BASIC SENSITIVITY ANALYSIS ====\n")

# --- Scenario 1: Shift pre-treatment period start by -1 month ---
# New pre-period starts one month earlier (e.g., 2003-04-01 instead of 2003-05-01)
cat("\n--- Scenario 1: Pre-period start shifted -1 month ---\n")
filtered_data_s1 <- all_data_df %>%
  filter(date >= as.Date("2003-04-01") & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= event_date, 1, 0)
  )
panel_data_s1 <- bind_rows(
  filtered_data_s1 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s1 %>% mutate(region = "Japan (Control)", treatment = 0, inflation_rate = japan_inflation)
) %>% na.omit()
s1_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s1)
s1_robust_se <- coeftest(s1_model, vcov = NeweyWest(s1_model, lag = L_nw, prewhite = FALSE))
print(s1_robust_se)


# --- Scenario 2: Shift pre-treatment period start by +1 month ---
# New pre-period starts one month later (e.g., 2003-06-01 instead of 2003-05-01)
cat("\n--- Scenario 2: Pre-period start shifted +1 month ---\n")
filtered_data_s2 <- all_data_df %>%
  filter(date >= as.Date("2003-06-01") & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= event_date, 1, 0)
  )
panel_data_s2 <- bind_rows(
  filtered_data_s2 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s2 %>% mutate(region = "Japan (Control)", treatment = 0, inflation_rate = japan_inflation)
) %>% na.omit()
s2_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s2)
s2_robust_se <- coeftest(s2_model, vcov = NeweyWest(s2_model, lag = L_nw, prewhite = FALSE))
print(s2_robust_se)


# --- Scenario 3: Shift event date by -1 month ---
# Treatment starts one month earlier (e.g., 2004-04-01 instead of 2004-05-01)
cat("\n--- Scenario 3: Event date shifted -1 month ---\n")
filtered_data_s3 <- all_data_df %>%
  filter(date >= start_analysis_date & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= as.Date("2004-04-01"), 1, 0) # Shifted event date
  )
panel_data_s3 <- bind_rows(
  filtered_data_s3 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s3 %>% mutate(region = "Japan (Control)", treatment = 0, inflation_rate = japan_inflation)
) %>% na.omit()
s3_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s3)
s3_robust_se <- coeftest(s3_model, vcov = NeweyWest(s3_model, lag = L_nw, prewhite = FALSE))
print(s3_robust_se)


# --- Scenario 4: Shift event date by +1 month ---
# Treatment starts one month later (e.g., 2004-06-01 instead of 2004-05-01)
cat("\n--- Scenario 4: Event date shifted +1 month ---\n")
filtered_data_s4 <- all_data_df %>%
  filter(date >= start_analysis_date & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= as.Date("2004-06-01"), 1, 0) # Shifted event date
  )
panel_data_s4 <- bind_rows(
  filtered_data_s4 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s4 %>% mutate(region = "Japan (Control)", treatment = 0, inflation_rate = japan_inflation)
) %>% na.omit()
s4_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s4)
s4_robust_se <- coeftest(s4_model, vcov = NeweyWest(s4_model, lag = L_nw, prewhite = FALSE))
print(s4_robust_se)

# --- Table for Sensitivity Analysis Comparison ---

# Function to extract relevant results from coeftest output
extract_did_results <- function(model_robust_se, scenario_name) {
  # The 'treatment:post' coefficient is typically the 4th row in the coeftest output matrix
  did_row <- model_robust_se[row.names(model_robust_se) == "treatment:post", ]
  
  # Extract values
  estimate <- did_row[1]
  std_error <- did_row[2]
  p_value <- did_row[4]
  
  # Determine significance stars
  stars <- ""
  if (p_value < 0.001) {
    stars <- "***"
  } else if (p_value < 0.01) {
    stars <- "**"
  } else if (p_value < 0.05) {
    stars <- "*"
  } else if (p_value < 0.1) {
    stars <- "."
  }
  
  tibble(
    Scenario = scenario_name,
    Estimate = round(estimate, 4),
    `Std. Error` = round(std_error, 4),
    `p-value` = round(p_value, 4),
    Significance = stars
  )
}

# Extract results for each model
results_original <- extract_did_results(robust_se_standard_did, "Original Model")
results_s1 <- extract_did_results(s1_robust_se, "Scenario 1: Pre-period start -1 month")
results_s2 <- extract_did_results(s2_robust_se, "Scenario 2: Pre-period start +1 month")
results_s3 <- extract_did_results(s3_robust_se, "Scenario 3: Event date shifted -1 month")
results_s4 <- extract_did_results(s4_robust_se, "Scenario 4: Event date shifted +1 month")

# Combine all results into a single table
sensitivity_table <- bind_rows(
  results_original,
  results_s1,
  results_s2,
  results_s3,
  results_s4
)

cat("\n\n==== SENSITIVITY ANALYSIS COMPARISON TABLE ====\n")
print(sensitivity_table)
cat("\nNote: Significance codes: '***' 0.001, '**' 0.01, '*' 0.05, '.' 0.1\n")
