# Load required packages
library(fredr)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(lmtest)      # For coeftest, dwtest, bptest functions
library(sandwich)    # For robust standard errors (NeweyWest)
library(tidyr)       # For pivot_wider function
library(tseries)     # For adf.test

# Set FRED API key
# Replace with your actual FRED API key
fredr_set_key("68930503314ecaad758dc1a52090ca2f")

# --- Data Download and Initial Preparation ---
start_data_cpi_rates <- as.Date("2015-01-01")
end_data_cpi_rates <- as.Date("2019-12-31")

us_cpi <- fredr(series_id = "CPIAUCSL", observation_start = start_data_cpi_rates)
german_cpi <- fredr(series_id = "DEUCPALTT01IXOBSAM", observation_start = start_data_cpi_rates)
fed_rate <- fredr(series_id = "FEDFUNDS", observation_start = start_data_cpi_rates)
euro_rate <- fredr(series_id = "ECBDFR", observation_start = start_data_cpi_rates, frequency = "m")

# Convert to time series format
us_cpi_ts <- ts(us_cpi$value, start = c(year(min(us_cpi$date)), month(min(us_cpi$date))), frequency = 12)
german_cpi_ts <- ts(german_cpi$value, start = c(year(min(german_cpi$date)), month(min(german_cpi$date))), frequency = 12)
fed_rate_ts <- ts(fed_rate$value, start = c(year(min(fed_rate$date)), month(min(fed_rate$date))), frequency = 12)
euro_rate_ts <- ts(euro_rate$value, start = c(year(min(euro_rate$date)), month(min(euro_rate$date))), frequency = 12)

# Compute Year-over-Year Inflation Rate
us_inflation <- 100 * diff(log(us_cpi_ts), lag = 12)
german_inflation <- 100 * diff(log(german_cpi_ts), lag = 12)

# Combine all series into a dataframe for consistent date ranges
all_data_ts_aligned <- ts.union(us_inflation, german_inflation, fed_rate_ts, euro_rate_ts)
all_data_df <- data.frame(
  date = as.Date(time(zoo::as.zoo(all_data_ts_aligned))),
  us_inflation = coredata(us_inflation)[match(time(all_data_ts_aligned), time(us_inflation))],
  german_inflation = coredata(german_inflation)[match(time(all_data_ts_aligned), time(german_inflation))],
  fed_rate = coredata(fed_rate_ts)[match(time(all_data_ts_aligned), time(fed_rate_ts))],
  euro_rate = coredata(euro_rate_ts)[match(time(all_data_ts_aligned), time(euro_rate_ts))]
) %>% na.omit()


# --- ANALYSIS: US (Treated) vs. Germany (Control) ---

# Define analysis period and event date
start_analysis_date <- as.Date("2016-01-01")
end_analysis_date <- as.Date("2018-12-31")
event_date <- as.Date("2017-01-01") # Start of Fed rate hike cycle

# Filter data for the analysis period and create DiD variables
filtered_data <- all_data_df %>%
  filter(date >= start_analysis_date & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= event_date, 1, 0),
    period = case_when(
      date < event_date ~ "pre",
      date < as.Date("2018-01-01") ~ "period1",
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
    region = "Germany (Control)",
    treatment = 0,
    inflation_rate = german_inflation
  )
) %>% na.omit()


# --- PART 1: Parallel Trends Assumption Check ---

# Calculate optimal lag for Newey-West HAC standard errors
T_obs <- nrow(panel_data)
L_nw <- floor(4 * (T_obs / 100)^(2/9))

# Visual Check for Parallel Trends in pre-period
p1_parallel_trends <- ggplot(filtered_data %>% filter(date < event_date), aes(x = date)) +
  geom_line(aes(y = us_inflation, color = "US Inflation")) +
  geom_line(aes(y = german_inflation, color = "German Inflation")) +
  geom_smooth(aes(y = us_inflation, color = "US Inflation"), method = "lm", se = FALSE, linetype = "dashed") +
  geom_smooth(aes(y = german_inflation, color = "German Inflation"), method = "lm", se = FALSE, linetype = "dashed") +
  labs(title = "Pre-Treatment Trends: US vs. German Inflation (2016)",
       subtitle = paste0("Data points in pre-period: ", nrow(filtered_data %>% filter(date < event_date))),
       x = "Date", y = "Inflation Rate (%)",
       color = "Series") +
  theme_minimal()
print(p1_parallel_trends)

# Statistical Check for Parallel Trends: Test for differential trend in pre-period
pre_treatment_data <- filtered_data %>% filter(date < event_date)
pre_diff_model <- lm(us_inflation - german_inflation ~ date, data = pre_treatment_data)
cat("--- Statistical Test for Differential Pre-Trend ---\n")
print(summary(pre_diff_model))


# --- PART 2: Standard DiD Analysis (Canonical Form) ---
# Estimates one average treatment effect for the entire post-period.

standard_did_canonical <- lm(inflation_rate ~ treatment * post, data = panel_data)

# Applying HAC standard errors using calculated L_nw
robust_se_standard_did <- coeftest(standard_did_canonical, vcov = NeweyWest(standard_did_canonical, lag = L_nw, prewhite = FALSE))
cat("\n==== STANDARD DiD MODEL (US vs. GERMANY) ====\n")
print(robust_se_standard_did)

# Diagnostic Checks for Standard DiD Model
cat("\n--- Diagnostic Checks for Standard DiD Model ---\n")
cat("Durbin-Watson Test for Autocorrelation:\n")
print(dwtest(standard_did_canonical))
cat("Breusch-Pagan Test for Heteroskedasticity:\n")
print(bptest(standard_did_canonical))


# --- PART 3: Sequential DiD Analysis (Event Study) ---
panel_data_seq_explicit <- panel_data %>%
  mutate(
    treat_period1 = ifelse(period == "period1" & treatment == 1, 1, 0),
    treat_period2 = ifelse(period == "period2" & treatment == 1, 1, 0)
  )
seq_did_model <- lm(inflation_rate ~ treatment + factor(period_num) + treat_period1 + treat_period2, data = panel_data_seq_explicit)
robust_se_seq <- coeftest(seq_did_model, vcov = NeweyWest(seq_did_model, lag = L_nw, prewhite = FALSE))
print(robust_se_seq)

# Diagnostic Checks for Sequential DiD Model
cat("\n--- Diagnostic Checks for Sequential DiD Model ---\n")
cat("Durbin-Watson Test for Autocorrelation:\n")
print(dwtest(seq_did_model))
cat("Breusch-Pagan Test for Heteroskedasticity:\n")
print(bptest(seq_did_model))


# --- PART 4: Visualization ---
# Filter first year of data for regression (e.g., 2016)
regression_data <- filtered_data[filtered_data$date < as.Date("2017-01-01"), ]

# Plot
p_combined_did <- ggplot() +
  # Shaded periods
  geom_rect(aes(xmin = as.Date("2017-01-01"), xmax = as.Date("2017-12-31"), ymin = -Inf, ymax = Inf),
            fill = "lightgray", alpha = 0.3) +
  geom_rect(aes(xmin = as.Date("2018-01-01"), xmax = max(filtered_data$date), ymin = -Inf, ymax = Inf),
            fill = "darkgray", alpha = 0.3) +
  
  # Regression lines for pre-treatment period
  geom_smooth(data = regression_data, aes(x = date, y = us_inflation),
              method = "lm", se = FALSE, color = "#FF9999", linetype = "dotted", size = 1.2, show.legend = FALSE) +
  geom_smooth(data = regression_data, aes(x = date, y = german_inflation),
              method = "lm", se = FALSE, color = "#9999FF", linetype = "dotted", size = 1.2, show.legend = FALSE) +
  
  # Main inflation lines
  geom_line(data = filtered_data, aes(x = date, y = us_inflation, color = "US Inflation", linetype = "US Inflation"), size = 1.2) +
  geom_line(data = filtered_data, aes(x = date, y = german_inflation, color = "German Inflation", linetype = "German Inflation"), size = 1.2) +
  
  # Vertical period markers
  geom_vline(xintercept = as.Date("2017-01-01"), linetype = "dashed", color = "black") +
  geom_vline(xintercept = as.Date("2018-01-01"), linetype = "dotted", color = "black") +
  
  # Period labels
  annotate("text", x = as.Date("2017-04-01"), y = 0.5,
           label = "Period 1", hjust = 0, size = 6) +
  annotate("text", x = as.Date("2018-04-01"), y = 0.5,
           label = "Period 2", hjust = 0, size = 6) +
  
  # Color/linetype/legend settings
  scale_color_manual(values = c("US Inflation" = "red", "German Inflation" = "blue")) +
  scale_linetype_manual(values = c("US Inflation" = "solid", "German Inflation" = "dashed")) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  
  # Theme
  labs(x = "Date", y = "Inflation Rate (%)", color = "Series", linetype = "Series") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_combined_did)


# --- PART 5: Basic Sensitivity Analysis ---

# This section explores the robustness of the Standard DiD model's results
# to minor changes in the definition of the pre-treatment period or event date.

cat("\n\n==== PART 5: BASIC SENSITIVITY ANALYSIS ====\n")

# --- Scenario 1: Shift pre-treatment period start by -1 month ---
# New pre-period starts one month earlier (e.g., 2015-12-01 instead of 2016-01-01)
cat("\n--- Scenario 1: Pre-period start shifted -1 month ---\n")
filtered_data_s1 <- all_data_df %>%
  filter(date >= as.Date("2015-12-01") & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= event_date, 1, 0)
  )
panel_data_s1 <- bind_rows(
  filtered_data_s1 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s1 %>% mutate(region = "Germany (Control)", treatment = 0, inflation_rate = german_inflation)
) %>% na.omit()
s1_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s1)
s1_robust_se <- coeftest(s1_model, vcov = NeweyWest(s1_model, lag = L_nw, prewhite = FALSE))
print(s1_robust_se)


# --- Scenario 2: Shift pre-treatment period start by +1 month ---
# New pre-period starts one month later (e.g., 2016-02-01 instead of 2016-01-01)
cat("\n--- Scenario 2: Pre-period start shifted +1 month ---\n")
filtered_data_s2 <- all_data_df %>%
  filter(date >= as.Date("2016-02-01") & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= event_date, 1, 0)
  )
panel_data_s2 <- bind_rows(
  filtered_data_s2 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s2 %>% mutate(region = "Germany (Control)", treatment = 0, inflation_rate = german_inflation)
) %>% na.omit()
s2_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s2)
s2_robust_se <- coeftest(s2_model, vcov = NeweyWest(s2_model, lag = L_nw, prewhite = FALSE))
print(s2_robust_se)


# --- Scenario 3: Shift event date by -1 month ---
# Treatment starts one month earlier (e.g., 2016-12-01 instead of 2017-01-01)
cat("\n--- Scenario 3: Event date shifted -1 month ---\n")
filtered_data_s3 <- all_data_df %>%
  filter(date >= start_analysis_date & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= as.Date("2016-12-01"), 1, 0) # Shifted event date
  )
panel_data_s3 <- bind_rows(
  filtered_data_s3 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s3 %>% mutate(region = "Germany (Control)", treatment = 0, inflation_rate = german_inflation)
) %>% na.omit()
s3_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s3)
s3_robust_se <- coeftest(s3_model, vcov = NeweyWest(s3_model, lag = L_nw, prewhite = FALSE))
print(s3_robust_se)


# --- Scenario 4: Shift event date by +1 month ---
# Treatment starts one month later (e.g., 2017-02-01 instead of 2017-01-01)
cat("\n--- Scenario 4: Event date shifted +1 month ---\n")
filtered_data_s4 <- all_data_df %>%
  filter(date >= start_analysis_date & date <= end_analysis_date) %>%
  mutate(
    post = ifelse(date >= as.Date("2017-02-01"), 1, 0) # Shifted event date
  )
panel_data_s4 <- bind_rows(
  filtered_data_s4 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s4 %>% mutate(region = "Germany (Control)", treatment = 0, inflation_rate = german_inflation)
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
