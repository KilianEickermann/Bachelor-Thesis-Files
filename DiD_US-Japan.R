
# --- Setup ---
# Load required R packages.
library(fredr)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(lmtest)
library(sandwich)
library(tidyr)
library(tseries)

# Set FRED API key for data access.
fredr_set_key("68930503314ecaad758dc1a52090ca2f")

# --- Data Preparation ---

# Define the data download period.
start_date <- as.Date("2000-01-01")
end_date <- as.Date("2021-06-01")

# Retrieve CPI and Federal Funds Rate data from FRED.
us_cpi <- fredr(series_id = "CPIAUCSL", observation_start = start_date)
japan_cpi <- fredr(series_id = "CPALTT01JPM661S", observation_start = start_date)
fed_rate <- fredr(series_id = "FEDFUNDS", observation_start = start_date)

# Load Japan interest rate data from CSV.
japan_rate <- read.csv("Interest_japan.csv")
japan_rate <- japan_rate[2:nrow(japan_rate), 1:2]
colnames(japan_rate) <- c("date", "interest_rate")

# Format Japan interest rate dates and values.
japan_rate$date <- as.Date(paste0(japan_rate$date, "/01"), format = "%Y/%m/%d")
rownames(japan_rate) <- NULL
japan_rate$interest_rate <- as.numeric(as.character(japan_rate$interest_rate))

# Convert data to time series objects.
us_cpi_ts <- ts(us_cpi$value, start = c(year(min(us_cpi$date)), month(min(us_cpi$date))), frequency = 12)
japan_cpi_ts <- ts(japan_cpi$value, start = c(year(min(japan_cpi$date)), month(min(japan_cpi$date))), frequency = 12)
fed_rate_ts <- ts(fed_rate$value, start = c(year(min(fed_rate$date)), month(min(fed_rate$date))), frequency = 12)
japan_rate_ts <- ts(japan_rate$interest_rate, start = c(year(min(japan_rate$date)), month(min(japan_rate$date))), frequency = 12)

# Compute Year-over-Year Inflation Rate.
us_inflation <- 100 * diff(log(us_cpi_ts), lag = 12)
japan_inflation <- 100 * diff(log(japan_cpi_ts), lag = 12)

# Combine and align all series into a single data frame.
all_data_ts_aligned <- ts.union(us_inflation, japan_inflation, fed_rate_ts, japan_rate_ts)
all_data_df <- data.frame(
  date = as.Date(time(zoo::as.zoo(all_data_ts_aligned))),
  us_inflation = coredata(us_inflation)[match(time(all_data_ts_aligned), time(us_inflation))],
  japan_inflation = coredata(japan_inflation)[match(time(all_data_ts_aligned), time(japan_inflation))],
  fed_rate = coredata(fed_rate_ts)[match(time(all_data_ts_aligned), time(fed_rate_ts))],
  japan_rate = coredata(japan_rate_ts)[match(time(all_data_ts_aligned), time(japan_rate_ts))]
) %>% na.omit()

# Define analysis period and event date for DiD.
start_analysis_date <- as.Date("2003-05-01")
end_analysis_date <- as.Date("2006-04-30")
event_date <- as.Date("2004-05-01")

# Filter data for the analysis period and create DiD variables.
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

# Create panel dataset for DiD analysis.
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

# --- Difference-in-Differences Analysis ---

# Calculate optimal lag for Newey-West HAC standard errors.
T_obs <- nrow(panel_data)
L_nw <- floor(4 * (T_obs / 100)^(2/9))

# Plot for Parallel Trends Assumption Check.
p1_parallel_trends <- ggplot(filtered_data %>% filter(date < event_date), aes(x = date)) +
  geom_line(aes(y = us_inflation, color = "US Inflation", linetype = "US Inflation"), size = 1.2) +
  geom_line(aes(y = japan_inflation, color = "Japan Inflation", linetype = "Japan Inflation"), size = 1.2) +
  geom_smooth(aes(y = us_inflation, color = "US Inflation"),
              method = "lm", se = FALSE, linetype = "dotted", size = 1.2, show.legend = FALSE) +
  geom_smooth(aes(y = japan_inflation, color = "Japan Inflation"),
              method = "lm", se = FALSE, color = "blue", linetype = "dotted", size = 1.2, show.legend = FALSE) +
  scale_color_manual(values = c("US Inflation" = "red", "Japan Inflation" = "blue")) +
  scale_linetype_manual(values = c("US Inflation" = "solid", "Japan Inflation" = "dashed")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m",
               limits = c(min(filtered_data$date[filtered_data$date < event_date]), max(filtered_data$date[filtered_data$date < event_date]))) +
  labs(x = "Date", y = "Inflation Rate (%)", color = "Series", linetype = "Series",
       title = "Parallel Trends Check: US vs. Japan Inflation (Pre-Treatment)",
       subtitle = "Solid lines: Actual values; Dotted lines: Linear trends") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5)
  )
print(p1_parallel_trends)

# Statistical Test for Parallel Trends (Pre-Trend Differential).
pre_treatment_data <- filtered_data %>% filter(date < event_date)
pre_diff_model <- lm(us_inflation - japan_inflation ~ date, data = pre_treatment_data)
cat("\n--- Statistical Test for Differential Pre-Trend ---\n")
print(summary(pre_diff_model))

# Standard DiD Analysis (Canonical Form).
standard_did_canonical <- lm(inflation_rate ~ treatment * post, data = panel_data)
robust_se_standard_did <- coeftest(standard_did_canonical, vcov = NeweyWest(standard_did_canonical, lag = L_nw, prewhite = FALSE))
cat("\n==== STANDARD DiD MODEL (US vs. JAPAN) ====\n")
print(robust_se_standard_did)

# Diagnostic Checks for Standard DiD Model.
cat("\n--- Diagnostic Checks for Standard DiD Model ---\n")
cat("Durbin-Watson Test for Autocorrelation:\n")
print(dwtest(standard_did_canonical))
cat("Breusch-Pagan Test for Heteroskedasticity:\n")
print(bptest(standard_did_canonical))

# Sequential DiD Analysis (Event Study).
panel_data_seq_explicit <- panel_data %>%
  mutate(
    treat_period1 = ifelse(period == "period1" & treatment == 1, 1, 0),
    treat_period2 = ifelse(period == "period2" & treatment == 1, 1, 0)
  )
seq_did_model <- lm(inflation_rate ~ treatment + factor(period_num) + treat_period1 + treat_period2, data = panel_data_seq_explicit)
robust_se_seq <- coeftest(seq_did_model, vcov = NeweyWest(seq_did_model, lag = L_nw, prewhite = FALSE))
cat("\n==== SEQUENTIAL DiD MODEL (Event Study - US vs. JAPAN) ====\n")
print(robust_se_seq)

# Diagnostic Checks for Sequential DiD Model.
cat("\n--- Diagnostic Checks for Sequential DiD Model ---\n")
cat("Durbin-Watson Test for Autocorrelation:\n")
print(dwtest(seq_did_model))
cat("Breusch-Pagan Test for Heteroskedasticity:\n")
print(bptest(seq_did_model))

# Combined DiD Visualization.
p_combined_did <- ggplot() +
  geom_rect(aes(xmin = as.Date("2004-05-01"), xmax = as.Date("2005-04-30"), ymin = -Inf, ymax = Inf),
            fill = "lightgray", alpha = 0.3) +
  geom_rect(aes(xmin = as.Date("2005-05-01"), xmax = max(filtered_data$date), ymin = -Inf, ymax = Inf),
            fill = "darkgray", alpha = 0.3) +
  geom_smooth(data = filtered_data %>% filter(date < as.Date("2004-05-01")), aes(x = date, y = us_inflation),
              method = "lm", se = FALSE, color = "#FF9999", linetype = "dotted", size = 1.2, show.legend = FALSE) +
  geom_smooth(data = filtered_data %>% filter(date < as.Date("2004-05-01")), aes(x = date, y = japan_inflation),
              method = "lm", se = FALSE, color = "#9999FF", linetype = "dotted", size = 1.2, show.legend = FALSE) +
  geom_line(data = filtered_data, aes(x = date, y = us_inflation, color = "US Inflation", linetype = "US Inflation"), size = 1.2) +
  geom_line(data = filtered_data, aes(x = date, y = japan_inflation, color = "Japan Inflation", linetype = "Japan Inflation"), size = 1.2) +
  geom_vline(xintercept = as.Date("2004-05-01"), linetype = "dashed", color = "black") +
  geom_vline(xintercept = as.Date("2005-05-01"), linetype = "dotted", color = "black") +
  annotate("text", x = as.Date("2004-08-01"), y = 1.5,
           label = "Period 1", hjust = 0, size = 6) +
  annotate("text", x = as.Date("2005-08-01"), y = 1.5,
           label = "Period 2", hjust = 0, size = 6) +
  scale_color_manual(values = c("US Inflation" = "red", "Japan Inflation" = "blue")) +
  scale_linetype_manual(values = c("US Inflation" = "solid", "Japan Inflation" = "dashed")) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  labs(x = "Date", y = "Inflation Rate (%)", color = "Series", linetype = "Series",
       title = "Inflation Trends: US (Treated) vs. Japan (Control)",
       subtitle = "Shaded regions indicate post-treatment periods") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5)
  )
print(p_combined_did)

# --- Sensitivity Analysis ---

cat("\n\n==== SENSITIVITY ANALYSIS ====\n")

# Scenario 1: Pre-period start shifted -1 month.
cat("\n--- Scenario 1: Pre-period start shifted -1 month ---\n")
filtered_data_s1 <- all_data_df %>%
  filter(date >= as.Date("2003-04-01") & date <= end_analysis_date) %>%
  mutate(post = ifelse(date >= event_date, 1, 0))
panel_data_s1 <- bind_rows(
  filtered_data_s1 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s1 %>% mutate(region = "Japan (Control)", treatment = 0, inflation_rate = japan_inflation)
) %>% na.omit()
s1_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s1)
s1_robust_se <- coeftest(s1_model, vcov = NeweyWest(s1_model, lag = L_nw, prewhite = FALSE))
print(s1_robust_se)

# Scenario 2: Pre-period start shifted +1 month.
cat("\n--- Scenario 2: Pre-period start shifted +1 month ---\n")
filtered_data_s2 <- all_data_df %>%
  filter(date >= as.Date("2003-06-01") & date <= end_analysis_date) %>%
  mutate(post = ifelse(date >= event_date, 1, 0))
panel_data_s2 <- bind_rows(
  filtered_data_s2 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s2 %>% mutate(region = "Japan (Control)", treatment = 0, inflation_rate = japan_inflation)
) %>% na.omit()
s2_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s2)
s2_robust_se <- coeftest(s2_model, vcov = NeweyWest(s2_model, lag = L_nw, prewhite = FALSE))
print(s2_robust_se)

# Scenario 3: Event date shifted -1 month.
cat("\n--- Scenario 3: Event date shifted -1 month ---\n")
filtered_data_s3 <- all_data_df %>%
  filter(date >= start_analysis_date & date <= end_analysis_date) %>%
  mutate(post = ifelse(date >= as.Date("2004-04-01"), 1, 0))
panel_data_s3 <- bind_rows(
  filtered_data_s3 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s3 %>% mutate(region = "Japan (Control)", treatment = 0, inflation_rate = japan_inflation)
) %>% na.omit()
s3_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s3)
s3_robust_se <- coeftest(s3_model, vcov = NeweyWest(s3_model, lag = L_nw, prewhite = FALSE))
print(s3_robust_se)

# Scenario 4: Event date shifted +1 month.
cat("\n--- Scenario 4: Event date shifted +1 month ---\n")
filtered_data_s4 <- all_data_df %>%
  filter(date >= start_analysis_date & date <= end_analysis_date) %>%
  mutate(post = ifelse(date >= as.Date("2004-06-01"), 1, 0))
panel_data_s4 <- bind_rows(
  filtered_data_s4 %>% mutate(region = "US (Treated)", treatment = 1, inflation_rate = us_inflation),
  filtered_data_s4 %>% mutate(region = "Japan (Control)", treatment = 0, inflation_rate = japan_inflation)
) %>% na.omit()
s4_model <- lm(inflation_rate ~ treatment * post, data = panel_data_s4)
s4_robust_se <- coeftest(s4_model, vcov = NeweyWest(s4_model, lag = L_nw, prewhite = FALSE))
print(s4_robust_se)

# Function to extract results for sensitivity table.
extract_did_results <- function(model_robust_se, scenario_name) {
  did_row <- model_robust_se[row.names(model_robust_se) == "treatment:post", ]
  estimate <- did_row[1]
  std_error <- did_row[2]
  p_value <- did_row[4]
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
    Estimate = format(round(estimate, 4), nsmall = 4),
    `Std. Error` = round(std_error, 4),
    `p-value` = round(p_value, 4),
    Significance = stars
  )
}

# Extract and combine results into a sensitivity table.
results_original <- extract_did_results(robust_se_standard_did, "Original Model")
results_s1 <- extract_did_results(s1_robust_se, "Scenario 1: Pre-period start -1 month")
results_s2 <- extract_did_results(s2_robust_se, "Scenario 2: Pre-period start +1 month")
results_s3 <- extract_did_results(s3_robust_se, "Scenario 3: Event date shifted -1 month")
results_s4 <- extract_did_results(s4_robust_se, "Scenario 4: Event date shifted +1 month")

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

