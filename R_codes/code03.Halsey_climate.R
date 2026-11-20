# Prepare a repeating five-year mean climate for a Halsey-style tick model.
# Run from the project root with:
#   source("R_codes/code03.Halsey_climate.R")
# Or run this file directly with Rscript. Only base R is required.
#
# This is the WEATHER AND CALENDAR module, not a tick population simulation.
# Deer movement remains at two-hour intervals. Daily weather is repeated across
# the day's 12 intervals; it is not interpolated into invented hourly weather.
# Average TMIN and TMAX separately by month/day across the selected five years,
# then repeat that seasonal cycle. No random weather is drawn between years.
# "Five-year climatology" here means a short-period average of this station's
# observations, not an official NOAA long-term climate-normal product.
#
# Source: Halsey & Miller (2018), Ecological Modelling 387:96-106,
# doi:10.1016/j.ecolmodel.2018.09.005, Eq. 1, sections 2.9.1-2.9.4,
# Table 1, and supplementary Table A.1.
# Important interpretation choices:
# * DD = (TMIN + TMAX)/2 - 6, exactly as printed (negative DD allowed).
# * CDD7 is a trailing SUM of seven daily DD values, including today.
#   The text also says "mean cumulative"; this implementation follows the
#   cumulative definition in Table A.1 and must be checked against source code.
# * The pre-oviposition temperature is a trailing 17-day MEAN of daily mean
#   temperatures, not a mean of DD. It is a temperature input, not the female's
#   separate 17-day pre-oviposition timer.
# * A zero-clipped DD alternative is exported separately for sensitivity checks.
# * Missing baseline temperatures stop climatology construction.
# * Feb 29 is interpolated from Feb 28/Mar 1 climatology (explicit policy).
# * Initial rolling windows use the preceding 16 days of the SAME repeated
#   climatology, so all climate indices are available on simulation day one.

# [1] Settings ---------------------------------------------------------------
input_filename <- "NOAA_ClimateData_2020-2025.csv"
expected_station <- "USW00004781"
expected_start <- as.Date("2020-01-01")
expected_end <- as.Date("2025-12-31")

# The CSV contains SIX years. Set to the latest FIVE complete years;
climatology_years <- 2021:2025
# Experiment duration: run code03 again after changing these settings, then code04.
burn_in_years <- 10
analysis_years <- 3
simulation_n_years <- burn_in_years + analysis_years

# Dummy calendar anchor used only for month/day arithmetic and dormancy rules.
# Simulation years have Gregorian month lengths, including leap days, but
# are not predictions for any particular historical/future calendar years.
simulation_start_year <- 2001

# This supplied NOAA export appears to use Celsius. CSV headers do not encode
# units, so verify the NOAA order's unit setting. Change to "F" for a Fahrenheit
# export; no automatic unit guessing or division by 10 is performed.
input_temperature_units <- "C"
base_temperature_C <- 6
movement_step_hours <- 2

project_root <- "D:/github/ABM_4-poster"
input_path <- file.path(project_root, input_filename)
output_dir <- file.path(project_root, "Results.climate")

# [2] Read and validate the station's daily observations -----------------------
# Read the CSV and convert its date/temperature columns.
raw <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE,
                na.strings = c("", "NA", "NaN"))
observations <- data.frame(
  date = as.Date(raw$DATE, format = "%Y-%m-%d"),
  tmin_C = suppressWarnings(as.numeric(raw$TMIN)),
  tmax_C = suppressWarnings(as.numeric(raw$TMAX))
)
# NOAA uses -9999 for missing temperatures. Missing baseline values are
# checked below before calculating the five-year means.
observations$tmin_C[observations$tmin_C == -9999] <- NA_real_
observations$tmax_C[observations$tmax_C == -9999] <- NA_real_
if (input_temperature_units == "F") {
  observations$tmin_C <- (observations$tmin_C - 32) * 5 / 9
  observations$tmax_C <- (observations$tmax_C - 32) * 5 / 9
}

# Complete the observation calendar so absent dates remain visible as NA.
weather_observed <- merge(
  data.frame(date = seq(expected_start, expected_end, by = "day")),
  observations, by = "date", all.x = TRUE, sort = TRUE
)
weather_observed$station <- expected_station
weather_observed$observation_present <- weather_observed$date %in% observations$date
weather_observed$temperature_complete <- complete.cases(
  weather_observed[c("tmin_C", "tmax_C")]
)

# [3] Build the five-year mean seasonal cycle ---------------------------------
if (length(climatology_years) != 5 || anyNA(climatology_years) || any(!is.finite(climatology_years)) ||
    any(climatology_years != floor(climatology_years)) || any(diff(sort(climatology_years)) != 1)) {
  stop("climatology_years must contain five consecutive, distinct years.")
}
climatology_years <- sort(climatology_years)
observed_year <- as.integer(format(weather_observed$date, "%Y"))
baseline <- weather_observed[observed_year %in% climatology_years, , drop = FALSE]
# Require one complete observation per year for each ordinary month/day.
# This avoids accidentally averaging only 3-4 years on dates with missing data.
template <- seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
month_days <- format(template, "%m-%d")
baseline_key <- format(baseline$date, "%m-%d")
records <- vector("list", length(month_days))
for (day_index in seq_along(month_days)) {
  key <- month_days[day_index]
  x <- baseline[baseline_key == key, , drop = FALSE]
  if (nrow(x) != 5 || any(!is.finite(x$tmin_C)) ||
      any(!is.finite(x$tmax_C))) {
    stop("Need five complete temperature pairs for ", key,
         " in the selected baseline years; no partial-year mean is used.")
  }
  records[[day_index]] <- data.frame(month_day = key, tmin_C = mean(x$tmin_C),
             tmax_C = mean(x$tmax_C), n_years = 5,
             method = "five_year_month_day_mean")
}
weather_climatology <- do.call(rbind, records)
# One leap-day observation in this five-year window would not be a five-year
# average. Instead derive Feb 29 from its two adjacent climatological days.
neighbors <- weather_climatology[weather_climatology$month_day %in% c("02-28", "03-01"), ]
weather_climatology <- rbind(weather_climatology, data.frame(
  month_day = "02-29", tmin_C = mean(neighbors$tmin_C),
  tmax_C = mean(neighbors$tmax_C), n_years = 0,
  method = "interpolated_Feb28_Mar01"
))
weather_climatology <- weather_climatology[order(weather_climatology$month_day), ]
rownames(weather_climatology) <- NULL
weather_climatology$tmean_C <- (weather_climatology$tmin_C + weather_climatology$tmax_C) / 2
weather_climatology$baseline_start_year <- min(climatology_years)
weather_climatology$baseline_end_year <- max(climatology_years)

# [4] Repeat the climatology and compute trailing-window indices --------------
# This helper is reused for three rolling indices.
# A window is valid only when ALL required calendar days have measurements.
# This avoids biased sums and averages caused by na.rm = TRUE.
trailing_complete <- function(x, width, statistic = c("sum", "mean")) {
  statistic <- match.arg(statistic)
  answer <- rep(NA_real_, length(x))
  if (length(x) < width) return(answer)
  for (i in seq.int(width, length(x))) {
    window <- x[seq.int(i - width + 1, i)]
    if (all(is.finite(window))) {
      answer[i] <- if (statistic == "sum") sum(window) else mean(window)
    }
  }
  answer
}

stopifnot(length(simulation_start_year) == 1, is.finite(simulation_start_year),
          simulation_start_year == floor(simulation_start_year), simulation_start_year >= 1,
          length(simulation_n_years) == 1, is.finite(simulation_n_years), simulation_n_years >= 1,
          simulation_n_years == floor(simulation_n_years))
simulation_start <- as.Date(sprintf("%04d-01-01", simulation_start_year))
simulation_end <- as.Date(sprintf("%04d-01-01", simulation_start_year + simulation_n_years)) - 1
# This lead-in represents the preceding cycle, not additional historical
# observations. It provides valid January windows even in simulation year 1.
dates <- seq(simulation_start - 16, simulation_end, by = "day")
indices <- match(format(dates, "%m-%d"), weather_climatology$month_day)
weather_daily <- data.frame(date = dates, station = expected_station,
                    month_day = format(dates, "%m-%d"),
                    tmin_C = weather_climatology$tmin_C[indices],
                    tmax_C = weather_climatology$tmax_C[indices],
                    climate_method = weather_climatology$method[indices])
# Average TMIN/TMAX first; calculate DD and its clipped alternative afterward.
# In particular, mean(max(0, historical DD)) is NOT the same calculation as
# max(0, DD from the mean climate), because clipping is nonlinear.
# NOAA TAVG is deliberately unused: Halsey's equation uses (TMIN+TMAX)/2.
weather_daily$tmean_C <- (weather_daily$tmin_C + weather_daily$tmax_C) / 2
weather_daily$dd6_raw_Cday <- weather_daily$tmean_C - base_temperature_C
weather_daily$cdd7_raw_Cday <- trailing_complete(weather_daily$dd6_raw_Cday, 7, "sum")
weather_daily$tmean17_C <- trailing_complete(weather_daily$tmean_C, 17, "mean")
weather_daily$dd6_nonnegative_Cday <- pmax(0, weather_daily$dd6_raw_Cday)
weather_daily$cdd7_nonnegative_Cday <- trailing_complete(
  weather_daily$dd6_nonnegative_Cday, 7, "sum"
)
weather_daily$climate_ready <- complete.cases(
  weather_daily[c("tmin_C", "tmax_C", "cdd7_raw_Cday", "tmean17_C")]
)
weather_daily <- weather_daily[weather_daily$date >= simulation_start, , drop = FALSE]
rownames(weather_daily) <- NULL
weather_daily$year <- as.integer(format(weather_daily$date, "%Y"))
weather_daily$simulation_year <- weather_daily$year - simulation_start_year + 1
weather_daily$day_of_year <- as.integer(format(weather_daily$date, "%j"))
weather_daily$simulation_day <- seq_len(nrow(weather_daily))
weather_daily$days_in_year <- as.integer(as.Date(sprintf("%04d-01-01", weather_daily$year + 1)) -
                                 as.Date(sprintf("%04d-01-01", weather_daily$year)))

weather_daily$phase <- ifelse(weather_daily$simulation_year <= burn_in_years,
                              "burn_in", "analysis")
weather_daily$analysis_year <- ifelse(weather_daily$phase == "analysis",
                                      weather_daily$simulation_year - burn_in_years, NA_integer_)

# [5] Construct the two-hour model calendar -----------------------------------
stopifnot(movement_step_hours > 0, movement_step_hours == as.integer(movement_step_hours),
          24 %% movement_step_hours == 0)
steps_per_day <- as.integer(24 / movement_step_hours)
rows <- rep(seq_len(nrow(weather_daily)), each = steps_per_day)
weather_2hour <- weather_daily[rows, , drop = FALSE]
rownames(weather_2hour) <- NULL
weather_2hour$step <- seq_len(nrow(weather_2hour))
weather_2hour$step_in_day <- rep(seq_len(steps_per_day), times = nrow(weather_daily))
weather_2hour$hour_start <- (weather_2hour$step_in_day - 1) * movement_step_hours
weather_2hour$dt_hours <- movement_step_hours
weather_2hour$dt_days <- movement_step_hours / 24
# This is a fixed model clock: 00,02,...,22 every calendar day. The string
# is deliberately timezone-free; daylight-saving time must not add/remove
# deer movement steps. Each row labels [start, start + two hours).
weather_2hour$model_datetime <- sprintf("Simulation year %02d %s %02d:00:00",
                                   weather_2hour$simulation_year, weather_2hour$month_day,
                                   weather_2hour$hour_start)
weather_2hour <- weather_2hour[, c("step", "model_datetime", "step_in_day", "hour_start",
              "dt_hours", "dt_days", names(weather_daily))]

# [6] Export model inputs and a concise QA record -----------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(weather_climatology, file.path(output_dir, "weather_5year_climatology.csv"),
          row.names = FALSE, na = "")
# Internal dates are retained only as explicitly named implementation fields.
# All interpretation/plotting uses relative simulation time.
daily_export <- weather_daily
hourly_export <- weather_2hour
daily_export$year <- NULL
hourly_export$year <- NULL
names(daily_export)[names(daily_export) == "date"] <- "internal_date"
names(hourly_export)[names(hourly_export) == "date"] <- "internal_date"
write.csv(daily_export, file.path(output_dir, "halsey_weather_daily.csv"),
          row.names = FALSE, na = "")
write.csv(hourly_export, file.path(output_dir, "halsey_weather_2hour.csv"),
          row.names = FALSE, na = "")

annual_records <- vector("list", simulation_n_years)
for (year_index in seq_len(simulation_n_years)) {
  x <- weather_daily[weather_daily$simulation_year == year_index, ]
  annual_records[[year_index]] <- data.frame(
    phase = x$phase[1], analysis_year = x$analysis_year[1], calendar_days = nrow(x),
    simulation_year = x$simulation_year[1],
    interpolated_leap_days = sum(x$month_day == "02-29"),
    climate_ready_days = sum(x$climate_ready),
    mean_temperature_C = mean(x$tmean_C)
  )
}
weather_annual <- do.call(rbind, annual_records)

rownames(weather_annual) <- NULL
write.csv(weather_annual, file.path(output_dir, "weather_annual_QA.csv"),
          row.names = FALSE, na = "")

ready_dates <- weather_daily$date[weather_daily$climate_ready]
qa_lines <- c(
  "NOAA / Halsey repeating five-year climatology", paste("Input:", input_path),
  paste("Station:", expected_station),
  paste("Configured input temperature units:", input_temperature_units),
  "CSV headers do not verify units; confirm Celsius in the NOAA download settings.",
  paste("Observed input calendar:", expected_start, "through", expected_end),
  paste("Five-year baseline:", paste(climatology_years, collapse = ", ")),
  "TMIN/TMAX are averaged separately by month/day, with no stochastic weather.",
  "Feb 29 is the mean of the Feb 28 and Mar 1 normal temperatures.",
  paste("Duration:", burn_in_years, "burn-in years +", analysis_years, "analysis years"),
  "internal_date is a dummy seasonal clock, not an observation/forecast date.",
  paste("Daily rows:", nrow(weather_daily)),
  paste("Absent input dates:", sum(!weather_observed$observation_present)),
  paste("Input days missing TMIN or TMAX:", sum(!weather_observed$temperature_complete)),
  paste("Two-hour rows:", nrow(weather_2hour)),
  paste("First simulation day with complete climate history:",
        if (length(ready_dates)) min(weather_daily$simulation_day[weather_daily$climate_ready]) else "none"),
  "Raw DD = (TMIN + TMAX)/2 - 6. Negative values are retained.",
  "CDD7 = trailing seven-calendar-day SUM, including the current date.",
  "Tmean17 = trailing 17-calendar-day MEAN of (TMIN + TMAX)/2.",
  "Zero-clipped DD and CDD7 are separate sensitivity-analysis columns.",
  "Windows use a 16-day lead-in from the preceding climatological cycle.",
  "Windows do not reset on January 1. Leap days can change nearby rolling indices.",
  "This is a five-year station climatology, not an official NOAA normal product.",
  "Daily weather repeats over 12 steps, not measured two-hour temperatures.",
  "No tick transitions are inferred solely from crossing a plotted threshold.",
  "NOAA observation quality flags are not included in this input CSV."
)
writeLines(qa_lines, file.path(output_dir, "weather_QA.txt"))

# [7] Diagnostic plots --------------------------------------------------------
png(file.path(output_dir, "weather_overview.png"), width = 1800, height = 1400, res = 160)
par(mfrow = c(3, 1), mar = c(4.5, 4.5, 2.6, 1), mgp = c(2.5, 0.7, 0))
model_time <- weather_daily$simulation_year + (weather_daily$day_of_year - 1) / weather_daily$days_in_year
temperature_limits <- range(c(weather_daily$tmin_C, weather_daily$tmax_C), na.rm = TRUE)
temperature_limits[2] <- temperature_limits[2] + 0.15 * diff(temperature_limits)
dd_limits <- range(c(weather_daily$cdd7_raw_Cday, weather_daily$cdd7_nonnegative_Cday), na.rm = TRUE)
dd_limits[2] <- dd_limits[2] + 0.15 * diff(dd_limits)
plot(model_time, weather_daily$tmax_C, type = "l", col = "firebrick", lwd = 0.8,
     ylim = temperature_limits,
     xlab = "", ylab = "Temperature (C)",
     main = "Repeated five-year mean climate: daily minimum and maximum")
lines(model_time, weather_daily$tmin_C, col = "steelblue", lwd = 0.8)
legend("topleft", c("TMAX", "TMIN"), col = c("firebrick", "steelblue"),
       lty = 1, bty = "n", ncol = 2)
plot(model_time, weather_daily$cdd7_raw_Cday, type = "l", col = "black",
     ylim = dd_limits,
     xlab = "", ylab = "7-day total (C-days)",
     main = "Trailing seven-day DD: printed equation versus zero-clipped DD")
lines(model_time, weather_daily$cdd7_nonnegative_Cday, col = "darkorange", lty = 2)
abline(h = 0, col = "grey70")
legend("topleft", c("Raw DD", "Zero-clipped DD"),
       col = c("black", "darkorange"), lty = c(1, 2), bty = "n", ncol = 2)
plot(model_time, weather_daily$tmean17_C, type = "l", col = "darkgreen",
     xlab = "Simulation year", ylab = "Temperature (C)",
     main = "Trailing 17-day mean temperature: pre-oviposition weather input")
dev.off()

# Lines show weather thresholds only. A tick also needs the correct stage,
# feeding history, age, and dormancy state before any transition can occur.
# Showing every year adds identical panels and becomes unreadable for long
# runs. Four years display the repeated pattern plus leap-day behavior.
years <- head(unique(weather_daily$simulation_year), 4)
png(file.path(output_dir, "weather_thresholds_by_year.png"), width = 1800, height = 1500, res = 160)
par(mfrow = c(ceiling(length(years) / 2), 2),
    mar = c(4.4, 4.5, 2.2, 1), oma = c(0, 0, 3.5, 0),
    mgp = c(2.5, 0.7, 0))
thresholds <- c(28, 58, 81, 110)
colors <- c("steelblue", "darkgreen", "purple", "firebrick")
limits <- range(c(weather_daily$cdd7_raw_Cday, thresholds), na.rm = TRUE)
for (year_value in years) {
  x <- weather_daily[weather_daily$simulation_year == year_value, ]
  plot(x$day_of_year, x$cdd7_raw_Cday, type = "l", ylim = limits,
       xlim = c(1, 366), xlab = "Day of year", ylab = "7-day DD (C-days)",
       main = paste("Simulation year", year_value))
  abline(h = thresholds, col = colors, lty = 2)
  if (year_value == years[1]) {
    legend("topleft", paste(thresholds, "C-days"), col = colors,
           lty = 2, bty = "n", cex = 0.75, ncol = 2)
  }
}
mtext("Repeated mean climate: Halsey thresholds (first four simulation years)",
      outer = TRUE, line = 1, cex = 1.05)
dev.off()

cat(paste(qa_lines[1:16], collapse = "\n"), "\n")
cat("Outputs written to:", output_dir, "\n")
print(weather_annual, row.names = FALSE)

# Objects weather_observed, weather_climatology, weather_daily and weather_2hour
# remain available after source(). Dates label simulation time, not forecasts.
# A later tick model can look up each date's weather and consume one row of the
# two-hour calendar per deer move. Require climate_ready before using the full
# set of indices; do not interpret NA weather history as zero tick activity.
