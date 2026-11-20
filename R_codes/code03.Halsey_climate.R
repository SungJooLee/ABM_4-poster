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
read_noaa_weather <- function(path, station, start, end, units = "C") {
  units <- match.arg(units, c("C", "F"))
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  na.strings = c("", "NA", "NaN"))
  required <- c("STATION", "DATE", "TMIN", "TMAX")
  if (!all(required %in% names(raw))) {
    stop("Missing required columns: ",
         paste(setdiff(required, names(raw)), collapse = ", "))
  }
  if (anyNA(raw$STATION) || any(raw$STATION != station)) {
    stop("Expected only station ", station, "; check the station selection.")
  }
  dates <- as.Date(raw$DATE, format = "%Y-%m-%d")
  if (anyNA(dates)) {
    stop("DATE must contain valid dates in YYYY-MM-DD format.")
  }
  if (anyDuplicated(dates)) stop("Duplicate station dates; resolve before modeling.")
  numeric_temperature <- function(x) {
    # Missing/invalid temperatures are caught when building the five-year mean.
    result <- suppressWarnings(as.numeric(x))
    result[result == -9999] <- NA_real_  # NOAA missing-temperature code
    if (units == "F") result <- (result - 32) * 5 / 9
    result
  }
  observations <- data.frame(
    date = dates,
    tmin_C = numeric_temperature(raw$TMIN),
    tmax_C = numeric_temperature(raw$TMAX)
  )
  if (any(observations$tmin_C > observations$tmax_C, na.rm = TRUE)) {
    stop("At least one day has TMIN > TMAX; inspect the observations.")
  }
  # Complete the calendar BEFORE taking rolling windows. An absent row must
  # become a missing day, not shorten a seven-day window to seven observations.
  daily <- merge(data.frame(date = seq(start, end, by = "day")),
                 observations, by = "date", all.x = TRUE, sort = TRUE)
  daily$station <- station
  daily$observation_present <- daily$date %in% dates
  daily$temperature_complete <- complete.cases(daily[c("tmin_C", "tmax_C")])
  daily
}

weather_observed <- read_noaa_weather(
  input_path, expected_station, expected_start, expected_end,
  input_temperature_units
)

# [3] Build the five-year mean seasonal cycle ---------------------------------
build_climatology <- function(observed, years) {
  if (length(years) != 5L || anyNA(years) || any(!is.finite(years)) ||
      any(years != floor(years)) || any(diff(sort(years)) != 1)) {
    stop("climatology_years must contain five consecutive, distinct years.")
  }
  years <- sort(years)
  observed_year <- as.integer(format(observed$date, "%Y"))
  baseline <- observed[observed_year %in% years, , drop = FALSE]
  # Require one complete observation per year for each ordinary month/day.
  # This avoids accidentally averaging only 3-4 years on dates with missing data.
  template <- seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  month_days <- format(template, "%m-%d")
  baseline_key <- format(baseline$date, "%m-%d")
  records <- lapply(month_days, function(key) {
    x <- baseline[baseline_key == key, , drop = FALSE]
    if (nrow(x) != 5L || any(!is.finite(x$tmin_C)) ||
        any(!is.finite(x$tmax_C))) {
      stop("Need five complete temperature pairs for ", key,
           " in the selected baseline years; no partial-year mean is used.")
    }
    data.frame(month_day = key, tmin_C = mean(x$tmin_C),
               tmax_C = mean(x$tmax_C), n_years = 5L,
               method = "five_year_month_day_mean")
  })
  normals <- do.call(rbind, records)
  # One leap-day observation in this five-year window would not be a five-year
  # average. Instead derive Feb 29 from its two adjacent climatological days.
  neighbors <- normals[normals$month_day %in% c("02-28", "03-01"), ]
  normals <- rbind(normals, data.frame(
    month_day = "02-29", tmin_C = mean(neighbors$tmin_C),
    tmax_C = mean(neighbors$tmax_C), n_years = 0L,
    method = "interpolated_Feb28_Mar01"
  ))
  normals <- normals[order(normals$month_day), ]
  rownames(normals) <- NULL
  normals$tmean_C <- (normals$tmin_C + normals$tmax_C) / 2
  normals$baseline_start_year <- min(years)
  normals$baseline_end_year <- max(years)
  normals
}
weather_climatology <- build_climatology(weather_observed, climatology_years)

# [4] Repeat the climatology and compute trailing-window indices --------------
# A window is valid only when ALL required calendar days have measurements.
# This avoids biased sums and averages caused by na.rm = TRUE.
trailing_complete <- function(x, width, statistic = c("sum", "mean")) {
  statistic <- match.arg(statistic)
  answer <- rep(NA_real_, length(x))
  if (length(x) < width) return(answer)
  for (i in seq.int(width, length(x))) {
    window <- x[seq.int(i - width + 1L, i)]
    if (all(is.finite(window))) {
      answer[i] <- if (statistic == "sum") sum(window) else mean(window)
    }
  }
  answer
}

repeat_climatology <- function(normals, start_year, n_years, station, base_C = 6) {
  stopifnot(length(start_year) == 1L, is.finite(start_year),
            start_year == floor(start_year), start_year >= 1,
            length(n_years) == 1L, is.finite(n_years), n_years >= 1,
            n_years == floor(n_years))
  start <- as.Date(sprintf("%04d-01-01", start_year))
  end <- as.Date(sprintf("%04d-01-01", start_year + n_years)) - 1L
  # This lead-in represents the preceding cycle, not additional historical
  # observations. It provides valid January windows even in simulation year 1.
  dates <- seq(start - 16L, end, by = "day")
  indices <- match(format(dates, "%m-%d"), normals$month_day)
  daily <- data.frame(date = dates, station = station,
                      month_day = format(dates, "%m-%d"),
                      tmin_C = normals$tmin_C[indices],
                      tmax_C = normals$tmax_C[indices],
                      climate_method = normals$method[indices])
  # Average TMIN/TMAX first; calculate DD and its clipped alternative afterward.
  # In particular, mean(max(0, historical DD)) is NOT the same calculation as
  # max(0, DD from the mean climate), because clipping is nonlinear.
  # NOAA TAVG is deliberately unused: Halsey's equation uses (TMIN+TMAX)/2.
  daily$tmean_C <- (daily$tmin_C + daily$tmax_C) / 2
  daily$dd6_raw_Cday <- daily$tmean_C - base_C
  daily$cdd7_raw_Cday <- trailing_complete(daily$dd6_raw_Cday, 7L, "sum")
  daily$tmean17_C <- trailing_complete(daily$tmean_C, 17L, "mean")
  daily$dd6_nonnegative_Cday <- pmax(0, daily$dd6_raw_Cday)
  daily$cdd7_nonnegative_Cday <- trailing_complete(
    daily$dd6_nonnegative_Cday, 7L, "sum"
  )
  daily$climate_ready <- complete.cases(
    daily[c("tmin_C", "tmax_C", "cdd7_raw_Cday", "tmean17_C")]
  )
  daily <- daily[daily$date >= start, , drop = FALSE]
  rownames(daily) <- NULL
  daily$year <- as.integer(format(daily$date, "%Y"))
  daily$simulation_year <- daily$year - start_year + 1L
  daily$day_of_year <- as.integer(format(daily$date, "%j"))
  daily$simulation_day <- seq_len(nrow(daily))
  daily$days_in_year <- as.integer(as.Date(sprintf("%04d-01-01", daily$year + 1L)) -
                                   as.Date(sprintf("%04d-01-01", daily$year)))
  daily
}
weather_daily <- repeat_climatology(
  weather_climatology, simulation_start_year, simulation_n_years,
  expected_station, base_temperature_C
)
weather_daily$phase <- ifelse(weather_daily$simulation_year <= burn_in_years,
                              "burn_in", "analysis")
weather_daily$analysis_year <- ifelse(weather_daily$phase == "analysis",
                                      weather_daily$simulation_year - burn_in_years, NA_integer_)

# [5] Construct the two-hour model calendar -----------------------------------
make_movement_calendar <- function(daily, step_hours = 2L) {
  stopifnot(step_hours > 0, step_hours == as.integer(step_hours),
            24 %% step_hours == 0)
  steps_per_day <- as.integer(24 / step_hours)
  rows <- rep(seq_len(nrow(daily)), each = steps_per_day)
  calendar <- daily[rows, , drop = FALSE]
  rownames(calendar) <- NULL
  calendar$step <- seq_len(nrow(calendar))
  calendar$step_in_day <- rep(seq_len(steps_per_day), times = nrow(daily))
  calendar$hour_start <- (calendar$step_in_day - 1L) * step_hours
  calendar$dt_hours <- step_hours
  calendar$dt_days <- step_hours / 24
  # This is a fixed model clock: 00,02,...,22 every calendar day. The string
  # is deliberately timezone-free; daylight-saving time must not add/remove
  # deer movement steps. Each row labels [start, start + two hours).
  calendar$model_datetime <- sprintf("Simulation year %02d %s %02d:00:00",
                                     calendar$simulation_year, calendar$month_day,
                                     calendar$hour_start)
  calendar[, c("step", "model_datetime", "step_in_day", "hour_start",
                "dt_hours", "dt_days", names(daily))]
}
weather_2hour <- make_movement_calendar(weather_daily, movement_step_hours)

# [6] Export model inputs and a concise QA record -----------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(weather_climatology, file.path(output_dir, "weather_5year_climatology.csv"),
          row.names = FALSE, na = "")
# Internal dates are retained only as explicitly named implementation fields.
# All interpretation/plotting uses relative simulation time.
climate_export <- function(x) {
  x$year <- NULL
  names(x)[names(x) == "date"] <- "internal_date"
  x
}
write.csv(climate_export(weather_daily), file.path(output_dir, "halsey_weather_daily.csv"),
          row.names = FALSE, na = "")
write.csv(climate_export(weather_2hour), file.path(output_dir, "halsey_weather_2hour.csv"),
          row.names = FALSE, na = "")

weather_annual <- do.call(rbind, lapply(split(weather_daily, weather_daily$year),
  function(x) data.frame(
    phase = x$phase[1], analysis_year = x$analysis_year[1], calendar_days = nrow(x),
    simulation_year = x$simulation_year[1],
    interpolated_leap_days = sum(x$month_day == "02-29"),
    climate_ready_days = sum(x$climate_ready),
    mean_temperature_C = mean(x$tmean_C)
  )
))
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
plot_weather_overview <- function(daily, path) {
  png(path, width = 1800, height = 1400, res = 160)
  on.exit(dev.off())
  par(mfrow = c(3, 1), mar = c(4.5, 4.5, 2.6, 1), mgp = c(2.5, 0.7, 0))
  model_time <- daily$simulation_year + (daily$day_of_year - 1) / daily$days_in_year
  temperature_limits <- range(c(daily$tmin_C, daily$tmax_C), na.rm = TRUE)
  temperature_limits[2] <- temperature_limits[2] + 0.15 * diff(temperature_limits)
  dd_limits <- range(c(daily$cdd7_raw_Cday, daily$cdd7_nonnegative_Cday), na.rm = TRUE)
  dd_limits[2] <- dd_limits[2] + 0.15 * diff(dd_limits)
  plot(model_time, daily$tmax_C, type = "l", col = "firebrick", lwd = 0.8,
       ylim = temperature_limits,
       xlab = "", ylab = "Temperature (C)",
       main = "Repeated five-year mean climate: daily minimum and maximum")
  lines(model_time, daily$tmin_C, col = "steelblue", lwd = 0.8)
  legend("topleft", c("TMAX", "TMIN"), col = c("firebrick", "steelblue"),
         lty = 1, bty = "n", ncol = 2)
  plot(model_time, daily$cdd7_raw_Cday, type = "l", col = "black",
       ylim = dd_limits,
       xlab = "", ylab = "7-day total (C-days)",
       main = "Trailing seven-day DD: printed equation versus zero-clipped DD")
  lines(model_time, daily$cdd7_nonnegative_Cday, col = "darkorange", lty = 2)
  abline(h = 0, col = "grey70")
  legend("topleft", c("Raw DD", "Zero-clipped DD"),
         col = c("black", "darkorange"), lty = c(1, 2), bty = "n", ncol = 2)
  plot(model_time, daily$tmean17_C, type = "l", col = "darkgreen",
       xlab = "Simulation year", ylab = "Temperature (C)",
       main = "Trailing 17-day mean temperature: pre-oviposition weather input")
}
plot_weather_overview(weather_daily, file.path(output_dir, "weather_overview.png"))

plot_annual_thresholds <- function(daily, path) {
  # Lines show weather thresholds only. A tick also needs the correct stage,
  # feeding history, age, and dormancy state before any transition can occur.
  # Showing every year adds identical panels and becomes unreadable for long
  # runs. Four years display the repeated pattern plus leap-day behavior.
  years <- head(unique(daily$simulation_year), 4L)
  png(path, width = 1800, height = 1500, res = 160)
  on.exit(dev.off())
  par(mfrow = c(ceiling(length(years) / 2), 2),
      mar = c(4.4, 4.5, 2.2, 1), oma = c(0, 0, 3.5, 0),
      mgp = c(2.5, 0.7, 0))
  thresholds <- c(28, 58, 81, 110)
  colors <- c("steelblue", "darkgreen", "purple", "firebrick")
  limits <- range(c(daily$cdd7_raw_Cday, thresholds), na.rm = TRUE)
  for (year_value in years) {
    x <- daily[daily$simulation_year == year_value, ]
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
}
plot_annual_thresholds(weather_daily,
                       file.path(output_dir, "weather_thresholds_by_year.png"))

cat(paste(qa_lines[1:16], collapse = "\n"), "\n")
cat("Outputs written to:", output_dir, "\n")
print(weather_annual, row.names = FALSE)

# Objects weather_observed, weather_climatology, weather_daily and weather_2hour
# remain available after source(). Dates label simulation time, not forecasts.
# A later tick model can look up each date's weather and consume one row of the
# two-hour calendar per deer move. Require climate_ready before using the full
# set of indices; do not interpret NA weather history as zero tick activity.
