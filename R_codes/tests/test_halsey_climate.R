# Rebuild climate outputs and check the fixed seasonal cycle end to end.
# Run from the project root: Rscript R_codes/tests/test_halsey_climate.R
source("R_codes/code03.Climate_data.R")
source("R_codes/tests/halsey_script_sections.R")

stopifnot(nrow(weather_climatology) == 365L,
          !any(weather_climatology$month_day == "02-29"),
          nrow(weather_daily) == 365L * simulation_n_years,
          nrow(weather_2hour) == 365L * simulation_n_years * 12L,
          all(weather_daily$climate_ready),
          all(table(weather_daily$simulation_year) == 365L),
          all(table(weather_2hour$simulation_day) == 12L),
          !any(c("date", "internal_date", "year") %in% names(weather_daily)))
tick_weather <- weather_daily
eval(code04_calendar_check)

# Every year's temperatures AND rolling indices must be identical, including
# January 1 and the Feb 28 -> Mar 1 transition (days 59 -> 60).
climate_columns <- c("tmin_C", "tmax_C", "dd6_raw_Cday", "cdd7_raw_Cday",
                     "tmean17_C", "cdd7_nonnegative_Cday")
for (column in climate_columns) {
  stopifnot(isTRUE(all.equal(weather_daily[[column]],
                            rep(weather_daily[[column]][1:365], simulation_n_years))))
}
for (day in c(1L, 60L, 366L)) {
  window7 <- (seq.int(day - 6L, day) - 1L) %% 365L + 1L
  window17 <- (seq.int(day - 16L, day) - 1L) %% 365L + 1L
  stopifnot(isTRUE(all.equal(weather_daily$cdd7_raw_Cday[day],
                            sum(weather_climatology$tmean_C[window7] - 6))),
            isTRUE(all.equal(weather_daily$tmean17_C[day],
                             mean(weather_climatology$tmean_C[window17]))))
}

# Verify the actual exported file is accepted by the tick-model consumer.
exported <- read.csv(file.path(output_dir, "halsey_weather_daily.csv"))
periods <- label_test_periods(exported, burn_in_years, analysis_years)
stopifnot(sum(periods$phase == "burn_in") == burn_in_years * 365L,
          sum(periods$phase == "analysis") == analysis_years * 365L)
probe <- run_test_cell(exported[350:380, ], initial = empty_cohorts(), verbose = FALSE)
stopifnot(nrow(probe$daily) == 31L, all(probe$daily$balance_error == 0))
cat("PASS: 365-day climate repetition, rolling windows, two-hour clock and tick input integration\n")
