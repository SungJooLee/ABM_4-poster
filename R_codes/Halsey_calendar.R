# Fixed 365-day seasonal calendar. Month/day labels do not carry a real year.
seasonal_month_days <- function() {
  month_lengths <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)
  sprintf("%02d-%02d", rep(seq_len(12L), month_lengths),
          unlist(lapply(month_lengths, seq_len)))
}

model_calendar <- function(simulation_day) {
  stopifnot(length(simulation_day) > 0L, all(is.finite(simulation_day)),
            all(simulation_day == floor(simulation_day)))
  # Nonpositive days are useful for the climate's preceding-cycle lead-in.
  day_of_year <- (simulation_day - 1L) %% 365L + 1L
  month_day <- seasonal_month_days()[day_of_year]
  data.frame(simulation_day = simulation_day,
             simulation_year = (simulation_day - 1L) %/% 365L + 1L,
             day_of_year = day_of_year,
             month = as.integer(substr(month_day, 1L, 2L)),
             month_day = month_day, days_in_year = 365L)
}
