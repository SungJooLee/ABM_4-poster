# STEP 2: run the tick life cycle in one homogeneous cell.
# Run: source("R_codes/code04.Halsey_tick_onecell.R") from the project root,
# or Rscript /full/path/to/R_codes/code04.Halsey_tick_onecell.R.
# Base R only. Biology functions are in halsey_tick_core.R.
#
# Read the daily repeating climate produced by code03. This script runs a
# stochastic, stage-structured tick model at TWO-HOUR intervals, with daily
# climate forcing. It does not yet call the deer ABM or move ticks among cells.
# Daily outputs are end-of-day stocks and sums of within-day events.
#
# Edit the settings at the top of code03 for burn-in/analysis duration, then run it.
# This script reads those periods from the exported climate calendar.
# Default: 10 burn-in years followed continuously by 3 analysis years. Burn-in
# ticks, ages, timers, dormancy and RNG carry into analysis without resetting.
# This script consumes that exported calendar as-is and checks its continuity.
# The default 900 m2 cell matches the ABM's 30 x 30 m grid. A small stochastic
# population can die out; this is a diagnostic trajectory, not a local forecast.

# [1] Locate inputs and load reusable functions -------------------------------
# Same explicit project path as code03. Update this if the project moves.
tick_project_root <- "D:/github/ABM_4-poster"

source(file.path(tick_project_root, "R_codes", "halsey_tick_core.R"), local = TRUE)
climate_path <- file.path(tick_project_root, "Results.climate", "halsey_weather_daily.csv")
if (!file.exists(climate_path)) stop("Run code03.Halsey_climate.R to create daily weather first.")
tick_weather <- read.csv(climate_path, stringsAsFactors = FALSE)
tick_weather$date <- as.Date(tick_weather$internal_date)
# Count the years assigned to each period by code03; no duplicate settings here.
burn_in_years <- length(unique(tick_weather$simulation_year[tick_weather$phase == "burn_in"]))
analysis_years <- length(unique(tick_weather$simulation_year[tick_weather$phase == "analysis"]))
simulation_start_year <- as.integer(format(tick_weather$date[1], "%Y"))
# Check duration before starting the expensive simulation. Internal dates
# supply seasons only; no forecast year is assigned to the experiment.
tick_weather <- label_tick_periods(tick_weather, burn_in_years, analysis_years)
tick_output_dir <- file.path(tick_project_root, "Results.ticks_onecell")

# [2] Experiment settings -----------------------------------------------------
tick_parameters <- halsey_parameters()
tick_seed <- 123
# All paper values and adaptation choices are annotated in halsey_tick_core.R.
# Optional overrides belong here; examples:
# tick_parameters$cell_area_m2 <- 10000  # one-hectare diagnostic, not a 30 m cell
# tick_parameters$contact_fraction["A"] <- 0.05  # uncalibrated deer exposure
# tick_parameters$hardening_hours["L"] <- 7 * 24 # main-text alternative
# Setting mice_per_ha or deer_equivalents_per_ha to 0 disables new attachments
# to that host. Hosts themselves have NO births, deaths, or state variables.

# [3] Run the model -----------------------------------------------------------
tick_elapsed <- system.time({
  tick_result <- run_tick_cell(tick_weather, tick_parameters, tick_seed)
})
# Reuse the calendar already labeled/checked above instead of checking it again.
reporting_columns <- c("simulation_day", "simulation_year", "day_of_year", "month_day",
                       "days_in_year", "phase", "analysis_year", "analysis_day")
tick_daily <- cbind(tick_result$daily, tick_weather[reporting_columns])
tick_analysis <- tick_daily[tick_daily$phase == "analysis", , drop = FALSE]
tick_burn_in <- tick_daily[tick_daily$phase == "burn_in", , drop = FALSE]
tick_result$daily <- tick_daily
tick_result$analysis_daily <- tick_analysis
tick_result$burn_in_years <- burn_in_years
tick_result$analysis_years <- analysis_years
# Beginning-of-analysis population is the preceding day's END stock.
tick_result$analysis_initial_count <- if (nrow(tick_burn_in)) {
  tail(tick_burn_in$total_ticks, 1)
} else {
  tick_result$initial_count
}

# [4] Export daily counts, final cohorts, parameters and checks ----------------
dir.create(tick_output_dir, recursive = TRUE, showWarnings = FALSE)
# Main results exclude burn-in. Keep the full run separately for diagnostics.
daily_tables <- list(tick_analysis, tick_daily, tick_burn_in)
daily_filenames <- c("ticks_daily.csv", "ticks_daily_full_run.csv", "ticks_daily_burn_in.csv")
for (table_index in seq_along(daily_tables)) {
  daily_export <- daily_tables[[table_index]]
  names(daily_export)[names(daily_export) == "date"] <- "internal_date"
  write.csv(daily_export, file.path(tick_output_dir, daily_filenames[table_index]),
            row.names = FALSE)
}

final_cohorts_export <- tick_result$final_cohorts
final_cohorts_export$eligible_simulation_year <- final_cohorts_export$eligible_year - simulation_start_year + 1
final_cohorts_export$eligible_year <- NULL
write.csv(final_cohorts_export, file.path(tick_output_dir, "ticks_final_cohorts.csv"),
          row.names = FALSE)
# RDS retains exact column types and settings for later inspection/restarts.
# A restart also needs RNG state and initial clock conventions; this file is
# an analysis artifact, not an automatic checkpoint/resume implementation.
saveRDS(tick_result, file.path(tick_output_dir, "ticks_onecell_result.rds"))
parameter_values <- unlist(tick_parameters)
write.csv(data.frame(parameter = names(parameter_values), value = unname(parameter_values)),
          file.path(tick_output_dir, "parameters.csv"), row.names = FALSE)

simulation_years <- sort(unique(tick_daily$simulation_year))
annual_records <- vector("list", length(simulation_years))
for (year_index in seq_along(simulation_years)) {
  x <- tick_daily[tick_daily$simulation_year == simulation_years[year_index], ]
  annual_records[[year_index]] <- data.frame(
    simulation_year = x$simulation_year[1], phase = x$phase[1],
    analysis_year = x$analysis_year[1],
    peak_questing_L = max(x$L_questing), peak_questing_N = max(x$N_questing),
    peak_questing_A_females = max(x$A_questing),
    attached_L = sum(x$attached_L), attached_N = sum(x$attached_N),
    attached_A = sum(x$attached_A), eggs_laid = sum(x$eggs_laid),
    year_end_ticks = tail(x$total_ticks, 1)
  )
}
annual_full <- do.call(rbind, annual_records)

rownames(annual_full) <- NULL
annual <- annual_full[annual_full$phase == "analysis", , drop = FALSE]
write.csv(annual, file.path(tick_output_dir, "ticks_annual_summary.csv"), row.names = FALSE)
write.csv(annual_full, file.path(tick_output_dir, "ticks_annual_summary_full_run.csv"), row.names = FALSE)

loss_columns <- c("deaths_natural", "deaths_grooming", "deaths_molting",
                  "males_removed", "mothers_removed")
total_losses <- sum(as.matrix(tick_daily[loss_columns]))
balance <- tail(tick_daily$total_ticks, 1) -
  (tick_result$initial_count + sum(tick_daily$eggs_laid) - total_losses)
analysis_balance <- tail(tick_analysis$total_ticks, 1) -
  (tick_result$analysis_initial_count + sum(tick_analysis$eggs_laid) -
     sum(as.matrix(tick_analysis[loss_columns])))
# Each two-hour step already checks population accounting in the core.
# These totals are reported below for inspection, without duplicate stops.
extinction <- tick_daily$simulation_day[tick_daily$total_ticks == 0]
qa <- c(
  "Halsey-style one-cell diagnostic: NOT a calibrated spatial tick forecast",
  paste("Climate input:", climate_path),
  paste("Duration:", burn_in_years, "burn-in years +", analysis_years, "analysis years"),
  paste("Seed:", tick_seed), paste("Cell area (m2):", tick_parameters$cell_area_m2),
  paste("Initial modeled ticks:", tick_result$initial_count),
  paste("Final modeled ticks:", tail(tick_daily$total_ticks, 1)),
  paste("Total eggs laid:", sum(tick_daily$eggs_laid)),
  paste("Population accounting error:", balance),
  paste("First extinction simulation day:", if (length(extinction)) min(extinction) else "none"),
  paste("Elapsed seconds:", round(tick_elapsed[["elapsed"]], 2)),
  paste("Beginning-of-analysis ticks (carried from burn-in):", tick_result$analysis_initial_count),
  paste("Analysis-only population accounting error:", analysis_balance),
  "Main daily/annual outputs and seasonal plots exclude burn-in; full-run files retain it.",
  "internal_date is a dummy Gregorian seasonal clock, not a forecast date.",
  "Burn-in and analysis form one continuous run; no population, timer, climate-history or RNG reset.",
  "Ten years of burn-in is a configured duration, not a demonstrated convergence criterion.",
  "Counts are integer stochastic cohorts, not expected fractional counts.",
  "A means reproductive female adults only; removed males are counted separately.",
  "Mice are a constant input. Feeding-on-mouse states describe ticks, not mice.",
  "Deer encounter opportunities are fixed; deer identities/transport come in the spatial step.",
  "Conditional hourly attachment probabilities are multiplied by provisional contact fractions.",
  "Reference densities/initial conditions are paper-based starting points, not local estimates.",
  "Initial overwintering questers start at stage age 0; early years include initialization effects.",
  "Feeding capacity = floor(fixed host abundance * per-host burden limit), by stage.",
  "Mortality and grooming use hourly Bernoulli trials within each two-hour interval.",
  "Transitions/new attachments occur at interval end; no same-interval transition cascades.",
  "Natural mortality applies before grooming; host slots free at detachment become usable next step.",
  "Winter dormancy pauses off-host development/preoviposition but not mortality or ongoing feeding.",
  "Autumn dormancy is latched until next calendar year; new nymphs also wait until next year.",
  "Larval hardening defaults to 21 days (supplement); main text says 7 days.",
  "No additional egg-development or engorged-nymph incubation duration was specified clearly.",
  "Adults can reproduce in warm autumn if gates/timers permit; no undocumented spring-only gate.",
  "Negative polynomial fecundity is clipped to 0; eggs per female are rounded to integers.",
  "Seven-day DD is a rolling sum; negative DD follows the printed equation in code03.",
  "All detachment/grooming survivors stay in the cell in this nonspatial diagnostic.",
  "Original NetLogo archive could not be retrieved; these timing choices need source-code review.",
  "Passing accounting/unit checks does not establish empirical phenology or abundance accuracy."
)
last_year <- tick_analysis[tick_analysis$analysis_year == analysis_years, ]
for (stage in c("L", "N", "A")) {
  values <- last_year[[paste0(stage, "_questing")]]
  qa <- c(qa, paste("Final-year questing peak", stage, ":",
                   if (max(values) > 0) paste("analysis year", analysis_years,
                                             last_year$month_day[which.max(values)])
                   else "none", "(count:", max(values), ")"))
}
qa <- c(qa, "Review seasonal peaks against field observations and source timing rules before spatial coupling.")
writeLines(qa, file.path(tick_output_dir, "model_assumptions_and_QA.txt"))

# [5] Plot population trajectories and seasonal questing ----------------------
# Reused for the analysis-only plot and the full-run plot.
plot_tick_population <- function(daily, path, full_run = FALSE) {
  png(path, width = 1800, height = 1300, res = 160)
  on.exit(dev.off())
  par(mfrow = c(2, 1), mar = c(4.5, 4.5, 3, 1), mgp = c(2.5, 0.7, 0))
  stocks <- daily[c("eggs", "L_total", "N_total", "A_total")]
  colors <- c("grey40", "darkorange", "steelblue", "firebrick")
  model_time <- (if (full_run) daily$simulation_year else daily$analysis_year) +
    (daily$day_of_year - 1) / daily$days_in_year
  year_ticks <- unique(if (full_run) daily$simulation_year else daily$analysis_year)
  time_label <- if (full_run) "Simulation year (dashed line: end of burn-in)" else "Analysis year (after burn-in)"
  # log10(count + 1) displays zero counts without inventing nonzero populations.
  matplot(model_time, log10(1 + as.matrix(stocks)), type = "l", lty = 1,
          ylim = c(0, max(1, log10(1 + as.matrix(stocks))) * 1.15),
          col = colors, xaxt = "n", xlab = "", ylab = "log10(count + 1)",
          main = "One-cell diagnostic: eggs, larvae, nymphs and adult females")
  axis(1, at = year_ticks, labels = year_ticks)
  if (full_run && burn_in_years > 0) abline(v = burn_in_years + 1, lty = 2, col = "grey40")
  legend("topleft", c("Eggs", "Larvae", "Nymphs", "Adult females"),
         col = colors, lty = 1, bty = "n", ncol = 4, cex = 0.85)
  feeding <- daily[c("L_feeding", "N_feeding", "A_feeding")]
  matplot(model_time, as.matrix(feeding), type = "l", lty = 1,
          ylim = c(0, max(1, as.matrix(feeding)) * 1.15),
          col = colors[-1], xaxt = "n", xlab = time_label, ylab = "Feeding ticks",
          main = "Ticks feeding on constant host availability")
  axis(1, at = year_ticks, labels = year_ticks)
  if (full_run && burn_in_years > 0) abline(v = burn_in_years + 1, lty = 2, col = "grey40")
  legend("topleft", c("Larvae on mice", "Nymphs on mice", "Females on deer"),
         col = colors[-1], lty = 1, bty = "n", ncol = 3, cex = 0.85)
}
plot_tick_population(tick_analysis, file.path(tick_output_dir, "tick_population.png"))
plot_tick_population(tick_daily, file.path(tick_output_dir, "tick_population_full_run.png"), full_run = TRUE)

# Seasonal comparison of the three analysis years.
png(file.path(tick_output_dir, "tick_seasonality.png"), width = 1800, height = 1400, res = 160)
par(mfrow = c(3, 1), mar = c(4, 4.5, 2.6, 1), mgp = c(2.5, 0.7, 0))
years <- unique(tick_analysis$analysis_year)
colors <- hcl.colors(length(years), "Dark 3")
# Use month/day on a leap-year reference axis so March is aligned across years.
seasonal_x <- as.Date(paste0("2000-", format(tick_analysis$date, "%m-%d")))
ticks <- as.Date(paste0("2000-", sprintf("%02d", 1:12), "-01"))
for (stage in c("L", "N", "A")) {
  values <- tick_analysis[[paste0(stage, "_questing")]]
  plot(range(seasonal_x), c(0, max(1, max(values)) * 1.15), type = "n", xaxt = "n",
       xlab = "", ylab = "Questing ticks", main = switch(stage,
         L = "Larvae", N = "Nymphs", A = "Adult females"))
  axis.Date(1, at = ticks, format = "%b")
  for (j in seq_along(years)) {
    idx <- tick_analysis$analysis_year == years[j]
    lines(seasonal_x[idx], values[idx], col = colors[j])
  }
  if (stage == "L") legend("topleft", paste("Analysis year", years), col = colors, lty = 1,
                             ncol = 5, cex = 0.75, bty = "n")
}
dev.off()

cat("\n", paste(qa[1:11], collapse = "\n"), "\n")
cat("Outputs:", tick_output_dir, "\n")
print(annual, row.names = FALSE)
