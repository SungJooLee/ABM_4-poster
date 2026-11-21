# Test the actual inline code04 sections without running exports and plots.
source("R_codes/Halsey_calendar.R")
source("R_codes/Halsey_tick_model_core.R")
script_lines <- readLines("R_codes/code04.Halsey_tick_onecell.R")
label_start <- grep("^# Label reporting periods once", script_lines)
label_end <- grep("^tick_output_dir <-", script_lines)
settings_start <- grep("^# \\[2\\]", script_lines)
run_start <- grep("^# \\[3\\]", script_lines)
timer_start <- grep("^tick_elapsed <-", script_lines)
run_end <- grep("^# Add the reporting labels", script_lines)
code04_labels <- parse(text = script_lines[label_start:(label_end - 1L)])
code04_calendar_check <- parse(text = script_lines[run_start:(timer_start - 1L)])
code04_simulation <- parse(text = script_lines[timer_start:(run_end - 1L)])
settings_env <- new.env()
eval(parse(text = script_lines[settings_start:(run_start - 1L)]), settings_env)
default_tick_parameters <- settings_env$tick_parameters

run_test_cell <- function(weather, p = default_tick_parameters, seed = 123,
                          initial = NULL, verbose = FALSE) {
  tick_weather <- weather
  tick_parameters <- p
  tick_seed <- seed
  tick_initial <- initial
  tick_verbose <- verbose
  eval(code04_calendar_check)
  eval(code04_simulation)
  tick_result
}

label_test_periods <- function(daily, burn_in_years, analysis_years) {
  tick_weather <- daily
  eval(code04_calendar_check)
  eval(code04_labels)
  tick_weather
}
