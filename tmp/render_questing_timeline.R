expressions <- parse("R_codes/code04.Halsey_tick_onecell.R")
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name("tick_elapsed"))) break
  eval(expression)
}
tick_result <- readRDS(file.path(tick_output_dir, "ticks_onecell_result.rds"))
tick_daily <- tick_result$daily
tick_analysis <- tick_result$analysis_daily
annual <- read.csv(file.path(tick_output_dir, "ticks_annual_summary.csv"))
tick_elapsed <- c(elapsed = NA_real_)
code <- readLines("R_codes/code04.Halsey_tick_onecell.R")
start <- grep("^# ---- Calculate the values used", code)
old_plots <- grep("^# \\[5\\]", code)
new_plot <- grep("^# \\[6\\]", code)
eval(parse(text = c(code[start:(old_plots - 1)], code[new_plot:length(code)])))
qa <- sub("^Elapsed seconds: NA$", "Elapsed seconds: unavailable (report regenerated from saved results)", qa)
writeLines(qa, file.path(tick_output_dir, "model_assumptions_and_QA.txt"))
stopifnot(nrow(tick_daily) == 13 * 365, balance == 0, analysis_balance == 0,
          tick_result$parameters$hardening_hours["L"] == 168,
          all(is.finite(as.matrix(tick_daily[c("L_questing", "N_questing", "A_questing")]))))
cat("PASS: 13 years, 7-day larval hardening, finite questing counts and zero accounting errors\n")