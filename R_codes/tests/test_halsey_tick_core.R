# Behavioral tests for tick accounting and biological timing; base R only.
# Run from the project root:
# Rscript R_codes/tests/test_halsey_tick_core.R
source("R_codes/tests/halsey_script_sections.R")

# Detailed programming checks belong here, outside the simulation flow.
# stopifnot() is intentional in this file: a failed expectation fails the test.
validate_cohorts <- function(x) {
  tick_stages <- c("E", "L", "N", "A")
  tick_states <- c("egg", "hardening", "questing", "feeding", "engorged", "preoviposition")
  stopifnot(all(x$stage %in% tick_stages), all(x$state %in% tick_states),
            all(is.finite(x$n)), all(x$n >= 0), all(x$n == floor(x$n)),
            all(is.finite(x$age_h)), all(x$age_h >= 0),
            all(is.finite(x$timer_h)), all(x$timer_h >= 0),
            !anyNA(x$dormant), !anyNA(x$eligible_year),
            all((x$stage == "E") == (x$state == "egg")),
            all(x$stage[x$state == "preoviposition"] == "A"))
  invisible(TRUE)
}

advance <- function(x, p, n = 1, year = 1L, month = 6L, cdd = 120, temp = 18) {
  result <- NULL
  for (i in seq_len(n)) {
    result <- advance_ticks(x, year, month, cdd, temp, p)
    x <- result$cohorts
    stopifnot(result$balance_error == 0)
    validate_cohorts(x)
  }
  result
}
passed <- function(label) cat("PASS:", label, "\n")

# Shared controlled settings: turn off random losses/attachments unless a test
# explicitly enables them. This setup runs once for the tests below.
p <- default_tick_parameters
for (key in c("death_unfed_young", "death_unfed_old", "death_replete",
               "grooming_hourly", "attachment_hourly")) p[[key]][] <- 0
p$death_egg_hourly <- 0
p$molt_success[] <- 1
p$female_fraction <- 1
for (stage in c("L", "N", "A")) {
  x <- new_cohort(stage, "feeding", 5, 1)
  before <- advance(x, p, p$feeding_hours[stage] / 2 - 1)
  stopifnot(all(before$cohorts$state == "feeding"))
  done <- advance(before$cohorts, p)
  stopifnot(all(done$cohorts$state == "engorged"), sum(done$cohorts$n) == 5,
            done$events[paste0("detached_", stage)] == 5)
}
passed("full 3/5/10-day blood meals; no early detachment")

x <- new_cohort("L", "engorged", 10, 1)
res <- advance(x, p)
stopifnot(all(res$cohorts$stage == "N"), all(res$cohorts$state == "hardening"),
          res$events["new_nymphs"] == 10, res$events["new_females"] == 0,
          all(res$cohorts$eligible_year == 2))
res <- advance(res$cohorts, p, p$hardening_hours["N"] / 2)
stopifnot(all(res$cohorts$dormant), all(res$cohorts$state == "questing"))
still_sleeping <- advance(res$cohorts, p, year = 1L, month = 8L)
stopifnot(all(still_sleeping$cohorts$dormant))
awake <- advance(res$cohorts, p, year = 2L, month = 5L, cdd = 58)
stopifnot(!any(awake$cohorts$dormant))
passed("one molt per step; new nymphs wait until the next year")

x <- new_cohort("N", "engorged", 10, 1)
blocked <- advance(x, p, cdd = 80)
stopifnot(all(blocked$cohorts$stage == "N"))
molted <- advance(x, p, cdd = 81)
stopifnot(all(molted$cohorts$stage == "A"), all(molted$cohorts$state == "hardening"),
          molted$events["eggs_laid"] == 0)
p_fail <- p
p_fail$molt_success["N"] <- 0
failed <- advance(x, p_fail)
stopifnot(nrow(failed$cohorts) == 0, failed$events["deaths_molting"] == 10)
p_male <- p
p_male$female_fraction <- 0
male <- advance(x, p_male)
stopifnot(nrow(male$cohorts) == 0, male$events["males_removed"] == 10)
passed("nymphal temperature gate, single molting-success draw, male accounting")

x <- new_cohort("E", "egg", 10, 1)
stopifnot(advance(x, p, cdd = 109)$cohorts$stage == "E")
hatched <- advance(x, p, cdd = 110)
stopifnot(hatched$cohorts$state == "hardening", hatched$cohorts$timer_h == 0)
before_hardened <- advance(hatched$cohorts, p, 7 * 24 / p$dt_hours - 1)
stopifnot(before_hardened$cohorts$state == "hardening")
hardened <- advance(before_hardened$cohorts, p)
stopifnot(hardened$cohorts$state == "questing", !hardened$cohorts$dormant)
passed("egg threshold and exactly 7 days of larval hardening")

x <- new_cohort("A", "preoviposition", 2, 1)
before <- advance(x, p, p$preov_hours / 2 - 1)
stopifnot(all(before$cohorts$state == "preoviposition"))
birth <- advance(before$cohorts, p)
expected_eggs <- 2 * round(-24.6 * 18^2 + 836 * 18 - 4106)
stopifnot(sum(birth$cohorts$n) == expected_eggs, all(birth$cohorts$stage == "E"),
          birth$events["mothers_removed"] == 2, birth$events["larvae_hatched"] == 0)
zero_birth <- advance(before$cohorts, p, temp = 0)
stopifnot(nrow(zero_birth$cohorts) == 0, zero_birth$events["eggs_laid"] == 0,
          zero_birth$events["fecundity_clipped_females"] == 2)
passed("17-day preoviposition clock, fecundity, mother removal, no instant hatching")

p_groom <- p
p_groom$grooming_hourly[] <- 1
p_groom$grooming_survival <- 1
p_groom$attachment_hourly[] <- 1
p_groom$contact_fraction[] <- 1
x <- new_cohort("L", "feeding", 10, 1, timer_h = 70)
groom <- advance(x, p_groom)
stopifnot(all(groom$cohorts$state == "questing"), all(groom$cohorts$timer_h == 0),
          groom$events["detached_L"] == 0, groom$events["attached_L"] == 0)
p_groom$grooming_survival <- 0
dead <- advance(x, p_groom)
stopifnot(nrow(dead$cohorts) == 0, dead$events["deaths_grooming"] == 10)
passed("grooming resets the meal and cannot reattach within the same interval")

p_dead <- p
p_dead$death_unfed_young[] <- 1
dead <- advance(x, p_dead)
stopifnot(nrow(dead$cohorts) == 0, dead$events["deaths_natural"] == 10,
          dead$events["detached_L"] == 0)
p_old <- p
p_old$death_unfed_old[] <- 1
old <- advance(new_cohort("N", "questing", 10, 1,
                           age_h = p$old_age_hours), p_old)
young <- advance(new_cohort("N", "questing", 10, 1,
                             age_h = p$old_age_hours - 2), p_old)
stopifnot(nrow(old$cohorts) == 0, sum(young$cohorts$n) == 10)
passed("death prevents detachment; 40-week age switch")

p_attach <- p
p_attach$attachment_hourly[] <- 1
p_attach$contact_fraction[] <- 1
p_attach$cell_area_m2 <- 10000
p_attach$mice_per_ha <- 1
p_attach$capacity_per_mouse[] <- 5
x <- rbind(new_cohort("L", "questing", 10, 1),
            new_cohort("L", "questing", 10, 1, age_h = 24))
limited <- advance(x, p_attach)
stopifnot(sum(limited$cohorts$n[limited$cohorts$state == "feeding"]) == 5,
          sum(limited$cohorts$n) == 20, limited$events["capacity_limited"] == 15)
p_attach$mice_per_ha <- 0
none <- advance(x, p_attach)
stopifnot(none$events["attached_L"] == 0, all(none$cohorts$state == "questing"))
p_attach$deer_equivalents_per_ha <- 0
none <- advance(new_cohort("A", "questing", 10, 1), p_attach)
stopifnot(none$events["attached_A"] == 0)
passed("capacity conserves questers; no host means no attachment")

x <- new_cohort("A", "preoviposition", 2, 1, timer_h = 12)
cold <- advance(x, p, year = 1L, month = 10L, cdd = 20)
stopifnot(cold$cohorts$dormant, cold$cohorts$state == "preoviposition",
          cold$cohorts$timer_h == 12, cold$cohorts$eligible_year == 2)
warm_fall <- advance(cold$cohorts, p, year = 1L, month = 11L, cdd = 120)
stopifnot(warm_fall$cohorts$dormant, warm_fall$cohorts$timer_h == 12)
spring <- advance(cold$cohorts, p, year = 2L, month = 4L, cdd = 30)
stopifnot(!spring$cohorts$dormant, spring$cohorts$timer_h == 14)
feeding_cold <- advance(new_cohort("A", "feeding", 2, 1), p,
                         year = 1L, month = 10L, cdd = 20)
stopifnot(feeding_cold$cohorts$timer_h == 2, !feeding_cold$cohorts$dormant)
passed("dormancy preserves state/timer; active feeding continues through cold weather")

weather <- data.frame(model_calendar(152:181),
                       cdd7_raw_Cday = 120, tmean17_C = 18)
r1 <- run_test_cell(weather, default_tick_parameters, seed = 321, verbose = FALSE)
r2 <- run_test_cell(weather, default_tick_parameters, seed = 321, verbose = FALSE)
stopifnot(identical(r1$daily, r2$daily), all(r1$daily$balance_error == 0),
          all(r1$daily$total_ticks == rowSums(r1$daily[c("eggs","L_total","N_total","A_total")])) )
for (stage in c("L", "N", "A")) {
  components <- paste0(stage, c("_questing", "_feeding", "_dormant", "_developing"))
  stopifnot(all(r1$daily[[paste0(stage,"_total")]] == rowSums(r1$daily[components])))
}
bad_weather <- weather[-2, ]
stopifnot(inherits(try(run_test_cell(bad_weather, verbose=FALSE), silent=TRUE), "try-error"))
passed("seed reproducibility, exclusive state counts, simulation-day continuity and mass balance")
# Reporting must exclude exactly ten complete years while keeping chronological
# continuity across fixed 365-day cycles and leaving simulated states untouched.
calendar <- model_calendar(seq_len(13L * 365L))
calendar$sentinel_state <- seq_len(nrow(calendar))
periods <- label_test_periods(calendar, 10, 3)
csv_calendar <- read.csv(text = paste(capture.output(write.csv(calendar, row.names = FALSE)), collapse = "\n"))
stopifnot(nrow(label_test_periods(csv_calendar, 10, 3)) == nrow(calendar))
burn <- periods[periods$phase == "burn_in", ]
analysis <- periods[periods$phase == "analysis", ]
stopifnot(nrow(burn) == 3650, nrow(analysis) == 1095,
          isTRUE(all.equal(unique(analysis$analysis_year), 1:3)),
          isTRUE(all.equal(unique(analysis$simulation_year), 11:13)),
          all(is.na(burn$analysis_year)), all(is.na(burn$analysis_day)),
          analysis$analysis_day[1] == 1, tail(analysis$analysis_day, 1) == 1095,
          analysis$simulation_day[1] == tail(burn$simulation_day, 1) + 1,
          analysis$sentinel_state[1] == tail(burn$sentinel_state, 1) + 1,
          identical(periods$sentinel_state, calendar$sentinel_state),
          !any(periods$month_day == "02-29"),
          all(table(periods$simulation_year) == 365L),
          all(periods$day_of_year[periods$month_day == "03-01"] == 60L),
          all(periods$days_in_year == 365L))
no_burn <- label_test_periods(calendar, 0, 13)
stopifnot(all(no_burn$phase == "analysis"), no_burn$analysis_year[1] == 1)
passed("inline burn-in labels, 365-day years and uninterrupted reporting boundary")

# Reject legacy leap-day calendars and inconsistent seasonal labels.
bad_calendar <- calendar
bad_calendar$month_day[60] <- "02-29"
stopifnot(inherits(try(label_test_periods(bad_calendar, 10, 3), silent = TRUE), "try-error"))
bad_calendar <- calendar
bad_calendar$simulation_year[366] <- 1L
stopifnot(inherits(try(label_test_periods(bad_calendar, 10, 3), silent = TRUE), "try-error"))

# A sleeping nymph wakes only after the model year changes, with no timer reset.
boundary_weather <- data.frame(model_calendar(365:366), cdd7_raw_Cday = 120, tmean17_C = 18)
boundary <- run_test_cell(boundary_weather, p, verbose = FALSE,
                          initial = new_cohort("N", "questing", 10, 2L, dormant = TRUE))
stopifnot(identical(boundary$daily$N_dormant, c(10, 0)),
          identical(boundary$daily$N_questing, c(0, 10)),
          boundary$final_cohorts$age_h == 48)
passed("seasonal calendar validation and dormancy across the 365-day year boundary")
cat("All tick-core behavioral tests passed. These are software checks, not field validation.\n")
