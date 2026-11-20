# One-cell, stochastic cohort implementation of Halsey-style tick biology.
# Source this file to define functions only. The runnable example is code04.
# Source: Halsey & Miller (2018), doi:10.1016/j.ecolmodel.2018.09.005,
# sections 2.6/2.9, Table 1 and supplementary Table A.1.
# This is an ADAPTATION, not a verified port of the unavailable NetLogo source.
# Counts are integers; ticks with identical state/history share one cohort.
# Mice have no agents or demographic compartments. Their fixed abundance only
# supplies encounter opportunities and an aggregate feeding capacity.

halsey_parameters <- function() {
  list(
    dt_hours = 2L,
    # Halsey Table 1. Interpret /h entries as hourly Bernoulli probabilities,
    # not instantaneous hazards. Hourly death/grooming are processed as two
    # successive substeps so a tick cannot both die and finish a blood meal.
    attachment_hourly = c(L = 0.103726853, N = 0.16788807, A = 0.011786875),
    death_egg_hourly = 0.00003017683,
    death_unfed_young = c(L = 0.000213472, N = 0.000003159, A = 0.000000183),
    death_unfed_old = c(L = 0.000517143, N = 0.000363855, A = 0.000360880),
    death_replete = c(L = 0.0001317776, N = 0.0000958566, A = 0.0000912404),
    old_age_hours = 40 * 7 * 24,
    feeding_hours = c(L = 72, N = 120, A = 240),
    grooming_hourly = c(L = 0.00734, N = 0.00734, A = 0.00495),
    grooming_survival = 0.2,
    molt_success = c(L = 0.415, N = 0.276),
    female_fraction = 0.5,
    # Table A.1 specifies a 3-week larval transition, whereas the main text
    # says 1 week. Use the more detailed appendix for new larvae; expose it
    # here for sensitivity analysis. Nymph/adult hardening: main text.
    hardening_hours = c(L = 21 * 24, N = 7 * 24, A = 14 * 24),
    awaken_cdd7 = c(L = 110, N = 58, A = 28),
    molt_cdd7 = c(L = 58, N = 81),
    egg_hatch_cdd7 = 110,
    preov_cdd7 = 28,
    preov_hours = 17 * 24,
    autumn_cdd7 = 58,
    autumn_start_month = 9L,
    # No extra egg-development or engorged-nymph delay is invented: the
    # paper does not specify their duration clearly. Temperature gates and
    # post-molt hardening are used. These are priorities for source-code review.
    fecundity_coefficients = c(quadratic = -24.6, linear = 836, constant = -4106),
    # PROVISIONAL exposure multipliers, NOT values estimated by Halsey.
    # Published attachment probabilities are conditional on nearby hosts.
    # Here these fractions approximate that proximity for a homogeneous cell.
    # Do not copy the deer fraction into the later passage-based ABM.
    contact_fraction = c(L = 0.02, N = 0.02, A = 0.10),
    # Paper's reference densities; held constant, not local abundance estimates.
    mice_per_ha = 50,
    deer_equivalents_per_ha = 2,
    cell_area_m2 = 900,
    # Implement the paper's juvenile limits separately for larvae and nymphs,
    # following the separate larvalthreshold/nymphalthreshold in Table A.1.
    capacity_per_mouse = c(L = 50, N = 50),
    capacity_per_deer = 100,
    initial_nymphs_per_ha = 1500,
    # Adult compartment contains FEMALES only. For this diagnostic interpret
    # the paper's initial 50 adults/ha as 50 reproductive females/ha. This
    # explicit initial-condition choice can be changed independently of the
    # 50:50 sex allocation applied to newly molted adults.
    initial_females_per_ha = 50,
    max_cohorts = 100000L,
    max_population = 1e8
  )
}

empty_cohorts <- function() {
  data.frame(stage = character(), state = character(), n = numeric(),
             age_h = numeric(), timer_h = numeric(), eligible_year = integer(),
             dormant = logical(), stringsAsFactors = FALSE)
}

new_cohort <- function(stage, state, n, year, dormant = FALSE,
                       age_h = 0, timer_h = 0) {
  data.frame(stage = stage, state = state, n = as.numeric(n), age_h = age_h,
             timer_h = timer_h, eligible_year = as.integer(year),
             dormant = dormant, stringsAsFactors = FALSE)
}

merge_cohorts <- function(x) {
  x <- x[x$n > 0, , drop = FALSE]
  if (!nrow(x)) return(empty_cohorts())
  # Age only matters until the 40-week switch (capped in advance_ticks).
  # Egg age is unused by this threshold-based hatch model.
  x$age_h[x$stage == "E"] <- 0
  keys <- paste(x$stage, x$state, x$age_h, x$timer_h, x$eligible_year, x$dormant)
  groups <- match(keys, unique(keys))
  out <- x[!duplicated(keys), , drop = FALSE]
  out$n <- as.numeric(rowsum(x$n, groups, reorder = FALSE))
  rownames(out) <- NULL
  out
}

validate_tick_parameters <- function(p) {
  # Check user-editable settings once, before running the model.
  stopifnot(p$dt_hours == 2L, p$cell_area_m2 > 0,
            p$mice_per_ha >= 0, p$deer_equivalents_per_ha >= 0,
            p$old_age_hours > 0, p$preov_hours > 0,
            p$initial_nymphs_per_ha >= 0, p$initial_females_per_ha >= 0)
  probability_names <- c("attachment_hourly", "death_egg_hourly", "death_unfed_young",
                         "death_unfed_old", "death_replete", "grooming_hourly",
                         "grooming_survival", "molt_success", "female_fraction",
                         "contact_fraction")
  for (name in probability_names) {
    value <- p[[name]]
    if (!length(value) || any(!is.finite(value)) || any(value < 0 | value > 1)) {
      stop("Invalid probability: ", name)
    }
  }
  stopifnot(all(p$feeding_hours > 0), all(p$hardening_hours >= 0),
            all(p$capacity_per_mouse >= 0), p$capacity_per_deer >= 0,
            all(p$feeding_hours %% p$dt_hours == 0),
            all(p$hardening_hours %% p$dt_hours == 0),
            p$preov_hours %% p$dt_hours == 0)
  invisible(TRUE)
}

initial_tick_cohorts <- function(p, start_year) {
  # Initialization assumption: overwintering questers begin with age_h = 0.
  # Their real prior age is unknown; treat early years as initialization effects.
  # Custom initial cohorts can instead be supplied to run_tick_cell(initial = ...).
  area_ha <- p$cell_area_m2 / 10000
  merge_cohorts(rbind(
    new_cohort("N", "questing", round(p$initial_nymphs_per_ha * area_ha),
               start_year, dormant = TRUE),
    new_cohort("A", "questing", round(p$initial_females_per_ha * area_ha),
               start_year, dormant = TRUE)
  ))
}

empty_tick_events <- function() {
  c(eggs_laid = 0, larvae_hatched = 0, new_nymphs = 0, new_females = 0,
    attached_L = 0, attached_N = 0, attached_A = 0,
    detached_L = 0, detached_N = 0, detached_A = 0,
    deaths_natural = 0, deaths_grooming = 0, deaths_molting = 0,
    males_removed = 0, mothers_removed = 0, groomed_survivors = 0,
    capacity_limited = 0, fecundity_clipped_females = 0)
}

# Select a uniform subset of successful attachments if host capacity is full.
# Multivariate hypergeometric allocation avoids favoring whichever cohort
# appears first in the table. No questing tick is killed by a capacity limit.
limit_attachments <- function(proposals, slots) {
  slots <- floor(max(0, slots))
  if (sum(proposals) <= slots) return(proposals)
  out <- numeric(length(proposals))
  remaining <- sum(proposals)
  for (j in seq_along(proposals)) {
    out[j] <- if (slots > 0) rhyper(1, proposals[j], remaining - proposals[j], slots) else 0
    slots <- slots - out[j]
    remaining <- remaining - proposals[j]
  }
  out
}

advance_ticks <- function(cohorts, date, cdd7, tmean17, p) {
  # Each call advances [time, time+2h). Entries created at the END of this
  # interval cannot feed, hatch, or molt again within the same interval.
  # Hourly mortality/grooming substeps preserve event ordering over two hours.
  events <- empty_tick_events()
  before <- sum(cohorts$n)
  if (!nrow(cohorts)) return(list(cohorts = cohorts, events = events, balance_error = 0))
  year <- as.integer(format(date, "%Y"))
  month <- as.integer(format(date, "%m"))
  x <- cohorts
  # Dormancy is an independent flag so overwintering never erases fed status
  # or accumulated developmental time. Dormant ticks still experience mortality.
  # Adaptation: already attached ticks finish feeding through cold weather.
  autumn <- month >= p$autumn_start_month && cdd7 < p$autumn_cdd7
  sleep <- x$stage != "E" & x$state != "feeding" & autumn
  x$dormant[sleep] <- TRUE
  x$eligible_year[sleep] <- pmax(x$eligible_year[sleep], year + 1L)
  threshold <- unname(p$awaken_cdd7[x$stage])
  wake <- x$stage != "E" & x$dormant & x$eligible_year <= year &
    !autumn & !is.na(threshold) & cdd7 >= threshold
  x$dormant[wake] <- FALSE
  # Canonicalize expired year restrictions to allow equivalent cohorts to merge.
  x$eligible_year[x$eligible_year < year] <- year
  old_state <- x$state
  old_stage <- x$stage
  active <- !x$dormant
  added <- list()

  for (hour in seq_len(p$dt_hours)) {
    # Age means time in current life stage, reset at a stage transition.
    # Feeding ticks use unfed-stage mortality until their meal is completed.
    rate <- unname(p$death_unfed_young[x$stage])
    old <- x$age_h >= p$old_age_hours
    rate[old] <- unname(p$death_unfed_old[x$stage[old]])
    replete <- x$state %in% c("engorged", "preoviposition")
    rate[replete] <- unname(p$death_replete[x$stage[replete]])
    rate[x$stage == "E"] <- p$death_egg_hourly
    died <- rbinom(nrow(x), x$n, rate)
    x$n <- x$n - died
    events["deaths_natural"] <- events["deaths_natural"] + sum(died)
    x$age_h <- pmin(x$age_h + 1, p$old_age_hours)

    feeding <- which(old_state == "feeding" & x$n > 0)
    if (length(feeding)) {
      removed <- rbinom(length(feeding), x$n[feeding], p$grooming_hourly[x$stage[feeding]])
      survivors <- rbinom(length(feeding), removed, p$grooming_survival)
      x$n[feeding] <- x$n[feeding] - removed
      events["deaths_grooming"] <- events["deaths_grooming"] + sum(removed - survivors)
      events["groomed_survivors"] <- events["groomed_survivors"] + sum(survivors)
      returned <- x[feeding[survivors > 0], , drop = FALSE]
      if (nrow(returned)) {
        returned$n <- survivors[survivors > 0]
        # Partial meals are discarded. Reattachment requires a new full meal.
        returned$state <- "questing"
        returned$timer_h <- 0
        returned$dormant <- FALSE
        # Surviving groomed ticks are returned at interval end. They cannot
        # reattach immediately; apply remaining natural mortality in this step.
        remaining_h <- p$dt_hours - hour
        if (remaining_h > 0) {
          return_rate <- ifelse(returned$age_h >= p$old_age_hours,
                                p$death_unfed_old[returned$stage],
                                p$death_unfed_young[returned$stage])
          return_dead <- rbinom(nrow(returned), returned$n,
                               1 - (1 - return_rate)^remaining_h)
          returned$n <- returned$n - return_dead
          events["deaths_natural"] <- events["deaths_natural"] + sum(return_dead)
          returned$age_h <- pmin(returned$age_h + remaining_h, p$old_age_hours)
        }
        added[[length(added) + 1L]] <- returned
      }
    }
  }

  # Timers only advance for cohorts already in the state at interval start.
  timed <- active & old_state %in% c("hardening", "feeding", "preoviposition")
  x$timer_h[timed] <- x$timer_h[timed] + p$dt_hours

  # Attachment is an interval-end event and applies ONLY to active questers.
  # Static host availability is not a host population simulation.
  area_ha <- p$cell_area_m2 / 10000
  capacities <- c(p$capacity_per_mouse * p$mice_per_ha * area_ha,
                  A = p$capacity_per_deer * p$deer_equivalents_per_ha * area_ha)
  host_present <- c(L = p$mice_per_ha > 0, N = p$mice_per_ha > 0,
                    A = p$deer_equivalents_per_ha > 0)
  for (stage in c("L", "N", "A")) {
    idx <- which(x$stage == stage & old_state == "questing" & active & x$n > 0)
    if (!length(idx) || !host_present[stage]) next
    effective_hourly <- p$attachment_hourly[stage] * p$contact_fraction[stage]
    proposed <- rbinom(length(idx), x$n[idx], 1 - (1 - effective_hourly)^p$dt_hours)
    occupied <- sum(x$n[x$stage == stage & old_state == "feeding"])
    accepted <- limit_attachments(proposed, capacities[stage] - occupied)
    events["capacity_limited"] <- events["capacity_limited"] + sum(proposed - accepted)
    events[paste0("attached_", stage)] <- sum(accepted)
    attached <- x[idx[accepted > 0], , drop = FALSE]
    if (nrow(attached)) {
      attached$n <- accepted[accepted > 0]
      attached$state <- "feeding"
      attached$timer_h <- 0
      x$n[idx] <- x$n[idx] - accepted
      added[[length(added) + 1L]] <- attached
    }
  }

  complete <- old_state == "feeding" & x$timer_h >= p$feeding_hours[x$stage] & x$n > 0
  complete[is.na(complete)] <- FALSE
  for (stage in c("L", "N", "A")) {
    events[paste0("detached_", stage)] <- sum(x$n[complete & x$stage == stage])
  }
  x$state[complete] <- "engorged"
  x$timer_h[complete] <- 0
  # All detachment is local in this one-cell check. The spatial extension must
  # add deer_id and deposit adult cohorts into that deer's current cell.

  hatched <- old_state == "egg" & cdd7 >= p$egg_hatch_cdd7 & x$n > 0
  events["larvae_hatched"] <- sum(x$n[hatched])
  x$stage[hatched] <- "L"
  x$state[hatched] <- "hardening"
  x$age_h[hatched] <- 0
  x$timer_h[hatched] <- 0
  x$dormant[hatched] <- FALSE

  # Molting success is a single transition trial, never a repeated hourly draw.
  for (stage in c("L", "N")) {
    idx <- which(old_stage == stage & old_state == "engorged" & active &
                   cdd7 >= p$molt_cdd7[stage] & x$n > 0)
    if (!length(idx)) next
    successful <- rbinom(length(idx), x$n[idx], p$molt_success[stage])
    events["deaths_molting"] <- events["deaths_molting"] + sum(x$n[idx] - successful)
    if (stage == "N") {
      females <- rbinom(length(idx), successful, p$female_fraction)
      events["males_removed"] <- events["males_removed"] + sum(successful - females)
      successful <- females
    }
    x$n[idx] <- successful
    x$stage[idx] <- if (stage == "L") "N" else "A"
    x$state[idx] <- "hardening"
    x$age_h[idx] <- 0
    x$timer_h[idx] <- 0
    if (stage == "L") {
      # New nymphs must not quest in the same year they fed as larvae.
      # They harden first, then wait until the following calendar year.
      x$eligible_year[idx] <- year + 1L
      events["new_nymphs"] <- sum(successful)
    } else {
      events["new_females"] <- sum(successful)
    }
  }

  hardened <- old_state == "hardening" & active &
    x$timer_h >= p$hardening_hours[x$stage] & x$n > 0
  hardened[is.na(hardened)] <- FALSE
  x$state[hardened] <- "questing"
  x$timer_h[hardened] <- 0
  x$dormant[hardened] <- x$eligible_year[hardened] > year

  start_preov <- old_state == "engorged" & x$stage == "A" & active &
    cdd7 >= p$preov_cdd7 & x$n > 0
  # Exclude nymphs that became adults during THIS step.
  start_preov <- start_preov & old_stage == "A"
  x$state[start_preov] <- "preoviposition"
  x$timer_h[start_preov] <- 0

  lay <- old_state == "preoviposition" & active & x$timer_h >= p$preov_hours & x$n > 0
  if (any(lay)) {
    cf <- p$fecundity_coefficients
    fecundity_raw <- cf["quadratic"] * tmean17^2 + cf["linear"] * tmean17 + cf["constant"]
    # Explicit safeguard: a polynomial outside its viable range cannot produce
    # negative eggs. Round once per female; do not treat fecundity as a rate.
    per_female <- round(max(0, fecundity_raw))
    mothers <- sum(x$n[lay])
    events["eggs_laid"] <- mothers * per_female
    events["mothers_removed"] <- mothers
    if (fecundity_raw < 0) events["fecundity_clipped_females"] <- mothers
    added[[length(added) + 1L]] <- new_cohort("E", "egg", mothers * per_female, year)
    x$n[lay] <- 0
  }

  out <- merge_cohorts(do.call(rbind, c(list(x), added)))
  losses <- sum(events[c("deaths_natural", "deaths_grooming", "deaths_molting",
                         "males_removed", "mothers_removed")])
  balance <- sum(out$n) - (before + events["eggs_laid"] - losses)
  if (balance != 0) stop("Tick births and losses do not balance in this step: ", balance)
  if (nrow(out) > p$max_cohorts || sum(out$n) > p$max_population) {
    stop("Diagnostic population/cohort limit exceeded; inspect parameters rather than truncate ticks.")
  }
  list(cohorts = out, events = events, balance_error = unname(balance))
}

tick_snapshot <- function(x) {
  result <- c(eggs = sum(x$n[x$stage == "E"]))
  for (stage in c("L", "N", "A")) {
    s <- x$stage == stage
    result[paste0(stage, "_total")] <- sum(x$n[s])
    result[paste0(stage, "_questing")] <- sum(x$n[s & x$state == "questing" & !x$dormant])
    result[paste0(stage, "_feeding")] <- sum(x$n[s & x$state == "feeding"])
    result[paste0(stage, "_dormant")] <- sum(x$n[s & x$dormant])
    result[paste0(stage, "_developing")] <- sum(x$n[s & !x$dormant &
      x$state %in% c("hardening", "engorged", "preoviposition")])
  }
  c(result, total_ticks = sum(x$n), cohorts = nrow(x))
}

# Assign reporting periods AFTER a single continuous biological simulation.
# Never restart ticks or RNG at the burn-in boundary. Require complete years
# so a partial climate file cannot silently shorten either requested period.
label_tick_periods <- function(daily, burn_in_years, analysis_years) {
  stopifnot(length(burn_in_years) == 1L, is.finite(burn_in_years),
            burn_in_years >= 0, burn_in_years == floor(burn_in_years),
            length(analysis_years) == 1L, is.finite(analysis_years),
            analysis_years >= 1, analysis_years == floor(analysis_years))
  dates <- as.Date(daily$date)
  if (!length(dates) || anyNA(dates)) stop("Missing simulation dates.")
  year0 <- as.integer(format(dates[1], "%Y"))
  expected <- seq(as.Date(sprintf("%04d-01-01", year0)),
                  as.Date(sprintf("%04d-01-01", year0 + burn_in_years + analysis_years)) - 1,
                  by = "day")
  # CSV-parsed Date values and seq.Date may use different numeric storage.
  if (length(dates) != length(expected) || any(dates != expected)) {
    stop("Need complete consecutive years matching burn-in + analysis duration.")
  }
  years <- as.integer(format(dates, "%Y"))
  daily$simulation_day <- seq_along(dates)
  daily$simulation_year <- years - year0 + 1L
  daily$day_of_year <- as.integer(format(dates, "%j"))
  daily$month_day <- format(dates, "%m-%d")
  daily$days_in_year <- as.integer(as.Date(sprintf("%04d-01-01", years + 1L)) -
                                   as.Date(sprintf("%04d-01-01", years)))
  daily$phase <- ifelse(daily$simulation_year <= burn_in_years, "burn_in", "analysis")
  daily$analysis_year <- ifelse(daily$phase == "analysis",
                                daily$simulation_year - burn_in_years, NA_integer_)
  daily$analysis_day <- NA_integer_
  daily$analysis_day[daily$phase == "analysis"] <- seq_len(sum(daily$phase == "analysis"))
  daily
}

run_tick_cell <- function(weather, p = halsey_parameters(), seed = 123L,
                           initial = NULL, verbose = TRUE) {
  validate_tick_parameters(p)
  required <- c("date", "cdd7_raw_Cday", "tmean17_C")
  if (!all(required %in% names(weather))) stop("Missing required daily climate columns.")
  weather$date <- as.Date(weather$date)
  if (!nrow(weather) || anyNA(weather$date) || any(diff(weather$date) != 1) ||
      any(!is.finite(weather$cdd7_raw_Cday)) || any(!is.finite(weather$tmean17_C))) {
    stop("Daily climate must be consecutive with complete rolling indices.")
  }
  set.seed(seed)
  year0 <- as.integer(format(weather$date[1], "%Y"))
  x <- if (is.null(initial)) initial_tick_cohorts(p, year0) else initial
  # Detailed state/age/timer checks live in tests/test_halsey_tick_core.R.
  # Retain the count check here because callers can supply custom initial ticks.
  if (any(!is.finite(x$n)) || any(x$n < 0 | x$n != floor(x$n))) {
    stop("Initial tick counts must be nonnegative whole numbers.")
  }
  initial_count <- sum(x$n)
  records <- vector("list", nrow(weather))
  for (day in seq_len(nrow(weather))) {
    event_sum <- empty_tick_events()
    balance_error <- 0
    for (step in seq_len(24L / p$dt_hours)) {
      result <- advance_ticks(x, weather$date[day], weather$cdd7_raw_Cday[day],
                               weather$tmean17_C[day], p)
      x <- result$cohorts
      event_sum <- event_sum + result$events
      balance_error <- max(balance_error, abs(result$balance_error))
    }
    records[[day]] <- c(tick_snapshot(x), event_sum, balance_error = balance_error)
    if (verbose && (format(weather$date[day], "%m-%d") == "12-31" || day == nrow(weather))) {
      cat("Simulation year", as.integer(format(weather$date[day], "%Y")) - year0 + 1L,
          format(weather$date[day], "%m-%d"), ":", sum(x$n), "ticks in", nrow(x), "cohorts\n")
    }
  }
  daily <- data.frame(date = weather$date, do.call(rbind, records), row.names = NULL)
  daily$cdd7 <- weather$cdd7_raw_Cday
  daily$tmean17_C <- weather$tmean17_C
  list(daily = daily, final_cohorts = x, initial_count = initial_count,
       parameters = p, seed = seed)
}
