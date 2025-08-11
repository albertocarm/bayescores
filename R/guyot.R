#' Helper to Place Censored Observations
#'
#' @description
#' Distributes a given number of censored observations within a time interval `[a, b]`.
#'
#' @param a (`numeric`)\cr Start of the interval.
#' @param b (`numeric`)\cr End of the interval.
#' @param ci (`integer`)\cr Number of censored observations to place.
#' @param eps (`numeric`)\cr A small epsilon value for numerical comparisons.
#' @param mode (`character`)\cr Method for placing censored observations, either "even"
#'   (evenly spaced) or "uniform" (from a uniform distribution).
#'
#' @return A numeric vector of censored times.
#' @keywords internal
.place_censors <- function(a, b, ci, eps = 1e-6, mode = c("even", "uniform")) {
  mode <- match.arg(mode)
  if (ci <= 0L) return(numeric(0))
  if (b <= a + eps) return(rep(a + eps, ci))
  if (mode == "even") {
    seq(a + (b - a) / (ci + 1), b - (b - a) / (ci + 1), length.out = ci)
  } else {
    sort(runif(ci, min = a, max = b))
  }
}


#' Reconstruct IPD for a Single Curve Using Guyot's Algorithm
#'
#' @description
#' This function takes digitized Kaplan-Meier (KM) curve data and a table of the
#' number of patients at risk at specific time points to reconstruct individual
#' patient data (IPD) for a single treatment arm or group. It implements the
#' algorithm described by Guyot et al. (2012).
#'
#' @param d_curve (`data.frame`)\cr The digitized KM curve points. Must contain 'time'
#'   and 'St' (survival probability) columns.
#' @param nr_curve (`data.frame`)\cr The number of patients at risk. Must contain
#'   'time_tick' and 'nrisk' columns.
#' @param totev (`integer` or `NULL`)\cr Optional total number of events. If
#'   provided, the algorithm will adjust the event count in the last interval to
#'   match this total.
#' @param eps (`numeric`)\cr A small epsilon value for numerical comparisons.
#' @param ensure_monotone (`logical`)\cr If `TRUE`, forces the survival curve to be
#'   monotonically decreasing.
#' @param censor_placement (`character`)\cr Method for placing censored observations,
#'   "even" or "uniform".
#' @param jitter_admin (`numeric`)\cr A small time value to add to the last time
#'   point for administrative censoring.
#' @param check (`logical`)\cr If `TRUE`, performs a check to verify if the
#'   reconstructed number at risk matches the published numbers.
#'
#' @return
#' A list containing: `ipd` (the reconstructed IPD `data.frame`), `drops`
#' (details of survival probability drops), `check` (the verification `data.frame`),
#' and `compat` (a logical indicating if the check was successful).
#'
#' @references
#' Guyot, P., Ades, A. E., Ouwens, M. J., & Welton, N. J. (2012).
#' Enhanced secondary analysis of survival data: reconstructing the data from
#' published Kaplan-Meier survival curves. *BMC Medical Research Methodology*, 12(1), 1-13.
#'
#' @keywords internal
recon_one_curve_guyot <- function(d_curve,
                                  nr_curve,
                                  totev = NULL,
                                  eps = 1e-6,
                                  ensure_monotone = TRUE,
                                  censor_placement = c("even", "uniform"),
                                  jitter_admin = 1e-6,
                                  check = TRUE) {
  censor_placement <- match.arg(censor_placement)
  stopifnot(all(c("time", "St") %in% names(d_curve)))
  stopifnot(all(c("time_tick", "nrisk") %in% names(nr_curve)))

  # --- Normalize/clean digitized KM data ---
  Smax <- max(d_curve$St, na.rm = TRUE)
  if (Smax > 1 + 1e-3) d_curve$St <- d_curve$St / 100
  d_curve$St <- pmin(pmax(d_curve$St, 0), 1)
  d_curve <- d_curve[order(d_curve$time, -d_curve$St), , drop = FALSE]
  d_curve <- stats::aggregate(St ~ time, data = d_curve, FUN = min)
  d_curve <- d_curve[order(d_curve$time), , drop = FALSE]
  if (nrow(d_curve) == 0 || d_curve$time[1] > 0 + eps) {
    d_curve <- rbind(data.frame(time = 0, St = 1), d_curve)
  } else {
    d_curve$St[1] <- 1
  }
  if (ensure_monotone) d_curve$St <- cummin(d_curve$St)

  # --- Published number at risk ---
  nr_curve <- nr_curve[order(nr_curve$time_tick), , drop = FALSE]
  Ttick <- nr_curve$time_tick
  Ntick <- as.integer(round(nr_curve$nrisk))
  stopifnot(length(Ttick) >= 2L)
  N0_pub <- Ntick[1L]
  n_last <- Ntick[length(Ntick)]
  nint <- length(Ttick) - 1L

  # --- Step functions ---
  step_S_right <- function(tt, tvec, Svec) {
    idx <- max(which(tvec <= tt))
    if (is.infinite(idx)) 1.0 else Svec[idx]
  }
  step_S_left <- function(tt, tvec, Svec) {
    idxs <- which(tvec < tt)
    if (!length(idxs)) 1.0 else Svec[max(idxs)]
  }

  # containers
  ipd_events <- list()
  ipd_cens <- list()
  out_drops <- list()
  tvec <- d_curve$time
  Svec <- d_curve$St

  # --- Guyot loop over intervals [t_i, t_{i+1}) ---
  for (i in seq_len(nint)) {
    t_lo <- Ttick[i]
    t_hi <- Ttick[i + 1]
    n_lo <- Ntick[i]
    n_hi <- Ntick[i + 1]
    if (t_hi <= t_lo + eps) next

    # S(t_i+), S(t_{i+1}-)
    S_lo <- step_S_right(t_lo, tvec, Svec)
    S_hi <- step_S_left(t_hi, tvec, Svec)
    S_ratio <- if (S_lo <= eps) 0 else pmin(pmax(S_hi / S_lo, 0), 1)

    # (1) censored observations in the interval from Guyot's formula
    c_i_real <- n_lo * S_ratio - n_hi
    c_i <- as.integer(round(c_i_real))
    c_i <- max(0L, min(c_i, max(0L, n_lo - n_hi))) # feasibility

    # (2) target events for the interval
    E_target <- n_lo - n_hi - c_i
    if (E_target < 0L) {
      c_i <- n_lo - n_hi
      E_target <- 0L
    } # safety bound

    # (3) internal drops within the interval
    dmask <- (d_curve$time > t_lo + eps) & (d_curve$time < t_hi - eps)
    d_int <- d_curve[dmask, , drop = FALSE]
    m <- nrow(d_int)

    if (m == 0L) {
      # no jumps: all c_i are distributed as censored
      if (c_i > 0L) {
        ct <- .place_censors(t_lo, t_hi, c_i, eps = eps, mode = censor_placement)
        ipd_cens[[length(ipd_cens) + 1L]] <- data.frame(time = ct, status = 0L)
      }
      out_drops[[length(out_drops) + 1L]] <- data.frame(
        time = numeric(0), S_prev = numeric(0), S_curr = numeric(0),
        n_i = n_lo, d_i = integer(0), c_i = as.integer(c_i)
      )
      next
    }

    # (4) sequential assignment of KM events with censored observations before each drop
    S_prev <- c(step_S_right(t_lo, tvec, Svec), head(d_int$St, -1))
    S_curr <- d_int$St
    timesJ <- d_int$time

    # distribute c_i into (m+1) sub-intervals
    seg_starts <- c(t_lo, timesJ)
    seg_ends <- c(timesJ, t_hi)
    lens <- pmax(seg_ends - seg_starts, 0)
    if (sum(lens) <= eps) {
      cens_sched <- integer(m + 1)
      cens_sched[1] <- c_i
    } else {
      frac <- lens / sum(lens)
      cand <- c_i * frac
      base <- floor(cand)
      rem <- c_i - sum(base)
      cens_sched <- base
      if (rem > 0L) {
        resid <- cand - base
        add_idx <- order(resid, decreasing = TRUE)[seq_len(rem)]
        cens_sched[add_idx] <- cens_sched[add_idx] + 1L
      }
      cens_sched <- as.integer(cens_sched)
    }

    n_now <- n_lo
    d_star <- numeric(m)
    d_ints <- integer(m)
    for (j in seq_len(m)) {
      n_now <- max(0L, n_now - cens_sched[j]) # censored before drop j
      ratio <- if (S_prev[j] <= eps) 0 else pmin(pmax(S_curr[j] / S_prev[j], 0), 1)
      d_star[j] <- n_now * (1 - ratio) # point-wise KM inversion
      d_ints[j] <- max(0L, min(as.integer(round(d_star[j])), n_now))
      n_now <- n_now - d_ints[j]
    }
    n_now <- max(0L, n_now - cens_sched[m + 1]) # censored after the last drop

    # minimal adjustment of events to match n_hi exactly (without touching c_i)
    D_sum <- sum(d_ints)
    need <- (n_now - n_hi) # >0 subjects left over => increase events; <0 subjects missing => decrease events
    if (need != 0L && m > 0L) {
      frac <- d_star - round(d_star)
      if (need > 0L) {
        # increase events (we prefer those with a larger fraction)
        ord <- order(frac, decreasing = TRUE)
        k <- 1L
        while (need > 0L && k <= m) {
          idx <- ord[k]
          if (d_ints[idx] + 1L <= (n_lo - sum(cens_sched[1:idx]) - sum(d_ints[seq_len(idx - 1)]))) {
            d_ints[idx] <- d_ints[idx] + 1L
            need <- need - 1L
          }
          k <- k + 1L
        }
      } else {
        # decrease events (we prefer those with a smaller fraction)
        ord <- order(frac, decreasing = FALSE)
        k <- 1L
        while (need < 0L && k <= m) {
          idx <- ord[k]
          if (d_ints[idx] > 0L) {
            d_ints[idx] <- d_ints[idx] - 1L
            need <- need + 1L
          }
          k <- k + 1L
        }
      }
    }

    # build IPD for this interval
    if (any(d_ints > 0L)) {
      ipd_events[[length(ipd_events) + 1L]] <- data.frame(time = rep(timesJ, d_ints), status = 1L)
    }
    for (seg in seq_len(m + 1)) {
      ci <- cens_sched[seg]
      if (ci <= 0L) next
      a <- seg_starts[seg]
      b <- seg_ends[seg]
      ct <- .place_censors(a, b, ci, eps = eps, mode = censor_placement)
      ipd_cens[[length(ipd_cens) + 1L]] <- data.frame(time = ct, status = 0L)
    }

    out_drops[[length(out_drops) + 1L]] <- data.frame(
      time = timesJ, S_prev = S_prev, S_curr = S_curr,
      n_i = n_lo, d_i = as.integer(d_ints), c_i = as.integer(c_i)
    )
  } # end intervals loop

  # (optional) force total number of events in the last interval (without touching censored)
  if (!is.null(totev) && length(out_drops) > 0L) {
    dtab <- do.call(rbind, out_drops)
    gap <- as.integer(totev) - sum(dtab$d_i)
    if (gap != 0L) {
      last_idx <- length(out_drops)
      d_last <- out_drops[[last_idx]]
      if (nrow(d_last) > 0L) {
        m <- nrow(d_last)
        ord <- if (gap > 0) rev(seq_len(m)) else seq_len(m)
        for (prop in rep(sign(gap), abs(gap))) {
          for (j in ord) {
            cand <- d_last$d_i[j] + prop
            if (cand >= 0L) {
              d_last$d_i[j] <- cand
              break
            }
          }
        }
        out_drops[[last_idx]] <- d_last
      }
    }
  }

  # --- combine intervals ---
  ipd_int <- rbind(
    if (length(ipd_events)) do.call(rbind, ipd_events) else NULL,
    if (length(ipd_cens)) do.call(rbind, ipd_cens) else NULL
  )
  if (is.null(ipd_int)) ipd_int <- data.frame(time = numeric(0), status = integer(0))

  # --- administrative censoring: add EXACTLY n_last at t_admin > last tick ---
  t_admin <- max(max(d_curve$time, na.rm = TRUE), Ttick[length(Ttick)]) + jitter_admin
  ipd_out <- if (n_last > 0L) rbind(ipd_int, data.frame(time = rep(t_admin, n_last), status = 0L)) else ipd_int

  # --- SORT and TRIM TAIL (censored only) to match N0 ---
  if (nrow(ipd_out)) ipd_out <- ipd_out[order(ipd_out$time, -ipd_out$status), , drop = FALSE]
  excess <- nrow(ipd_out) - N0_pub
  if (excess > 0L) {
    tiny <- 1e-9
    # 1) remove from administrative block (t >= t_admin - tiny)
    idx_admin <- which(ipd_out$status == 0L & ipd_out$time >= t_admin - tiny)
    if (length(idx_admin) > 0L) {
      rm1 <- min(length(idx_admin), excess)
      ipd_out <- ipd_out[-tail(idx_admin, rm1), , drop = FALSE]
      excess <- nrow(ipd_out) - N0_pub
    }
    # 2) if still in excess, remove censored from the last interval backwards
    if (excess > 0L) {
      for (ii in seq(from = length(Ttick) - 1L, to = 1L, by = -1L)) {
        a <- Ttick[ii]
        b <- Ttick[ii + 1]
        idx_i <- which(ipd_out$status == 0L & ipd_out$time >= a - tiny & ipd_out$time < b - tiny)
        if (!length(idx_i)) next
        # remove those with the latest times first
        idx_i <- idx_i[order(ipd_out$time[idx_i], decreasing = TRUE)]
        rm2 <- min(length(idx_i), excess)
        ipd_out <- ipd_out[-idx_i[seq_len(rm2)], , drop = FALSE]
        excess <- nrow(ipd_out) - N0_pub
        if (excess <= 0L) break
      }
    }
    if (nrow(ipd_out) != N0_pub) warning("Could not trim exactly to N0; delta = ", nrow(ipd_out) - N0_pub)
  }

  # drops output
  drops_out <- if (length(out_drops)) {
    do.call(rbind, out_drops)
  } else {
    data.frame(
      time = numeric(0), S_prev = numeric(0), S_curr = numeric(0),
      n_i = integer(0), d_i = integer(0), c_i = integer(0)
    )
  }

  # verification
  nrisk_check <- NULL
  compat_ok <- NA
  if (check) {
    ticks <- Ttick
    n_at_risk <- function(ipd, tt) sum(ipd$time >= tt)
    rec <- vapply(ticks, function(tt) n_at_risk(ipd_out, tt), FUN.VALUE = integer(1))
    nrisk_check <- data.frame(
      time_tick = ticks,
      nrisk_published = Ntick,
      nrisk_reconstructed = as.integer(rec)
    )
    compat_ok <- all(nrisk_check$nrisk_published == nrisk_check$nrisk_reconstructed)
  }

  list(ipd = ipd_out, drops = drops_out, check = nrisk_check, compat = compat_ok)
}


#' Reconstruct Individual Patient Data from Multiple Kaplan-Meier Curves
#'
#' @description
#' A wrapper function that applies the Guyot reconstruction algorithm to multiple
#' curves (e.g., different treatment arms) provided in a single data frame. It
#' splits the data by curve, applies `recon_one_curve_guyot` to each, and then
#' combines the results.
#'
#' @param plot_km (`data.frame`)\cr The digitized KM curve points for one or more
#'   curves. Must contain 'time', 'St', and an optional 'curve' identifier column.
#'   If 'curve' is missing, all data is assumed to belong to a single curve.
#' @param nrisk_tbl (`data.frame`)\cr The number of patients at risk for one or more
#'   curves. Must contain 'time_tick', 'nrisk', and an optional 'curve' identifier
#'   column.
#' @param totev (`integer`, named `list`, or `NULL`)\cr Optional total number of
#'   events for each curve. If a single integer, it's applied to all curves. If a
#'   named list, names should match the 'curve' identifiers.
#' @param ensure_monotone (`logical`)\cr If `TRUE`, forces the survival curves to be
#'   monotonically decreasing.
#' @param censor_placement (`character`)\cr Method for placing censored observations,
#'   "even" or "uniform".
#' @param jitter_admin (`numeric`)\cr A small time value to add to the last time
#'   point for administrative censoring.
#' @param verbose (`logical`)\cr If `TRUE`, prints a summary of the verification
#'   check for each curve.
#'
#' @return
#' A list containing the combined results: `ipd` (`data.frame`), `drops` (`data.frame`),
#' `check` (`data.frame`), and `compat` (`logical`).
#'
#' @importFrom dplyr bind_rows group_by summarise
#' @importFrom stats aggregate
#' @importFrom utils head tail
#' @export
reconstruct_ipd <- function(plot_km, nrisk_tbl, totev = NULL,
                                        ensure_monotone = TRUE,
                                        censor_placement = c("even", "uniform"),
                                        jitter_admin = 1e-6,
                                        verbose = TRUE) {
  censor_placement <- match.arg(censor_placement)
  stopifnot(all(c("time", "St") %in% names(plot_km)))
  if (!("curve" %in% names(plot_km))) plot_km$curve <- 1L
  plot_km$curve <- as.integer(plot_km$curve)

  stopifnot(all(c("time_tick", "nrisk") %in% names(nrisk_tbl)))
  if (!("curve" %in% names(nrisk_tbl))) nrisk_tbl$curve <- 1L
  nrisk_tbl$curve <- as.integer(nrisk_tbl$curve)

  curves <- sort(unique(plot_km$curve))
  ipd_all <- list()
  drops_all <- list()
  checks <- list()
  compats <- c()

  for (k in curves) {
    d_k <- subset(plot_km, curve == k, select = c("time", "St"))
    nr_k <- subset(nrisk_tbl, curve == k, select = c("time_tick", "nrisk"))
    te_k <- if (is.null(totev)) {
      NULL
    } else {
      if (length(totev) == 1L) {
        totev
      } else {
        nm <- names(totev)
        if (length(nm) && as.character(k) %in% nm) totev[[as.character(k)]] else NULL
      }
    }

    out <- recon_one_curve_guyot(d_k, nr_k,
                                 totev = te_k,
                                 ensure_monotone = ensure_monotone,
                                 censor_placement = censor_placement,
                                 jitter_admin = jitter_admin,
                                 check = TRUE
    )

    if (nrow(out$ipd)) {
      out$ipd$curve <- k
      ipd_all[[length(ipd_all) + 1]] <- out$ipd
    }
    if (nrow(out$drops)) {
      out$drops$curve <- k
      drops_all[[length(drops_all) + 1]] <- out$drops
    }
    if (!is.null(out$check)) checks[[length(checks) + 1]] <- transform(out$check, curve = k)
    compats <- c(compats, isTRUE(out$compat))
  }

  res <- list(
    ipd = if (length(ipd_all)) dplyr::bind_rows(ipd_all) else data.frame(time = numeric(0), status = integer(0), curve = integer(0)),
    drops = if (length(drops_all)) dplyr::bind_rows(drops_all) else data.frame(time = numeric(0), S_prev = numeric(0), S_curr = numeric(0), n_i = integer(0), d_i = integer(0), c_i = integer(0), curve = integer(0)),
    check = if (length(checks)) dplyr::bind_rows(checks) else NULL,
    compat = all(compats)
  )
  if (verbose && !is.null(res$check)) {
    print(
      res$check |>
        dplyr::group_by(curve) |>
        dplyr::summarise(ok = all(nrisk_published == nrisk_reconstructed))
    )
  }
  res
}
