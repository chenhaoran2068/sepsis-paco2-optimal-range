# Daily TWA-PaCO2 utilities for the ICU sepsis PaCO2 multicohort study.
#
# Main exposure rule:
#   - ICU day 1 is [ICU admission, ICU admission + 24 h); day 2 is
#     [ICU admission + 24 h, ICU admission + 48 h), and so on through day 7.
#   - The day-window end is truncated at ICU discharge and, when available,
#     death time.
#   - Daily TWA-PaCO2 uses only PaCO2 measurements inside the same day-window.
#   - No cross-day carry-forward is used in the primary analysis.
#   - Within a day-window, values are treated as a step function: the first
#     in-window value extends backward to the window start; each value is
#     carried forward until the next measurement; the last value extends to
#     the window end.

requireNamespace("dplyr")

.require_columns <- function(data, cols, data_name = "data") {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      data_name,
      " is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

.infer_tz <- function(x, default_tz = "UTC") {
  tz <- attr(x, "tzone")
  if (is.null(tz) || length(tz) == 0 || is.na(tz[1]) || identical(tz[1], "")) {
    default_tz
  } else {
    tz[1]
  }
}

.as_posix <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXct")) {
    as.POSIXct(x, tz = tz)
  } else {
    as.POSIXct(x, tz = tz)
  }
}

.pmin_posix_ignore_na <- function(..., tz = "UTC") {
  time_list <- list(...)
  numeric_list <- lapply(time_list, as.numeric)
  numeric_matrix <- do.call(cbind, numeric_list)
  min_numeric <- apply(numeric_matrix, 1, function(row) {
    if (all(is.na(row))) {
      NA_real_
    } else {
      min(row, na.rm = TRUE)
    }
  })
  as.POSIXct(min_numeric, origin = "1970-01-01", tz = tz)
}

make_icu_day_windows <- function(
    stays,
    id_col = "stay_id",
    intime_col = "icu_intime",
    outtime_col = "icu_outtime",
    death_time_col = NULL,
    max_day = 7,
    timezone = NULL) {
  .require_columns(stays, c(id_col, intime_col, outtime_col), "stays")
  if (!is.null(death_time_col)) {
    .require_columns(stays, death_time_col, "stays")
  }
  if (max_day < 1) {
    stop("max_day must be >= 1.", call. = FALSE)
  }

  if (is.null(timezone)) {
    timezone <- .infer_tz(stays[[intime_col]])
  }

  base <- data.frame(
    stay_id = stays[[id_col]],
    icu_intime = .as_posix(stays[[intime_col]], tz = timezone),
    icu_outtime = .as_posix(stays[[outtime_col]], tz = timezone),
    stringsAsFactors = FALSE
  )

  if (!is.null(death_time_col)) {
    base$death_time <- .as_posix(stays[[death_time_col]], tz = timezone)
  } else {
    base$death_time <- as.POSIXct(NA_real_, origin = "1970-01-01", tz = timezone)
  }

  day_index <- seq_len(max_day)
  expanded <- base[rep(seq_len(nrow(base)), each = max_day), , drop = FALSE]
  expanded$icu_day <- rep(day_index, times = nrow(base))

  expanded$window_start <- expanded$icu_intime + (expanded$icu_day - 1L) * 86400
  nominal_end <- expanded$icu_intime + expanded$icu_day * 86400
  expanded$window_end <- .pmin_posix_ignore_na(
    nominal_end,
    expanded$icu_outtime,
    expanded$death_time,
    tz = timezone
  )
  expanded$window_duration_min <- as.numeric(
    difftime(expanded$window_end, expanded$window_start, units = "mins")
  )
  expanded$observable <- !is.na(expanded$window_duration_min) &
    expanded$window_duration_min > 0

  dplyr::as_tibble(expanded[, c(
    "stay_id",
    "icu_day",
    "window_start",
    "window_end",
    "window_duration_min",
    "observable"
  )])
}

compute_daily_twa_paco2 <- function(
    measurements,
    windows,
    id_col = "stay_id",
    time_col = "charttime",
    value_col = "paco2",
    day_col = "icu_day",
    window_start_col = "window_start",
    window_end_col = "window_end",
    min_paco2 = 1,
    max_paco2 = 200,
    duplicate_method = c("mean", "last", "first"),
    timezone = NULL) {
  duplicate_method <- match.arg(duplicate_method)

  .require_columns(measurements, c(id_col, time_col, value_col), "measurements")
  .require_columns(
    windows,
    c(id_col, day_col, window_start_col, window_end_col),
    "windows"
  )

  if (is.null(timezone)) {
    timezone <- .infer_tz(windows[[window_start_col]])
  }

  m <- data.frame(
    stay_id_key = as.character(measurements[[id_col]]),
    charttime = .as_posix(measurements[[time_col]], tz = timezone),
    paco2 = suppressWarnings(as.numeric(measurements[[value_col]])),
    row_order = seq_len(nrow(measurements)),
    stringsAsFactors = FALSE
  )

  w <- data.frame(
    stay_id = windows[[id_col]],
    stay_id_key = as.character(windows[[id_col]]),
    icu_day = as.integer(windows[[day_col]]),
    window_start = .as_posix(windows[[window_start_col]], tz = timezone),
    window_end = .as_posix(windows[[window_end_col]], tz = timezone),
    stringsAsFactors = FALSE
  )

  m <- m |>
    dplyr::filter(
      !is.na(.data$stay_id_key),
      !is.na(.data$charttime),
      !is.na(.data$paco2),
      .data$paco2 >= min_paco2,
      .data$paco2 <= max_paco2
    ) |>
    dplyr::arrange(.data$stay_id_key, .data$charttime, .data$row_order)

  if (nrow(m) > 0) {
    if (duplicate_method == "mean") {
      m <- m |>
        dplyr::group_by(.data$stay_id_key, .data$charttime) |>
        dplyr::summarise(
          paco2 = mean(.data$paco2),
          row_order = min(.data$row_order),
          .groups = "drop"
        )
    } else if (duplicate_method == "last") {
      m <- m |>
        dplyr::group_by(.data$stay_id_key, .data$charttime) |>
        dplyr::slice_tail(n = 1) |>
        dplyr::ungroup()
    } else {
      m <- m |>
        dplyr::group_by(.data$stay_id_key, .data$charttime) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup()
    }
    m <- m |>
      dplyr::arrange(.data$stay_id_key, .data$charttime, .data$row_order)
  }

  measurements_by_stay <- split(m, m$stay_id_key, drop = TRUE)

  compute_one_window <- function(i) {
    wrow <- w[i, , drop = FALSE]
    window_duration_min <- as.numeric(
      difftime(wrow$window_end, wrow$window_start, units = "mins")
    )
    observable <- !is.na(window_duration_min) && window_duration_min > 0

    empty_result <- data.frame(
      stay_id = wrow$stay_id,
      icu_day = wrow$icu_day,
      window_start = wrow$window_start,
      window_end = wrow$window_end,
      window_duration_min = window_duration_min,
      observable = observable,
      n_paco2 = 0L,
      first_paco2_time = as.POSIXct(NA_real_, origin = "1970-01-01", tz = timezone),
      last_paco2_time = as.POSIXct(NA_real_, origin = "1970-01-01", tz = timezone),
      twa_paco2 = NA_real_
    )

    if (!observable) {
      return(empty_result)
    }

    sub_m <- measurements_by_stay[[wrow$stay_id_key]]
    if (is.null(sub_m) || nrow(sub_m) == 0) {
      return(empty_result)
    }

    sub_m <- sub_m[
      sub_m$charttime >= wrow$window_start &
        sub_m$charttime <= wrow$window_end,
      ,
      drop = FALSE
    ]
    if (nrow(sub_m) == 0) {
      return(empty_result)
    }

    sub_m <- sub_m[order(sub_m$charttime, sub_m$row_order), , drop = FALSE]
    time_points <- c(wrow$window_start, sub_m$charttime, wrow$window_end)
    interval_values <- c(sub_m$paco2[1], sub_m$paco2)
    durations_min <- as.numeric(
      difftime(time_points[-1], time_points[-length(time_points)], units = "mins")
    )
    positive_duration <- durations_min > 0

    if (!any(positive_duration)) {
      return(empty_result)
    }

    empty_result$n_paco2 <- nrow(sub_m)
    empty_result$first_paco2_time <- sub_m$charttime[1]
    empty_result$last_paco2_time <- sub_m$charttime[nrow(sub_m)]
    empty_result$twa_paco2 <- sum(
      durations_min[positive_duration] * interval_values[positive_duration]
    ) / sum(durations_min[positive_duration])
    empty_result
  }

  if (nrow(w) == 0) {
    return(dplyr::as_tibble(data.frame(
      stay_id = windows[[id_col]][0],
      icu_day = integer(),
      window_start = as.POSIXct(character()),
      window_end = as.POSIXct(character()),
      window_duration_min = numeric(),
      observable = logical(),
      n_paco2 = integer(),
      first_paco2_time = as.POSIXct(character()),
      last_paco2_time = as.POSIXct(character()),
      twa_paco2 = numeric()
    )))
  }

  dplyr::bind_rows(lapply(seq_len(nrow(w)), compute_one_window))
}

summarise_twa_paco2_coverage <- function(
    twa_paco2,
    id_col = "stay_id") {
  .require_columns(
    twa_paco2,
    c(id_col, "observable", "n_paco2", "twa_paco2"),
    "twa_paco2"
  )

  twa_paco2 |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(
      n_observable_days = sum(.data$observable, na.rm = TRUE),
      n_days_with_paco2 = sum(.data$observable & .data$n_paco2 > 0, na.rm = TRUE),
      has_missing_observable_paco2_day = any(
        .data$observable & .data$n_paco2 == 0,
        na.rm = TRUE
      ),
      all_observable_days_have_paco2 = !.data$has_missing_observable_paco2_day,
      .groups = "drop"
    )
}
