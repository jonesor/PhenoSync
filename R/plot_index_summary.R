#' Plot summary estimates from multiple phenological synchrony indices
#'
#' Creates a base R two-panel dot plot from the summary table returned by
#' `calculate_indices()`. Indices bounded between 0 and 1 are shown in the left
#' panel with a fixed x-axis from 0 to 1. Indices on other scales are shown in
#' the right panel with a data-driven x-axis.
#'
#' @param summary_df A summary data frame returned by `calculate_indices()`.
#' @param fixed_left_xlim Numeric vector of length two giving the x-axis limits
#'   for the left panel containing indices bounded between 0 and 1.
#' @param main Character. Main title for the complete plot.
#' @param point_cex Numeric. Point size for plotted index estimates.
#' @param label_cex Numeric. Axis-label size for index names.
#'
#' @return Invisibly returns the input summary data frame. The plot is drawn as a
#'   base R plotting side effect.
#'
#' @examples
#' x <- make_pheno_dummy(format = "wide", n_species = 5, n_days = 30, seed = 1)
#' res <- calculate_indices(x, columns_are = "species", verbose = FALSE, make_plot = FALSE)
#' plot_index_summary(res$summary)
#'
#' @export
plot_index_summary <- function(summary_df,
                               fixed_left_xlim = c(0, 1),
                               main = "Phenological synchrony index summary",
                               point_cex = 1.1,
                               label_cex = 0.8) {
  if (!is.data.frame(summary_df) || nrow(summary_df) == 0) {
    graphics::plot.new()
    graphics::title(main = main)
    graphics::text(0.5, 0.5, "No index estimates available")
    return(invisible(summary_df))
  }

  required_cols <- c("index", "estimate")
  if (!all(required_cols %in% names(summary_df))) {
    stop("summary_df must contain at least 'index' and 'estimate' columns.", call. = FALSE)
  }

  df <- summary_df
  if (!("range_class" %in% names(df))) {
    df$range_class <- "other"
  }
  if (!("range_panel" %in% names(df))) {
    df$range_panel <- ifelse(
      df$range_class == "bounded_0_1",
      "Indices bounded 0 to 1",
      "Indices on other scales"
    )
  }

  df$estimate <- suppressWarnings(as.numeric(df$estimate))
  df <- df[is.finite(df$estimate), , drop = FALSE]

  if (nrow(df) == 0) {
    graphics::plot.new()
    graphics::title(main = main)
    graphics::text(0.5, 0.5, "No finite index estimates available")
    return(invisible(summary_df))
  }

  bounded_df <- df[df$range_class == "bounded_0_1", , drop = FALSE]
  other_df <- df[df$range_class != "bounded_0_1", , drop = FALSE]

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mfrow = c(1, 2), mar = c(5, 11, 4, 2), oma = c(0, 0, 3, 0))

  draw_panel <- function(panel_df, panel_title, xlim, fixed_axis = FALSE) {
    if (nrow(panel_df) == 0) {
      graphics::plot.new()
      graphics::title(main = panel_title)
      graphics::text(0.5, 0.5, "No indices in this group")
      return(invisible(NULL))
    }

    panel_df <- panel_df[order(panel_df$estimate, decreasing = FALSE), , drop = FALSE]
    y <- seq_len(nrow(panel_df))

    if (is.null(xlim)) {
      xr <- range(panel_df$estimate, finite = TRUE, na.rm = TRUE)
      if (!all(is.finite(xr))) xr <- c(0, 1)
      if (diff(xr) == 0) {
        pad <- ifelse(abs(xr[1]) < 1, 0.1, abs(xr[1]) * 0.1)
        xr <- xr + c(-pad, pad)
      } else {
        pad <- diff(xr) * 0.08
        xr <- xr + c(-pad, pad)
      }
      xlim <- xr
    }

    graphics::plot(
      panel_df$estimate,
      y,
      xlim = xlim,
      ylim = c(0.5, length(y) + 0.5),
      yaxt = "n",
      xlab = "Index value",
      ylab = "",
      main = panel_title,
      pch = 19,
      cex = point_cex
    )

    graphics::axis(2, at = y, labels = panel_df$index, las = 1, cex.axis = label_cex)
    graphics::abline(h = y, col = "grey90", lty = "dotted")
    graphics::points(panel_df$estimate, y, pch = 19, cex = point_cex)

    if (isTRUE(fixed_axis)) {
      graphics::axis(1, at = seq(0, 1, by = 0.25))
      graphics::abline(v = c(0, 1), col = "grey80", lty = "dashed")
    }

    invisible(NULL)
  }

  draw_panel(
    bounded_df,
    panel_title = "Indices bounded 0 to 1",
    xlim = fixed_left_xlim,
    fixed_axis = TRUE
  )

  draw_panel(
    other_df,
    panel_title = "Indices on other scales",
    xlim = NULL,
    fixed_axis = FALSE
  )

  graphics::mtext(main, outer = TRUE, cex = 1.1, font = 2)

  invisible(summary_df)
}
