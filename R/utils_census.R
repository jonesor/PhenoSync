# Internal helpers for preserving census identifiers in wide phenology matrices.

resolve_column_position <- function(data, col, arg_name = "census_col") {
  if (is.null(col)) return(NULL)

  nms <- names(data)

  if (is.character(col) && length(col) == 1 && !is.na(col)) {
    if (is.null(nms) || !(col %in% nms)) {
      stop(arg_name, " = '", col, "' was not found in the data columns.", call. = FALSE)
    }
    return(match(col, nms))
  }

  if (is.numeric(col) && length(col) == 1 && !is.na(col) &&
      col == as.integer(col) && col >= 1 && col <= ncol(data)) {
    return(as.integer(col))
  }

  stop(arg_name, " must be NULL, one column name, or one valid column index.", call. = FALSE)
}

prepare_census_matrix <- function(data,
                                  census_col = NULL,
                                  data_arg = "data",
                                  default_col_prefix = "col") {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop(data_arg, " must be a matrix or data frame.", call. = FALSE)
  }

  if (!is.null(census_col)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
    census_pos <- resolve_column_position(df, census_col, "census_col")
    census_values <- as.character(df[[census_pos]])
    census_values[is.na(census_values) | census_values == ""] <- paste0("census", seq_along(census_values))
    df <- df[, -census_pos, drop = FALSE]

    if (ncol(df) < 1) {
      stop("After removing census_col, at least one phenology column is required.", call. = FALSE)
    }

    converted <- suppressWarnings(
      as.data.frame(lapply(df, function(z) as.numeric(as.character(z))), check.names = FALSE)
    )

    bad_cols <- names(df)[vapply(seq_along(df), function(j) {
      original <- as.character(df[[j]])
      original[trimws(original) == ""] <- NA_character_
      converted_j <- suppressWarnings(as.numeric(original))
      any(is.na(converted_j) & !is.na(original))
    }, logical(1))]

    if (length(bad_cols) > 0) {
      stop(
        "All phenology columns after census_col must be numeric. Non-numeric columns: ",
        paste(bad_cols, collapse = ", "),
        call. = FALSE
      )
    }

    m <- as.matrix(converted)
    storage.mode(m) <- "numeric"
    rownames(m) <- make.unique(census_values)
    if (is.null(colnames(m))) {
      colnames(m) <- paste0(default_col_prefix, seq_len(ncol(m)))
    }
    attr(m, "census_col") <- names(as.data.frame(data, check.names = FALSE))[census_pos]
    attr(m, "census") <- rownames(m)
    return(m)
  }

  m <- suppressWarnings(as.matrix(data))

  if (!is.numeric(m)) {
    stop(
      data_arg, " must be numeric after removing the census column. ",
      "If your first column contains census names or dates, specify it with census_col, ",
      "for example census_col = 1 or census_col = 'Date'.",
      call. = FALSE
    )
  }

  storage.mode(m) <- "numeric"

  if (is.null(rownames(m))) {
    rownames(m) <- paste0("census", seq_len(nrow(m)))
  }
  if (is.null(colnames(m))) {
    colnames(m) <- paste0(default_col_prefix, seq_len(ncol(m)))
  }

  attr(m, "census") <- rownames(m)
  m
}


.prepare_pheno_matrix <- function(data,
                                  census_col = NULL,
                                  first_row_colnames = TRUE,
                                  first_col_rownames = FALSE) {
  prepare_wide_pheno_data(
    data = data,
    census_col = census_col,
    first_row_colnames = first_row_colnames,
    first_col_rownames = first_col_rownames
  )
}
