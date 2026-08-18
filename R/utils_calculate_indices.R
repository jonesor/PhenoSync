# Internal metadata used by calculate_indices().
#
# The goal is not to forbid users from requesting an index explicitly, but to
# keep indices = "auto" biologically conservative.

pheno_index_rules <- function() {
  data.frame(
    index = c(
      "albert_Si",
      "augspurger_Xi",
      "bolmgren_S",
      "dainou_SI",
      "freitas_bolmgren_S",
      "koenig_rp",
      "mahoro_Si",
      "mahoro_Ri",
      "marquis_SM",
      "michalski_durka_CV",
      "pianka_overlap",
      "poole_rathcke",
      "primack_Z",
      "schoener_SI",
      "shrotri_SI_temp",
      "shrotri_SI_abd",
      "wadgymar_Si",
      "weis_sync"
    ),
    individual_default = c(
      TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE
    ),
    species_default = c(
      TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE,
      TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE
    ),
    binary_compatible = c(
      TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE,
      TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE
    ),
    numeric_compatible = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
    ),
    range_class = c(
      "bounded_0_1",
      "bounded_0_1",
      "other",
      "bounded_0_1",
      "bounded_0_1",
      "other",
      "bounded_0_1",
      "other",
      "bounded_0_1",
      "other",
      "bounded_0_1",
      "other",
      "bounded_0_1",
      "bounded_0_1",
      "bounded_0_1",
      "bounded_0_1",
      "bounded_0_1",
      "other"
    ),
    stringsAsFactors = FALSE
  )
}

wide_candidate_indices <- function() {
  c(
    "albert_Si",
    "augspurger_Xi",
    "freitas_bolmgren_S",
    "koenig_rp",
    "mahoro_Si",
    "mahoro_Ri",
    "marquis_SM",
    "michalski_durka_CV",
    "pianka_overlap",
    "poole_rathcke",
    "primack_Z",
    "shrotri_SI_temp",
    "shrotri_SI_abd",
    "wadgymar_Si"
  )
}

filter_indices_by_context <- function(candidate_indices, columns_are, binary) {
  rules <- pheno_index_rules()
  rules <- rules[rules$index %in% candidate_indices, , drop = FALSE]

  if (identical(columns_are, "species")) {
    keep <- rules$species_default
  } else {
    keep <- rules$individual_default
  }

  if (isTRUE(binary)) {
    keep <- keep & rules$binary_compatible
  } else if (identical(binary, FALSE)) {
    keep <- keep & rules$numeric_compatible
  }

  rules$index[keep]
}

validate_binary_argument <- function(binary) {
  if (is.null(binary)) return(invisible(TRUE))
  if (!is.logical(binary) || length(binary) != 1 || is.na(binary)) {
    stop("binary must be TRUE, FALSE, or NULL for automatic detection.", call. = FALSE)
  }
  invisible(TRUE)
}

resolve_columns_are <- function(data, columns_are) {
  if (!identical(columns_are, "auto")) return(columns_are)

  cn <- colnames(data)
  if (!is.null(cn)) {
    cn_lower <- tolower(cn)
    if (any(grepl("^(sp|species|taxon)", cn_lower))) return("species")
    if (any(grepl("^(ind|individual|plant|id)", cn_lower))) return("individuals")
  }

  "individuals"
}

detect_binary_values <- function(data, binary = NULL) {
  validate_binary_argument(binary)
  if (!is.null(binary)) return(binary)

  x <- suppressWarnings(as.numeric(as.matrix(data)))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  all(x %in% c(0, 1))
}

binary_matrix <- function(data) {
  m <- suppressWarnings(as.matrix(data))
  storage.mode(m) <- "numeric"
  m <- ifelse(is.na(m), NA_real_, ifelse(m > 0, 1, 0))
  colnames(m) <- colnames(data)
  rownames(m) <- rownames(data)
  m
}

prepare_wide_pheno_data <- function(data,
                                    census_col = NULL,
                                    first_row_colnames = TRUE,
                                    first_col_rownames = FALSE) {
  if (!is.matrix(data) && !is.data.frame(data)) return(data)

  # If the user explicitly identifies the census column, remove it and preserve
  # its values as row names. This is now the preferred workflow for ordinary
  # wide tables imported from CSV/Excel.
  if (!is.null(census_col)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)

    names_are_generic <- function(x) {
      nms <- names(x)
      is.null(nms) || all(grepl("^(V|X)[0-9]+$", nms))
    }

    numeric_like <- function(x) {
      x <- trimws(as.character(x))
      x[x == ""] <- NA_character_
      !is.na(suppressWarnings(as.numeric(x)))
    }

    # Raw tables imported with header = FALSE often have the true column names
    # in the first row. Promote that row to names before resolving census_col.
    if (isTRUE(first_row_colnames) && nrow(df) >= 2 && names_are_generic(df)) {
      header_values <- unlist(df[1, , drop = FALSE], use.names = FALSE)
      header_is_text <- any(!numeric_like(header_values) & !is.na(as.character(header_values)))

      if (isTRUE(header_is_text)) {
        new_names <- as.character(header_values)
        new_names[is.na(new_names) | new_names == ""] <- paste0("col", seq_along(new_names))
        names(df) <- make.unique(new_names)
        df <- df[-1, , drop = FALSE]
        rownames(df) <- NULL
      }
    }

    return(prepare_census_matrix(df, census_col = census_col, data_arg = "data"))
  }

  # Backwards-compatible option. This is no longer the preferred workflow;
  # use census_col instead when the census labels are stored in a data column.
  if (isTRUE(first_col_rownames)) {
    return(prepare_census_matrix(data, census_col = 1, data_arg = "data"))
  }

  # If the input is already a numeric matrix/data frame, keep row names and
  # column names exactly as supplied.
  m0 <- suppressWarnings(as.matrix(data))
  if (is.numeric(m0)) {
    return(prepare_census_matrix(data, census_col = NULL, data_arg = "data"))
  }

  # Raw header-row handling without a census column. This can only work when all
  # remaining rows/columns are numeric after the first row has been promoted.
  df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) < 2 || ncol(df) < 2) return(data)

  names_are_generic <- function(x) {
    nms <- names(x)
    is.null(nms) || all(grepl("^(V|X)[0-9]+$", nms))
  }

  numeric_like <- function(x) {
    x <- trimws(as.character(x))
    x[x == ""] <- NA_character_
    !is.na(suppressWarnings(as.numeric(x)))
  }

  if (isTRUE(first_row_colnames) && nrow(df) >= 2 && names_are_generic(df)) {
    header_values <- unlist(df[1, , drop = FALSE], use.names = FALSE)
    header_is_text <- any(!numeric_like(header_values) & !is.na(as.character(header_values)))

    if (isTRUE(header_is_text)) {
      new_names <- as.character(header_values)
      new_names[is.na(new_names) | new_names == ""] <- paste0("col", seq_along(new_names))
      names(df) <- make.unique(new_names)
      df <- df[-1, , drop = FALSE]
      rownames(df) <- NULL
    }
  }

  prepare_census_matrix(df, census_col = NULL, data_arg = "data")
}
