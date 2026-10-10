# Internal helpers -----------------------------------------------------------
#
# None of the functions in this file are exported. They take care of input
# validation, fast grouping of identical rows, and other small utilities shared
# by the user-facing functions.

#' Coerce user input to a numeric (dense or sparse) matrix and validate it
#'
#' @param x A matrix, data frame or `Matrix` object.
#' @param arg Name of the argument, used in error messages.
#' @return A base numeric matrix, or a `Matrix::sparseMatrix` (column
#'   compressed, `dgCMatrix`) when the input was sparse.
#' @noRd
as_composition_matrix <- function(x, arg = "meta_composition") {
  if (inherits(x, "Matrix")) {
    if (!inherits(x, "sparseMatrix")) {
      x <- as.matrix(x)
    } else {
      x <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix")
      if (!is.numeric(x@x) && !is.logical(x@x)) {
        stop("`", arg, "` must contain numeric values.", call. = FALSE)
      }
      x <- methods::as(x, "dMatrix")
      check_finite(x@x, arg)
      return(x)
    }
  }
  if (is.data.frame(x)) {
    ok <- vapply(x, function(col) is.numeric(col) || is.logical(col), logical(1))
    if (!all(ok)) {
      stop("All columns of `", arg, "` must be numeric; non-numeric column(s): ",
        paste(names(x)[!ok], collapse = ", "), ".",
        call. = FALSE
      )
    }
    x <- as.matrix(x)
  }
  if (is.vector(x) && !is.list(x)) {
    x <- matrix(x, nrow = 1L)
  }
  if (!is.matrix(x)) {
    stop("`", arg, "` must be a matrix, data frame or sparse `Matrix`.", call. = FALSE)
  }
  if (!is.numeric(x) && !is.logical(x)) {
    stop("`", arg, "` must contain numeric values.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  check_finite(x, arg)
  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop("`", arg, "` must have at least one row and one column.", call. = FALSE)
  }
  x
}

check_finite <- function(values, arg) {
  if (anyNA(values)) {
    stop("`", arg, "` contains missing values (NA/NaN).", call. = FALSE)
  }
  if (any(!is.finite(values))) {
    stop("`", arg, "` contains infinite values.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Character keys identifying identical rows of a matrix
#'
#' Equivalent to the keys used internally by `base::duplicated.matrix()`, but
#' vectorised over columns instead of looping over rows with `apply()`, which
#' makes it orders of magnitude faster on tall matrices. Works for dense and
#' sparse matrices.
#' @noRd
row_keys <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    tx <- methods::as(x, "TsparseMatrix")
    i <- tx@i + 1L
    o <- order(i, tx@j)
    entry <- paste(tx@j[o], tx@x[o], sep = ":")
    keys <- character(nrow(x))
    if (length(entry) > 0L) {
      pasted <- vapply(split(entry, factor(i[o], levels = seq_len(nrow(x)))),
        paste,
        character(1),
        collapse = ","
      )
      keys <- unname(pasted)
    }
    return(keys)
  }
  if (ncol(x) == 1L) {
    return(as.character(x[, 1L]))
  }
  do.call(paste, c(lapply(seq_len(ncol(x)), function(j) x[, j]), sep = "\r"))
}

#' Group identical rows
#'
#' @return A list with `group` (integer id of the unique composition of each
#'   row, in order of first appearance), `counts` (multiplicity of each unique
#'   composition) and `first` (index of the first row of each group).
#' @noRd
row_groups <- function(x) {
  keys <- row_keys(x)
  first <- which(!duplicated(keys))
  group <- match(keys, keys[first])
  list(group = group, counts = tabulate(group, nbins = length(first)), first = first)
}

is_sparse <- function(x) inherits(x, "sparseMatrix")

transpose <- function(x) if (is_sparse(x)) Matrix::t(x) else t(x)

#' Resolve the hypervolume method
#' @noRd
resolve_method <- function(hypervolume_method, d, dim_threshold) {
  if (is.null(hypervolume_method)) {
    return(if (d > dim_threshold) "hyper_normal" else "deterministic")
  }
  if (!is.character(hypervolume_method) || length(hypervolume_method) != 1L) {
    stop("`hypervolume_method` must be a single string.", call. = FALSE)
  }
  method <- match.arg(hypervolume_method, c("deterministic", "exact", "hyper_normal"))
  if (method == "exact") method <- "deterministic"
  method
}

check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}
