#' @importFrom stats is.empty.model model.matrix model.response terms na.pass
parseFormula <- function(formula, data, intercept = TRUE, justX = FALSE) {

    # extract Y, X, and variable names for model formula and frame
    mt <- terms(formula, data = data)
    if (missing(data))
        data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.pass                      # keep all rows; missing values are checked explicitly
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    if (!intercept) {
        attributes(mt)$intercept <- 0
    }

    # null model support
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)
    X <- as.matrix(X)  # X matrix
    xvars <- dimnames(X)[[2L]]  # X variable names
    xobs <- dimnames(X)[[1L]]  # X observation names
    if (justX) {
        Y <- NULL
    } else {
        Y <- as.matrix(model.response(mf, "numeric"))  # Y matrix
    }

    return(list(Y, X, xvars, xobs))

}

#' @importFrom stats is.empty.model model.matrix model.response terms model.frame as.formula
parseFormula2 <- function(formula, data, intercept = TRUE, justX = FALSE) {

  # extract Y, X, and variable names for model formula and frame
  mt <- terms(formula, data = data)
  if (missing(data))
    data <- sys.frame(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  mf$intercept <- mf$justX <- NULL
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass                        # keep all rows; missing values are checked explicitly
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  if (!intercept) {
    attributes(mt)$intercept <- 0
  }

  # null model support
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf)
  X <- as.matrix(X)  # X matrix
  xvars <- dimnames(X)[[2L]]  # X variable names
  xobs <- dimnames(X)[[1L]]  # X observation names
  if (justX) {
    Y <- NULL
  } else {
    Y <- as.matrix(model.response(mf, "numeric"))  # Y matrix
  }

  # Parse for variables in parentheses
  # Match the part of the formula inside parentheses
  matches <- regmatches(as.character(formula)[3], gregexpr("\\((.*?)\\)",
                        as.character(formula)[3]))
  if (length(matches[[1]]) > 0) {
    # Filter out invalid cases like (0)
    inner_parts <- matches[[1]]
    inner_parts <- inner_parts[inner_parts != "(0)"]  # Remove "(0)"

    if (length(inner_parts) > 0) {
      # Construct a new formula for valid parentheses content
      inner_formula <- as.formula(paste("~", paste(inner_parts, collapse = "+")))
      inner_terms <- terms(inner_formula, data = data)
      mf_inner <- model.frame(inner_terms, data, drop.unused.levels = TRUE)
      X_tilde <- as.matrix(model.matrix(inner_terms, mf_inner))
    } else {
      # No valid content in parentheses
      X_tilde <- NULL
    }
  } else {
    # No parentheses, X_tilde is NULL
    X_tilde <- NULL
  }
  X_tilde_vars <- dimnames(X_tilde)[[2L]]

  return(list(Y, X, xvars, xobs, X_tilde, X_tilde_vars))

}

# internal function: checks if input is integer
is_integer <- function(x) {
    is.numeric(x) && (floor(x) == x)
}

# internal function: pretty prints a matrix with colnames and rownames
pretty_print_matrix <- function(mat, heading = NULL){

  # Check if the input is a matrix
  if (!is.matrix(mat)) {
    stop("Input must be a matrix.")
  }

  # Get row and column names
  row_names <- rownames(mat)
  col_names <- colnames(mat)

  # If no row names or column names, create defaults
  if(is.null(row_names)){
    row_names <- as.character(1:nrow(mat))
  }

  if(is.null(col_names)){
    col_names <- as.character(1:ncol(mat))
  }

  # Function to determine the maximum length and format the column
  format_column <- function(col){
    # Convert column to character (missing values printed as "NA")
    col_char <- as.character(col)
    col_char[is.na(col_char)] <- "NA"

    # Find the maximum length of entries including the decimal point
    max_length <- max(nchar(col_char))

    # Function to format each number to match the maximum length
    format_to_max_length <- function(x){
      x_char <- as.character(x)
      if(x_char == "NA"){
        return(x_char)
      }
      # Split into integer and decimal parts
      if(grepl("\\.", x_char)){
        integer_part_length <- nchar(sub("\\..*", "", x_char))
        decimal_part_length <- nchar(sub("^[^.]*\\.", "", x_char))
      }else{
        integer_part_length <- nchar(x_char)
        decimal_part_length <- 0
      }

      # Calculate the total length of the number with decimal point
      total_length <- integer_part_length + decimal_part_length + 1
      required_length <- max_length

      # Determine padding
      if(total_length < required_length){
        num_trailing_zeros <- required_length - total_length
        formatted_x <- sprintf(paste0("%.",
            decimal_part_length + num_trailing_zeros, "f"), as.numeric(x))
      }else{
        formatted_x <- x_char
      }

      return(formatted_x)
    }

    # Apply formatting function to each element in the column
    padded_col <- sapply(col_char, format_to_max_length)
    return(padded_col)
  }

  # Apply formatting to each column
  formatted_mat <- apply(mat, 2, format_column)

  # Determine the width for each column
  col_widths <- apply(formatted_mat, 2, function(col) max(nchar(col)))
  col_widths <- pmax(nchar(col_names), col_widths)

  # Print heading is not NULL
  if(!is.null(heading)){
    cat(paste("\n", as.character(heading), "\n\n", sep = ""))
  }

  # Calculate total width of the table, including the vertical lines and spacing
  # +1 for the | separator
  row_name_width <- max(nchar(row_names)) + 1
  # Add row_name_width
  total_width <- sum(col_widths) + length(col_widths) + row_name_width

  # Function to create a separator line
  create_separator <- function(col_widths, row_name_width) {
    separator_parts <- sapply(col_widths,
        function(width) paste0(rep("-", width + 2), collapse = ""))
    separator_line <- paste0("+", paste(rep("-", row_name_width + 1),
        collapse = ""), "+", paste(separator_parts, collapse = "+"), "+\n")
    return(separator_line)
  }

  # Print the header with shifted column names
  # Empty space for row names
  cat(sprintf(" %-*s", row_name_width - 1, ""), " ")
  for (j in seq_along(col_names)) {
    cat(sprintf("| %-*s", col_widths[j] + 1, col_names[j]))
  }
  cat("|\n")

  # Print the separator line
  separator <- create_separator(col_widths, row_name_width)
  cat(separator)

  # Print each row with row names and vertical lines
  for (i in seq_along(row_names)) {
    row_str <- sprintf("| %-*s", row_name_width, row_names[i])
    for (j in seq_along(col_names)) {
      if(j == length(col_names)){
        row_str <- paste0(row_str, sprintf("| %-*s", col_widths[j],
                        formatted_mat[i, j]))
      }else{
        row_str <- paste0(row_str, sprintf("| %*s", col_widths[j] + 1,
                        formatted_mat[i, j]))
      }
    }
    cat(row_str, "|\n")
  }

  # Print the bottom border
  cat(separator)
  cat("\n")
}

# internal function: inverse logit transformation
ilogit <- function(x){
  return(1.0 / (1.0 + exp(- x)))
}

#' Create a collection of candidate models for stacking
#'
#' @description Creates an object of class \code{'candidateModels'} that
#' contains a list of candidate models for stacking. The function takes a list
#' of candidate values for each model parameter and returns a list of possible
#' combinations of these values based on either simple aggregation or Cartesian
#' product of indivdual candidate values.
#' @param params_list a list of candidate values for each model parameter. See
#' examples for details.
#' @param aggregation a character string specifying the type of aggregation to
#' be used. Options are \code{'simple'} and \code{'cartesian'}. Default is
#' \code{'simple'}.
#' @return an object of class \code{'candidateModels'}
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @seealso [stvcGLMstack()]
#' @examples
#' m1 <- candidateModels(list(phi_s = c(1, 1), phi_t = c(1, 2)), "simple")
#' m1
#' m2 <- candidateModels(list(phi_s = c(1, 1), phi_t = c(1, 2)), "cartesian")
#' m2
#' m3 <- candidateModels(list(phi_s = list(c(1, 1), c(1, 2)),
#'                           phi_t = list(c(1, 3), c(2, 3)),
#'                           boundary = c(0.5, 0.75)),
#'                       "simple")
#' @export
candidateModels <- function(params_list, aggregation = "simple"){

    if(!is.list(params_list)){
      stop("params_list must be a list")
    }

    if(!is.character(aggregation)){
      stop("aggregation must be a character string")
    }

    aggregation <- tolower(aggregation)

    if(!aggregation %in% c("simple", "cartesian")){
      stop("aggregation must be either 'simple' or 'cartesian'")
    }

    if(aggregation == "simple"){

      if(!all(sapply(params_list, is.vector))){
        stop("All elements of params_list must be vectors")
      }

      if(length(unique(sapply(params_list, length))) > 1){
        stop("All vectors in params_list must be of the same length")
      }

      models_list <- do.call("cbind", params_list)
      models_list <- apply(models_list, 1, function(x) as.vector(x, mode = "list"))

    }else{

      models_list <- expand.grid(params_list)
      models_list <- apply(models_list, 1, function(x) as.vector(x, mode = "list"))

    }

    class(models_list) <- "candidateModels"
    return(models_list)

}

# internal function: find duplicate rows in two matrices
find_approx_matches <- function(pred, obs, tol = 1e-6) {
  if (ncol(pred) != ncol(obs)) {
    stop("Prediction and observed matrices must have the same number of columns.")
  }

  # rows of pred within tol of some row of obs in every coordinate (L-infinity
  # distance <= tol). The coordinates are binned on a grid of width tol: a match
  # can only lie in the same or an adjacent bin, so each row of pred is compared
  # with the obs rows in its 3^d neighbouring bins (hashed keys) instead of with
  # all of obs. O((n_pred + n_obs) 3^d) instead of O(n_pred n_obs).
  d <- ncol(pred)
  bin_obs <- floor(obs / tol)
  bin_pred <- floor(pred / tol)
  key <- function(b) do.call("paste", c(lapply(seq_len(d), function(j) sprintf("%.0f", b[, j])), sep = ":"))
  obs_keys <- key(bin_obs)
  obs_rows <- split(seq_len(nrow(obs)), obs_keys)
  offsets <- as.matrix(expand.grid(rep(list(-1:1), d)))

  is_match <- logical(nrow(pred))
  for (o in seq_len(nrow(offsets))) {
    cand <- obs_rows[key(sweep(bin_pred, 2, offsets[o, ], "+"))]
    todo <- which(!is_match & lengths(cand) > 0)
    for (i in todo) {
      dif <- abs(sweep(obs[cand[[i]], , drop = FALSE], 2, pred[i, ]))
      if (any(apply(dif, 1, max) <= tol)) {
        is_match[i] <- TRUE
      }
    }
  }
  match_idx <- which(is_match)

  matched <- pred[match_idx, , drop = FALSE]

  list(
    any_match = length(match_idx) > 0,
    num_matches = length(match_idx),
    matched_rows = matched
  )
}

# internal function: sample size specific threshold for the Pareto k diagnostic
# of PSIS, min(1 - 1/log10(S), 0.7) (Vehtari et al. 2024; loo::ps_khat_threshold)
psis_khat_threshold <- function(S){
  min(1 - 1 / log10(S), 0.7)
}

# TRUE if x is a single, finite, whole number (e.g., a 1-based index)
is_whole_scalar <- function(x){
  is.numeric(x) && length(x) == 1 && is.finite(x) && x == round(x)
}

# internal function: stops if any row of coords is duplicated. The spatial and
# spatial-temporal process models assume distinct locations: repeated rows make
# the correlation matrix singular. For spatial-temporal models, coords is the
# matrix of spatial and temporal coordinates, so a row is a duplicate only if
# it coincides in space and in time.
check_distinct_coords <- function(coords, what = "spatial locations",
                                  hint = "Average the observations at a common location, or model them as spatial-temporal data."){

  dup <- which(duplicated(coords))
  if(length(dup) > 0){
    first <- match(data.frame(t(coords[dup, , drop = FALSE])),
                   data.frame(t(coords)))
    n_show <- min(6, length(dup))
    preview <- vapply(seq_len(n_show), function(k){
      paste0("  row ", dup[k], " duplicates row ", first[k], ": (",
             paste(format(coords[dup[k], ], digits = 6), collapse = ", "), ")")
    }, character(1))
    stop(length(dup), " duplicated ", what, " found; the model requires ",
         "distinct ", what, ". ", hint, "\nFirst ", n_show,
         " duplicate(s):\n", paste(preview, collapse = "\n"), call. = FALSE)
  }

  invisible(TRUE)

}

# Is R linked to an optimized BLAS? Judged from the BLAS and LAPACK library
# paths reported by R (see sessionInfo()): OpenBLAS, Intel MKL, BLIS, Apple
# Accelerate (vecLib) or ATLAS. Returns FALSE when unsure, e.g., for R's
# reference BLAS, or an optimized BLAS installed under the reference name
# (common on Windows, where Rblas.dll is replaced).
blas_is_optimized <- function(){

  libs <- tryCatch(c(extSoftVersion()[["BLAS"]], La_library()),
                   error = function(e) character(0))
  libs <- tolower(paste(libs, collapse = " "))
  grepl("openblas|mkl|blis|accelerate|veclib|atlas", libs)

}

# Pre-processing of the K-fold cross-validation subsets: "update" (deletion
# updates of the full-data factors, scalar loops, fastest with the reference
# BLAS) or "direct" (recomputation on each subset, level-3 BLAS, fastest with an
# optimized BLAS). "auto" picks "direct" if blas_is_optimized(), else "update".
# Both give the same results up to floating-point rounding.
resolve_CV_update <- function(CV.update = "auto"){

  if(!is.character(CV.update) || length(CV.update) != 1){
    stop("error: CV.update in loopd.controls must be one of 'auto', 'update' or 'direct'.")
  }
  CV.update <- tolower(CV.update)
  if(!CV.update %in% c("auto", "update", "direct")){
    stop("CV.update = '", CV.update, "' is not a valid option; choose from c('auto', 'update', 'direct').")
  }
  if(CV.update == "auto"){
    CV.update <- if(blas_is_optimized()) "direct" else "update"
  }

  CV.update

}

# Fit diagnostics. Each C++ fit returns c(min.pivot, min.cor, max.cor): the smallest
# relative Cholesky pivot of the n x n factorizations (its inverse bounds the
# condition number below), and the smallest and largest off-diagonal entry of the
# correlation matrix (the correlations of the farthest-apart and of the closest
# locations). Computed from quantities the fit already forms. The 'diagnostics'
# element of a fit is a list with
#   numerical: data frame of these values, one row per model;
#   pareto:    (PSIS only) list(k, threshold, n.high) of Pareto k diagnostics;
#   solver:    (stacking only) details of the optimization for the weights.
diagnostics_thresholds <- list(min.pivot = 1e-8, min.cor = 0.95, max.cor = 0.05)

# data frame of the numerical diagnostics of a list of fits, one row per fit, or,
# for fits with several correlation matrices (stvc 'independent'), one row per
# fit and process, with columns 'model' and 'process'
collect_diagnostics <- function(fits, model.names = NULL){

  parts <- lapply(seq_along(fits), function(m){
    x <- fits[[m]][["diagnostics"]]
    if(!is.matrix(x)){
      x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
    }
    data.frame(model = m, process = seq_len(nrow(x)), x, check.names = FALSE)
  })
  d <- do.call("rbind", parts)
  multi <- any(d$process > 1)
  labels <- if(is.null(model.names)) rep("", nrow(d)) else model.names[d$model]
  if(multi){
    labels <- paste0(labels, ifelse(nchar(labels) > 0, ", ", ""), "process ", d$process)
  }else{
    d$model <- NULL
    d$process <- NULL
  }
  if(any(nchar(labels) > 0)){
    rownames(d) <- labels
  }else{
    rownames(d) <- NULL
  }
  d

}

# Pareto k diagnostics: k is a vector (one model) or a list of vectors (one per
# model); returns list(k, threshold, n.high)
pareto_diagnostics <- function(k, n.samples){

  threshold <- psis_khat_threshold(n.samples)
  n.high <- if(is.list(k)) vapply(k, function(x) sum(x > threshold), integer(1)) else sum(k > threshold)
  list(k = k, threshold = threshold, n.high = n.high)

}

# issues flagged by the numerical diagnostics: a list with one character vector
# per row. pivot.hint: likely causes of a small pivot, model specific
diagnostics_issues <- function(d, pivot.hint = "nearly coincident locations, or a very small decay parameter"){

  th <- diagnostics_thresholds
  lapply(seq_len(nrow(d)), function(i){
    out <- character(0)
    if(isTRUE(d$min.pivot[i] < th$min.pivot)){
      out <- c(out, paste0("smallest relative Cholesky pivot is ",
                           format(signif(d$min.pivot[i], 2)), ": the covariance",
                           " matrix is nearly singular and results may lose",
                           " accuracy (", pivot.hint, ")."))
    }
    if(isTRUE(d$min.cor[i] > th$min.cor)){
      out <- c(out, paste0("correlation between the two farthest-apart",
                           " locations is ", format(round(d$min.cor[i], 3)),
                           ": the effective range far exceeds the extent of the",
                           " data; consider a larger decay parameter."))
    }
    if(isTRUE(d$max.cor[i] < th$max.cor)){
      out <- c(out, paste0("correlation between the two closest locations is ",
                           format(signif(d$max.cor[i], 2)), ": the locations",
                           " are nearly uncorrelated and the spatial effect is",
                           " hard to separate from noise; consider a smaller",
                           " decay parameter."))
    }
    out
  })

}

# issue flagged by the Pareto k diagnostics of model i (character(0) if none)
pareto_issue <- function(pareto, i = 1){

  n.high <- pareto$n.high[i]
  if(n.high == 0){
    return(character(0))
  }
  n.obs <- length(if(is.list(pareto$k)) pareto$k[[i]] else pareto$k)
  paste0(n.high, " of ", n.obs, " Pareto k diagnostic values exceed ",
         format(round(pareto$threshold, 2)), ": the PSIS estimates of the",
         " corresponding leave-one-out predictive densities may be unreliable;",
         " consider loopd.method = 'exact'.")

}

# issues of the optimization for the stacking weights (character(0) if none)
solver_issues <- function(solver){

  out <- character(0)
  if(length(solver$missing.requested) > 0){
    out <- c(out, paste0("requested solver(s) not installed: ",
                         paste(solver$missing.requested, collapse = ", "),
                         "; the default order was used."))
  }
  if(identical(solver$used, "none")){
    out <- c(out, paste0("no solver reached an optimal solution and the",
                         " fallback needs the 'loo' package, which is not",
                         " installed; the stacking weights are NA (see the",
                         " message above for how to compute them)."))
  }else if(isTRUE(solver$fallback)){
    tried <- solver$attempts
    tried_txt <- if(nrow(tried) > 0) paste0(" (", paste(tried$solver, tried$status, sep = ": ", collapse = "; "), ")") else ""
    out <- c(out, paste0("the CVXR solvers failed or did not reach optimality",
                         tried_txt, "; the weights are from loo::stacking_weights()."))
  }else if(identical(solver$status, "optimal_inaccurate")){
    out <- c(out, paste0("solver ", solver$used, " returned status",
                         " 'optimal_inaccurate'."))
  }
  out

}

# Prints a "Diagnostics" section if there is any issue. With stacking weights,
# numerical issues are detailed only for models with weight above weight.min
# (extreme candidates with negligible weight are expected in a stacking grid and
# are counted); Pareto k issues are shown for every model, since they affect the
# weights themselves. ... is passed to diagnostics_issues().
print_diagnostics <- function(diag, weights = NULL, weight.min = 0.05, ...){

  d <- diag$numerical
  M <- nrow(d)
  midx <- if(is.null(d$model)) seq_len(M) else d$model            # model of each row
  num_issues <- diagnostics_issues(d, ...)
  first_row <- !duplicated(midx)                                  # Pareto k values are per model, not per process
  par_issues <- lapply(seq_len(M), function(i){
    if(is.null(diag$pareto) || !first_row[i]) character(0) else pareto_issue(diag$pareto, midx[i])
  })
  show_num <- rep(TRUE, M)
  if(!is.null(weights)){
    show_num <- is.na(weights[midx]) | weights[midx] > weight.min
  }
  blocks <- lapply(seq_len(M), function(i){
    c(if(show_num[i]) num_issues[[i]] else character(0), par_issues[[i]])
  })
  n_hidden <- length(unique(midx[lengths(num_issues) > 0 & !show_num]))
  sol_issues <- if(is.null(diag$solver)) character(0) else solver_issues(diag$solver)

  if(all(lengths(blocks) == 0) && n_hidden == 0 && length(sol_issues) == 0){
    return(invisible(NULL))
  }

  bullet <- function(msg){
    cat(paste(strwrap(msg, width = 76, initial = "  - ", prefix = "    "),
              collapse = "\n"), "\n", sep = "")
  }

  cat("----------------------------------------\n")
  cat("\tDiagnostics\n")
  cat("----------------------------------------\n")
  for(i in which(lengths(blocks) > 0)){
    if(M > 1 || !is.null(weights)){
      lab <- rownames(d)[i]
      if(!is.null(weights)){
        lab <- paste0(lab, " (stacking weight ", format(round(weights[midx[i]], 3)), ")")
      }
      cat(lab, ":\n", sep = "")
    }
    for(msg in blocks[[i]]){
      bullet(msg)
    }
  }
  if(n_hidden > 0){
    cat(paste(strwrap(paste0(n_hidden, " other candidate model(s) with stacking",
                             " weight at most ", weight.min, " have numerical",
                             " flags; see the 'diagnostics' element."), width = 76),
              collapse = "\n"), "\n", sep = "")
  }
  if(length(sol_issues) > 0){
    cat("Stacking weights:\n")
    for(msg in sol_issues){
      bullet(msg)
    }
  }
  cat("----------------------------------------\n")

  invisible(NULL)

}

# message when no solver reached an optimal solution and loo is not installed:
# the fit is returned with NA weights, and this gives the code that computes
# them. A function cannot see the name its output is assigned to, so the code
# uses 'fit' as a placeholder.
loo_fallback_message <- function(){

  paste0("None of the CVXR solvers reached an optimal solution, and the fallback",
         " loo::stacking_weights() needs the 'loo' package, which is not",
         " installed. The stacking weights are set to NA; the fitted models are",
         " kept. To compute the weights, run the following, with 'fit' replaced",
         " by the name the output was saved as:\n\n",
         "  install.packages(\"loo\")\n",
         "  w <- as.numeric(loo::stacking_weights(do.call(\"cbind\", fit$loopd)))\n",
         "  fit$stacking.weights <- w\n",
         "  fit$stacking.summary[, \"weight\"] <- round(w, 3)\n")

}

# Parameters of a fitted model, as a named list read by name: for a stacked fit
# (spLMstack, spGLMstack, stvcGLMstack and their posteriorPredict outputs) those
# of candidate model i, from 'model.params' (one named list per candidate); for
# an exact fit, its 'model.params'. For the spatial models nu is NA when the
# correlation function is exponential (it is not used).
model_params <- function(fit, i = 1){

  mp <- fit$model.params
  if(is.null(mp)){
    stop("the fitted object has no 'model.params'; it was made by an earlier version of spStack, refit it with this version.")
  }
  pars <- if(length(mp) > 0 && is.list(mp[[1]])) mp[[i]] else mp
  if(!is.null(pars[["phi"]]) && is.null(pars[["nu"]])){
    pars[["nu"]] <- NA_real_
  }
  pars

}

# warns about arguments passed through ... that the function does not use
check_dots <- function(...){

  k <- ...length()
  if(k > 0){
    nms <- ...names()
    if(is.null(nms)){
      nms <- rep("", k)
    }
    nms[is.na(nms) | nms == ""] <- "<unnamed>"
    for(nm in nms){
      warning("'", nm, "' is not an argument", call. = FALSE)
    }
  }
  invisible(NULL)

}

# stops if any of the named inputs has missing values
check_no_missing <- function(...){

  args <- list(...)
  bad <- names(args)[vapply(args, anyNA, logical(1))]
  if(length(bad) > 0){
    stop("missing values (NA) in ", paste(bad, collapse = ", "), "; remove the",
         " incomplete observations (and the corresponding rows of the",
         " coordinates) before fitting.", call. = FALSE)
  }
  invisible(TRUE)

}

# stops unless V is a symmetric positive definite k x k matrix
check_cov_matrix <- function(V, k, name){

  V <- matrix(as.numeric(V), k, k)
  if(!isSymmetric(V, tol = 100 * .Machine$double.eps * max(1, max(abs(V))))){
    stop(name, " must be a symmetric matrix.", call. = FALSE)
  }
  if(inherits(tryCatch(chol(V), error = function(e) e), "error")){
    stop(name, " must be positive definite.", call. = FALSE)
  }
  invisible(TRUE)

}

# priors of the conjugate Gaussian models: "flat" assigns p(beta, sigma.sq)
# proportional to 1/sigma.sq; a list with tags 'beta.norm' and/or 'sigma.sq.ig'
# assigns N(mu, sigma.sq*V) to beta and/or IG(a, b) to sigma.sq, a component not
# supplied receiving its flat prior. Returns the arguments of the C++ routines
# (beta.prior, beta.Norm, sigma.sq.IG) and their record for the output (out).
parse_lm_priors <- function(priors, p){

  msg <- "priors must be either 'flat' or a named list with tags 'beta.norm' and/or 'sigma.sq.ig'."
  beta.prior <- "flat"
  beta.Norm <- 0
  sigma.sq.prior <- "flat"
  sigma.sq.IG <- c(0.0, 0.0)

  if(is.character(priors)){
    if(length(priors) != 1 || tolower(priors) != "flat"){
      stop(msg)
    }
  }else if(is.list(priors)){
    if(is.null(names(priors))){
      stop(msg)
    }
    names(priors) <- tolower(names(priors))
    if(any(!names(priors) %in% c("beta.norm", "sigma.sq.ig"))){
      stop("invalid tag(s) in priors: '",
           paste(setdiff(names(priors), c("beta.norm", "sigma.sq.ig")), collapse = "', '"),
           "'. Valid tags are 'beta.norm' and 'sigma.sq.ig'.")
    }
    if("beta.norm" %in% names(priors)){
      beta.Norm <- priors[["beta.norm"]]
      if(!is.list(beta.Norm) || length(beta.Norm) != 2){
        stop("priors[['beta.norm']] must be a list of length 2.")
      }
      if(!is.numeric(beta.Norm[[1]]) || length(beta.Norm[[1]]) != p){
        stop("priors[['beta.norm']][[1]] must be a numeric vector of length ", p, ".")
      }
      if(!is.numeric(beta.Norm[[2]]) || length(beta.Norm[[2]]) != p^2){
        stop("priors[['beta.norm']][[2]] must be a ", p, "x", p, " covariance matrix.")
      }
      check_cov_matrix(beta.Norm[[2]], p, "the prior covariance of beta (beta.norm[[2]])")
      beta.Norm <- list(as.double(beta.Norm[[1]]), matrix(as.double(beta.Norm[[2]]), p, p))
      beta.prior <- "normal"
    }
    if("sigma.sq.ig" %in% names(priors)){
      sigma.sq.IG <- priors[["sigma.sq.ig"]]
      if(!is.numeric(sigma.sq.IG) || length(sigma.sq.IG) != 2 || any(!is.finite(sigma.sq.IG)) ||
         any(sigma.sq.IG <= 0)){
        stop("priors[['sigma.sq.ig']] must be a positive numeric vector of length 2.")
      }
      sigma.sq.prior <- "ig"
    }
  }else{
    stop(msg)
  }
  storage.mode(sigma.sq.IG) <- "double"

  out <- list(beta.Norm = if(beta.prior == "normal") list(mu = beta.Norm[[1]], V = beta.Norm[[2]]) else "flat",
              sigma.sq.IG = if(sigma.sq.prior == "ig") sigma.sq.IG else "flat")
  list(beta.prior = beta.prior, beta.Norm = beta.Norm, sigma.sq.IG = sigma.sq.IG, out = out)

}

# response, design matrices and space-time coordinates of the spatially-temporally
# varying coefficients models (Gaussian response), with input checks
stvc_lm_data <- function(formula, data, sp_coords, time_coords){

  if(missing(formula) || !inherits(formula, "formula")){
    stop("formula must be specified as a formula, e.g. y ~ x1 + (x1).")
  }
  holder <- parseFormula2(formula, data)
  if(ncol(holder[[1L]]) != 1){
    stop("the response must be a single numeric variable.")
  }
  y <- as.numeric(holder[[1L]])
  X <- as.matrix(holder[[2L]])
  X_tilde <- holder[[5L]]
  if(is.null(X_tilde)){
    stop("formula does not indicate varying coefficient terms; put them in parentheses, e.g. y ~ x1 + (x1).")
  }
  X_tilde <- as.matrix(X_tilde)
  n <- nrow(X)

  if(!is.matrix(sp_coords) || ncol(sp_coords) != 2 || nrow(sp_coords) != n){
    stop("sp_coords must be an n x 2 matrix of spatial coordinates, with n = ", n,
         " the number of observations in the model formula.")
  }
  if(is.data.frame(time_coords)){
    time_coords <- as.matrix(time_coords)
  }
  if(is.vector(time_coords) && is.numeric(time_coords)){
    time_coords <- matrix(time_coords, ncol = 1)
  }
  if(!is.matrix(time_coords) || ncol(time_coords) != 1 || nrow(time_coords) != n){
    stop("time_coords must be an n x 1 matrix (or a vector of length n) of temporal coordinates, with n = ", n,
         " the number of observations in the model formula.")
  }
  if(nrow(X_tilde) != n){
    stop("the varying coefficient terms and the model formula have different numbers of observations.")
  }

  check_no_missing(y = y, X = X, X_tilde = X_tilde, sp_coords = sp_coords, time_coords = time_coords)
  check_distinct_coords(cbind(sp_coords, time_coords),
                        what = "spatial-temporal coordinates",
                        hint = "Average the observations that share both location and time.")

  storage.mode(y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(X_tilde) <- "double"
  storage.mode(sp_coords) <- "double"
  storage.mode(time_coords) <- "double"

  list(y = y, X = X, X.names = holder[[3L]], X_tilde = X_tilde, X_tilde.names = holder[[6L]],
       sp_coords = sp_coords, time_coords = time_coords,
       n = as.integer(n), p = as.integer(ncol(X)), r = as.integer(ncol(X_tilde)))

}

# process type of the spatially-temporally varying coefficients linear models
check_stvc_lm_process_type <- function(process.type){

  if(missing(process.type)){
    stop("process.type must be specified. Choose from c('independent', 'independent.shared').")
  }
  if(!is.character(process.type) || length(process.type) != 1){
    stop("process.type must be one of 'independent' or 'independent.shared'.")
  }
  if(process.type == "multivariate"){
    stop("process.type = 'multivariate' is not available for the Gaussian model; choose 'independent' or 'independent.shared'.")
  }
  if(!process.type %in% c("independent", "independent.shared")){
    stop("Invalid process.type. Choose from c('independent', 'independent.shared').")
  }
  process.type

}
