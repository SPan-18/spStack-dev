#' @importFrom stats is.empty.model model.matrix model.response terms
parseFormula <- function(formula, data, intercept = TRUE, justX = FALSE) {

    # extract Y, X, and variable names for model formula and frame
    mt <- terms(formula, data = data)
    if (missing(data))
        data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
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
    # Convert column to character
    col_char <- as.character(col)

    # Find the maximum length of entries including the decimal point
    max_length <- max(nchar(col_char))

    # Function to format each number to match the maximum length
    format_to_max_length <- function(x){
      x_char <- as.character(x)
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

# internal function: input a list of candidate values of all model parameters
# expands list of vectors, called inside spLMstack and similar functions
candidate_models <- function(params_list){

    models <- expand.grid(params_list)
    models_list <- apply(models, 1, function(x) as.vector(x, mode = "list"))

    return(models_list)

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

  match_idx <- integer(0)

  for (i in seq_len(nrow(pred))) {
    row_pred <- pred[i, ]
    diffs <- abs(sweep(obs, 2, row_pred))      # matrix of |obs - row_pred|
    rowwise_max_diff <- apply(diffs, 1, max)   # L∞ norm (max elementwise difference)

    if (any(rowwise_max_diff <= tol)) {
      match_idx <- c(match_idx, i)
    }
  }

  matched <- pred[match_idx, , drop = FALSE]

  list(
    any_match = length(match_idx) > 0,
    num_matches = length(match_idx),
    matched_rows = matched
  )
}


