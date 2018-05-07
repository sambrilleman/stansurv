#' Calculate basis terms for restricted cubic splines
#'
#' The \code{rcs} function calculates the basis terms for the restricted
#' cubic splines as defined in Royston and Parmar (2002). The \code{drcs}
#' function calculates the corresponding first derivative of the basis terms.
#'
#' @export
#'
#' @param x A numeric vector for which to evaluate the spline basis terms.
#' @param df Integer specifying the degrees of freedom for the splines. If
#'   the \code{iknots} argument is \code{NULL}, then \code{df} - 1
#'   internal knots are placed at equally spaced percentiles of \code{x}.
#' @param iknots An optional vector of internal knot positions. Note that the
#'   \code{df} argument is ignored if \code{iknots} is specified.
#' @param bknots A numeric vector of length 2 specifying the position
#'   of the lower and upper boundary knots, beyond which a linearity
#'   restriction is enforced. If \code{NULL}, then the boundary knots are set
#'   to the minimum and maximum of \code{x}.
#'
#' @return A matrix, with dimension \code{length(x)} by \code{df}, where
#'   \code{df} is one greater than the number of internal knots.
#'
#' @references Royston P, Parmar MKB. Flexible parametric proportional-hazards
#'   and proportional-odds models for censored survival data, with application
#'   to prognostic modelling and estimation of treatment effects.
#'   \emph{Statistics in Medicine}. 2002;21:2175--2197. \doi{10.1002/sim.1203}
#'
#' @examples
#' times <- runif(40,0,10)
#'
#' rcs1 <- rcs(times, df = 2) # one internal knot at median
#' rcs2 <- rcs(times, iknots = c(2,5,8))
#' rcs3 <- rcs(times, iknots = c(2,5), bknots = c(1,8))
#'
#' drcs1 <- drcs(times, df = 2)
#' drcs2 <- drcs(times, iknots = c(2,5,8))
#'
rcs <- function(x, df = 3, iknots = NULL, bknots = NULL) {

  knots <- get_knots(x, df = df, iknots = iknots, bknots = bknots)

  iknots <- knots$iknots # internal knot locations
  bknots <- knots$bknots # boundary knot locations
  df     <- knots$df     # validated degrees of freedom

  lambdas <- (bknots[2] - iknots) / (bknots[2] - bknots[1])

  rcs <- uapply(1:df, function(i) {
    if (i == 1) {
      out <- x
    } else {
      j <- i - 1
      x1 <- pmax(x - iknots[j], 0)
      x2 <- pmax(x - bknots[1], 0)
      x3 <- pmax(x - bknots[2], 0)
      out <- I(x1) ^ 3 - lambdas[j] * I(x2) ^ 3 - (1 - lambdas[j]) * I(x3) ^ 3
    }
    out
  })

  matrix(rcs, ncol = df)
}

#' @export
#' @rdname rcs
drcs <- function(x, df = 3, iknots, bknots = NULL) {

  knots <- get_knots(x, df = df, iknots = iknots, bknots = bknots)

  iknots <- knots$iknots # internal knot locations
  bknots <- knots$bknots # boundary knot locations
  df     <- knots$df

  lambdas <- (bknots[2] - iknots) / (bknots[2] - bknots[1])

  drcs <- uapply(1:df, function(i) {
    if (i == 1) {
      out <- rep(1, length(x))
    } else {
      j <- i - 1
      x1 <- pmax(x - iknots[j], 0)
      x2 <- pmax(x - bknots[1], 0)
      x3 <- pmax(x - bknots[2], 0)
      out <- 3 * I(x1) ^ 2 - 3 * lambdas[j] * I(x2) ^ 2 - 3 * (1 - lambdas[j]) * I(x3) ^ 2
    }
    out
  })

  matrix(drcs, ncol = df)
}


#----- internal

# Get the internal and boundary knot locations from a numeric vector
get_knots <- function(x, df = 3, iknots = NULL, bknots = NULL) {

  if (is.null(iknots)) {
    # user did not specify internal knot locations
    nk <- df - 1
    df <- validate_positive_scalar(df)
    iknots <- qtile(x, nq = df)
  } else {
    # user did specify internal knot locations
    nk <- length(iknots)
    df <- validate_positive_scalar(nk + 1)
  }

  if (is.null(bknots)) {
    # user did not specify boundary knot locations
    bknots <- c(min(x), max(x))
  }

  validate_knots(iknots = iknots, bknots = bknots)

  nlist(iknots, bknots, df = length(iknots) + 1)
}

# Check the knot locations are valid
#
# @param iknots The vector of internal knot locations.
# @param bknots The vector (length 2) of boundary knot locations.
validate_knots <- function(iknots, bknots) {

  check1 <- !is.numeric(bknots)
  check2 <- !length(bknots) == 2L
  if (any(check1, check2)) {
    stop2("'bknots' should a length 2 numeric vector.")
  }

  if (!is.null(iknots)) {
    check3 <- !is.numeric(iknots)
    if (check3) {
      stop2("'iknots' must be numeric.")
    }

    check4 <- any(iknots < bknots[1L])
    check5 <- any(iknots > bknots[2L])
    if (any(check3, check4, check5)) {
      stop2("'iknots' cannot be outside the boundary knots.")
    }
  }
}

# Return the cutpoints for a specified number of (evenly spaced) quantiles
#
# @param x A numeric vector.
# @param nq Integer specifying the desired number of quantiles.
qtile <- function(x, nq = 2) {
  if (nq > 1) {
    probs <- seq(1, nq - 1) / nq
    return(quantile(x, probs = probs))
  } else if (nq == 1) {
    return(NULL)
  } else {
    stop("'nq' must be >= 1.")
  }
}

# Check the number of df is valid for the type of spline basis
#
# @param df The degrees of freedom for the splines.
validate_df <- function(df, spline_type = c("i", "b")) {
  type <- match.arg(spline_type)
  if (type == "i") { # I-splines
    if (df < 5L) {
      stop2("'df' must be >= 5 for I-splines.")
    }
  } else if (type == "b") {
    if (df < 1L) {
      stop2("'df' must be >= 1 for B-splines.")
    }
  }
  df
}

# Return error if not a positive scalar
validate_positive_scalar <- function(x) {
  if (!x > 0) {
    stop2("'x' should be a positive scalar.")
  }
  x
}
