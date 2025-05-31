#' Calculate the loss function of the A- or D-optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#' @param FUN The function to calculate the derivative of the given model.
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param A The calculated covariance matrix
#' @param criterion The criterion to be used for the design, either "D" for D-optimality or "A" for A-optimality. Default is "D".
#'
#' @details This function calculates the loss function of the design problem under the A- or D-optimality. The loss functions under A-, or D-optimality are defined as the trace and log determinant of the inverse of the Fisher information matrix
#'
#' @import CVXR
#'
#' @return The loss of the model at each design points
#'
#' @examples
#' my_design <- data.frame(location = c(0, 180), weight = c(1/2, 1/2))
#' theta <- c(0.05, 0.5)
#' peleg <- function(xi, theta){
#'    deno <- (theta[1] + xi * theta[2])^2
#'    rbind(-xi/deno, -xi^2/deno)
#' }
#' A <- matrix(c(1, 0, 0, 0, 0.2116, 1.3116, 0, 1.3116, 15.462521), byrow = TRUE, ncol = 3)
#' res <- calc_phi(my_design, theta, peleg, 0, A, criterion = "A")
#' res
#' @export

calc_phi <- function(design, theta, FUN, tt, A, criterion = "D") {
  u <- design$location
  N <- length(u)
  q <- length(theta)
  phi <- numeric(N)

  if (criterion == "A") {
    BI <- solve(A)
    C <- diag(c(0, rep(1, q)))
  }

  for (i in 1:N) {
    f <- FUN(u[i], theta)
    I <- rbind(cbind(1, sqrt(tt) * t(f)),
               cbind(sqrt(tt) * f, f %*% t(f)))

    if (criterion == "A") {
      phi[i] <- sum(diag(I %*% BI %*% t(C) %*% C %*% BI))
    } else if (criterion == "D") {
      phi[i] <- sum(diag(solve(A, I)))
    }
  }
  return(phi)
}
