#' Calculate the loss function of the D-optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#' @param FUN The function to calculate the derivative of the given model.
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param A The calculated covariance matrix
#'
#' @details This function calculates the loss function of the design problem under the D-optimality. The loss function under D-optimality is defined as the log determinant of the inverse of the Fisher information matrix
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
#' res <- calc_phiA(my_design, theta, peleg, 0, A)
#' res
#' @export

calc_phiD <- function(design, theta, FUN, tt, A){
  u <- design$location
  w_hat <- design$weight
  N <- length(u)
  q <- length(theta)
  g1 <- matrix(0, q, 1)
  G2 <- matrix(0, q, q)
  phi_D <- rep(0, N)

  for( i in 1:N){
    f <- FUN(u[i], theta)
    I <- rbind(cbind(1, sqrt(tt) * t(f)),
               cbind(sqrt(tt) * f, f %*% t(f)))
    phi_D[i] = sum(diag(solve(A, I)));
  }
  phi_D
}
