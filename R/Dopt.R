#' Calculate the D-optimal design under the SLSE
#'
#' @param FUN The function to calculate the derivative of the given model.
#' @param N The number of sample points in the design space.
#' @param u The discretized design space.
#' @param tt The level of skewness. When tt=0, it is equivalent to compute the D-optimal design under the ordinary least squares estimator.
#' @param theta The parameter value of the model.
#' @param num_iter Maximum number of iteration.
#'
#' @details This function calculates the D-optimal design and the loss function under the D-optimality. The loss function under D-optimality is defined as the log determinant of the inverse of the Fisher information matrix.
#'
#' @import CVXR
#'
#' @return A list that contains 1. Value of the objective function at solution. 2. Status. 3. Optimal design
#'
#' @examples
#' poly3 <- function(xi, theta){
#'   matrix(c(1, xi, xi^2, xi^3), ncol = 1)
#' }
#' Npt <- 101
#' my_design <- Dopt(N = Npt, u = seq(-1, +1, length.out = Npt),
#'    tt = 0, FUN = poly3, theta = rep(0,4), num_iter = 2000)
#' round(my_design$design, 3)
#' my_design$val
#'
#' @export

Dopt <- function(N, u, tt, FUN, theta, num_iter = 1000){
  n <- length(theta)

  w <- CVXR::Variable(N)
  multi_f <- sapply(u, FUN, theta)
  g1 <- multi_f %*% w
  G2 <- multi_f %*% CVXR::diag(w) %*% t(multi_f)# Set up constraints --------------------------------------------------------------------------

  my_constraints <- list(w >= 0, sum(w) == 1)

  B <- rbind(cbind(1, sqrt(tt) * t(g1)),
             cbind(sqrt(tt) * g1, G2))

  # Solve
  objective <- -CVXR::log_det(B)
  problem <- CVXR::Problem(CVXR::Minimize(objective),
                           constraints = my_constraints)
  res <- CVXR::solve(problem, num_iter = num_iter,
                     ignore_dcp = TRUE)

  # figure out the location of the design points
  tb <- data.frame(location = u,
                   weight = c(res$getValue(w)))
  tb <- tb[tb$weight > 1E-2, ]
  list(val = res$value, status = res$status, design = tb)
}

