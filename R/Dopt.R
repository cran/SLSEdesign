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
#' \donttest{
#' poly3 <- function(xi, theta){
#'   matrix(c(1, xi, xi^2, xi^3), ncol = 1)
#' }
#' my_design <- Dopt(N = 11, u = seq(-1, +1, length.out = 11),
#'    tt = 0, FUN = poly3, theta = rep(0,4), num_iter = 500)
#' my_design$design
#' my_design$val
#' }
#' @export

Dopt <- function(N, u, tt, FUN, theta, num_iter = 1000){
  n <- length(theta)
  g1 <- matrix(0, n, 1)
  G2 <- matrix(0, n, n)
  obj_val <- 0
  # Set up constraints --------------------------------------------------------------------------
  w <- CVXR::Variable(N)
  my_constraints <- list(w >= 0, sum(w) == 1)
  for (i in 1:N) {
    f <- FUN(u[i], theta)
    g1 <- g1 + w[i] * f
    G2 <- G2 + w[i] * f %*% t(f)
  }

  B <- rbind(cbind(1, sqrt(tt) * t(g1)),
             cbind(sqrt(tt) * g1, G2))

  # Solve
  objective <- -CVXR::log_det(B)
  problem <- CVXR::Problem(CVXR::Minimize(objective),
                           constraints = my_constraints)
  res <- CVXR::solve(problem, num_iter = num_iter)
  res$getValue(w)
  # figure out the location of the design points
  tb <- data.frame(location = u,
                   weight = c(res$getValue(w)))
  tb <- tb[tb$weight > 1E-2, ]
  list(val = res$value, status = res$status, design = tb)
}

