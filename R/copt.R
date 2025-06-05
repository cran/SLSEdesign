#' Calculate the c-optimal design under the SLSE with the given combination of the parameters
#'
#' @param FUN The function to calculate the derivative of the given model.
#' @param N The number of sample points in the design space.
#' @param u The discretized design space.
#' @param tt The level of skewness. When tt=0, it is equivalent to compute the c-optimal design under the ordinary least squares estimator.
#' @param theta The parameter value of the model.
#' @param num_iter Maximum number of iteration.
#' @param cVec c vector used to determine the combination of the parameters
#'
#' @details This function calculates the c-optimal design and the loss function under the c-optimality. The loss function under c-optimality is defined as the log determinant of the inverse of the Fisher information matrix.
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
#' my_design <- copt(N = Npt, u = seq(-1, +1, length.out = Npt),
#'    tt = 0, FUN = poly3, theta = rep(0,4), num_iter = 2000,
#'    cVec = c(0,1,1,1))
#' round(my_design$design, 3)
#' my_design$val
#'
#' @export

copt <- function(N, u, tt, FUN, theta, num_iter = 1000, cVec){
  n <- length(theta)
  w <- CVXR::Variable(N)
  multi_f <- sapply(u, FUN, theta)
  g1 <- multi_f %*% w # # g = sum w_i f(u_i)
  G2 <- multi_f %*% CVXR::diag(w) %*% t(multi_f) # G = sum w_i f(u_i) f(u_i)^T


  # Construct B matrix
  B <- rbind(cbind(1, sqrt(tt) * t(g1)),
             cbind(sqrt(tt) * g1, G2))

  # Pad cVec with 0 at the top to match dimension of B
  cVec_aug <- c(0, cVec)  # size: (n+1) x 1

  # Objective: c^T B^{-1} c
  objective <- CVXR::matrix_frac(cVec_aug, B)

  # Constraints
  my_constraints <- list(w >= 0, sum(w) == 1)

  # Solve the problem
  problem <- CVXR::Problem(Minimize(objective), constraints = my_constraints)
  # res <- solve(problem, num_iter = num_iter)
  res <- solve(problem, eps = 1e-8, max_iters = num_iter)

  # Extract non-negligible weights
  tb <- data.frame(location = u,
                   weight = as.vector(res$getValue(w)))
  tb <- tb[tb$weight > 1E-2, ]

  # normalize the weights
  tb[, "weight"] <- tb[, "weight"]/sum(tb[, "weight"])

  list(val = res$value, status = res$status, design = tb)
}
