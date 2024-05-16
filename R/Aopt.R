#' Calculate the A-optimal design under the second-order Least squares estimator
#'
#' @param FUN The function to calculate the derivative of the given model.
#' @param N The number of sample points in the design space
#' @param u The discretized design space
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param num_iter Maximum number of iteration
#'
#' @details This function calculates the loss function of the design problem under the A-optimality. The loss function under A-optimality is defined as the trace of the inverse of the Fisher information matrix
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#' @importFrom tibble tibble
#'
#' @return A list that contains 1. Value of the objective function at solution. 2. Status. 3. Optimal design
#'
#' @examples
#' \donttest{
#' poly1 <- function(xi, theta){
#'   matrix(c(1, xi), ncol = 1)
#' }
#' my_design <- Aopt(N = 11, u = seq(-1, +1, length.out = 11),
#'    tt = 0, FUN = poly1, theta = rep(0,2), num_iter = 50)
#' my_design$design
#' my_design$val
#' }
#' @export

Aopt <- function(N, u, tt, FUN, theta, num_iter = 1000){

  n <- length(theta)
  g1 <- matrix(0, n, 1)
  G2 <- matrix(0, n, n)
  obj_val <- 0
  C <- rbind(0, diag(1, n, n))

  w <- CVXR::Variable(N)

  # Set up constraints
  for (i in 1:N) {
    f <- FUN(u[i], theta)
    g1 <- g1 + w[i] * f
    G2 <- G2 + w[i] * f %*% t(f)
  }

  B <- rbind(cbind(1, sqrt(tt) * t(g1)),
             cbind(sqrt(tt) * g1, G2))

  C <- pracma::blkdiag(matrix(0), diag(1, n))

  for(k in 1:n){
    obj_val <- obj_val + CVXR::matrix_frac(C[, k], B)
  }

  my_constraints <- list(w >=0, sum(w) == 1)

  # Solve the optimization problem
  problem <- CVXR::Problem(CVXR::Minimize(obj_val),
                           constraints = my_constraints)
  res <- CVXR::solve(problem, num_iter = num_iter)

  # figure out the location of the design points
  tb <- tibble(location = u,
                   weight = c(res$getValue(w)))
  tb <- tb[tb$weight > 1E-2, ]
  list(val = res$value, status = res$status, design = tb)
}
