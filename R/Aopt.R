#' Calculate the A-optimal design under the second-order Least squares estimator
#'
#' @param FUN The function to calculate the derivative of the given model.
#' @param N The number of sample points in the design space.
#' @param u The discretized design space.
#' @param tt The level of skewness between 0 to 1 (inclusive). When tt=0, it is equivalent to compute the A-optimal design under the ordinary least squares estimator.
#' @param theta The parameter value of the model.
#' @param show_cvxr_status A boolean variable to indicate whether to show the status of the CVXR optimization. By default, it is set to FALSE.
#'
#' @details This function calculates the A-optimal design and the loss function under the A-optimality. The loss function under A-optimality is defined as the trace of the inverse of the Fisher information matrix
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
#' my_design <- Aopt(N = Npt, u = seq(-1, +1, length.out = Npt),
#'    tt = 0, FUN = poly3, theta = rep(0,4))
#' round(my_design$design, 3)
#' my_design$val
#' @export

Aopt <- function(N, u, tt, FUN, theta, show_cvxr_status = FALSE){

  n <- length(theta)
  C <-  rbind(0, diag(1, n))
  # Set up constraints
  w <- CVXR::Variable(N)
  multi_f <- sapply(u, FUN, theta)
  g1 <- multi_f %*% w
  G2 <- multi_f %*% CVXR::DiagVec(w) %*% t(multi_f)


  B <- CVXR::bmat(list(
    list(1, sqrt(tt)*t(g1)),
    list(sqrt(tt)*g1, G2)
  ))
  obj_val <- CVXR::matrix_frac(C, B)

  my_constraints <- list(w >= 0, sum(w) == 1)

  # Solve the optimization problem
  problem <- CVXR::Problem(CVXR::Minimize(obj_val),
                           constraints = my_constraints)
  res <- CVXR::psolve(problem,
                      solver = "SCS",
                      verbose = show_cvxr_status,
                      reltol = 1e-6,
                      abstol = 1e-6)

  # figure out the location of the design points
  tb <- data.frame(location = u,
                   weight = c(CVXR::value(w)))
  tb <- tb[tb$weight > 1E-3, ]
  # normalize the weights
  tb[, "weight"] <- tb[, "weight"]/sum(tb[, "weight"])
  list(val = res, status = CVXR::status(problem), design = tb)
}
