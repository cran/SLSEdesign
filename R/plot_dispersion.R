#' Verify the optimality condition for an optimal design (A- or D-optimality)
#'
#' @param u The discretized design points
#' @param design The optimal design containing the design points and the associated weights
#' @param tt The level of skewness
#' @param FUN The function to calculate the derivative of the given model
#' @param theta The parameter value of the model
#' @param criterion The optimality criterion: one of "A" or "D"
#'
#' @details This function visualizes the directional derivative under A- or D-optimality using the general equivalence theorem. For an optimal design, the directional derivative should not exceed the reference threshold (0 for A-optimality, q+1 for D-optimality).
#'
#' @return A plot verifying the general equivalence condition for the specified optimal design
#'
#' @examples
#' poly3 <- function(xi, theta){
#'   matrix(c(1, xi, xi^2, xi^3), ncol = 1)
#' }
#' design_A <- data.frame(location = c(-1, -0.464, 0.464, 1),
#'                        weight = c(0.151, 0.349, 0.349, 0.151))
#' design_D = data.frame(location = c(-1, -0.447, 0.447, 1),
#'                       weight = rep(0.25, 4))
#' u <- seq(-1, 1, length.out = 201)
#' par(mfrow = c(2,2))
#' plot_dispersion(u, design_A, tt = 0, FUN = poly3, theta = rep(0, 4), criterion = "A")
#' plot_dispersion(u, design_A, tt = 0, FUN = poly3, theta = rep(0, 4), criterion = "D")
#'
#' plot_dispersion(u, design_D, tt = 0, FUN = poly3, theta = rep(0, 4), criterion = "A")
#' plot_dispersion(u, design_D, tt = 0, FUN = poly3, theta = rep(0, 4), criterion = "D")
#'
#' @import CVXR
#' @importFrom graphics abline legend lines points
#' @export

plot_dispersion <- function(u, design, tt, FUN, theta, criterion = "D") {
  N <- length(u)
  q <- length(theta)
  sqt <- sqrt(tt)
  g1 <- matrix(0, ncol = 1, nrow = q)
  G2 <- matrix(0, nrow = q, ncol = q)

  # Compute moments from design
  for (j in 1:nrow(design)) {
    uj <- design$location[j]
    wj <- design$weight[j]
    f <- FUN(uj, theta)
    g1 <- g1 + wj * f
    G2 <- G2 + wj * tcrossprod(f)
  }

  # Compute Fisher information matrix and its inverse
  B <- rbind(cbind(1, sqt * t(g1)),
             cbind(sqt * g1, G2))
  BI <- solve(B)

  y <- numeric(N)
  if (criterion == "A") {
    C <- diag(c(0, rep(1, q)))
    term <- sum(diag(C %*% BI %*% t(C)))  # constant to subtract
    for (i in 1:N) {
      f <- FUN(u[i], theta)
      M <- rbind(cbind(1, sqt * t(f)),
                 cbind(sqt * f, tcrossprod(f)))
      y[i] <- sum(diag(M %*% BI %*% t(C) %*% C %*% BI)) - term
    }
  } else if (criterion == "D") {
    for (i in 1:N) {
      f <- FUN(u[i], theta)
      M <- rbind(cbind(1, sqt * t(f)),
                 cbind(sqt * f, tcrossprod(f)))
      y[i] <- sum(diag(BI %*% M)) - (q + 1)
    }
  }

  # Plot
  plot(u, y, col = "blue", type = "l",
       xlab = "Design space", ylab = expression(d(x, theta)),
       main = paste(criterion, "-optimality directional derivative"), ylim = range(y))
  abline(h = 0, col = "black")
  points(design$location, rep(0, nrow(design)), col = "red", cex = 2)
  legend("bottomright",
         legend = c("Reference Line", expression(d(x, theta)), "Discretized point", "Support point"),
         col = c("black", "blue", "blue", "red"),
         lty = c(1, 1, NA, NA), pch = c(NA, NA, 1, 1), cex = 0.8)
}
