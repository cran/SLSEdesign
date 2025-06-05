#' Verify the optimality condition for an optimal design (A-, c- or D-optimality)
#'
#' @param u The discretized design points
#' @param design The optimal design containing the design points and the associated weights
#' @param tt The level of skewness
#' @param FUN The function to calculate the derivative of the given model
#' @param theta The parameter value of the model
#' @param criterion The optimality criterion: one of "A", "c", or "D"
#' @param cVec c vector used to determine the combination of the parameters. This is only used in c-optimality
#'
#' @details This function visualizes the directional derivative under A-, c-, or D-optimality using the general equivalence theorem. For an optimal design, the directional derivative should not exceed the reference threshold
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

plot_dispersion <- function(u, design, tt, FUN, theta, criterion = "D", cVec = rep(0, length(theta))) {
  N <- length(u)
  q <- length(theta)
  sqt <- sqrt(tt)
  x_star <- design$location
  w_star <- design$weight
  n <- length(theta)

  if (criterion == "A") {
    C <-  rbind(0, diag(1, n))
  } else if (criterion == "c") {
    cVec_arg <- c(0, cVec)
  }


  multi_f <- sapply(x_star, FUN, theta)
  g1 <- multi_f %*% w_star
  G2 <- multi_f %*% diag(w_star) %*% t(multi_f)

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
  } else if (criterion == "c") {
    for (i in 1:N) {
      f <- FUN(u[i], theta)
      M <- rbind(cbind(1, sqt * t(f)),
                 cbind(sqt * f, tcrossprod(f)))
      y[i] <- t(cVec_arg) %*% BI %*% M %*% BI %*% cVec_arg - t(cVec_arg) %*% BI %*% cVec_arg
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
