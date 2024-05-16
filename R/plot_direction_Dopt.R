#' Calculate the loss function of the D-optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#' @param FUN The function to calculate the derivative of the given model.
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param u The discretized design points
#'
#' @details This function produces the figure for the directional derivative of the given D-optimal design of the compact supports. According to the general equivalence theorem, for an optimal design, all the directional derivative should be below zero line.
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#'
#' @examples
#' poly3 <- function(xi, theta){
#'     matrix(c(1, xi, xi^2, xi^3), ncol = 1)
#' }
#' design = data.frame(location = c(-1, -0.447, 0.447, 1),
#'  weight = rep(0.25, 4))
#' u = seq(-1, 1, length.out = 201)
#' plot_direction_Dopt(u, design, tt=0, FUN = poly3,
#'   theta = rep(0, 4))
#'
#'
#' @return The plot of the directional derivative of a D-optimal design
#'
#' @export

plot_direction_Dopt <- function(u, design, tt, FUN, theta){
  N <- length(u)
  S <- u[c(1, N)]
  q <- length(theta)
  sqt <- sqrt(tt)
  g1 <- matrix(0, ncol = 1, nrow = q)
  G2 <- matrix(0, nrow = q, ncol = q)
  for (j in 1:nrow(design)) {
    uj <- design$location[j]
    wj <- design$weight[j]
    f <- FUN(uj, theta)
    g1 <- g1 + wj * f
    G2 <- G2 + wj * tcrossprod(f)
  }
  B <- rbind(cbind(1, sqt*t(g1)),
             cbind(sqt*g1, G2))
  BI <- solve(B)
  phiD <- rep(0, N)
  for(i in 1:N){
    f <- FUN(u[i], theta)
    M <- rbind(cbind(1, sqt * t(f)),
               cbind(sqt*f, tcrossprod(f)))
    phiD[i] <-  sum(diag(BI %*% M))
  }
  # error = max(phiD - (q+1))
  plot(x = u, phiD - (q+1), col = "blue")
  lines(x = u, phiD - (q+1), col = "blue")
  abline(h=0, col="black")
  points(design$location, rep(0,nrow(design)), col = "red", cex = 2)
  legend("bottomright",
         legend=c("Reference Line", expression(d(x, theta)), "discretized point", "Support point"),
         col = c("black", "blue", "blue", "red"),
         lty = c(1, 1, NA, NA), pch = c(NA, NA, 1, 1), cex = 0.8)
}
