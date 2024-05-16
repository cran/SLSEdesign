#' Calculate the loss function of the A-optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#' @param FUN The function to calculate the derivative of the given model.
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param u The discretized design points
#'
#' @details This function produces the figure for the directional derivative of the given A-optimal design of the compact supports. According to the general equivalence theorem, for an optimal design, all the directional derivative should be below zero line.
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#'
#' @return The plot of the directional derivative of a A-optimal design
#'
#' @examples
#' poly3 <- function(xi, theta){
#'   matrix(c(1, xi, xi^2, xi^3), ncol = 1)
#' }
#' design = data.frame(location = c(-1, -0.464, 0.464, 1),
#'                       weight = c(0.151, 0.349, 0.349, 0.151))
#' u = seq(-1, 1, length.out = 201)
#' plot_direction_Aopt(u=u, design=design, tt=0, FUN = poly3, theta = rep(0,4))
#'
#' @export

plot_direction_Aopt <- function(u, design, tt, FUN, theta){
  q <- length(theta)
  N <- length(u)
  S <- u[c(1, N)]
  sqt <- sqrt(tt)
  g1 <- matrix(0, ncol = 1, nrow = q)
  G2 <- matrix(0, nrow = q, ncol = q)
  # Compute the B matrix

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

  C <- pracma::blkdiag(matrix(0), diag(1, q))
  term <- sum(diag(C * BI * t(C)))
  # mini <- min(phi - q)
  ff <- function(x){
    f <- FUN(x, theta)
    M <- rbind(cbind(1, sqt * t(f)),
               cbind(sqt * f, tcrossprod(f)))
    sum(diag(M %*% BI %*% t(C) %*% C %*% BI)) - term
  }

  y <- rep(0, N)
  for(i in 1:N){
    y[i] <- ff(u[i])
  }


  plot(u, y, col = "blue", xlim  = S,
       xlab = "Design space", ylab = expression(d(x, theta)))
  points(u, y, type = "l", xlim  = S, col = "blue")
  abline(h=0, col="black")
  points(design$location, rep(0,nrow(design)), col = "red", cex = 2)
  legend("bottomright",
         legend=c("Reference Line", expression(d(x, theta)), "discretized point", "Support point"),
         col = c("black", "blue", "blue", "red"),
         lty = c(1, 1, NA, NA), pch = c(NA, NA, 1, 1), cex = 0.8)
}
