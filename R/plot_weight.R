#' Plot the weight distribution of the optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#'
#' @details This functions produce a figure that contains the location and their associated weights of the resulted optimal design measures.
#'
#' @importFrom graphics abline legend lines points
#'
#' @return The plot that shows the given optimal design
#'
#' @examples
#' Des = list(location = c(-1, +1), weight = c(0.5, 0.5))
#' plot_weight(Des)
#'
#' @export

plot_weight <- function(design){
  u <- design$location
  w_hat <- design$weight
  S <- u[c(1, length(u))]
  plot(u, w_hat, type = "p", xlim = S, ylim = c(0, 1),
       xlab  = "Design space", ylab = "weight")
  lines(u, w_hat, type = "l")
  points(u[which(w_hat > 1E-4)], w_hat[which(w_hat > 1E-4)],
         col = "red", cex = 2)
}
