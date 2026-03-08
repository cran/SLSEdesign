## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)
original <- options(digits = 3)

## ----load package-------------------------------------------------------------
# required dependencies
require(SLSEdesign)
require(CVXR)

## ----Define inputs, cache = TRUE----------------------------------------------
N <- 101
S <- c(-1, 1)
tt <- 0
theta <- rep(1, 4)

poly3 <- function(xi,theta){
    matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}

u <- seq(from = S[1], to = S[2], length.out = N)

res <- Aopt(N = N, u = u, tt = tt, FUN = poly3, 
            theta = theta)

## ----output-------------------------------------------------------------------
res$val
res$status
round(res$design, 4)

## ----weight-------------------------------------------------------------------
plot_weight(res$design)

## ----D-optimality-------------------------------------------------------------
poly3 <- function(xi,theta){
    matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}
design <- data.frame(location = c(-1, -0.447, 0.447, 1),
 weight = rep(0.25, 4))
u = seq(-1, 1, length.out = 201)
plot_dispersion(u, design, tt = 0, FUN = poly3,
  theta = rep(0, 4), criterion = "D")

## ----A-optimality-------------------------------------------------------------
poly3 <- function(xi, theta){
  matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}
design <- data.frame(location = c(-1, -0.464, 0.464, 1),
                     weight = c(0.151, 0.349, 0.349, 0.151))
u = seq(-1, 1, length.out = 201)
plot_dispersion(u, design, tt = 0, 
                FUN = poly3, theta = rep(0,4), criterion = "A")

## ----c-optimality-------------------------------------------------------------
my_peleg <- function(xi, theta) {
  deno <- (theta[1] + theta[2]*xi)
  matrix(c(-xi/deno^2, -xi^2/deno^2), ncol = 1)
}
Npt <- 1001
my_u <- seq(0, 100, length.out = Npt)
my_theta <- c(0.5, 0.05)
my_cVec <- c(1, 1)
my_design <- copt(
  N = Npt, u = my_u,
  tt = 0, FUN = my_peleg, theta = my_theta,
  cVec = my_cVec
)

plot_dispersion(my_u, my_design$design, tt = 0, 
                FUN = my_peleg, theta = my_theta, 
                criterion = "c", cVec = my_cVec)

## ----include = FALSE----------------------------------------------------------
options(original) # reset to old settings

