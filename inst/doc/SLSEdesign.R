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
require(tibble)
require(pracma)
require(gridExtra)

## ----Define inputs, cache = TRUE----------------------------------------------
N <- 21
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
res$design

## ----weight-------------------------------------------------------------------
plot_weight(res$design)

## -----------------------------------------------------------------------------
poly3 <- function(xi,theta){
    matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}
design = data.frame(location = c(-1, -0.447, 0.447, 1),
 weight = rep(0.25, 4))
u = seq(-1, 1, length.out = 201)
plot_direction_Dopt(u, design, tt=0, FUN = poly3,
  theta = rep(0, 4))

## -----------------------------------------------------------------------------
poly3 <- function(xi, theta){
  matrix(c(1, xi, xi^2, xi^3), ncol = 1)
}
design = data.frame(location = c(-1, -0.464, 0.464, 1),
                    weight = c(0.151, 0.349, 0.349, 0.151))
u = seq(-1, 1, length.out = 201)
plot_direction_Aopt(u, design, tt=0, FUN = poly3, theta = rep(0,4))

## ----include = FALSE----------------------------------------------------------
options(original) # reset to old settings

