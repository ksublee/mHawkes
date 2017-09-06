## ---- eval=FALSE---------------------------------------------------------
#  # install.packages("devtools")  #if devtools is not installed

## ---- eval=FALSE---------------------------------------------------------
#  #devtools::install_github("ksublee/mHawkes", build_vignettes=FALSE, force=TRUE)

## ------------------------------------------------------------------------
library("mHawkes")

## ------------------------------------------------------------------------
MU1 <- 0.3
ALPHA1 <- 1.5
BETA1 <- 2
mHSpec1 <- new("mHSpec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1)
show(mHSpec1)

## ------------------------------------------------------------------------
res1 <- mHSim(mHSpec1,  n=20)

## ------------------------------------------------------------------------
res1

## ------------------------------------------------------------------------
# plot(res1$arrival, res1$N[,'N1'], 's', xlab="t", ylab="N")

## ------------------------------------------------------------------------
# plotlambda(res1$arrival, res1$N[,'N1'], BETA1)

## ------------------------------------------------------------------------
inter_arrival1 <- res1$inter_arrival
logLik(mHSpec1, inter_arrival = inter_arrival1)

## ---- warning=FALSE------------------------------------------------------
mHSpec0 <- mHSpec1
mle <- mHFit(mHSpec0, inter_arrival = inter_arrival1)
summary(mle)

## ------------------------------------------------------------------------
ETA1 <- 0.15
JUMP1 <- function(n,...) rgeom(n, 0.65) + 1

mHSpec1 <- new("mHSpec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1, ETA=ETA1, Jump =JUMP1)

## ------------------------------------------------------------------------
res1 <- mHSim(mHSpec1,  n=10)

## ------------------------------------------------------------------------
# plot(res1$arrival, res1$N[,'N1'], 's', xlab="t", ylab="N")

## ------------------------------------------------------------------------
MU2 <- matrix(c(0.2), nrow = 2)
ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
BETA2 <- matrix(c(2.25, 2.25, 2.25, 2.25), nrow = 2, byrow=TRUE)
ETA2 <- matrix(c(0.19, 0.19, 0.19, 0.19), nrow = 2, byrow=TRUE)
JUMP2 <- function(n,...) rgeom(n, 0.65) + 1
LAMBDA0 <- matrix(c(0.1, 0.1, 0.1, 0.1), nrow = 2, byrow=TRUE)
mHSpec2 <- new("mHSpec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2, Jump =JUMP2)

## ------------------------------------------------------------------------
mHSpec2

## ------------------------------------------------------------------------
res2 <- mHSim(mHSpec2,  n=100)
summary(res2)

## ------------------------------------------------------------------------
# plot(res2$arrival[1:10], res2$N[1:10,1], 's')

## ------------------------------------------------------------------------
# plot(res2)

## ------------------------------------------------------------------------
# plotlambda(res2$arrival[1:10], res2$lambda[1:10,1], BETA2[1,1])

## ------------------------------------------------------------------------
inter_arrival1 <- res1$inter_arrival

## ---- warning=FALSE------------------------------------------------------
mle <- mHFit(mHSpec1, inter_arrival = inter_arrival1)

## ------------------------------------------------------------------------
inter_arrival2 <- res2$inter_arrival
mark2 <- res2$mark
jump_type2 <- res2$jump_type

## ------------------------------------------------------------------------
logLik(mHSpec2,  inter_arrival = inter_arrival2, jump_type = jump_type2, mark = mark2)

## ---- warning=FALSE, error=FALSE-----------------------------------------
mHSpec0 <- mHSpec2
mle <- mHFit(mHSpec0, inter_arrival = inter_arrival2, jump_type = jump_type2, mark = mark2)
summary(mle)

## ------------------------------------------------------------------------
MU0 <- matrix(c(0.2, 0.21), nrow = 2)
ALPHA0 <- matrix(c(0.75, 0.75, 0.75, 0.75), nrow = 2, byrow=TRUE)
BETA0 <- matrix(c(2.25, 2.251, 2.251, 2.25), nrow = 2, byrow=TRUE)
ETA0 <- matrix(c(0.19, 0.19, 0.19, 0.19), nrow = 2, byrow=TRUE)
JUMP0 <- function(n,...) rgeom(n, 0.65) + 1

mHSpec0 <- new("mHSpec", MU=MU0, ALPHA=ALPHA0, BETA=BETA0, ETA=ETA0, Jump =JUMP0)

## ---- warning=FALSE, error=FALSE-----------------------------------------

mle <- mHFit(mHSpec0, inter_arrival = inter_arrival2, jump_type = jump_type2, mark = mark2)
summary(mle)

