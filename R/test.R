
G <- distr::Dirac(location = 1)

MU2 <- matrix(c(0.2, 0.2), nrow=2)
ALPHA2 <- matrix(c(0.95, 0.82, 0.82, 0.95), nrow=2, byrow=TRUE)
BETA2 <- matrix(rep(2.25, 4), nrow=2, byrow=TRUE)
ETA2 <- matrix(rep(0.0, 4), nrow=2)


mHSpec2 <- new("mHSpec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA <- ETA2)
L0 <- matrix(c(0.03,0.03,0.03,0.03), nrow=2, byrow=TRUE)
#mHP <- new("mHProcess", spec=mHSpec1, LAMBDA0=LAMBDA0)
res <- mHSim(mHSpec2, LAMBDA0=L0, n=1000)

plot(res$arrival, res$N$N1 - res$N$N2, 's')

plot(res$arrival, res$lambda$lambda1, 'l')

N1 <- res$N$N1
N2 <- res$N$N2

Nu <- c( N1[2:length((N1))]-N1[1:(length(N1)-1)])
Nd <- c( N2[2:length((N1))]-N2[1:(length(N1)-1)])
direction <- Nu + 2*Nd
interarrival <- res$arrival[2:length(res$arrival)] - res$arrival[1:(length(res$arrival)-1)]
mark_size <- rep(1, length(direction))

loglikelihood(mHSpec2, c(0,interarrival), c(0,direction), c(1,mark_size), LAMBDA0 = L0)

mle <- mHFit(mHSpec2, inter_arrival = c(0,interarrival), jump_type = c(0,direction),  LAMBDA0 = L0)
print(summary(mle))

logLikFun(c(mu = 0.02, alpha_s = 0.95, alpha_c = 0.82, beta = 2.25, eta =0))



mle<-maxLik::maxLik(logLik=logLikFun, start=c(mu = 0.02, alpha_s = 0.95, alpha_c = 0.82, beta = 2.25, eta =0))
print(summary(mle))

templlh <- function(param){

  MU <- matrix(rep(param[[1]], 2), nrow=2)
  ALPHA <- matrix(c(param[[2]], param[[3]], param[[3]], param[[2]]), nrow=2, byrow=TRUE)
  BETA <- matrix(rep(param[[4]], 4), nrow=2, byrow=TRUE)
  ETA <- matrix(rep(param[[5]], 4), nrow=2)
  G <- distr::Dirac(location = 1)

  mHSpec1 <- new("mHSpec", MU=MU, ALPHA=ALPHA, BETA=BETA, ETA=ETA, Jump=G)

  LAMBDA0 <- matrix(c(0.03,0.03,0.03,0.03), nrow=2, byrow=TRUE)
  mHP <- new("mHProcess", spec=mHSpec1, LAMBDA0=LAMBDA0)

  llh <- loglikelihood(mHP, c(0,interarrival), c(0,direction), c(1,mark_size), LAMBDA0)
  print(param)
  return(llh)

}

templlh(c(mu = 0.02, alpha_s = 0.95, alpha_c = 0.82, beta = 2.25, eta =0))

#constraint
A <- matrix(c(0, -1, -1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0), 5, 5, byrow = TRUE)
B <- matrix(c(0,0,0,0,0),5,1)

mle<-maxLik::maxLik(logLik=templlh, start=c(mu = 0.02, alpha_s = 0.95, alpha_c = 0.82, beta = 2.25, eta =0), constraint=list(ineqA=A, ineqB=B))

print(summary(mle))

