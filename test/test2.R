
MU2 <- matrix(c(0.2, 0.2), nrow=2)
ALPHA2 <- matrix(c(0.95, 0.82, 0.82, 0.95), nrow=2, byrow=TRUE)
BETA2 <- matrix(rep(2.25, 4), nrow=2, byrow=TRUE)
ETA2 <- matrix(rep(0.0, 4), nrow=2)


mHSpec2 <- new("mHSpec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA <- ETA2)
L0 <- matrix(c(0.03,0.03,0.03,0.03), nrow=2, byrow=TRUE)

res2 <- mHSim(mHSpec2, LAMBDA0=L0, n=1000)

inter_arrival <- res2$inter_arrival
mark <- res2$mark
jump_type <- res2$jump_type


loglikelihood(mHSpec2, inter_arrival=inter_arrival, jump_type = jump_type, LAMBDA0 = L0)
mle <- mHFit(mHSpec2, inter_arrival = inter_arrival, jump_type = jump_type,  LAMBDA0 = L0)
print(summary(mle))
