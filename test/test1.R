MU1 <- 0.2
ALPHA1 <- 2.2
BETA1 <- 2.4
mHSpec1 <- new("mHSpec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1)

L0 <- 0.2*BETA1/(BETA1-ALPHA1)

res1 <- mHSim(mHSpec1, LAMBDA0=L0, n=400)

inter_arrival <- res1$inter_arrival
mark <- res1$mark
jump_type <- res1$jump_type



loglikelihood(mHSpec1, inter_arrival = inter_arrival, jump_type = jump_type,  mark = mark, LAMBDA0 = L0)

mle <- mHFit(mHSpec1, inter_arrival = inter_arrival)
print(summary(mle))
