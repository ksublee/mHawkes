# This code implements the markes Hawkes process simulation.

setGeneric("mHSim", function(object, ...) standardGeneric("mHSim"))

#' Simulate an n-dimensional marked Hawkes process
#'
#' This is a generic function.
#' The mark (jump) structure may or may not be included.
#' It returns an object of class 'mHreal'.
#'
#' @param object mHSpec
#' @param LAMBDA0 The starting values of lambda. Must have the same dimensional matrix (n by n) with mHSpec.
#' @param n The number of observations.
#'
#' @examples
#' # example for one dimensional Hawkes process (without mark)
#' # Define the model.
#' MU1 <- 0.3; ALPHA1 <- 1.5; BETA1 <- 2
#' mHSpec1 <- new("mHSpec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1)
#' # Simulate with mHSim funciton.
#' res1 <- mHSim(mHSpec1,  n=100)
#' summary(res1)
#' as.matrix(res1)
#'
#' # example for two dimensional Hawkes process
#' # Define the model.
#' MU2 <- matrix(c(0.2), nrow = 2)
#' ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
#' BETA2 <- matrix(c(2.25, 2.25, 2.25, 2.25), nrow = 2, byrow=TRUE)
#' ETA2 <- matrix(c(0.19, 0.19, 0.19, 0.19), nrow = 2, byrow=TRUE)
#' JUMP2 <- function(n,...) rgeom(n, 0.65) + 1   # mark size follows a geometric distribution
#' mHSpec2 <- new("mHSpec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2, Jump =JUMP2)
#' # Simulate with mHSim function.
#' LAMBDA0 <- matrix(c(0.1, 0.1, 0.1, 0.1), nrow = 2, byrow=TRUE)
#' res2 <- mHSim(mHSpec2, LAMBDA0 = LAMBDA, n = 100)
#' summary(res2)
#' as.matrix(res2)
setMethod(
  f="mHSim",
  signature(object = "mHSpec"),
  definition = function(object, LAMBDA0 = NULL, n = 1000){

    # dimension of Hawkes process
    dimens <- length(object@MU)

    if ( dimens >= 2 ){

      if( is.matrix(LAMBDA0) && nrow(LAMBDA0) != ncol(LAMBDA0) )
        stop("LAMBDA0 matrix should be n by n.")

      if( is.matrix(LAMBDA0) && dimens != nrow(LAMBDA0) )
        stop("Check the dimension of LAMBDA0 matrix.")

    }

    # parameter setting

    MU <- object@MU
    ALPHA <- object@ALPHA
    BETA <- object@BETA
    ETA <- object@ETA

    # default LAMBDA0
    if(is.null(LAMBDA0)) {

      if (dimens == 1){
        LAMBDA0 <- (MU * BETA / (BETA - ALPHA) - MU) / 2
      } else if (dimens == 2) {

        #default LAMBDA0 with dimesion 2
        LAMBDA0 <- matrix(rep(0, dimens^2), nrow=dimens)

        H <- ALPHA[1,1]*BETA[1,2]*BETA[2,1]*(ALPHA[2,2] - BETA[2,2]) - BETA[1,1]*(ALPHA[2,2]*BETA[1,2]*BETA[2,1] + ALPHA[1,2]*ALPHA[2,1]*BETA[2,2] - BETA[1,2]*BETA[2,1]*BETA[2,2])

        LAMBDA0[1, 1] <- ALPHA[1,1]*BETA[2,1]* ((BETA[2,2] - ALPHA[2,2])*BETA[1,2]*MU[1] + ALPHA[1,2]*BETA[2,2]*MU[2]) / H
        LAMBDA0[1, 2] <- ALPHA[1,2]*BETA[2,2]* ((BETA[1,1] - ALPHA[1,1])*BETA[2,1]*MU[2] + ALPHA[2,1]*BETA[1,1]*MU[1]) / H
        LAMBDA0[2, 1] <- ALPHA[2,1]*BETA[1,1]* ((BETA[2,2] - ALPHA[2,2])*BETA[1,2]*MU[1] + ALPHA[1,2]*BETA[2,2]*MU[2]) / H
        LAMBDA0[2, 2] <- ALPHA[2,2]*BETA[1,2]* ((BETA[1,1] - ALPHA[1,1])*BETA[2,1]*MU[2] + ALPHA[2,1]*BETA[1,1]*MU[1]) / H


      } else {
        # for higher dimension, default LAMBDA0 will be ...
        lamb0 <- (MU[1]*BETA[1]/(BETA[1]- sum(ALPHA[1:dimens])) - MU[1])/2
        LAMBDA0 <- matrix(rep(lamb0, dimens^2), nrow=dimens, byrow=TRUE)
      }
    }

    # We need an appropriate algorithm for buffer size n.

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(LAMBDA0, c, numeric(length = n - 1)), ncol = dimens^2)
    rowSums_LAMBDA0 <- rowSums(matrix(LAMBDA0, nrow=dimens))

    lambda_process   <- matrix(sapply(MU + rowSums_LAMBDA0, c, numeric(length = n - 1)), ncol = dimens)

    Ng <- matrix(numeric(length = dimens * n), ncol = dimens)
    N  <- matrix(numeric(length = dimens * n), ncol = dimens)

    jump_type <- numeric(length = n)
    inter_arrival <- numeric(length = n)
    mark <- numeric(length = n)


    # Exact method
    for (k in 1:(n-1)) {

      # Generate candidate arrivals
      # arrival due to mu
      candidate_arrival <- rexp(dimens, rate = MU)
      current_LAMBDA <- matrix(as.numeric(lambda_component[k, ]), nrow = dimens, byrow = TRUE)

      # arrival due to components

      matrixD <- 1 + BETA * log(runif(dimens^2)) / current_LAMBDA
      candidate_arrival <- cbind(candidate_arrival, -1 / BETA * log(pmax(matrixD, 0)))

      # The minimum is inter arrival time
      inter_arrival[k+1] <- min(candidate_arrival)
      minIndex <- which(candidate_arrival == inter_arrival[k+1], arr.ind = TRUE) #row and col

      jumpType <- minIndex[1]  # row
      jump_type[k+1] <- jumpType

      mark[k+1] <- object@Jump(n = 1)  # generate one random number

      Ng[k+1, ] <- Ng[k, ]
      Ng[k+1, jumpType] <- Ng[k, jumpType] + 1
      N[k+1, ] <- N[k, ]
      N[k+1, jumpType] <- N[k, jumpType] + mark[k+1]

      # update lambda
      if (dimens == 1) {
        Impact <- ALPHA * (1 + (mark[k+1] - 1) * ETA )
      } else {
        Impact <- matrix(rep(0, dimens^2), nrow = dimens)
        Impact[ , jumpType] <- ALPHA[ , jumpType] * (1 + (mark[k+1] - 1) * ETA[ , jumpType])
      }

      new_LAMBDA <- current_LAMBDA * exp(-BETA * inter_arrival[k+1]) + Impact

      # new_LAMBDA = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      # lambda_component = {"lambda11", "lambda12", ..., "lambda21", "lambda22", ...}
      lambda_component[k+1, ] <- t(new_LAMBDA)
      lambda_process[k+1, ] <- MU + rowSums(new_LAMBDA)
    }


    # Set column names
    colnames(lambda_process) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))
    colnames(N)  <- paste0("N", 1:dimens)
    colnames(Ng) <- paste0("Ng", 1:dimens)

    realization <- list(object, inter_arrival, cumsum(inter_arrival), jump_type, mark, N, Ng, lambda_process, lambda_component)
    names(realization) <- c("mHSpec", "inter_arrival", "arrival", "jump_type", "mark", "N", "Ng", "lambda", "lambda_component")
    class(realization) <- c("mHreal")

    return(realization)
  }
)
