# This code implements the markes Hawkes process simulation.

setGeneric("mhsim", function(object, ...) standardGeneric("mhsim"))

#' Simulate an n-dimensional (marked) Hawkes process
#'
#' The method simulate n-dimensional marked Hawkes processes.
#' The object \code{\link{mhspec-class}} contains the parameter values such as alpha, beta, eta.
#' The mark (jump) structure may or may not be included.
#' It returns an object of class 'mhreal'.
#'
#' @param object \code{\link{mhspec-class}}. This object includes the parameter values.
#' @param LAMBDA0 the starting values of lambda (intensity process). Must have the same dimensional matrix (n by n) with the parameters in mhspec.
#' @param n the number of observations.
#'
#' @examples
#' # Example for one dimensional Hawkes process (without mark)
#' # Simple simulation for example
#' mhsim()
#' mhsim(dimens = 2)
#'
#' # Define a one-dimensional model.
#' MU1 <- 0.3; ALPHA1 <- 1.5; BETA1 <- 2
#' mhspec1 <- new("mhspec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1)
#' # Simulate with mhsim funciton.
#' res1 <- mhsim(mhspec1,  n=100)
#' summary(res1)
#' as.matrix(res1)
#'
#' # Example for two dimensional Hawkes process
#' # Define a two-dimensional model.
#' MU2 <- matrix(c(0.2), nrow = 2)
#' ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
#' BETA2 <- matrix(c(2.25, 2.25, 2.25, 2.25), nrow = 2, byrow=TRUE)
#' ETA2 <- matrix(c(0.19, 0.19, 0.19, 0.19), nrow = 2, byrow=TRUE)
#' mark2 <- function(n,...) rgeom(n, 0.65) + 1   # mark size follows a geometric distribution
#' mhspec2 <- new("mhspec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2, mark =mark2)
#' # Simulate with mhsim function.
#' LAMBDA0 <- matrix(c(0.1, 0.1, 0.1, 0.1), nrow = 2, byrow=TRUE)
#' res2 <- mhsim(mhspec2, LAMBDA0 = LAMBDA0, n = 100)
#' class(res2)
#' summary(res2)
#' as.matrix(res2)
setMethod(
  f="mhsim",
  signature(object = "mhspec"),
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
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")

      if (dimens == 1){
        #default LAMBDA0 with dimesion 1
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

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(LAMBDA0, c, numeric(length = n - 1)), ncol = dimens^2)
    rowSums_LAMBDA0 <- rowSums(matrix(LAMBDA0, nrow=dimens))

    lambda   <- matrix(sapply(MU + rowSums_LAMBDA0, c, numeric(length = n - 1)), ncol = dimens)

    Ng <- matrix(numeric(length = dimens * n), ncol = dimens)
    N  <- matrix(numeric(length = dimens * n), ncol = dimens)

    mark_type <- numeric(length = n)
    inter_arrival <- numeric(length = n)
    mark <- numeric(length = n)

    # Set column names
    colnames(lambda) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    colnames(lambda_component) <- paste0("lambda", indxM, t(indxM))
    colnames(N)  <- paste0("N", 1:dimens)
    colnames(Ng) <- paste0("Ng", 1:dimens)

    # Exact method
    for (i in 2:n) {

      # Generate candidate arrivals
      # arrival due to mu
      candidate_arrival <- stats::rexp(dimens, rate = MU)
      current_LAMBDA <- matrix(as.numeric(lambda_component[i-1, ]), nrow = dimens, byrow = TRUE)

      # arrival due to components

      matrixD <- 1 + BETA * log(stats::runif(dimens^2)) / current_LAMBDA
      candidate_arrival <- cbind(candidate_arrival, -1 / BETA * log(pmax(matrixD, 0)))

      # The minimum is inter arrival time
      inter_arrival[i] <- min(candidate_arrival)
      minIndex <- which(candidate_arrival == inter_arrival[i], arr.ind = TRUE) #row and col

      mark_type[i] <- minIndex[1]  # row

      # lambda decayed due to time, impact due to mark is not added yet
      decayled_lambda <- current_LAMBDA * exp(-BETA * inter_arrival[i])
      lambda_component[i, ] <- t(decayled_lambda)
      lambda[i, ] <- MU + rowSums(decayled_lambda)

      # generate one random number
      mark[i] <- object@mark(n = 1, i = i, N = N, Ng = Ng,
                             lambda = lambda, lambda_component = lambda_component,
                             mark_type = mark_type)

      Ng[i, ] <- Ng[i-1, ]
      Ng[i, mark_type[i]] <- Ng[i-1, mark_type[i]] + 1
      N[i, ] <- N[i-1, ]
      N[i, mark_type[i]] <- N[i-1, mark_type[i]] + mark[i]

      # update lambda
      if (dimens == 1) {
        Impact <- ALPHA * (1 + (mark[i] - 1) * ETA )
      } else {
        Impact <- matrix(rep(0, dimens^2), nrow = dimens)
        Impact[ , mark_type[i]] <- ALPHA[ , mark_type[i]] * (1 + (mark[i] - 1) * ETA[ , mark_type[i]])
      }

      # new_LAMBDA = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      # lambda_component = {"lambda11", "lambda12", ..., "lambda21", "lambda22", ...}
      #
      # Impact is added.
      new_lambda <- decayled_lambda + Impact
      lambda_component[i, ] <- t(new_lambda)
      lambda[i, ] <- MU + rowSums(new_lambda)
    }


    realization <- list(object, inter_arrival, cumsum(inter_arrival), mark_type, mark, N, Ng, lambda, lambda_component)
    names(realization) <- c("mhspec", "inter_arrival", "arrival", "mark_type", "mark", "N", "Ng", "lambda", "lambda_component")
    class(realization) <- c("mhreal")

    return(realization)
  }
)

setMethod(
  "mhsim",
  signature("missing"),
  function(object, dimens = 1, n = 1000) {
    # default values

      if (dimens == 1){
      MU1 <- 0.2
      ALPHA1 <- 1.5
      BETA1 <- 2

      mhspec1 <- methods::new("mhspec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1)

      mhsim(mhspec1, n = n)

    } else if (dimens == 2){

      MU2 <- matrix(c(0.2), nrow = 2)
      ALPHA2 <- matrix(c(0.7, 0.9, 0.9, 0.7), nrow = 2, byrow=TRUE)
      BETA2 <- matrix(c(2, 2, 2, 2), nrow = 2, byrow=TRUE)

      mhspec2 <- methods::new("mhspec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2)
      mhsim(mhspec2, n = n)

    } else {

      stop("One or two dimesinoal models is supported for default simulation.")
    }


  }

)

