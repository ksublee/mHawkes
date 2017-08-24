#' Compute the loglikelihood function
#'
#' This is a generic function.
#' The loglikelihood of the ground process of the Hawkes model.
#' (The estimation for jump distribution is not provided.)
#'
#' @param inter_arrival Inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param jump_type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param LAMBDA0 The starting values of lambda. Must have the same dimensional matrix (n by n) with \code{mHSpec}.
#'
#' @examples
#' # construct a mHSpec
#' MU1 <- 0.2; ALPHA1 <- 1.0; BETA1 <- 2; ETA1 <- 0.2
#' JUMP1 <- function(n,...) rgeom(n, 0.7) + 1
#' mHSpec1 <- new("mHSpec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1, ETA=ETA1, Jump =JUMP1)
#' # simualte a path
#' res1 <- mHSim(mHSpec1,  LAMBDA0 = MU1, n=1000)
#' inter_arrival <- res1$inter_arrival
#' mark <- res1$mark
#' # compute a loglikelihood function with parameter values in mHSpec1
#' # LAMBDA0 = MU1 is a naive way of starting point choice.
#' logLik(mHSpec1, LAMBDA0 = MU1, inter_arrival = inter_arrival, mark = mark)
#'
#' @seealso \code{\link{mHSpec-class}}, \code{\link{mHFit,mHSpec-method}}
setMethod(
  f="logLik",
  signature(object="mHSpec"),
  function(object, inter_arrival, jump_type=NULL, mark=NULL, LAMBDA0=NULL){

    # When the mark sizes are not provided, the jumps are all unit jumps
    if(is.null(mark)) {
      mark <- rep(1, length(inter_arrival))
    }

    # dimension of Hawkes process
    dimens <- length(object@MU)

    # if dimens == 1 and jump_type is not provided, then all jump_type is 1.
    if(dimens==1 & is.null(jump_type)) {
      jump_type <- rep(1, length(inter_arrival))
    } else if (dimens != 1 & is.null(jump_type)) {
      stop("The argument jump_type should be provided.")
    }

    # parameter setting

    MU <- object@MU
    ALPHA <- object@ALPHA
    BETA <- object@BETA
    ETA <- object@ETA

    # default LAMBDA0
    if(is.null(LAMBDA0)) {

      if (dimens == 1){
        lamb0 <- (MU[1]*BETA[1]/(BETA[1]- ALPHA[1]) - MU[1])/2
        LAMBDA0 <- matrix(rep(lamb0, dimens^2), nrow=dimens, byrow=TRUE)

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


    # n is length(inter_arrival) - 1
    n <- length(inter_arrival) - 1

    #if (dimens==1) rowSums_LAMBDA0 <- LAMBDA0
    #else rowSums_LAMBDA0 <- rowSums(LAMBDA0)
    rowSums_LAMBDA0 <- rowSums(matrix(LAMBDA0, nrow=dimens))

    sum_log_lambda <- 0
    sum_integrated_lambda_component <- 0

    for (k in 1:n) {
      # current_LAMBDA <- matrix(lambda_component[k, ], nrow = dimens, byrow = TRUE)
      if (k == 1) current_LAMBDA <- LAMBDA0
      else current_LAMBDA <- new_LAMBDA  # LAMBDA determined in the previous loop

      # update lambda
      if (dimens == 1) {
        Impact <- ALPHA * (1 + (mark[k+1] - 1) * ETA )
      } else {
        Impact <- matrix(rep(0, dimens^2), nrow = dimens)
        Impact[ , jump_type[k+1]] <- ALPHA[ , jump_type[k+1]] * (1 + (mark[k+1] - 1) * ETA[ , jump_type[k+1]])
      }


      decayed <- exp(-BETA * inter_arrival[k+1])
      decayed_LAMBDA <- current_LAMBDA * decayed
      new_LAMBDA <- decayed_LAMBDA + Impact

      # sum of integrated_lambda_component
      sum_integrated_lambda_component <- sum_integrated_lambda_component + sum(current_LAMBDA / BETA * ( 1 - decayed ))

      # sum of log lambda when jump occurs
      if (dimens == 1) lambda_lc <- MU + decayed_LAMBDA
      else lambda_lc <- MU + rowSums(decayed_LAMBDA)
      sum_log_lambda <- sum_log_lambda + log(lambda_lc[jump_type[k+1]])

    }

    # log likelihood for ground process
    sum_log_lambda - sum(MU*sum(inter_arrival)) - sum_integrated_lambda_component

  }
)


setGeneric("mHFit", function(object, ...) standardGeneric("mHFit"))

#' Perform a maximum likelihood estimation
#'
#' This function uses \code{\link[maxLik]{maxLik}} for the optimizer.
#'
#'
#' @param object mHSpec, or can be omitted.
#' @param inter_arrival Inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param jump_type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param LAMBDA0 The starting values of lambda. Must have the same dimensional matrix (n by n) with mHSpec.
#' @param llh_fun user provided log-likelihood function.
#' @param constraint constraint matrix. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param method method for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#'
#' @examples
#' # Generate sample path
#' MU1 <- 0.3; ALPHA1 <- 1.5; BETA1 <- 2
#' mHSpec1 <- new("mHSpec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1)
#' res1 <- mHSim(mHSpec1,  n=5000)
#'
#' # Perform maximum likelihood estimation with a starting point defined by mHSpec0.
#' mHSpec0 <- new("mHSpec", MU=0.2, ALPHA=1.2, BETA=1.8)
#' mle <- mHFit(mHSpec0, inter_arrival = res1$inter_arrival)
#' summary(mle)
#'
#' MU2 <- matrix(c(0.2), nrow = 2)
#' ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
#' BETA2 <- matrix(c(2.90, 2.90, 2.90, 2.90), nrow = 2, byrow=TRUE)
#' ETA2 <- matrix(c(0.19, 0.19, 0.19, 0.19), nrow = 2, byrow=TRUE)
#' JUMP2 <- function(n,...) rgeom(n, 0.65) + 1
#' mHSpec2 <- new("mHSpec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2, Jump =JUMP2)
#' res2 <- mHSim(mHSpec2)
#'
#' # Perform maximum likelihood estimation
#' summary(mHFit(mHSpec2, arrival = res2$arrival, N = res2$N))
#' summary(mHFit(mHSpec2, inter_arrival = res2$inter_arrival,
#'         mark = res2$mark2, jump_type = res2$jump_type))
#' summary(mHFit(arrival = res2$arrival, N = res2$N))
#'
#' @seealso \code{\link{mHSpec-class}}, \code{\link{mHSim,mHSpec-method}}
setMethod(
  f="mHFit",
  signature(object="mHSpec"),
  function(object, arrival = NULL, inter_arrival = NULL, N = NULL,
           jump_type = NULL, mark = NULL, LAMBDA0 = NULL, llh_fun = NULL,
          grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",  ...){

    # dimension of Hawkes process
    dimens <- length(object@MU)

    # argument check
    if(is.null(arrival) & is.null(inter_arrival)){
      stop("One of arrival and inter_arrival should be provided.")
    } else if(is.null(inter_arrival)){
      inter_arrival <- c(0, arrival[-1] - arrival[-length(arrival)])
    }

    # argument check
    if(dimens != 1 & is.null(N) & is.null(jump_type)){
      stop("One of N and jump_type should be provided.")
    } else if(dimens != 1 & is.null(jump_type)){

      if(!is.matrix(N)) N <- matrix(N, nrow = length(N))

      jump_type <- numeric(nrow(N))
      mark <- numeric(nrow(N))

      for (i in 2:nrow(N)){
        jump_type[i] <- which(N[i,] != N[i-1,])
        mark[i] <- N[i, jump_type[i]] - N[i-1, jump_type[i]]
      }

    }


    # When the mark sizes are not provided or max(mark) == 1, the jumps are all unit jumps.
    unit <- FALSE
    if(is.null(mark) | max(mark) == 1) {
      mark <- c(0, rep(1, length(inter_arrival)-1))
      unit <- TRUE
    }



    # parameter setting
    MU <- matrix(object@MU, nrow=dimens)
    ALPHA <- matrix(object@ALPHA, nrow=dimens)
    BETA <- matrix(object@BETA, nrow=dimens)
    ETA <- matrix(object@ETA, nrow=dimens)


    ref_mu <- name_unique_coef_mtrx(MU, "mu")
    unique_mus <- unique(as.vector(MU))
    names(unique_mus) <- unique(ref_mu)

    # ref_alpha looks like ["alpha11", "alpha12", "alpha12", "alpha11"] when alpha11==alpha22, alpha12==alpha21
    ref_alpha <- name_unique_coef_mtrx(ALPHA, "alpha")
    unique_alphas <- unique(as.vector(t(ALPHA)))
    names(unique_alphas) <-  unique(ref_alpha)


    ref_beta <- name_unique_coef_mtrx(BETA, "beta")
    unique_betas <- unique(as.vector(t(BETA)))
    names(unique_betas) <-  unique(ref_beta)

    # constant unit jump or not
    if (unit) starting_point <- c(unique_mus, unique_alphas, unique_betas)
    else {
      ref_eta <- name_unique_coef_mtrx(ETA, "eta")
      unique_etas <- unique(as.vector(t(ETA)))
      names(unique_etas) <-  unique(ref_eta)
      starting_point <- c(unique_mus, unique_alphas, unique_betas, unique_etas)
    }

    len_mu <- length(unique_mus)
    len_alpha <- length(unique_alphas)
    len_beta <- length(unique_betas)

    if (!unit) len_eta <- length(unique_etas)
    else len_eta <- 0


    # constraint matrix
    # mu, alpha, beta should be larger than zero
    if (unit) A <- diag(1, nrow = length(starting_point) - len_eta)
    else A <- cbind(diag(1, nrow = length(starting_point) - length(unique_etas)), rep(0, length(starting_point) - length(unique_etas)))



    # constraint : sum of alpha < beta
    A <- rbind(A, c(0, rep(-1, len_alpha), 1, rep(0, len_eta)))
    B <- rep(0, nrow(A))



    # loglikelihood function for maxLik
    if (is.null(llh_fun)) {
      llh_fun <- function(param){

        # redefine unique vectors from param
        unique_mus <- param[1:len_mu]
        unique_alphas <- param[(len_mu + 1):(len_mu + len_alpha)]
        unique_betas <- param[(len_mu + len_alpha + 1):(len_mu + len_alpha + len_beta)]
        if (!unit) unique_etas <- param[(len_mu + len_alpha + len_beta + 1):(len_mu + len_alpha + len_beta + len_eta)]

        # retreive MU, ALPHA, BETA, ETA matrix
        MU <- matrix( rep(0, dimens))
        k <- 1
        for  (i in 1:dimens){
          MU[i] <- unique_mus[ref_mu[k]]
          k <- k + 1
        }

        ALPHA <- matrix( rep(0, dimens^2), nrow=dimens)
        k <- 1
        for  (i in 1:dimens){
          for (j in 1:dimens) {
            ALPHA[i,j] <- unique_alphas[ref_alpha[k]]
            k <- k + 1
          }
        }

        BETA <- matrix( rep(0, dimens^2), nrow=dimens)
        k <- 1
        for  (i in 1:dimens){
          for (j in 1:dimens) {
            BETA[i,j] <- unique_betas[ref_beta[k]]
            k <- k + 1
          }
        }

        if (unit) ETA <- matrix(rep(0, dimens^2), nrow=dimens)
        else{
          ETA <- matrix( rep(0, dimens^2), nrow=dimens)
          k <- 1
          for  (i in 1:dimens){
            for (j in 1:dimens) {
              ETA[i,j] <- unique_etas[ref_eta[k]]
              k <- k + 1
            }
          }
        }

        mHSpec0 <- new("mHSpec", MU=MU, ALPHA=ALPHA, BETA=BETA, ETA=ETA, Jump=object@Jump)


        llh <- logLik(mHSpec0, inter_arrival = inter_arrival, jump_type = jump_type, mark = mark, LAMBDA0)
        return(llh)

      }
    }

    maxLik::maxLik(logLik=llh_fun,
                    start=starting_point, grad, hess, constraint, method = method)

  }
)


setMethod(
  f="mHFit",
  signature(object="missing"),
  function(object, arrival = NULL, inter_arrival = NULL, N = NULL,
           jump_type = NULL, mark = NULL, ...){


    # argument check
    if(is.null(arrival) & is.null(inter_arrival)){
      stop("One of arrival and inter_arrival should be provided.")
    } else if(is.null(inter_arrival)){
      inter_arrival <- c(0, arrival[-1] - arrival[-length(arrival)])
    }



    # argument check
    if(is.null(N) & is.null(jump_type)){

      #assuming one dimensional model
      dimens <- 1


    } else if(is.null(jump_type)){

      if(!is.matrix(N)) N <- matrix(N, nrow = length(N))

      jump_type <- numeric(nrow(N))
      mark <- numeric(nrow(N))

      for (i in 2:nrow(N)){
        jump_type[i] <- which(N[i,] != N[i-1,])
        mark[i] <- N[i, jump_type[i]] - N[i-1, jump_type[i]]
      }

      dimens <- ncol(N)

    } else if(is.null(N)){

      dimens <- max(jump_type)
    }

    if(dimens != 1 & dimens !=2 ){
      stop("One or two dimesinoal models is supported for default estimation.")
    }


    # set default mHSpec0
    if (dimens ==1 ){

      MU1 <- 0.2
      ALPHA1 <- 1.0
      BETA1 <- 2
      ETA1 <- 0

      mHSpec0 <- new("mHSpec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1, ETA=ETA1)

    } else if(dimens == 2){

      MU2 <- matrix(c(0.2), nrow = 2)
      ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
      BETA2 <- matrix(c(2.90, 2.90, 2.90, 2.90), nrow = 2, byrow=TRUE)
      ETA2 <- matrix(c(0, 0, 0, 0), nrow = 2, byrow=TRUE)

      mHSpec0 <- new("mHSpec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2)
    }

    if (is.null(mark)){
      mark <- rep(1, length(inter_arrival))
      mark[1] <- 0
    }

    mHFit(mHSpec0, inter_arrival = inter_arrival, jump_type = jump_type, mark = mark)

  }
)
