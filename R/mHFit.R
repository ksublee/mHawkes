#' Compute the loglikelihood function of the ground process
#'
#' @param inter_arrival
#' @param jump_type
#' @param mark
#' @param LAMBDA0
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
      Impact <- matrix(rep(0, dimens^2), nrow = dimens)
      Impact[ , jump_type[k+1]] <- ALPHA[ , jump_type[k+1]] * (1 + (mark[k+1] - 1) * ETA[ , jump_type[k+1]])

      decayed <- exp(-BETA * inter_arrival[k+1])
      decayed_LAMBDA <- current_LAMBDA * decayed
      new_LAMBDA <- decayed_LAMBDA + Impact

      # sum of integrated_lambda_component
      sum_integrated_lambda_component <- sum_integrated_lambda_component + sum(current_LAMBDA / BETA * ( 1 - decayed ))

      # sum of log lambda when jump occurs
      lambda_lc <- MU + rowSums(decayed_LAMBDA)
      sum_log_lambda <- sum_log_lambda + log(lambda_lc[jump_type[k+1]])

    }

    # log likelihood for ground process
    sum_log_lambda - sum(MU*sum(inter_arrival)) - sum_integrated_lambda_component

  }
)

setGeneric("mHFit", function(object, ...) standardGeneric("mHFit"))

setMethod(
  f="mHFit",
  signature(object="mHSpec"),
  function(object, inter_arrival, jump_type=NULL, mark=NULL, LAMBDA0=NULL, model="symmetric", constraint=FALSE, method="BFGS"){

    cat("starting optimization...\n")
    # When the mark sizes are not provided, the jumps are all unit jumps.
    unit <- FALSE
    if(is.null(mark)) {
      mark <- c(0, rep(1, length(inter_arrival)-1))
      unit <- TRUE
    }

    # dimension of Hawkes process
    dimens <- length(object@MU)

    # parameter setting
    MU <- matrix(object@MU, nrow=dimens)
    ALPHA <- matrix(object@ALPHA, nrow=dimens)
    BETA <- matrix(object@BETA, nrow=dimens)
    ETA <- matrix(object@ETA, nrow=dimens)

    # MU <- object@MU
    # ALPHA <- object@ALPHA
    # BETA <- object@BETA
    # ETA <- object@ETA

    # simple reference function to find unique value
    find_ref <- function(M, notation){
      reference <- character(length(M))

      if (ncol(M) == 1){
        k <- 1
        for (i in 1:nrow(M)){
          if (reference[k] == "")
            reference[which(M == M[i])] <- paste0(notation, toString(i))
          k <- k + 1
        }
      } else {
        k <- 1
        for  (i in 1:nrow(M)){
          for (j in 1:ncol(M)) {
            if (reference[k] == "")
              reference[which( t(M) == M[i,j])] <- paste0(notation, toString(i), toString(j))
            k <- k + 1
          }
        }
      }
      reference
    }

    ref_mu <- find_ref(MU, "mu")
    unique_mus <- unique(as.vector(MU))
    names(unique_mus) <- unique(ref_mu)

    ref_alpha <- find_ref(ALPHA, "alpha")
    unique_alphas <- unique(as.vector(t(ALPHA)))
    names(unique_alphas) <-  unique(ref_alpha)


    ref_beta <- find_ref(BETA, "beta")
    unique_betas <- unique(as.vector(t(BETA)))
    names(unique_betas) <-  unique(ref_beta)

    # constant unit jump or not
    if (unit) starting_point <- c(unique_mus, unique_alphas, unique_betas)
    else {
      ref_eta <- find_ref(ETA, "eta")
      unique_etas <- unique(as.vector(t(ETA)))
      names(unique_etas) <-  unique(ref_eta)
      starting_point <- c(unique_mus, unique_alphas, unique_betas, unique_etas)
    }

    len_mu <- length(unique_mus)
    len_alpha <- length(unique_alphas)
    len_beta <- length(unique_betas)
    len_eta <- length(unique_etas)

    # constraint matrix
    # mu, alpha, beta should be larger than zero
    if (unit) A <- diag(1, nrow = length(starting_point) - length(unique_etas))
    else A <- cbind(diag(1, nrow = length(starting_point) - length(unique_etas)), rep(0, length(starting_point) - length(unique_etas)))
    # constraint : sum of alpha < beta
    A <- rbind(A, c(0, rep(-1, length(unique_alphas)), 1, rep(0, length(unique_etas))))

    B <- rep(0, nrow(A))

    # loglikelihood function for maxLik
    llh_maxLik <- function(param){

      cat("#")
      # redefine unique vectors from param
      unique_mus <- param[1:len_mu]
      unique_alphas <- param[(len_mu + 1):(len_mu + len_alpha)]
      unique_betas <- param[(len_mu + len_alpha + 1):(len_mu + len_alpha + len_beta)]
      unique_etas <- param[(len_mu + len_alpha + len_beta + 1):(len_mu + len_alpha + len_beta + len_eta)]

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

    if (constraint)
      mle<-maxLik::maxLik(logLik=llh_maxLik,
                          start=starting_point, constraint=list(ineqA=A, ineqB=B),
                          method=method)
    else
      mle<-maxLik::maxLik(logLik=llh_maxLik,
                          start=starting_point,
                          method=method)

    cat("\n")
    cat("ending procedure.")
    mle

  }
)


setGeneric("mHFit_old", function(object, ...) standardGeneric("mHFit_old"))

setMethod(
  f="mHFit_old",
  signature(object="mHSpec"),
  function(object, inter_arrival, jump_type=NULL, mark=NULL, LAMBDA0=NULL, model="symmetric", constraint=TRUE, method="BFGS"){

    # When the mark sizes are not provided, the jumps are all unit jumps.
    unit <- FALSE
    if(is.null(mark)) {
      mark <- c(0, rep(1, length(inter_arrival)-1))
      unit <- TRUE
    }

    # dimension of Hawkes process
    dimens <- length(object@MU)

    # parameter setting
    MU <- matrix(object@MU, nrow=dimens)
    ALPHA <- matrix(object@ALPHA, nrow=dimens)
    BETA <- matrix(object@BETA, nrow=dimens)
    ETA <- matrix(object@ETA, nrow=dimens)


    # symmetric for alpha and the same beta
    mu <- MU[1]

    # alphas = [alpha11, alpha12, alpha13, ..., alpha1n, alpha23, ..., alpha2n, alpha34, ..., alpha3n , ...]
    alphas <- numeric(length = dimens*(dimens-1)/2 + 1)
    alphas[1] <- ALPHA[1,1]
    alpha_names <- "alpha11"
    k <- 2
    for (i in 1:dimens) {
      j <- i + 1
      while (j <=  dimens) {
        alphas[k] <- ALPHA[i, j]
        alpha_names <- c(alpha_names, paste0("alpha",i,j))
        j <- j + 1
        k <- k + 1
      }
    }

    names(alphas) <- alpha_names

    beta <- BETA[1]
    eta <- numeric()

    # constant unit jump or not
    if (unit) starting_point <- c(mu = mu, alphas, beta = beta)
    else {
      eta <- object@ETA[1]
      starting_point <- c(mu = mu, alphas, beta = beta, eta = eta)
    }

    # constraint matrix
    # mu, alpha, beta should be larger than zero
    if (unit) A <- diag(1, nrow = length(starting_point) - length(eta))
    else A <- cbind(diag(1, nrow = length(starting_point) - length(eta)), rep(0, length(starting_point) - length(eta)))
    # constraint : sum of alpha < beta
    A <- rbind(A, c(0, rep(-1, length(alphas)), 1, rep(0, length(eta))))

    B <- rep(0, nrow(A))

    # loglikelihood function for maxLik
    llh_maxLik <- function(param){

      MU <- matrix(rep(param[[1]], dimens), nrow=dimens)
      # diagonal part of ALPHA
      ALPHA <- diag(param[[2]], nrow=dimens)
      # upper tri part of ALPHA
      k <- 3
      for (i in 1:dimens) {
        j <- i + 1
        while (j <=  dimens) {
          ALPHA[i, j] <- param[k]
          j <- j + 1
          k <- k + 1
        }
      }
      # symmetric ALPHA
      ALPHA[lower.tri(ALPHA)] = t(ALPHA)[lower.tri(ALPHA)]

      # ALPHA <- matrix(c(param[[2]], param[[3]], param[[3]], param[[2]]), nrow=2, byrow=TRUE)
      BETA <- matrix(rep(param[[1 + length(alphas) + 1]], dimens^2), nrow=dimens, byrow=TRUE)

      if (unit) ETA <- matrix(rep(0, dimens^2), nrow=dimens)
      else ETA <- matrix(rep(param[[1 + length(alphas) + 1 + 1]], dimens^2), nrow=dimens)

      # starting point
      mHSpec0 <- new("mHSpec", MU=MU, ALPHA=ALPHA, BETA=BETA, ETA=ETA, Jump=object@Jump)


      llh <- logLik(mHSpec0, inter_arrival = inter_arrival, jump_type = jump_type, mark = mark, LAMBDA0)
      return(llh)

    }

    if (constraint)
      mle<-maxLik::maxLik(logLik=llh_maxLik,
                          start=starting_point, constraint=list(ineqA=A, ineqB=B),
                          method=method)
    else
      mle<-maxLik::maxLik(logLik=llh_maxLik,
                          start=starting_point,
                          method=method)

  }
)

