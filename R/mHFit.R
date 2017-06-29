setGeneric("loglikelihood", function(object, ...) standardGeneric("loglikelihood"))

#' Compute the loglikelihood function of the ground process
#'
#' @param inter_arrival
#' @param jump_type
#' @param mark
#' @param LAMBDA0
setMethod(
  f="loglikelihood",
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
    } else if (dimens!=1 & is.null(jump_type)) {
      stop("The argument jump_type should be provided.")
    }

    # parameter setting
    MU <- matrix(object@MU, nrow=dimens)
    ALPHA <- matrix(object@ALPHA, nrow=dimens)
    BETA <- matrix(object@BETA, nrow=dimens)
    ETA <- matrix(object@ETA, nrow=dimens)

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


    # n is length(inter_arrival) then the length of lambda is n+1
    n <- length(inter_arrival)

    #if (dimens==1) rowSums_LAMBDA0 <- LAMBDA0
    #else rowSums_LAMBDA0 <- rowSums(LAMBDA0)
    rowSums_LAMBDA0 <- rowSums(matrix(LAMBDA0, nrow=dimens))

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(LAMBDA0, c, numeric(length = n)), ncol=dimens^2)
    integrated_lambda_component <- matrix(numeric(length=dimens^2 * (n+1)), ncol = dimens^2)

    lambda_lc_process   <- matrix(sapply(MU + rowSums_LAMBDA0, c, numeric(length = n)), ncol = dimens)
    lambda_lc_component <- matrix(sapply(LAMBDA0, c, numeric(length = n)), ncol = dimens^2)


    sum_log_lambda <- 0

    for (k in 1:n) {

      current_LAMBDA <- matrix(lambda_component[k, ], nrow = dimens, byrow = TRUE)

      # update lambda
      Impact <- matrix(rep(0, dimens^2), nrow = dimens)
      Impact[ , jump_type[k]] <- ALPHA[ , jump_type[k]] * (1 + (mark[k] - 1) * ETA[ , jump_type[k]])

      decayed <- exp(-BETA * inter_arrival[k])
      decayed_LAMBDA <- current_LAMBDA * decayed
      new_LAMBDA <- decayed_LAMBDA + Impact


      # new_LAMBDA = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      # lambda_component = ["lambda11", "lambda12", ..., "lambda21", "lambda22", ...]
      # lambda_process = [lambda1, lambda2, ...]
      # integrated_lambda_component
      lambda_component[k+1, ] <- t(new_LAMBDA)
      lambda_lc_component[k+1, ] <- t(decayed_LAMBDA)
      integrated_lambda_component[k+1, ] <- current_LAMBDA/BETA*( 1- decayed)

      #lambda_process[k+1, ] <- MU + rowSums(new_LAMBDA)
      lambda_lc_process[k+1, ] <- MU + rowSums(decayed_LAMBDA)

      sum_log_lambda <- sum_log_lambda + log(lambda_lc_process[k+1, jump_type[k]])
    }

    # log likelihood for ground process
    sum_log_lambda - sum(MU*sum(inter_arrival)) - sum(integrated_lambda_component)

  }
)


setGeneric("mHFit", function(object, ...) standardGeneric("mHFit"))

setMethod(
  f="mHFit",
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

      mHSpec0 <- new("mHSpec", MU=MU, ALPHA=ALPHA, BETA=BETA, ETA=ETA, Jump=object@Jump)


      llh <- loglikelihood(mHSpec0, inter_arrival = inter_arrival, jump_type = jump_type, mark = mark, LAMBDA0)
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
