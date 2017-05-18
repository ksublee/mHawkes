

setClassUnion("matrixORvector", c("matrix", "vector"))

#' check the validity of parameter matrix dimesions
check_dimension <- function(object) {

  errors <- character()
  print(length(object@MU))
  print(length(object@ALPHA))
  print(length(object@BETA))
  print(length(object@ETA))

  if ( length(object@MU)!=1 | length(object@ALPHA)!=1 | length(object@BETA)!=1 | length(object@ETA)!=1 ) {

    # The dimensions of parameter matrices should be consistent.
    if ( !(nrow(object@MU) == nrow(object@ALPHA) & nrow(object@ALPHA) == nrow(object@BETA) & nrow(object@BETA) == nrow(object@ETA)) ) {
      msg <- "The number of rows of parameter matrix should be equal"
      errors <- c(errors, msg)
    }

    # The dimension check for each parameter matrix.
    if ( ncol(object@MU) > 1 ) {
      msg <- "The number of columns of MU matrix should be one."
      errors <- c(errors, msg)
    }

    if ( ncol(object@ALPHA) != nrow(object@ALPHA) ) {
      msg <- "ALPHA matrix should be n by n matrix."
      errors <- c(errors, msg)
    }

    if ( ncol(object@BETA) != nrow(object@BETA) ) {
      msg <- "BETA matrix should be n by n matrix."
      errors <- c(errors, msg)
    }

    if ( ncol(object@ETA) != nrow(object@ETA) ) {
      msg <- "ETA matrix should be n by n matrix."
      errors <- c(errors, msg)
    }

    if (length(errors) == 0) TRUE else errors
  }
}


#' An S4 calss to represent the specification of a marked Hawkes model
#'
#' @slot MU
setClass(
  "mHSpec",
  # slots = list(MU = "matrix", ALPHA = "matrix", BETA = "matrix", ETA = "matrix", Jump = "distr::Distribution"),
  representation(
    MU = "matrixORvector",
    ALPHA = "matrixORvector",
    BETA = "matrixORvector",
    ETA = "matrixORvector",
    Jump = "Distribution"
    )

  #prototype = list(
  #  Jump = distr::Dirac(location = 1)
  #)

  #validity = check_dimension
)


setMethod(
  "initialize",
  "mHSpec",
  function(.Object, MU, ALPHA, BETA, ETA=NULL, Jump=NULL){

    .Object@MU <- MU
    .Object@ALPHA <- ALPHA
    .Object@BETA <- BETA

    if (is.null(ETA)) {
      if (length(.Object@MU) == 1) .Object@ETA <- 0
      else .Object@ETA <- matrix(rep(0, length(.Object@BETA), nrow = length(.Object@BETA)/2))
    } else {
      .Object@ETA <- ETA
    }

    if (is.null(Jump)) .Object@Jump <- distr::Dirac(location = 1)
    else .Object@Jump <- Jump

    .Object
  }
)

setGeneric("mHSim", function(object, ...) standardGeneric("mHSim"))

setMethod(
  f="mHSim",
  signature(object = "mHSpec"),
  definition = function(object, LAMBDA0, n=1000){

    # dimension of Hawkes process
    dimens <- length(object@MU)

    # parameter setting
    MU <- matrix(object@MU, nrow=dimens)
    ALPHA <- matrix(object@ALPHA, nrow=dimens)
    BETA <- matrix(object@BETA, nrow=dimens)
    ETA <- matrix(object@ETA, nrow=dimens)

    # We need an appropriate algorithm for buffer size n.

    # Preallocation for lambdas and Ns and set initial values for lambdas
    lambda_component <- matrix(sapply(LAMBDA0, c, numeric(length = n - 1)), ncol = dimens^2)
    rowSums_LAMBDA0 <- rowSums(matrix(LAMBDA0, nrow=dimens))

    lambda_process   <- matrix(sapply(MU + rowSums_LAMBDA0, c, numeric(length = n - 1)), ncol = dimens)

    Ng <- matrix(numeric(length = dimens * n), ncol = dimens)
    N  <- matrix(numeric(length = dimens * n), ncol = dimens)

    jump_type <- numeric(length = n-1)
    inter_arrival <- numeric(length = n-1)
    mark <- numeric(length = n-1)

    names(lambda_process) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    names(lambda_component) <- paste0("lambda", indxM, t(indxM))
    names(N)  <- paste0("N", 1:dimens)
    names(Ng) <- paste0("Ng", 1:dimens)

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
      inter_arrival[k] <- min(candidate_arrival)
      minIndex <- which(candidate_arrival == inter_arrival[k], arr.ind = TRUE) #row and col

      jumpType <- minIndex[1]  # row
      jump_type[k] <- jumpType

      mark[k] <- distr::r(object@Jump)(1) # generate one random number

      Ng[k+1, ] <- Ng[k, ]
      Ng[k+1, jumpType] <- Ng[k, jumpType] + 1
      N[k+1, ] <- N[k,]
      N[k+1, jumpType] <- N[k, jumpType] + mark[k]

      # update lambda
      Impact <- matrix(rep(0, dimens^2), nrow = dimens)
      Impact[ , jumpType] <- ALPHA[ , jumpType] + (1 + (mark[k] - 1)) * ETA[ , jumpType]

      new_LAMBDA <- current_LAMBDA * exp(-BETA * inter_arrival[k]) + Impact

      # new_LAMBDA = [[lambda11, lambda12, ...], [lambda21, lambda22, ...], ...]
      # lambda_component = {"lambda11", "lambda12", ..., "lambda21", "lambda22", ...}
      lambda_component[k+1, ] <- t(new_LAMBDA)
      lambda_process[k+1, ] <- MU + rowSums(new_LAMBDA)
    }

    # convert to data frame
    lambda_component <- data.frame(lambda_component)
    lambda_process   <- data.frame(lambda_process)

    Ng <- data.frame(Ng)
    N  <- data.frame(N)

    # naming the dataframes
    names(lambda_process) <- paste0("lambda", 1:dimens)
    indxM <- matrix(rep(1:dimens, dimens), byrow = TRUE, nrow = dimens)
    names(lambda_component) <- paste0("lambda", indxM, t(indxM))
    names(N)  <- paste0("N", 1:dimens)
    names(Ng) <- paste0("Ng", 1:dimens)

    realization <- list(inter_arrival, cumsum(inter_arrival), jump_type, mark, N, Ng, lambda_process, lambda_component)
    names(realization) <- c("inter_arrival", "arrival", "jump_type", "mark", "N", "Ng", "lambda", "lambda_component")

    return(realization)
  }
)

setGeneric("loglikelihood", function(object, ...) standardGeneric("loglikelihood"))

#' Estimate a marked Hawkes process
#' @param inter_arrival
#' @param jump_type
#' @param mark
setMethod(
  f="loglikelihood",
  signature(object="mHSpec"),
  function(object, inter_arrival, jump_type=NULL, mark=NULL, LAMBDA0){

    # When the mark sizes are not provided, the jumps are all unit jumps
    if(is.null(mark)) {
      mark <- rep(1, length(inter_arrival))
    }

    # dimension of Hawkes process
    dimens <- length(object@MU)

    # if dimes == 1 and jump_type is not provided, then all jump_type is 1.
    if(dimens==1 & is.null(jump_type)) {
      jump_type <- rep(1, length(inter_arrival))
    } else if (dimens!=1 & is.null(jump_type)) {
      stop("The argument mark should be provided.")
    }

    # parameter setting
    MU <- matrix(object@MU, nrow=dimens)
    ALPHA <- matrix(object@ALPHA, nrow=dimens)
    BETA <- matrix(object@BETA, nrow=dimens)
    ETA <- matrix(object@ETA, nrow=dimens)

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

      current_LAMBDA <- matrix(as.numeric(lambda_component[k, ]), nrow = dimens, byrow = TRUE)

      # update lambda
      Impact <- matrix(rep(0, dimens^2), nrow = dimens)
      Impact[ , jump_type[k]] <- ALPHA[ , jump_type[k]] + (1 + (mark[k] - 1)) * ETA[ , jump_type[k]]

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
  function(object, inter_arrival, jump_type=NULL, mark=NULL, LAMBDA0, model="symmetric"){

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

    # unit jump or not
    if (unit) starting_point <- c(mu = mu, alphas, beta = beta)
    else {
      eta <- object@ETA[1]
      starting_point <- c(mu = mu, alphas, beta = beta, eta = eta)
    }

    G <- object@Jump


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
      else ETA <- matrix(rep(param[[1 + length(alphas) + 1 + 1]], 4), nrow=dimens)

      mHSpec1 <- new("mHSpec", MU=MU, ALPHA=ALPHA, BETA=BETA, ETA=ETA, Jump=G)
      #LAMBDA0 <- matrix(c(0.03,0.03,0.03,0.03), nrow=2, byrow=TRUE)

      llh <- loglikelihood(mHSpec1, inter_arrival = inter_arrival, jump_type = jump_type, mark, LAMBDA0)
      return(llh)

    }

    mle<-maxLik::maxLik(logLik=llh_maxLik,
                        start=starting_point, constraint=list(ineqA=A, ineqB=B),
                        method="BFGS")


  }
)
