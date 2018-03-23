#' Compute the loglikelihood function
#'
#' This is a generic function.
#' The loglikelihood of the ground process of the Hawkes model.
#' (The estimation for jump distribution is not provided.)
#'
#' @param object \code{\link{mhspec-class}}. The parameter values in the object are used to compute the log-likelihood.
#' @param inter_arrival Inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param mark_type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param LAMBDA0 The starting values of lambda. Must have the same dimensional matrix (n by n) with \code{mhspec}.
#'
#' @examples
#' # construct a mhspec
#' MU1 <- 0.2; ALPHA1 <- 1.0; BETA1 <- 2; ETA1 <- 0.2
#' mark1 <- function(n,...) rgeom(n, 0.7) + 1
#' mhspec1 <- new("mhspec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1, ETA=ETA1, mark =mark1)
#' # simualte a path
#' res1 <- mhsim(mhspec1,  LAMBDA0 = MU1, n=1000)
#' inter_arrival <- res1$inter_arrival
#' mark <- res1$mark
#' # compute a loglikelihood function with parameter values in mhspec1
#' # LAMBDA0 = MU1 is a naive way of starting point choice.
#' logLik(mhspec1, LAMBDA0 = MU1, inter_arrival = inter_arrival, mark = mark)
#'
#' @seealso \code{\link{mhspec-class}}, \code{\link{mhfit,mhspec-method}}
setMethod(
  f="logLik",
  signature(object="mhspec"),
  function(object, inter_arrival, mark_type=NULL, mark=NULL, LAMBDA0=NULL){

    # When the mark sizes are not provided, the jumps are all unit jumps
    if(is.null(mark)) {
      mark <- rep(1, length(inter_arrival))
    }

    # dimension of Hawkes process
    dimens <- length(object@MU)

    # if dimens == 1 and mark_type is not provided, then all mark_type is 1.
    if(dimens==1 & is.null(mark_type)) {
      mark_type <- rep(1, length(inter_arrival))
    } else if (dimens != 1 & is.null(mark_type)) {
      stop("The argument mark_type should be provided.")
    }

    # parameter setting

    MU <- object@MU
    ALPHA <- object@ALPHA
    BETA <- object@BETA
    ETA <- object@ETA

    # default LAMBDA0
    if(is.null(LAMBDA0)) {
       LAMBDA0 <- get_lambda0(object)
    }

    # n is length(inter_arrival) - 1
    n <- length(inter_arrival) - 1

    #if (dimens==1) rowSums_LAMBDA0 <- LAMBDA0
    #else rowSums_LAMBDA0 <- rowSums(LAMBDA0)
    rowSums_LAMBDA0 <- rowSums(matrix(LAMBDA0, nrow=dimens))

    sum_log_lambda <- 0
    sum_integrated_lambda_component <- 0

    for (i in 1:n) {
      # current_LAMBDA <- matrix(lambda_component[i, ], nrow = dimens, byrow = TRUE)
      if (i == 1) current_LAMBDA <- LAMBDA0
      else current_LAMBDA <- new_LAMBDA  # LAMBDA determined in the previous loop

      # update lambda
      if (dimens == 1) {
        Impact <- ALPHA * (1 + (mark[i+1] - 1) * ETA )
      } else {
        Impact <- matrix(rep(0, dimens^2), nrow = dimens)
        Impact[ , mark_type[i+1]] <- ALPHA[ , mark_type[i+1]] * (1 + (mark[i+1] - 1) * ETA[ , mark_type[i+1]])
      }


      decayed <- exp(-BETA * inter_arrival[i+1])
      decayed_LAMBDA <- current_LAMBDA * decayed
      new_LAMBDA <- decayed_LAMBDA + Impact

      # sum of integrated_lambda_component
      sum_integrated_lambda_component <- sum_integrated_lambda_component + sum(current_LAMBDA / BETA * ( 1 - decayed ))

      # sum of log lambda when jump occurs
      if (dimens == 1) lambda_lc <- MU + decayed_LAMBDA
      else lambda_lc <- MU + rowSums(decayed_LAMBDA)
      sum_log_lambda <- sum_log_lambda + log(lambda_lc[mark_type[i+1]])

    }

    # log likelihood for ground process
    sum_log_lambda - sum(MU*sum(inter_arrival)) - sum_integrated_lambda_component

  }
)


setGeneric("mhfit", function(object, ...) standardGeneric("mhfit"))

#' Perform a maximum likelihood estimation
#'
#' This function uses \code{\link[maxLik]{maxLik}} for the optimizer.
#'
#'
#' @param object mhspec, or can be omitted.
#' @param arrival arrival times of events and hence monotonically increases. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param inter_arrival inter-arrival times of events. Includes inter-arrival for events that occur in all dimensions. Start with zero.
#' @param mark_type a vector of dimensions. Distinguished by numbers, 1, 2, 3, and so on. Start with zero.
#' @param mark a vector of mark (jump) sizes. Start with zero.
#' @param N a realization of n-dimensional Hawkes process.
#' @param LAMBDA0 the starting values of lambda. Must have the same dimensional matrix (n by n) with mhspec.
#' @param constraint constraint matrix. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param method method for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param grad gradient matrix for the likelihood function. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param hess Hessian matrix for the likelihood function. For more information, see \code{\link[maxLik]{maxLik}}.
#' @param ... other parameters for optimization. For more information, see \code{\link[maxLik]{maxLik}}.
#'
#' @examples
#' # Generate sample path
#' MU1 <- 0.3; ALPHA1 <- 1.5; BETA1 <- 2
#' mhspec1 <- new("mhspec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1)
#' res1 <- mhsim(mhspec1,  n=5000)
#'
#' # Perform maximum likelihood estimation with a starting point defined by mhspec0.
#' mhspec0 <- new("mhspec", MU=0.2, ALPHA=1.2, BETA=1.8)
#' mle <- mhfit(mhspec0, inter_arrival = res1$inter_arrival)
#' summary(mle)
#'
#' MU2 <- matrix(c(0.2), nrow = 2)
#' ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
#' BETA2 <- matrix(c(2.90, 2.90, 2.90, 2.90), nrow = 2, byrow=TRUE)
#' ETA2 <- matrix(c(0.19, 0.19, 0.19, 0.19), nrow = 2, byrow=TRUE)
#' mark2 <- function(n,...) rgeom(n, 0.65) + 1
#' mhspec2 <- new("mhspec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2, mark =mark2)
#' res2 <- mhsim(mhspec2)
#'
#' # Perform maximum likelihood estimation
#' summary(mhfit(mhspec2, arrival = res2$arrival, N = res2$N))
#' summary(mhfit(mhspec2, inter_arrival = res2$inter_arrival,
#'         mark = res2$mark, mark_type = res2$mark_type))
#' summary(mhfit(arrival = res2$arrival, N = res2$N))
#'
#' @seealso \code{\link{mhspec-class}}, \code{\link{mhsim,mhspec-method}}
setMethod(
  f="mhfit",
  signature(object="mhspec"),
  function(object, arrival = NULL, inter_arrival = NULL, N = NULL,
           mark_type = NULL, mark = NULL, LAMBDA0 = NULL,
           grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",  ...){

    # dimension of Hawkes process
    dimens <- length(object@MU)


    # argument check : one of arrival or inter_arrival is needed
    if(is.null(arrival) & is.null(inter_arrival)){
      stop("One of arrival or inter_arrival should be provided.")
    } else if(is.null(inter_arrival)){
      inter_arrival <- c(0, arrival[-1] - arrival[-length(arrival)])
    }

    # argument check for two or higher dimension
    if(dimens != 1){

      if(!is.null(N)){
        # when N is provided as matrix, it's ok to proceed.
        if(!is.matrix(N)) stop("N should be a matrix.")

        # extract mark_type from N, when mark_type is null
        if(is.null(mark_type)){

          mark_type <- numeric(nrow(N))

          for (i in 2:nrow(N)){
            mark_type[i] <- which(N[i,] != N[i-1,])
          }
        }


        # extrat mark from N, when mark is null

        if(is.null(mark)){
          mark <- numeric(nrow(N))

          for (i in 2:nrow(N)){
            mark[i] <- N[i, mark_type[i]] - N[i-1, mark_type[i]]
          }
        }



      } else {
        # when N is not provided as matrix, at least we need mark_type

        if(is.null(mark_type)){

          stop("One of N or mark_type should be provided.")

        } else {

          if(is.null(mark)){

            # if mark is not provided, default values with 1 are used.
            warning("Mark is not provided. Default values are used.")

            mark <- c(0, rep(1, length(mark_type-1)))

          }

        }

      }

    }


    if(is.null(LAMBDA0)){
      warning("The initial values for intensity processes are not provided. Internally determined initial values are used.\n")
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

    llh_function <- function(param){

      # redefine unique vectors from param
      unique_mus <- param[1:len_mu]
      unique_alphas <- param[(len_mu + 1):(len_mu + len_alpha)]
      unique_betas <- param[(len_mu + len_alpha + 1):(len_mu + len_alpha + len_beta)]
      if (!unit) unique_etas <- param[(len_mu + len_alpha + len_beta + 1):(len_mu + len_alpha + len_beta + len_eta)]

      # retreive MU, ALPHA, BETA, ETA matrix
      MU <- matrix( rep(0, dimens))
      i <- 1
      for  (m in 1:dimens){
        MU[m] <- unique_mus[ref_mu[i]]
        i <- i + 1
      }

      ALPHA <- matrix( rep(0, dimens^2), nrow=dimens)
      i <- 1
      for  (m in 1:dimens){
        for (n in 1:dimens) {
          ALPHA[m,n] <- unique_alphas[ref_alpha[i]]
          i <- i + 1
        }
      }

      BETA <- matrix( rep(0, dimens^2), nrow=dimens)
      i <- 1
      for  (m in 1:dimens){
        for (n in 1:dimens) {
          BETA[m,n] <- unique_betas[ref_beta[i]]
          i <- i + 1
        }
      }

      if (unit) ETA <- matrix(rep(0, dimens^2), nrow=dimens)
      else{
        ETA <- matrix( rep(0, dimens^2), nrow=dimens)
        i <- 1
        for  (m in 1:dimens){
          for (n in 1:dimens) {
            ETA[m,n] <- unique_etas[ref_eta[i]]
            i <- i + 1
          }
        }
      }

      mhspec0 <- methods::new("mhspec", MU=MU, ALPHA=ALPHA, BETA=BETA, ETA=ETA, mark=object@mark)

      llh <- logLik(mhspec0, inter_arrival = inter_arrival, mark_type = mark_type, mark = mark, LAMBDA0)
      return(llh)

    }

    maxLik::maxLik(logLik=llh_function,
                    start=starting_point, grad, hess, constraint, method = method)

  }
)


setMethod(
  f="mhfit",
  signature(object="missing"),
  function(object, arrival = NULL, inter_arrival = NULL, N = NULL,
           mark_type = NULL, mark = NULL,
           grad = NULL, hess = NULL, constraint = NULL, method = "BFGS",  ...){


    # argument check
    if(is.null(arrival) & is.null(inter_arrival)){
      stop("One of arrival and inter_arrival should be provided.")
    } else if(is.null(inter_arrival)){
      inter_arrival <- c(0, arrival[-1] - arrival[-length(arrival)])
    }


    # argument check
    if(is.null(N) & is.null(mark_type)){

      #assuming one dimensional model
      dimens <- 1


    } else if(is.null(mark_type)){

      if(!is.matrix(N)) N <- matrix(N, nrow = length(N))

      mark_type <- numeric(nrow(N))
      mark <- numeric(nrow(N))

      for (i in 2:nrow(N)){
        mark_type[i] <- which(N[i,] != N[i-1,])
        mark[i] <- N[i, mark_type[i]] - N[i-1, mark_type[i]]
      }

      dimens <- ncol(N)

    } else if(is.null(N)){

      dimens <- max(mark_type)
    }

    if(dimens != 1 & dimens !=2 ){
      stop("One or two dimesinoal models is supported for default estimation.")
    }


    # set default mhspec0
    if (dimens ==1 ){

      MU1 <- 0.2
      ALPHA1 <- 1.0
      BETA1 <- 2
      ETA1 <- 0

      mhspec0 <- methods::new("mhspec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1, ETA=ETA1)

    } else if(dimens == 2){

      MU2 <- matrix(c(0.2), nrow = 2)
      ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
      BETA2 <- matrix(c(2.90, 2.90, 2.90, 2.90), nrow = 2, byrow=TRUE)
      ETA2 <- matrix(c(0, 0, 0, 0), nrow = 2, byrow=TRUE)

      mhspec0 <- methods::new("mhspec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2)
    }

    if (is.null(mark)){
      mark <- rep(1, length(inter_arrival))
      mark[1] <- 0
    }

    mhfit(mhspec0, inter_arrival = inter_arrival, mark_type = mark_type, mark = mark,
          grad=grad, hess=hess, constraint=constraint, method = method)


  }
)
