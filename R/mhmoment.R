setGeneric("mean_lambda", function(object, ...) standardGeneric("mean_lambda"))

#' Compute the long-run mean of lambda
#'
#' This method only works for a one-dimensional or two-dimensional symmetric model.
#'
#' @param object mhspec, a one-dimensional or two-dimensional symmetricl model
#' @param mean_jump the mean of jump distribution. If not specified, simulated mean value will be used.
#' @param sample_size the number of simulation to compute the mean of jump and squared jump
#' @param seed_value seed for random number generation.
setMethod(
  f = "mean_lambda",
  signature(object = "mhspec"),
  definition = function(object, mean_jump = NULL, sample_size = 10000, seed_value = 903){

    dimens <- length(object@MU)

    if (dimens == 1){

      mu <- object@MU
      alpha <- object@ALPHA
      beta <- object@BETA
      eta <- object@ETA

      set.seed(seed_value)

      if(is.null(mean_jump)){
        K <- mean(object@mark(sample_size))
      } else {
        K <- mean_jump
      }

      mu*beta/(beta - alpha*(1 + (K - 1)*eta))

    } else if (dimens == 2){

      mu <- object@MU[1]
      alpha_s <- object@ALPHA[1,1]
      alpha_c <- object@ALPHA[1,2]
      beta <- object@BETA[1,1]
      eta <- object@ETA[1,1]

      if(is.null(mean_jump)){
        K <- mean(object@mark(sample_size))
      } else {
        K <- mean_jump
      }


      mu*beta/(beta - (alpha_s + alpha_c)*(1 + (K - 1)*eta))

    }

  }
)



setGeneric("var_diff", function(object, ...) standardGeneric("var_diff"))

#' Compute the variance of Hawkes difference process
#'
#' This method computes Var(N1(t) - N2(t)).
#' This method only works for a two-dimensional symmetric model.
#' Assume that a two-dimensional Hawkes process describes the tick dynamics of financial price process.
#' One of the two Hawkes processes is responsible for the upward movement and the other for the downward movement.
#' Therefore, the difference between two process, N1 - N2, describes the tick price process.
#' This method is useful in quantitative finance and financial econometics.
#'
#'
#' @param object mhspec a two dimensional symmetricl model
#' @param time_length time horizon
#' @param mean_jump the mean of jump distribution. If not specified, simulated mean value will be used.
#' @param mean_jump_square the mean of the square of jump distribution. If not specified, simulated value will be used.
#' @param sample_size the number of simulation to compute the mean of jump and squared jump
#' @param seed_value seed for random number generation.
#'
#' @examples
#' # two dimensional symmetric Hawkes model
#' MU2 <- matrix(c(0.2), nrow = 2)
#' ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
#' BETA2 <- matrix(c(2.90, 2.90, 2.90, 2.90), nrow = 2, byrow=TRUE)
#' ETA2 <- matrix(c(0.19, 0.19, 0.19, 0.19), nrow = 2, byrow=TRUE)
#' mark2 <- function(n,...) rgeom(n, 0.7) + 1
#' mhspec2 <- new("mhspec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2, mark =mark2)
#' var_diff(mhspec2, 1)
setMethod(
  f = "var_diff",
  signature(object = "mhspec"),
  definition = function(object, time_length, mean_jump = NULL, mean_jump_square = NULL, sample_size = 10000, seed_value = 903){

    dimens <- length(object@MU)
    if(dimens != 2) stop("This method only works for a two-dimensional symmetric model with i.i.d. jump distribution.\n")

    if( object@MU[1] != object@MU[2] | object@ALPHA[1,1] != object@ALPHA[2,2] | object@ALPHA[1,2] != object@ALPHA[2,1] |
        max(object@BETA) != min(object@BETA) | max(object@ETA) != min(object@ETA) ) {

      warning("This method only works for a two-dimensional symmetric model with i.i.d. jump distribution.\n")

    }

    set.seed(seed_value)
    mu <- object@MU[1]
    alpha_s <- object@ALPHA[1,1]
    alpha_c <- object@ALPHA[1,2]
    beta <- object@BETA[1,1]
    eta <- object@ETA[1,1]

    if(is.null(mean_jump)){
      K <- mean(object@mark(sample_size))
    } else {
      K <- mean_jump
    }

    if(is.null(mean_jump_square)){
      K2 <- mean(object@mark(sample_size)^2)
    } else {
      K2 <- mean_jump_square
    }


    K_bb <- 1 + 2*(K - 1)*eta + (K2 - 2*K +1)*eta^2
    K_b <- K + (K2 - K)*eta

    alpha_s_ <- alpha_s*(1 + (K - 1)*eta)
    alpha_c_ <- alpha_s*(1 + (K - 1)*eta)

    E_lambda <- mu*beta / (beta - (alpha_s + alpha_c)*(1 + (K - 1)*eta))

    C1 <- K*K_bb*(alpha_s - alpha_c)^2/(beta - alpha_s_ + alpha_c_)^2
    C2 <- 2*(alpha_s - alpha_c)*K_b/(beta - alpha_s_ + alpha_c_)
    C3 <- K^2/K

    2*K*E_lambda*(C1 + C2 + C3)*time_length

  }
)
