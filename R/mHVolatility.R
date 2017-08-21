
setGeneric("mHVar", function(object, ...) standardGeneric("mHVar"))

#' Compute the Hawkes variance
#'
#' This method is useful in quantitative finance or financial econometics.
#' This method only works for a two-dimensional symmetric model.
#' Assume that a two-dimensional Hawkes process describes the tick dynamics of financial price process.
#' One of the two Hawkes processes is responsible for the upward movement and the other for the downward movement.
#' Therefore, the difference between two process, N1 - N2, describes the tick price process.
#' This method computes Var(N1(t) - N2(t)).
#'
#'
#' @param object
#' @param tick_ratio
#' @time_length
#' @mean_jump
#' @mean_jump_square
#'
setMethod(
  f="mHVar",
  signature(object = "mHSpec"),
  definition = function(object, time_length, mean_jump = NULL, mean_jump_square = NULL, sample_size = 10000){

    cat("This method only works for a two-dimensional symmetric model with i.i.d. jump distribution.\n")

    mu <- object@MU[1]
    alpha_s <- object@ALPHA[1,1]
    alpha_c <- object@ALPHA[1,2]
    beta <- object@BETA[1,1]
    eta <- object@ETA[1,1]

    if(is.null(mean_jump)){
      K <- mean(object@Jump(sample_size))
    } else {
      K <- mean_jump
    }

    if(is.null(mean_jump_square)){
      K2 <- mean(object@Jump(sample_size)^2)
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

    variance <- 2*K*E_lambda*(C1 + C2 + C3)*time_length

  }
)
