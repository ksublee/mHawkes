#'
setClassUnion("matrixORnumeric", c("matrix", "numeric"))

validmHSpec <- function(object) {


  if(!is.matrix(object@MU) | !is.matrix(object@ALPHA) | !is.matrix(object@BETA) | (!is.null(object@ETA) & !is.matrix(object@ETA))){

    # If one of the parameter is not matrix, then all parameters should not be matrix.

    if(!(!is.matrix(object@MU) & !is.matrix(object@ALPHA) & !is.matrix(object@BETA) & (!is.null(object@ETA) & !is.matrix(object@ETA)))){
      show(object)
      return("If one of the parameter is not matrix, then all parameters should not be matrix.")
    }

    # If one of the parameter is atomic, then all parameters should be atomic.

    len_mu <- length(object@MU)
    len_alpha <- length(object@ALPHA)
    len_beta <- length(object@BETA)
    if (!is.null(object@ETA))
      len_eta <- length(object@ETA)
    else
      len_eta <- len_mu
    if (max(c(len_mu, len_alpha, len_beta, len_eta)) != 1 ){
      show(object)
      return("If one of the parameter is atomic, then all parameters should be atomic.")
    }


  } else{

    # If all parameters are matrix, then the nrow should be equal.

    dim_mu <- nrow(object@MU)
    dim_alpha <- nrow(object@ALPHA)
    dim_beta <- nrow(object@BETA)
    if (!is.null(object@ETA))
      dim_eta <- nrow(object@ETA)
    else
      dim_eta <- dim_mu

    if (dim_mu > 9)
      return("The dimension of the model is too large.")

    if( max(c(dim_mu, dim_alpha, dim_beta, dim_eta)) != min(c(dim_mu, dim_alpha, dim_beta, dim_eta)) ){
      show(object)
      return("The number of rows of parameter matrix should be equal.")
    }

    if ( ncol(object@MU) > 1 ){
      show(object)
      return("The number of columns of MU matrix should be one.")
    }

    if ( ncol(object@ALPHA) != nrow(object@ALPHA) ){
      show(object)
      return("ALPHA matrix should be n by n matrix.")
    }

    if ( ncol(object@BETA) != nrow(object@BETA) ){
      show(object)
      return("BETA matrix should be n by n matrix.")
    }

    if ( !is.null(object@ETA) )
      if (ncol(object@ETA) != nrow(object@ETA) ){
        show(object)
        return("ETA matrix should be n by n matrix.")
      }
  }

  return(TRUE)

}


#' An S4 calss to represent the specification of a marked Hawkes model with exponential kernel.
#'
#'
#' @slot MU
#' @slot ALPHA
#' @slot BETA
#' @slot ETA
#' @slot Jump
setClass(
  "mHSpec",
  slots = list(
    MU = "matrixORnumeric",
    ALPHA = "matrixORnumeric",
    BETA = "matrixORnumeric",
    ETA = "matrixORnumeric",
    Jump = "function"
    ),
  validity = validmHSpec
)


setMethod(
  "initialize",
  "mHSpec",
  function(.Object, MU, ALPHA, BETA, ETA=NULL, Jump=NULL){

    # If Jump is not provided, then Jump is constant 1.
    if (is.null(Jump)) Jump <- function(n,...) rep(1,n)

    # If ETA is not provided, then ETA = 0 or zero matrix with the same dimension of BETA
    if (is.null(ETA)) {
      if (length(MU) == 1) ETA <- 0
      else {
        if (is.atomic(MU))
          ETA <- 0
        else
          ETA <- matrix(rep(0, length(BETA)), nrow = length(MU))
      }
    }

    .Object@MU <- MU
    .Object@ALPHA <- ALPHA
    .Object@BETA <- BETA
    .Object@ETA <- ETA
    .Object@Jump <- Jump

    # Check spectral radius
    if ( max(abs(eigen(ALPHA/BETA)$values)) >= 1)
      warning("This model does not satisfy the stability condition.")

    callNextMethod()

    .Object


  }
)



setMethod(
  "show",
  "mHSpec",
  function(object){
    cat("Marked Hawkes model with linear impact function\n")
    cat("The intensity process is defined by\n\n")
    cat("LAMBDA(t) = MU + int ALPHA/BETA (1+(k-1)*ETA) exp(-BETA(t-u)) d N(t)\n" )
    cat("\n")
    cat("The parameters are : \n\n")
    callNextMethod(object)
  }
)
