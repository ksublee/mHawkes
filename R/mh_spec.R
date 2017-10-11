setClassUnion("matrixORnumeric", c("matrix", "numeric"))

#' Check the validity of mhspec
#'
#' This function checks the validity of mhspec.
#' For one dimensional case, if one of the parameter is not matrix, then all parameters should not be matrix.
#' If one of the parameter is atomic, then all parameters should be atomic.
#' If all parameters are matrix, then the nrow should be equal.
#' The number of columns of MU matrix should be one.
#' ALPHA, BETA, ETA matrices should be sqaure matrices.
#' The dimension of the model is less than 10.
#'
#' @param object S4-class of mhspec
valid_mhspec <- function(object) {

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

#' An S4 class to represent a marked Hawkes model
#'
#' This class represents a specification of a marked Hawkes model with exponential kernel.
#' The intensity of the ground process is defined by:
#' \deqn{\lambda(t) = \mu + \int \alpha / \beta (1 + (k - 1) \eta) exp( -\beta (t-u)) d N(t)}
#'
#' @slot MU numeric value or matrix. Shoud be a matrix for two or larger dimensional model.
#' @slot ALPHA numeric value or matrix. Shoud be a matrix for two or larger dimensional model.
#' @slot BETA numeric value or matrix. Shoud be a matrix for two or larger dimensional model.
#' @slot ETA numeric value or matrix. Shoud be a matrix for two or larger dimensional model.
#' @slot mark a function that generate a random number with specific distribution.
#' @slot impact a function that describe the mark impact
#'
#' @examples
#' MU2 <- matrix(c(0.2), nrow = 2)
#' ALPHA2 <- matrix(c(0.75, 0.92, 0.92, 0.75), nrow = 2, byrow=TRUE)
#' BETA2 <- matrix(c(2.25, 2.25, 2.25, 2.25), nrow = 2, byrow=TRUE)
#' ETA2 <- matrix(c(0.19, 0.19, 0.19, 0.19), nrow = 2, byrow=TRUE)
#' mark2 <- function(n,...) rgeom(n, 0.65) + 1
#' mhspec2 <- new("mhspec", MU=MU2, ALPHA=ALPHA2, BETA=BETA2, ETA=ETA2, mark = mark2)
setClass(
  "mhspec",
  slots = list(
    MU = "matrixORnumeric",
    ALPHA = "matrixORnumeric",
    BETA = "matrixORnumeric",
    ETA = "matrixORnumeric",
    mark = "function"
    #impact = "function"
  ),
  validity = valid_mhspec
)

setMethod(
  "initialize",
  "mhspec",
  function(.Object, MU, ALPHA, BETA, ETA=NULL, mark=NULL, impact=NULL, stability_check=FALSE){

    # If mark is not provided, then mark is constant 1.
    if (is.null(mark)) mark <- function(n,...) rep(1,n)

    # If ETA is not provided, then ETA = 0 or zero matrix with the same dimension of BETA
    if (is.null(ETA)) {
      if (length(MU) == 1){
        ETA <- 0
      } else {
        ETA <- matrix(rep(0, length(BETA)), nrow = length(MU))
      }
    }

    .Object@MU <- MU
    .Object@ALPHA <- ALPHA
    .Object@BETA <- BETA
    .Object@ETA <- ETA

    # check the number of arguments
    if(length(formals(mark)) == 1){
      .Object@mark <- function(n, ...) mark(n)
    } else {
      .Object@mark <- mark
    }

    # Check spectral radius
    if ( stability_check==TRUE && max(abs(eigen(ALPHA/BETA)$values)) >= 1)
      warning("This model does not satisfy the stability condition.")

    callNextMethod()

    .Object


  }
)


setMethod(
  "show",
  "mhspec",
  function(object){

    dimens <- length(object@MU)
    cat(paste0(toString(dimens), "-dimensional (marked) Hawkes model with linear impact function.\n"))
    cat("The intensity process is defined by\n")
    cat("LAMBDA(t) = MU + int ALPHA %/% BETA (1+(k-1)ETA) %*% exp(-BETA(t-u)) d N(t)\n" )
    cat("\n")
    cat("Parameters: \n")

    MU <- object@MU
    cat("MU: \n")
    print(MU)

    ALPHA <- object@ALPHA
    cat("ALPHA: \n")
    print(ALPHA)

    BETA <- object@BETA
    cat("BETA: \n")
    print(BETA)

    ETA <- object@ETA
    cat("ETA: \n")
    print(ETA)

    cat("Mark distribution: \n")
    print(object@mark)

    #if(!is.null(impact)){
    #  cat("Impact function: \n")
    #  print(object@impact)
    #}

    cat("------------------------------------------\n")



  }
)


name_unique_coef_mtrx <- function(M, notation){
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


setMethod(
  "coef",
  "mhspec",
  function(object, uniqueness = FALSE){

    dimens <- length(object@MU)

    coeff <- numeric()
    name <- character()

    if (uniqueness == FALSE | dimens == 1){

      #MU
      name <- c(name, sapply(seq(1, dimens), function(x) paste0("mu", x)) )

      #ALPHA
      for (i in 1:dimens){
        for (j in 1:dimens){
          name <- c(name, paste0("alpha", i, j))
        }
      }

      #BETA
      for (i in 1:dimens){
        for (j in 1:dimens){
          name <- c(name, paste0("beta", i, j))
        }
      }

      #ETA
      for (i in 1:dimens){
        for (j in 1:dimens){
          name <- c(name, paste0("eta", i, j))
        }
      }

      coeff <- c(as.vector(object@MU), as.vector(t(object@ALPHA)), as.vector(t(object@BETA)), as.vector(t(object@ETA)))


    } else {

      #MU
      unique_mus <- unique(as.vector(object@MU))
      name <- c(name, unique(name_unique_coef_mtrx(object@MU, "mu")))

      #ALPHA
      unique_alphas <- unique(as.vector(t(object@ALPHA)))
      name <- c(name, unique(name_unique_coef_mtrx(object@ALPHA, "alpha")))

      #BETA
      unique_betas <- unique(as.vector(t(object@BETA)))
      name <- c(name, unique(name_unique_coef_mtrx(object@BETA, "beta")))

      #ETA
      unique_etas <- unique(as.vector(t(object@ETA)))
      name <- c(name, unique(name_unique_coef_mtrx(object@ETA, "eta")))

      coeff <- c(unique_mus, unique_alphas, unique_betas, unique_etas)

    }

    names(coeff) <- name
    coeff

  }
)


