

setClassUnion("matrixORvector", c("matrix", "vector"))

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
    MU = "matrixORvector",
    ALPHA = "matrixORvector",
    BETA = "matrixORvector",
    ETA = "matrixORvector",
    Jump = "function"
    )
)


setMethod(
  "initialize",
  "mHSpec",
  function(.Object, MU, ALPHA, BETA, ETA=NULL, Jump=NULL){

    # If Jump is not provided, then Jump is constant 1.
    if (is.null(Jump)) Jump <- function(n,...) rep(1,n)

    # length of each parameter should be 1 or even
    if ( length(ALPHA) != 1 & length(ALPHA) %% 2 != 0 ) stop("The length of ALPHA must be 1 or even")
    if ( length(BETA) != 1 & length(BETA) %% 2 != 0 ) stop("The length of BETA must be 1 or even")

    # If ETA is not provided, then ETA = 0 or zero matrix with the same dimension of BETA
    if (is.null(ETA)) {
      if (length(MU) == 1) ETA <- 0
      else ETA <- matrix(rep(0, length(BETA)), nrow = length(BETA)/2)
    }
    if ( length(ETA) != 1 & length(ETA)%%2 != 0 ) stop("The length of ETA must be 1 or even")



    # When one of MU, ALPHA, BETA, ETA is not one dimensional
    if ( length(MU)!=1 | length(ALPHA)!=1 | length(BETA)!=1 | length(ETA)!=1 ) {


      if ( !is.matrix(MU) ) MU <- matrix(MU)
      if ( !is.matrix(ALPHA) ) ALPHA <- matrix(ALPHA, nrow = length(ALPHA)/2)
      if ( !is.matrix(BETA) ) BETA <- matrix(BETA, nrow = length(BETA)/2)
      if ( !is.matrix(ETA) ) ETA <- matrix(ETA, nrow = length(ETA)/2)

      # The dimensions of parameter matrices should be consistent.
      if ( !( nrow(MU) == nrow(ALPHA) & nrow(ALPHA) == nrow(BETA) & nrow(BETA) == nrow(ETA)) )
        stop("The number of rows of parameter matrix should be equal")

      # The dimension check for each parameter matrix.
      if ( ncol(MU) > 1 )
        stop("The number of columns of MU matrix should be one.")

      if ( ncol(ALPHA) != nrow(ALPHA) )
        stop("ALPHA matrix should be n by n matrix.")

      if ( ncol(BETA) != nrow(BETA) )
        stop("BETA matrix should be n by n matrix.")

      if ( ncol(ETA) != nrow(ETA) )
        stop("ETA matrix should be n by n matrix.")

    }

    # Check spectral radius
    if ( max(abs(eigen(ALPHA/BETA)$values)) >= 1) warning("This model does not satisfy the stability condition.")

    .Object@MU <- MU
    .Object@ALPHA <- ALPHA
    .Object@BETA <- BETA
    .Object@ETA <- ETA
    .Object@Jump <- Jump

    .Object
  }
)


