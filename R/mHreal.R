print.mHreal <- function(res, n=10){
  cat(".........\n")
  print.default(res)
}


as.matrix.mHreal <- function(res){

  mtrx <- numeric()
  for (i in 2:length(res)){
    mtrx <- cbind(mtrx, res[[i]])
    if(is.vector(res[[i]])){
      colnames(mtrx)[i-1] <- names(res)[i]
    }
  }
  mtrx
}

as.data.frame.mHreal <- function(res){
  as.data.frame(as.matrix(res))
}

summary.mHreal <- function(res, n=20){

  cat("------------------------------------------\n")
  cat("Simulation result of marked Hawkes model.\n")
  cat("Realized path:\n")
  mtrx <- as.matrix(res)
  dimens <- length(res$mHSpec@MU)
  name_N  <- paste0("N", 1:dimens)
  name_lambda  <- paste0("lambda", 1:dimens)

  len <- min(n, length(mtrx[,"arrival"]))

  print(mtrx[1:len, c("arrival", name_N, name_lambda)])
  if ( length(mtrx[,"arrival"]) > len){

    remaning <- length(mtrx[,"arrival"]) - len

    cat("... with ")
    cat(remaning)
    cat(" more rows \n")
  }

  cat("------------------------------------------\n")
}

#' Get left continuous version of lambda process
#'
#' The realized version of the lambda process in \code{mHreal} is right continuous version.
#' If the left continuous version is needed, this function is applied.
#'
#' @param res \code{mHreal} an S3 class contains the realized lambda processes.
#'
#'
#' @examples
get_lc_lambda <- function(res){

  dimens <- length(res$mHSpec@MU)
  lc_lambda_component <- res$lambda_component

  for (i in 2:nrow(lc_lambda_component)) {

    if (dimens == 1) {
      impact <- res$mHSpec@ALPHA * ( 1 + (res$mark[i] - 1 ) * res$mHSpec@ETA )
      print(impact)
      lc_lambda_component[i] <- lc_lambda_component[i] - impact

    } else {

      col_indx <- seq(res$jump_type[i], dimens^2, dimens)
      lc_lambda_component[i, col_indx] <- lc_lambda_component[i, col_indx] -
        res$mHSpec@ALPHA[, res$jump_type[i]] * ( 1 + (res$mark[i] - 1 ) * res$mHSpec@ETA[, res$jump_type[i]])

    }

  }


  lc_lambda_component
}


plot.mHreal <- function(res, ...){

  dimens <- ncol(res$N)
  par(mfrow=c(dimens, 1))

  n <- length(res$arrival)

  for (i in 1:dimens) {
    plot(res$arrival[1:n], res$N[ ,i][1:n], 's', xlab='t', ylab=colnames(res$N)[i])
    #points(res$arrival[1:n], res$N[ ,i][1:n])
  }

}



#' Plot exponentially decaying lambda process
#'
#' This plot method describes the exponentially decaying lambda (intensity) process.
#'
#'
#' @param arrival a vector of arrival times.
#' @param lambda a vector of lambda processs.
#' @param beta a decaying parameter lambda.
#' @param dt a small step size on time horizon.
#'
#'
#' @examples
#' MU1 <- 0.3; ALPHA1 <- 1.5; BETA1 <- 2
#' mHSpec1 <- new("mHSpec", MU=MU1, ALPHA=ALPHA1, BETA=BETA1)
#' # Simulate with mHSim funciton.
#' res1 <- mHSim(mHSpec1,  n=100)
#' plotlambda(res1$arrival, res1$lambda, BETA1)
plotlambda <- function(arrival, lambda, beta, dt = NULL, ...){
  maxT <- tail(arrival, n=1)

  if (is.null(dt)){
    inter_arrival <- arrival[-1] - arrival[-length(arrival)]
    dt <- mean(inter_arrival)/100
  }

  time_vector <- seq(0, maxT, dt)

  lambda_vector <- numeric(length=length(time_vector))

  lambda_vector[1] <- lambda[1]
  j <- 2
  for (i in 2:length(time_vector)){
    next_arrival <- arrival[j]
    next_lambda <- lambda[j]
    if (time_vector[i] <  next_arrival){
      lambda_vector[i] <- lambda_vector[i-1]*exp(-beta*dt)
    } else{
      lambda_vector[i] <- next_lambda*exp(-beta * (time_vector[i] - next_arrival))
      j <- j + 1
    }
  }
  plot(time_vector, lambda_vector, 'l', xlab='t', ylab='lambda')

}

