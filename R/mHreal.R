print.mHreal <- function(res, n=10){
  cat(".........\n")
  print.default(res)
  # for (i in 1:length(res)){
  #
  #   if(is.vector(res[i])){
  #     print(res[[i]][1:n])
  #   } else{
  #     print(res[i])
  #   }
  # }
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

summary.mHreal <- function(res){

  cat("----------------------------------------\n")
  cat("Simulation result of marked Hawkes model\n")
  cat("realized path:\n")
  mtrx <- as.matrix(res)
  dimens <- length(res$mHSpec@MU)
  name_N  <- paste0("N", 1:dimens)

  len <- min(10, length(mtrx[,"arrival"]))

  print(mtrx[1:len, c("arrival",name_N)])
  if ( length(mtrx[,"arrival"]) > len){

    remaning <- length(mtrx[,"arrival"]) - len

    cat("... with ")
    cat(remaning)
    cat(" more rows \n")
  }

  cat("----------------------------------------\n")
}


get_lc_lambda <- function(res){

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

ggplot.mHreal <- function(res){

  dfr <- as.data.frame(res)

  dfr2 <- melt(dfr[,c("arrival", "N1", "N2")], id.vars="arrival", value.name="N", variable.name="Ns")

  ggplot(data=dfr2, aes(x=arrival, y=N, group = Ns, colour = Ns)) +
    geom_step()
}

plotlambda <- function(arrival, lambda, beta, dt = NULL, ...){
  maxT <- tail(arrival, n=1)

  if (is.null(dt)){
    inter_arrival <- arrival[2:length(arrival)] - arrival[1:(length(arrival)-1)]
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

