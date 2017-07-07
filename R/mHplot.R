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
