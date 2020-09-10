# MAXIMUM LIKELIHOOD ESTIMATION: NEWTON-RAPHSON ALGORITHM
# EXPONENTIAL DISTRIBUTION

MLE_NR_Exponential <- function(y, maxiter = 100, epsilon = 0.000000001, stop_criteria = 10000){

  n <- length(y)
  ysum <- sum(y)

  U <- function(n, ysum, beta){
    -(n/beta) + (1/beta^2) * ysum
  }

  Uprime <- function(n, ysum, beta){
    (n/beta^2) - (2/beta^3) * ysum
  }

  # first guess:
  beta0 <- min(y)
  Estimator <- matrix(NA, ncol = 1)
  Estimator[1, 1] <- beta0

  iter <- 1
  while( (stop_criteria > epsilon) & (iter <= maxiter) ){

    num <- U(n, ysum, Estimator[iter, 1])
    den <- Uprime(n, ysum, Estimator[iter, 1])

    UPD <- -(num/den)
    Estimator_iter <- as.matrix(Estimator[iter, 1]) + UPD
    Estimator <- rbind(Estimator, Estimator_iter)
    stop_criteria <- UPD^2
    iter <- iter + 1
  }

  Likelihood <- matrix(NA, ncol = 1)
  LogLikelihood <- matrix(NA, ncol = 1)

  for(i in 1:length(Estimator)){
    Likelihood[i] <- prod( (1/Estimator[i,1]) * exp(-(y/Estimator[i,1])) )
  }

  for(i in 1:length(Estimator)){
    LogLikelihood[i] <- sum( -log(Estimator[i,1]) - (y/Estimator[i,1]) )
  }

  results <- cbind(format(Estimator, nsmall = 6), Likelihood, LogLikelihood)
  colnames(results) <- c("ML Estimator", "Likelihood", "Log-Likelihood")
  return(results)

}

