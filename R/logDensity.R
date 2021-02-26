#' This function is in development stage, but used to generate the posterior predictive distribution score, or log score, from a fitted model object. Currently, it's only applicable to the gamma, lognormal, or inverse Gaussian likelihoods.
#'
#' @param obj A fitted modle object
#' @return avg.bin The log score for the presence-absence model
#'
logDensity <- function(obj) {
  #################################
  # To do: ECEs
  #################################
  attach.jags(obj) # attach object to workspace

  # Divide data into a zero and non-zero component
  # For each data point (i = 1:N), calculate the appropriate
  # likelihood, and avg over all mcmc draws. Then multiply these
  # values OR sum in log-space over data points.

  strata <- Data$strata
  year <- Data$year
  strataYear <- Data$strataYear
  vesselYear <- Data$vesselYear
  effort <- Data$effort
  effort2 <- Data$effort2
  logeffort <- Data$logeffort
  logeffort2 <- Data$logeffort2

  n.mcmc <- dim(Sdev)[1]
  # calculate the predicted components for the positive model
  avg.pos <- 0
  # calculate the predicted components for the positive model
  for (j in 1:length(nonZeros)) { # loop over records
    # for the jth data point, over all mcmc draws calculate expected value
    exp.value <- Sdev[, as.numeric(strata[nonZeros[j]])] + Ydev[, as.numeric(year[nonZeros[j]])] + SYdev[, as.numeric(strataYear[nonZeros[j]])] + VYdev[, as.numeric(vesselYear[nonZeros[j]])] + logeffort[nonZeros[j]] * B.pos[, 1] + logeffort2[nonZeros[j]] * B.pos[, 2]
    if (obj$likelihood == "gamma") {
      gamma.a <- (1 / CV[, 1])^2 # gamma.a is constant over data pts
      gamma.b <- gamma.a / exp(exp.value) # b = a/E[x]
      avg.pos[j] <- log(mean(dgamma(y[nonZeros[j]], shape = gamma.a, rate = gamma.b)))
    }
    if (obj$likelihood == "lognormal") {
      sigma <- sqrt(log(CV[, 1]^2 + 1))
      avg.pos[j] <- log(mean(dlnorm(y[nonZeros[j]], meanlog = exp.value, sdlog = sigma)))
    }
    if (obj$likelihood == "invGaussian") {
      oneOverCV2 <- (1 / CV[, 1])^2
      u <- exp(exp.value)
      lambda <- u / oneOverCV2
      avg.pos[j] <- log(mean(dinvgauss(rep(y[nonZeros[j]], n.mcmc), u, lambda)))
    }
  }

  # do the same thing with the binomial model
  avg.bin <- 0
  for (j in 1:dim(Data)[1]) {
    # for this data point calculate posterior distribution of expectation
    exp.value <- plogis(pSdev[, as.numeric(strata[j])] + pYdev[, as.numeric(year[j])] + pSYdev[, as.numeric(strataYear[j])] + pVYdev[, as.numeric(vesselYear[j])] + effort[j] * B.zero[, 1] + effort2[j] * B.zero[, 2])
    # calculate log of the integral over parameters, P(y|y,theta)
    avg.bin[j] <- log(mean(dbinom(rep(Data$isNonZeroTrawl[j], length = dim(Data)[1]), prob = exp.value, size = 1)))
  }
  # return(sum(avg.bin)+sum(avg.pos))
  return(list("Presence.score" = avg.bin, "Positive.score" = avg.pos))
}
