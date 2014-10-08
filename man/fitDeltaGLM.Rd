\name{fitDeltaGLM}
\alias{fitDeltaGLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Delta-Generalized Linear Mixed Models
}
\description{
\code{fitDeltaGLM} is the core estimation function that constructs the statistical delta GLMM model from user input, and performs Bayesian estimation using JAGS. 
}
\usage{
fitDeltaGLM(modelStructure = list(StrataYear.positiveTows = "random", VesselYear.positiveTows = "random", StrataYear.zeroTows = "random", VesselYear.zeroTows = "random", Catchability.positiveTows = "one", Catchability.zeroTows = "zero", year.deviations = "uncorrelated", strata.deviations = "uncorrelated"), covariates = list(positive = FALSE, binomial = FALSE), likelihood = "gamma", model.name = "deltaGLM.txt", fit.model = TRUE, mcmc.control = list(chains = 5, thin = 1, burn = 5000, iterToSave = 2000), Parallel = TRUE, Species = "NULL", logitBounds = c(-20, 20), logBounds = c(-20, 20), prior.scale = rep(25, 4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelStructure}{List specifying the model structure, including strata, year, and vessel effects
%%     ~~Describe \code{modelStructure} here~~
\itemize{
    \item \code{StrataYear.positiveTows} Character string that is one of the following: "zero", "fixed", random", "randomExpanded", or "correlated". The "random" and "randomExpanded" treatments of normally distributed random effects differ in that "random" assigns an Inverse Gamma prior on the random effects variance, while "randomExpanded" assigns the folded non-central t-distribution as a prior (Gelman 2006, Gelman et al. 2007). Specifying the effects as "correlated" also models deviations as random effects, but treats them as multivariate normal, with a correlation between the presence-absence model estimated. Defaults to "random".
    \item \code{VesselYear.positiveTows} Character string that is one of the following: "zero", "fixed", random", "randomExpanded", or "correlated". The "random" and "randomExpanded" treatments of normally distributed random effects differ in that "random" assigns an Inverse Gamma prior on the random effects variance, while "randomExpanded" assigns the folded non-central t-distribution as a prior (Gelman 2006, Gelman et al. 2007). Specifying the effects as "correlated" also models deviations as random effects, but treats them as multivariate normal, with a correlation between the presence-absence model estimated. Defaults to "random".
    \item \code{StrataYear.zeroTows} Character string that is one of the following: "zero", "fixed", random", "randomExpanded", or "correlated". The "random" and "randomExpanded" treatments of normally distributed random effects differ in that "random" assigns an Inverse Gamma prior on the random effects variance, while "randomExpanded" assigns the folded non-central t-distribution as a prior (Gelman 2006, Gelman et al. 2007). Specifying the effects as "correlated" also models deviations as random effects, but treats them as multivariate normal, with a correlation between the positive model estimated. Defaults to "random".
    \item \code{VesselYear.zeroTows} Character string that is one of the following: "zero", "fixed", random", "randomExpanded", or "correlated". The "random" and "randomExpanded" treatments of normally distributed random effects differ in that "random" assigns an Inverse Gamma prior on the random effects variance, while "randomExpanded" assigns the folded non-central t-distribution as a prior (Gelman 2006, Gelman et al. 2007). Specifying the effects as "correlated" also models deviations as random effects, but treats them as multivariate normal, with a correlation between the positive model estimated. Defaults to "random". 
    \item \code{Catchability.positiveTows}Specify offset for the positive model. Can be fixed ("one"), or estimated ("linear", "quadratic"). If "E" represents effort, the offset is included in the following form, ln(u) = ... + B1*ln(E) + B2*(ln(E)^2). Defaults to "one".
    \item \code{Catchability.zeroTows}Specify offset for the presence-absence model. Can be fixed ("zero","one"), or estimated ("linear", "quadratic"). If "E" represents effort, the offset is included in the following form, logit(p) = ... + B1*E + B2*(E^2). Defaults to "zero".
    \item \code{year.deviations} Argument ("correlated","uncorrelated") specifying whether year deviations should be estimated as correlated random effects (correlation between presence-absence and positive model deviates estimated). By default, this is "uncorrelated", meaning year deviations are included as fixed effects. 
    \item \code{strata.deviations} Argument ("correlated","uncorrelated") specifying whether strata deviations should be estimated as correlated random effects (correlation between presence-absence and positive model deviates estimated). By default, this is "uncorrelated", meaning strata deviations are included as fixed effects.
    }
}
  \item{covariates}{List specifying whether covariates should be included in either model
    \itemize{  
    \item \code{positive} Whether to include covariates in positive model (defaults to FALSE)
    \item \code{binomial} Whether to include covariates in presence-absence model (defaults to FALSE)
    }
  }
  \item{likelihood}{Character string specifying the form of the positive model. Can be one of the following: "gamma" (or "gammaFixedCV"), "lognormal" (or "lognormalFixedCV"), "invGaussian" (or "invGaussianFixedCV"), "lognormalECE", "gammaECE", "poisson", "zt.poisson", or "negbin". Defaults to "gamma". The forms of the model are as follows:
  \itemize{  
  \item \code{gamma} Models the response as continuous from a Gamma distribution, with a log-link. The form of the Gamma distribution used is Y ~ Gamma(shape = a, rate = b), where E[Y] = a / b, and CV(Y) = 1 / sqrt(a). For consistency with the \code{lognormal}, the CV^2 is assigned an Inverse Gamma (0.001,0.001) prior. If the distribution is specified as "gammaFixedCV", the CV = 1.
  \item \code{lognormal} Models the response as continuous from a lognormal distribution, with a log-link. The CV^2 is assigned an Inverse Gamma (0.001,0.001) prior, and the variance of the lognormal distribution is calculated as var = (log(CV^2)+1). If the distribution is specified as "lognormalFixedCV", the CV = 1.
  \item \code{invGaussian} Models the response as continuous from an inverse Gaussian (Wald) distribution, with a log-link. The parameterization is in terms of the mean (u) and variance (u^3 / lambda). This model is considerably slower than the \code{lognormal} or \code{gamma} because the "zeros" trick needs to be implemented. For consistency with the other distributions, the CV^2 is assigned an Inverse Gamma (0.001,0.001) prior.  If the distribution is specified as "invGaussianFixedCV", the CV = 1.
  \item \code{gammaECE} Extends the \code{gamma} distribution to model the positive distribution as a 2-component mixture (normal and "extreme catch events", Thorson 2011). The probability of membership in each class ("normal", "extreme") is estimated. Extreme catch events are allowed to have a separate estimated variance parameter (or CV). The mean of extreme catch events is also estimated; for identifiability, this is estimated as an offset from normal catch events in link-space, e.g. log(u_extreme) = log(u_normal) + logratio, where the "logratio" parameter is assigned a log-uniform(0,5) prior distribution.  
  \item \code{lognormalECE} Extends the \code{lognormal} distribution to model the positive distribution as a 2-component mixture (normal and "extreme catch events", Thorson 2011). The probability of membership in each class ("normal", "extreme") is estimated. Extreme catch events are allowed to have a separate estimated variance parameter (or CV). The mean of extreme catch events is also estimated; for identifiability, this is estimated as an offset from normal catch events in link-space, e.g. log(u_extreme) = log(u_normal) + logratio, where the "logratio" parameter is assigned a log-uniform(0,5) prior distribution.    
  \item \code{poisson} Models the response as discrete counts from the Poisson distribution, parameterized using the log-link function, Y ~ Poisson(u), with log(u) = B0 + B1*X...Because this is a delta-GLM model, applying the Poisson distribution isn't really appropriate when counts are small (zeros are a possibility). The \code{zt.poisson} distribution is more appropriate for these settings.
  \item \code{zt.poisson} Models the response as discrete counts from the zero-truncated Poisson distribution, parameterized using the log-link function, Y ~ Poisson(u), with log(u) = B0 + B1*X...This is more appropriate for small count data than the \code{poisson}. Because this implements the "zeros" trick to calculate the likelihood, estimation is slower than the traditional Poisson.
  \item \code{negbin} Models the response as discrete counts from the Negative binomial distribution, parameterized using the log-link function, Y ~ NegBin(u, r), with log(u) = B0 + B1*X... This parameterization 
  }
}
  \item{model.name}{Character string specifying the name of the JAGS txt file that is written to the working directory. Defaults "to deltaGLM.txt"
}
  \item{fit.model}{Boolean, specifying whether the model should be fit or not. Defaults to TRUE; if FALSE, just the JAGS script file is written.
}
  \item{mcmc.control}{List, used to specify the MCMC parameters. These are 
    \itemize{  
    \item \code{chains} Number of MCMC chains. More is slower (defaults to 5)
    \item \code{thin} Thinning rate of saved samples (defaults to 1)
    \item \code{burn} Length of burn-in period (defaults to 5000)
    \item \code{iterToSave} Number of samples to save after the burn-in (defaults to 2000)  
    }
}
  \item{Parallel}{Whether to conduct estimation using parallel computing (faster), using the function jags.parallel(). Defaults to TRUE, but is currently only possible on Windows machines. If specified as TRUE on another operating system, estimation will be done, but not in parallel. 
}
  \item{Species}{Character string specifying the name of the species analyzed. This is mostly to prevent confusion with output, but must be entered as something other than "NULL".
}
  \item{logitBounds}{2-element vector specifying bounds in the logit-link space, defaults to (-20, 20)
}
  \item{logBounds}{2-element vector specifying bounds in the log-link space, defaults to (-20, 20)
}
  \item{prior.scale}{Scale parameter for the variance parameter of the "randomExpanded" random effects. This is a 4 element vector, defaults to (25, 25, 25, 25), where the elements are (strataYear.positive, strataYear.presenceAbsence, vesselYear.positive, vesselYear.presenceAbsence)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An list object with the following components:
  %\item{modelFit}{ The fitted JAGS object using R2jags. To access elements, use attach.jags(modelFit). }
    %\itemize {
    \item{modelFit}{ The fitted JAGS object using R2jags. To access elements, use attach.jags(modelFit).}
    \item{functionCall}{ A list containing the following elements. Many recycled from the \code{fitDeltaGLM} function arguments, but designed to attach all important components to the fitted model object.
    \itemize{
    \item{modelStructure}{ The arguments used to specify the form of the model, set by user}
    \item{model.name}{ The .txt file specified by the user}
    \item{fit.model}{ Whether the model was fit}
    \item{mcmc.control}{ The list of MCMC parameters}
    \item{Parallel}{ Whether parallel computing was used}
    \item{likelihood}{ The distribution for the positive model}
    \item{covariates}{ Whether covariates were included}
    }
    }
    \item{estimatedParameters}{ A data frame containing the names of estimated parameters.}
    \item{Species}{ The species name inputted by the user.}
    \item{Data}{ The data frame (raw) used to estimate all parameters.}
    %}
}
\references{
%% ~put references to the literature/web site here ~
Gelman, A. 2006. Prior distributions for variance parameters in hierarchical models. Bayesian analysis, 1(3):515-533.

Gelman, A., D.A. Van Dyk, Z. Huang, and W.J. Boscardin. 2007. Using redundant parameterizations to fit hierarchical models. Journal of Computational and Graphical Statistics, 17(1), 95-122.

Plummer, M. 2012. JAGS Version 3.3.0 user manual. http://mcmc-jags.sourceforge.net/

Thorson, J.T., I. Stewart, and A.E. Punt.  2011. Accounting for fish shoals in single- and multi-species survey data using mixture distribution models.  Can. J. Fish. Aquat. Sci., 68(9): 1681-1693.

Thorson, J.T. and E.J. Ward. 2013. Accounting for space-time interactions in index standardization models. Fisheries Research, 147: 426-433.
}
\author{
%%  ~~who you are~~
Eric Ward and Jim Thorson, with help from the FRAM assessment team at the Northwest Fisheries Science Center. 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
