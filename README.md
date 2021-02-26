
Description
=============

nwfscDeltaGLM
* Is an R package for implementing a Bayesian stratified delta-GLMM for use when standardizing fishery-independent index data for U.S. West Coast surveys.
* Has built in diagnostic functions and model-comparison tools
* Is intended to improve analysis speed, replicability, peer-review, and interpretation of index standardization methods


Instructions
=============

First, please install JAGS (Just Another Gibbs Sampler) here: http://mcmc-jags.sourceforge.net/

Next, please use R version >=3.1.1 and install the package:

    # Install package
    install.packages("devtools")
    library("devtools")
    # Current releases can be installed with
    install_github("nwfsc-assess/nwfscDeltaGLM")
    # Specific version can be installed as
    # install_github("nwfsc-assess/nwfscDeltaGLM", ref="1.0.0") 
    # Load package
    library(nwfscDeltaGLM) 

Please see vignettes folder (i.e., https://github.com/nwfsc-assess/nwfscDeltaGLM/blob/master/vignettes/example.Rmd) for an example of how to run the model. Or https://github.com/nwfsc-assess/nwfscDeltaGLM/blob/master/vignettes/simulation_example.Rmd for an example of how to simulate and fit some different data structures. 

Known installation/usage issues
=============
Please read the following list of known problems and solutions
* Users sometimes have trouble loading the package, with R throwing an error about package "rjags".  Please reinstall the latest version of JAGS and then re-install nwfscDeltaGLM

Further reading
=============

For more details regarding development and testing of this delta-GLMM software please see:
* Shelton, A. O., J. T. Thorson, E. J. Ward, and B. E. Feist. 2014. Spatial semiparametric models improve estimates of species abundance and distribution. Can. J. Fish. Aquat. Sci. 71, 1655–1666. doi:10.1139/cjfas-2013-0508. URL: http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0508
* Thorson, J. T., and E. J. Ward. 2014. Accounting for vessel effects when standardizing catch rates from cooperative surveys. Fisheries Research 155:168–176.
* Thorson, J. T., and E. Ward. 2013. Accounting for space-time interactions in index standardization models. Fisheries Research 147:426–433.
* Thorson, J. T., I. J. Stewart, and A. E. Punt. 2012. Development and application of an agent-based model to evaluate methods for estimating relative abundance indices for shoaling fish such as Pacific rockfish (Sebastes spp.). ICES Journal of Marine Science 69:635–647.
* Thorson, J. T., I. Stewart, and A. Punt. 2011. Accounting for fish shoals in single- and multi-species survey data using mixture distribution models. Canadian Journal of Fisheries and Aquatic Sciences 68:1681–1693.
* Helser, T. E., A. E. Punt, and R. D. Methot. 2004. A generalized linear mixed model analysis of a multi-vessel fishery resource survey. Fisheries Research 70:251–264.


