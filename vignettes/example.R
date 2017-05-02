## ----install-------------------------------------------------------------
# Install package
install.packages("devtools", repos = "http://cran.us.r-project.org")
library("devtools")
install_github("nwfsc-assess/nwfscDeltaGLM")
  
# Load package
library(nwfscDeltaGLM)

## ------------------------------------------------------------------------
data(SA3)

## ------------------------------------------------------------------------
data(Example_Species)
names(Example_Species)

## ----load----------------------------------------------------------------
# Load data and strata
data(Example_Species)
masterDat = Example_Species

strata.limits = data.frame("STRATA" = c("First", "Second", "Mid", "D", "South"),
  NLat = c(49, 49, 49, 36, 34.5), 
  SLat = c(36, 36, 34.5, 34.5, 34.5),
  MinDepth = c(55, 183, 549, 55, 55),
  MaxDepth = c(183, 549, 700, 549, 549), stringsAsFactors = FALSE)

## ------------------------------------------------------------------------
names(masterDat)

species = "Example_Species"
names(masterDat)[which(names(masterDat)=="EXPANDED_WT_KG")] = "Example_Species"

