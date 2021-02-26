#' Process raw haul data
#'
#' @param masterDat Data frame with raw haul data
#' @param Truncate Lower bound of hauls, defaults to 0. Values < Truncate will be turned to 0
#' @param strata.limits Data frame of strata limits
#' @param X.pos Optional matrix of covariates for the positive model, dimensioned as same rows as raw data
#' @param X.bin Optional matrix of covariates for the presence-absence model, dimensioned as same rows as raw data
#' @return DataList Data frame of processed data to fit
#'
#' @import grDevices
#' @export
#'
processData <- function(masterDat = NULL, Truncate = 0, strata.limits = NULL, X.pos = NULL, X.bin = NULL) {
  Letters <- apply(MARGIN = 1, FUN = paste, collapse = "", expand.grid(letters, letters))

  print("Necessary column names for masterDat:")
  print("1. BEST_DEPTH_M -> tow depth in meters")
  print("2. BEST_LAT_DD -> tow latitude in degrees ")
  print("3. [insert species name] -> tow catch in kilograms ")
  print("4. YEAR -> calendar year, or time-strata")
  print("5. AREA_SWEPT_MSQ -> area-swept in square meters, or effort offset ")
  print("5. VESSEL -> vessel ID ")
  print("Please ensure that latitude and depth in strata.limits match the following boundaries:")
  print("Latitude: 42-49 in 0.5 increments")
  print("Depth (meters): 55, 75, 100, 125, 155, 183, 200, 250, 300, 350, 400, 450, 500, 549, 600, 700, 800, 900, 1000, 1100, 1200, 1280")

  # Give information about necessary headers
  if (!all(c("BEST_DEPTH_M", "BEST_LAT_DD", species, "YEAR", "AREA_SWEPT_MSQ", "VESSEL") %in% colnames(masterDat))) {
    print("Warning: processData() terminated unsuccessfully.")
    print("Please ensure that masterDat has appropriate column names")
    stop()
  }

  # Start tagged DataList
  DataList <- list()

  # set up the generic data frame for this species
  Data <- data.frame(
    "PROJECT_CYCLE" = masterDat[, "YEAR"],
    "BEST_DEPTH_M" = masterDat[, "BEST_DEPTH_M"],
    "BEST_LAT_DD" = masterDat[, "BEST_LAT_DD"],
    "HAUL_WT_KG" = masterDat[, which(dimnames(masterDat)[[2]] == species)],
    "year" = as.factor(masterDat[, "YEAR"]),
    "effort" = masterDat[, "AREA_SWEPT_MSQ"] * 0.0001,
    "VESSEL" = masterDat[, "VESSEL"]
  )

  # Error checking for covariates
  if (!is.null(X.pos)) {
    if (!is.matrix(X.pos)) X.pos <- matrix(X.pos, ncol = 1)
    if (nrow(X.pos) != nrow(masterDat)) {
      print("Error: covariate matrix for positive model input, but the rows do not match the rows of masterDat")
      stop()
    }
  }
  if (!is.null(X.bin)) {
    if (!is.matrix(X.bin)) X.bin <- matrix(X.bin, ncol = 1)
    if (nrow(X.bin) != nrow(masterDat)) {
      print("Error: covariate matrix for binomial model input, but the rows do not match the rows of masterDat")
      stop()
    }
  }

  if (Truncate > 0) {
    print(paste("Changing any observation with less than ", Truncate, " kilograms to 0 kilograms", sep = ""))
    Data[, "HAUL_WT_KG"] <- ifelse(Data[, "HAUL_WT_KG"] < Truncate, 0, Data[, "HAUL_WT_KG"])
  }
  Data <- cbind(Data, "y" = Data[, "HAUL_WT_KG"])
  Data <- cbind(Data, "strata" = apply(masterDat, 1, strata.fn, Strata.df = strata.limits))
  Data <- cbind(Data, "isNonZeroTrawl" = ifelse(Data[, "y"] > 0, 1, 0))
  Data <- cbind(Data, "ones.vec" = rep(0, length(Data[, "y"]))) # this is just for the 'ones-trick', inv Gaussian
  Data <- cbind(Data, "logy3" = log(pi * 2 * (Data[, "y"]^3))) # this is a constant for the invGaussian
  Data <- cbind(Data, "logy" = log((Data[, "y"]))) # this is a constant for the LognormalECE2 model
  Data <- cbind(Data, "logeffort" = log(Data[, "effort"]))
  Data <- cbind(Data, "effort2" = Data[, "effort"]^2)
  Data <- cbind(Data, "logeffort2" = Data[, "logeffort"]^2)
  Data <- cbind(Data, "lfacty" = lfactorial(Data[, "y"]))
  # Exclude tows from strata that are not included
  Exclude_NoStratum <- which(is.na(Data$strata))
  print(paste("Excluded ", length(Exclude_NoStratum), " observations that were not assigned to any strata", sep = ""))
  if (length(Exclude_NoStratum) > 0) {
    print(paste("Observations that were not assigned to any strata are shown in 'Tows_outside_strata.csv'", sep = ""))
    write.table(Data[Exclude_NoStratum, ], "Tows_outside_strata.csv", row.names = F, col.names = T, sep = ",")
    Data <- Data[-Exclude_NoStratum, ]
  }

  # Exclude tows with some missing entry
  Exclude_Missing <- which(apply(Data, MARGIN = 1, FUN = function(Vec) {
    any(is.na(Vec))
  }))
  print(paste("Excluded ", length(Exclude_Missing), " additional observations that had some missing data", sep = ""))
  write.table(Data[Exclude_Missing, ], "Tows_with_missing_data.csv", row.names = F, col.names = T, sep = ",")
  if (length(Exclude_Missing) < 10 & length(Exclude_Missing) > 0) print(Data[Exclude_Missing, ])
  if (length(Exclude_Missing) >= 10) print("Entries are not printed to the screen due to having 10 or more")
  if (length(Exclude_Missing) > 0) Data <- Data[-Exclude_Missing, ]

  # Count number of tows per strata and year
  TowsPerStrataYear <- table(Data[, "strata"], Data[, "year"])
  write.table(TowsPerStrataYear, "Tows_Per_StrataYear.csv", row.names = F, col.names = T, sep = ",")
  print(paste("Tows per strata and year are displayed below"))
  print(TowsPerStrataYear)

  # Count number of positive catches per strata and year
  EncountersPerStrataYear <- table(Data[, "strata"], Data[, "year"], Data[, "isNonZeroTrawl"])[, , "1"]
  write.table(EncountersPerStrataYear, "Encounters_Per_StrataYear.csv", row.names = F, col.names = T, sep = ",")
  print(paste("Encounters per strata and year are displayed below"))
  print(EncountersPerStrataYear)

  # Redefine variables that depend on year
  Data[, "year"] <- factor(as.numeric(as.character(Data[, "year"]))) # Record year factor in case years have been eliminated due to missing observations
  Data <- data.frame(Data, "vessel" = Letters[as.numeric(as.factor(as.character(Data[, "VESSEL"])))]) # Record year factor in case years have been eliminated due to missing observations
  # Define derived variables involving year
  Data <- cbind(Data, "strataYear" = factor(paste(Data[, "strata"], ":", Data[, "year"], sep = ""), levels = as.vector(outer(sort(unique(Data[, "strata"])), sort(unique(Data[, "year"])), FUN = paste, sep = ":"))))
  Data <- cbind(Data, "vesselYear" = factor(paste(Data[, "vessel"], ":", Data[, "year"], sep = "")))
  # Attach and calculate other values that aren't in the data.frame
  DataList[["Data"]] <- Data
  DataList[["nonZeros"]] <- which(Data[, "isNonZeroTrawl"] == TRUE)
  # Number of elements
  DataList[["nS"]] <- length(unique(Data[, "strata"]))
  DataList[["nY"]] <- length(unique(Data[, "year"]))
  DataList[["nNonZeros"]] <- length(DataList[["nonZeros"]])
  DataList[["n"]] <- length(Data[, "y"])
  DataList[["nSY"]] <- nlevels(Data[, "strataYear"])
  DataList[["nVY"]] <- length(unique(Data[, "vesselYear"]))
  DataList[["nV"]] <- length(unique(Data[, "vessel"]))
  # Diagonal matrix for the wishart / correlation model
  DataList[["R"]] <- diag(2)

  # If the covariates aren't in the R environment, create them
  if (is.null(X.bin)) {
    DataList[["X.bin"]] <- matrix(NA, ncol = 0, nrow = nrow(Data))
  } else {
    DataList[["X.bin"]] <- X.bin
    if (length(Exclude_NoStratum) > 0) DataList[["X.bin"]] <- DataList[["X.bin"]][-Exclude_NoStratum]
    if (length(Exclude_Missing) > 0) DataList[["X.bin"]] <- DataList[["X.bin"]][-Exclude_Missing]
  }
  if (is.null(X.pos)) {
    DataList[["X.pos"]] <- matrix(NA, ncol = 0, nrow = nrow(Data))
  } else {
    DataList[["X.pos"]] <- X.pos
    if (length(Exclude_NoStratum) > 0) DataList[["X.pos"]] <- DataList[["X.pos"]][-Exclude_NoStratum]
    if (length(Exclude_Missing) > 0) DataList[["X.pos"]] <- DataList[["X.pos"]][-Exclude_Missing]
    DataList[["X.pos"]] <- DataList[["X.pos"]][which(Data$isNonZeroTrawl == 1), ]
    if (!is.matrix(DataList[["X.pos"]])) DataList[["X.pos"]] <- as.matrix(DataList[["X.pos"]], ncol = 1)
  }
  DataList[["nX.binomial"]] <- ncol(DataList[["X.bin"]])
  DataList[["nX.pos"]] <- ncol(DataList[["X.pos"]])

  DataList[["log2pi"]] <- log(2 * pi)

  # return
  return(DataList)
}
