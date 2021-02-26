#' Plot data
#'
#' @param Data The data frame storing raw tow data
#' @param FileName The name of the output file -- this will be species and model specific. Defaults to 'TowMap.jpg'
#' @param Folder Where the output is stored, defaults to NA
#'
#' @import grDevices
#' @export
#'
PlotData <- function(Data, FileName, Folder = NA) {
  if (is.na(Folder)) Folder <- getwd()

  Data <- cbind(Data, "Pres" = ifelse(Data[, "HAUL_WT_KG"] > 0, 1, 0))
  Pos <- Data[which(Data[, "HAUL_WT_KG"] > 0), ]
  Nyears <- length(unique(Data[, "PROJECT_CYCLE"]))

  # Histogram of positive catch | year
  Ncol <- ceiling(sqrt(Nyears))
  Nrow <- ceiling(Nyears / Ncol)
  jpeg(filename = paste(Folder, FileName, "Positive catch BY Year.jpg", sep = ""), width = Ncol * 3, height = Nrow * 3, units = "in", res = 200)
  par(mfrow = c(Nrow, Ncol), mar = c(2, 2, 2, 0), mgp = c(1.25, 0.25, 0), tck = -0.02, oma = c(4, 4, 0, 0))
  for (YearI in 1:Nyears) {
    Which <- which(Pos[, "PROJECT_CYCLE"] == unique(Pos[, "PROJECT_CYCLE"])[YearI])
    hist(Pos[Which, "HAUL_WT_KG"], main = unique(Pos[, "PROJECT_CYCLE"])[YearI], xlab = "", ylab = "", breaks = seq(0, max(Pos[, "HAUL_WT_KG"]) + 1, length = 100), cex.main = 1.5)
  }
  mtext("Positive catch rates", outer = TRUE, line = 2, side = 1, cex = 2)
  mtext("Frequency", outer = TRUE, line = 2, side = 2, cex = 2)
  dev.off()

  # Scatterplot of positive catch by depth
  Ncol <- Nrow <- 1
  jpeg(filename = paste(Folder, FileName, "Positive catch and depth.jpg", sep = ""), width = Ncol * 4, height = Nrow * 4, units = "in", res = 200)
  par(mar = c(2, 2, 2, 0), mgp = c(1.25, 0.25, 0), tck = -0.02)
  Y <- Pos[, "HAUL_WT_KG"]
  plot(x = Pos[, "BEST_DEPTH_M"], y = ifelse(Y == 0, NA, Y), log = "y", main = "Positive catch by depth", xlab = "Depth", ylab = "Positive catch rates", pch = 20, col = rgb(0, 0, 0, alpha = 0.2))
  lines(lowess(x = Pos[, "BEST_DEPTH_M"], y = ifelse(Y == 0, NA, Y)), lwd = 2)
  dev.off()

  # Scatterplot of positive catch by depth | Year
  Ncol <- ceiling(sqrt(Nyears))
  Nrow <- ceiling(Nyears / Ncol)
  jpeg(filename = paste(Folder, FileName, "Positive catch and depth BY Year.jpg", sep = ""), width = Ncol * 3, height = Nrow * 3, units = "in", res = 200)
  par(mfrow = c(Nrow, Ncol), mar = c(2, 2, 2, 0), mgp = c(1.25, 0.25, 0), tck = -0.02, oma = c(4, 4, 0, 0))
  for (YearI in 1:Nyears) {
    Which <- which(Pos[, "PROJECT_CYCLE"] == unique(Pos[, "PROJECT_CYCLE"])[YearI])
    Y <- Pos[Which, "HAUL_WT_KG"]
    plot(x = Pos[Which, "BEST_DEPTH_M"], y = ifelse(Y == 0, NA, Y), log = "y", main = unique(Pos[, "PROJECT_CYCLE"])[YearI], xlab = "", ylab = "", pch = 20, col = rgb(0, 0, 0, alpha = 0.2))
    lines(lowess(x = Pos[Which, "BEST_DEPTH_M"], y = ifelse(Y == 0, NA, Y)), lwd = 2)
  }
  mtext("Depth", outer = TRUE, line = 2, side = 1, cex = 2)
  mtext("Positive catch rates", outer = TRUE, line = 2, side = 2, cex = 2)
  dev.off()
  # Line 405 from John Wallace's "Survey.Biomass.GlmmBUGS.ver.3.00.R"
  # xyplot(ifelse(HAUL_WT_KG==0,min(HAUL_WT_KG[HAUL_WT_KG!=0])/2,HAUL_WT_KG) ~ BEST_DEPTH_M | factor(PROJECT_CYCLE), data=Data, ylab = "Log of Weight (kg)", xlab="Depth (m)")

  # Scatterplot of positive catch by lattitude
  Ncol <- Nrow <- 1
  jpeg(filename = paste(Folder, FileName, "Positive catch and latitude.jpg", sep = ""), width = Ncol * 4, height = Nrow * 4, units = "in", res = 200)
  par(mar = c(2, 2, 2, 0), mgp = c(1.25, 0.25, 0), tck = -0.02)
  Y <- Pos[, "HAUL_WT_KG"]
  plot(x = Pos[, "BEST_LAT_DD"], y = ifelse(Y == 0, NA, Y), log = "y", main = "Positive catch by Latitude", xlab = "Latitude", ylab = "Positive catch rates", pch = 20, col = rgb(0, 0, 0, alpha = 0.2))
  lines(lowess(x = Pos[, "BEST_LAT_DD"], y = ifelse(Y == 0, NA, Y)), lwd = 2)
  dev.off()

  # Scatterplot of positive catch by lattitude | Year
  Ncol <- ceiling(sqrt(Nyears))
  Nrow <- ceiling(Nyears / Ncol)
  jpeg(filename = paste(Folder, FileName, "Positive catch and latitude BY Year.jpg", sep = ""), width = Ncol * 3, height = Nrow * 3, units = "in", res = 200)
  par(mfrow = c(Nrow, Ncol), mar = c(2, 2, 2, 0), mgp = c(1.25, 0.25, 0), tck = -0.02, oma = c(4, 4, 0, 0))
  for (YearI in 1:Nyears) {
    Which <- which(Pos[, "PROJECT_CYCLE"] == unique(Pos[, "PROJECT_CYCLE"])[YearI])
    Y <- Pos[Which, "HAUL_WT_KG"]
    plot(x = Pos[Which, "BEST_LAT_DD"], y = ifelse(Y == 0, NA, Y), log = "y", main = unique(Pos[, "PROJECT_CYCLE"])[YearI], xlab = "", ylab = "", pch = 20, col = rgb(0, 0, 0, alpha = 0.2))
    lines(lowess(x = Pos[Which, "BEST_LAT_DD"], y = ifelse(Y == 0, NA, Y)), lwd = 2)
  }
  mtext("Latitude", outer = TRUE, line = 2, side = 1, cex = 2)
  mtext("Positive catch rates", outer = TRUE, line = 2, side = 2, cex = 2)
  dev.off()

  # Bar graph of proportion positive by year
  Ncol <- Nrow <- 1
  jpeg(filename = paste(Folder, FileName, "Presence and year.jpg", sep = ""), width = Ncol * 4, height = Nrow * 4, units = "in", res = 200)
  par(mar = c(2, 2, 2, 0), mgp = c(1.25, 0.25, 0), tck = -0.02)
  SumPres <- tapply(Data[, "Pres"], INDEX = Data[, "PROJECT_CYCLE"], FUN = sum, na.rm = TRUE)
  SumAbs <- tapply(1 - Data[, "Pres"], INDEX = Data[, "PROJECT_CYCLE"], FUN = sum, na.rm = TRUE)
  Prop <- tapply(Data[, "Pres"], INDEX = Data[, "PROJECT_CYCLE"], FUN = mean, na.rm = TRUE)
  # plot(x=unique(Data[,'PROJECT_CYCLE']), y=Prop, main="Presence/absence by year", xlab="Year", ylab="Proportion positive", pch=20, type="l", ylim=c(0,1))
  barplot(Prop, , main = "Presence/absence by year", xlab = "Year", ylab = "Proportion positive", ylim = c(0, 1))
  # barplot(-SumAbs, ylim=c(0,max(SumPres+SumAbs)),add=TRUE)
  dev.off()

  # Bar graph of proportion positive by 25 meter depth bins | Year
  Ncol <- ceiling(sqrt(Nyears))
  Nrow <- ceiling(Nyears / Ncol)
  BinWidth <- 50
  jpeg(filename = paste(Folder, FileName, "Presence and depth BY year.jpg", sep = ""), width = Ncol * 3, height = Nrow * 3, units = "in", res = 200)
  par(mfrow = c(Nrow, Ncol), mar = c(2, 2, 2, 0), mgp = c(1.25, 0.25, 0), tck = -0.02, oma = c(4, 4, 0, 0))
  for (YearI in 1:Nyears) {
    Which <- which(Data[, "PROJECT_CYCLE"] == unique(Data[, "PROJECT_CYCLE"])[YearI])
    X <- unique(BinWidth * floor(Data[Which, "BEST_DEPTH_M"] / BinWidth))
    Order <- order(X)
    Prop <- tapply(Data[Which, "Pres"], INDEX = BinWidth * floor(Data[Which, "BEST_DEPTH_M"] / BinWidth), FUN = mean, na.rm = TRUE)[Order]
    # plot(x=X[Order], y=Prop, main=unique(Data[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", pch=20, type="l", ylim=c(0,1))
    barplot(Prop, main = unique(Data[, "PROJECT_CYCLE"])[YearI], xlab = "", ylab = "", ylim = c(0, 1))
  }
  mtext("Depth bin", outer = TRUE, line = 2, side = 1, cex = 2)
  mtext("Proportion positive", outer = TRUE, line = 2, side = 2, cex = 2)
  dev.off()

  # Bar graph of proportion positive by 1 degree latitude bins | Year
  Ncol <- ceiling(sqrt(Nyears))
  Nrow <- ceiling(Nyears / Ncol)
  BinWidth <- 1
  jpeg(filename = paste(Folder, FileName, "Presence and latitude BY year.jpg", sep = ""), width = Ncol * 3, height = Nrow * 3, units = "in", res = 200)
  par(mfrow = c(Nrow, Ncol), mar = c(2, 2, 2, 0), mgp = c(1.25, 0.25, 0), tck = -0.02, oma = c(4, 4, 0, 0))
  for (YearI in 1:Nyears) {
    Which <- which(Data[, "PROJECT_CYCLE"] == unique(Data[, "PROJECT_CYCLE"])[YearI])
    X <- unique(BinWidth * floor(Data[Which, "BEST_LAT_DD"] / BinWidth))
    Order <- order(X)
    Prop <- tapply(Data[Which, "Pres"], INDEX = BinWidth * floor(Data[Which, "BEST_LAT_DD"] / BinWidth), FUN = mean, na.rm = TRUE)[Order]
    # plot(x=X[Order], y=Prop, main=unique(Data[,'PROJECT_CYCLE'])[YearI], xlab="", ylab="", pch=20, type="l", ylim=c(0,1))
    barplot(Prop, main = unique(Data[, "PROJECT_CYCLE"])[YearI], xlab = "", ylab = "", ylim = c(0, 1))
  }
  mtext("Latitude bin", outer = TRUE, line = 2, side = 1, cex = 2)
  mtext("Proportion positive", outer = TRUE, line = 2, side = 2, cex = 2)
  dev.off()
}
