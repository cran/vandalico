#' @title Calibration graph
#'
#' @description A function to plot a calibration graph.
#' @param mat A matrix with two columns. The first column must contain the
#' suitability values (i.e., the classification rule); the second column must
#' contain the presences and absences.
#' @param by Size of the suitability intervals (bins). By default,
#' \code{by} = 0.1.
#' @details Dots for bins with 15 or more cases are shown in solid black; dots
#' for bins with less than 15 cases are shown empty (see Jiménez-Valverde et
#' al. 2013). This way, by plotting the calibration graph before running
#' \code{\link{AUCuniform}}, one can get a glimpse of how reliable \emph{uAUC}
#' or \emph{uSe*} can be expected to be.
#' @return This function returns a calibration plot
#' @examples
#' suit<-rbeta(100, 2, 2) #Generate suitability values
#' random<-runif(100)
#' sp<-ifelse(random < suit,1 , 0) #Generate presence-absence data
#' CALplot(cbind(suit, sp))
#' @encoding UTF-8
#' @references Jiménez-Valverde, A., Acevedo, P., Barbosa, A. M., Lobo, J. M. &
#'     Real, R. (2013).  Discrimination capacity in species distribution models
#'     depends on the representativeness of the environmental domain.
#'     \emph{Global Ecol. Biogeogr.}, 22, 508-516.
#' @export

CALplot <- function(mat, by = 0.1) {
    bins <- cut(mat[, 1], seq(0, 1, by), include.lowest = T, right = F)
    tableBins <- table(mat[, 2], bins)
    prevalBins <- as.matrix(tableBins[2, ]/colSums(tableBins))
    prevalBins.2 <- as.matrix(tapply(mat[, 1], bins, mean))
    colDots <- ifelse(colSums(tableBins) < 14, 21, 19)
    plot(prevalBins.2, prevalBins, main = "Calibration plot", pch = colDots, ylab = "observed probability", xlab = "predicted probability", 
        xlim = c(0, 1), ylim = c(0, 1), yaxt = "n", cex.lab = 1.3, cex = 1.1)
    axis(side = 2, las = 2, mgp = c(3, 0.75, 0))
    abline(a = 0, b = 1, col = "black", lty = 2)
}
