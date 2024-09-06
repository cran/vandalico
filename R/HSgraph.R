#' @title Suitability values distribution graph
#'
#' @description A function to visualize the distribution of the suitability
#' values associated to presences, absences, and all cases together.
#' @param mat A matrix with two columns. The first column must contain the
#' suitability values (i.e., the classification rule); the second column must
#' contain the presences and absences.
#' @param breaks Number of cells for the total histogram. By default,
#' \code{breaks} = 10.
#' @param hist.total Logical. Indicates whether or not the distribution of
#' suitability values for all the cases together is graphed.
#' @details In blue, the distribution of the suitability values associated to
#' presences. In red, the distribution of the suitability values associated to
#' absences. This graph helps to understand why the \emph{AUC} (or \emph{Se*})
#' is greater, equal to, or less than the \emph{uAUC} (or \emph{uSe*}) (see
#' Jiménez-Valverde 2022).
#' @return This function returns a multiple histogram.
#' @examples
#' suit<-rbeta(100, 2, 2) #Generate suitability values
#' random<-runif(100)
#' sp<-ifelse(random < suit, 1 , 0) #Generate presence-absence data
#' HSgraph(cbind(suit, sp), breaks = 20, hist.total = TRUE)
#' @encoding UTF-8
#' @references Jiménez-Valverde, A. (2022). The uniform AUC: dealing with the
#'     representativeness effect in presence-absence models. \emph{Methods Ecol.
#'     Evol.}, 13, 1224-1236.
#' @importFrom graphics hist
#' @importFrom grDevices rgb
#' @importFrom graphics lines
#' @export

HSgraph <- function(mat, breaks = 10, hist.total = TRUE) {
    breaks <- seq(min(mat[, 1]), max(mat[, 1]), length.out = breaks + 1)
    histpresence <- hist(mat[, 1][mat[, 2] == 1], breaks = breaks)
    histabsence <- hist(mat[, 1][mat[, 2] == 0], breaks = breaks)
    histtotal <- hist(mat[, 1], breaks = breaks)
    plot(histpresence, col = rgb(0, 0, 1, 1 / 4), xlim = c(min(mat[, 1]),
        max(mat[, 1])), ylim = c(0, max(histpresence$counts, histabsence$counts,
        histtotal$counts)), border = rgb(0, 0, 1, 0), main = "",
        xlab = "Suitability", ylab = "Frequency", las = 1)
    plot(histabsence, col = rgb(1, 0, 0, 1 / 4), xlim = c(min(mat[, 1]),
        max(mat[, 1])), border = rgb(1, 0, 0, 0), add = T)
    if (hist.total == TRUE) {
        lines(c(0, histtotal$breaks), c(0, histtotal$counts, 0), type = "s",
        lwd = 2)
    }
}
