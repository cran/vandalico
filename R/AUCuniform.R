#' @title Calculates the uniform \emph{AUC} and uniform \emph{Se*}
#'
#' @description This function computes the uniform \emph{AUC} (\emph{uAUC}) and
#' uniform \emph{Se*} (\emph{uSe*}) following Jiménez-Valverde (2022).
#' @param mat A matrix with two columns. The first column must contain the
#' suitability values (i.e., the classification rule); the second column must
#' contain the presences and absences.
#' @param rep Number of sampling replications. By default, \code{rep} = 100.
#' @param by Size of the suitability intervals (i.e., bins). By default,
#' \code{by} = 0.1.
#' @param deleteBins A vector (e.g., from 1 to 10 if \code{by} = 0.1) with the
#' bins that have to be excluded (1 for [0,0.1), 10 for [0.9,1]) from the
#' resampling procedure (trimming); \code{NULL} by default.
#' @param plot Logical. Indicates whether or not the observed ROC curve is
#' plotted.
#' @param plot.adds Logical. Indicates whether or not the negative diagonal and
#' the point of equivalence are added to the observed ROC plot.
#' @details This function performs the stratified weighted bootstrap to
#' calculate the uniform \emph{AUC} (\emph{uAUC}) and uniform \emph{Se*}
#' (\emph{uSe*}) as suggested in Jiménez-Valverde (2022).  A warning message
#' will be shown if the sample size of any bin is zero. Another warning message
#' will be shown if the sample size of any bin is lower than 15. In such case,
#' trimming should be considered. The \emph{AUC} (non-uniform) is estimated
#' non-parametrically (Bamber 1975). \emph{Se*} is calculated by selecting the
#' point that minimizes the absolute difference between sensitivity and
#' specificity and by doing the mean of those values (Jiménez-Valverde 2020).
#' @return A list with the following elements:
#' @return \code{AUC}: the \emph{AUC} value (non-uniform), a numeric value 
#' between 0 and 1.
#' @return \code{Se}: the \emph{Se*} value (non-uniform), a numeric value 
#' between 0 and 1. 
#' @return \code{bins}: a table with the sample size of each bin.
#' @return \code{suit.sim}: a matrix with the bootstrapped suitability values.
#' @return \code{sp.sim}: a matrix with the bootstrapped presence-absence data.
#' @return \code{uAUC}: a numeric vector with the (\emph{uAUC}) values for each
#' replication.
#' @return \code{uAUC.95CI}: a numeric vector with the sample (\emph{uAUC}) 
#' quantiles corresponding to the probabilities 0.025, 0.5 and 0.975.
#' @return \code{uSe}: a numeric vector with the (\emph{uSe*}) values for each
#' replication.
#' @return \code{uSe.95CI}: a numeric vector with the sample (\emph{uSe*}) 
#' quantiles corresponding to the probabilities 0.025, 0.5 and 0.975.
#' @examples
#' suit<-rbeta(100, 2, 2) #Generate suitability values
#' random<-runif(100)
#' sp<-ifelse(random < suit, 1, 0) #Generate presence-absence data
#' result<-AUCuniform(cbind(suit, sp), plot = TRUE, plot.adds = TRUE)
#' result$uAUC.95CI[2] #Get the uAUC
#' @encoding UTF-8
#' @references Bamber, D. (1975). The Area above the Ordinal Dominance Graph and
#'     the Area below the Receiver Operating Characteristic Graph.
#'     \emph{J. Math. Psychol}., 12, 387-415.
#'
#'     Jiménez-Valverde, A. (2020). Sample size for the evaluation of
#'     presence-absence models. \emph{Ecol. Indic}., 114, 106289.
#'
#'     Jiménez-Valverde, A. (2022). The uniform AUC: dealing with the
#'     representativeness effect in presence-absence models. \emph{Methods Ecol.
#'     Evol.}, 13, 1224-1236.
#' @importFrom stats quantile wilcox.test
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics points
#' @export

AUCuniform <- function(mat, rep = 100, by = 0.1, deleteBins = NULL, plot = FALSE, plot.adds = FALSE) {
    # Non-uniform:
    Wilcox <- wilcox.test(mat[, 1] ~ mat[, 2], exact = FALSE)
    AUC <- 1 - (Wilcox["statistic"]$statistic/((length(mat[, 2]) - (sum(mat[, 2]))) * sum(mat[, 
        2])))
    cosa <- ROCR::prediction(mat[, 1], mat[, 2])
    temporal <- ROCR::performance(cosa, "sens", "spec")
    diferencia <- abs(temporal@x.values[[1]] - temporal@y.values[[1]])
    SE <- (temporal@y.values[[1]][which.min(diferencia)] + temporal@x.values[[1]][which.min(diferencia)])/2
    if (plot == TRUE) {
        plot(1 - temporal@x.values[[1]], temporal@y.values[[1]], pch = 16, xlab = "false positive rate", 
            ylab = "sensitivity", main = "ROC curve", yaxt = "n", cex.lab = 1.3, cex.axis = 1)
        axis(side = 2, las = 2, mgp = c(3, 0.75, 0))
        abline(a = 0, b = 1, lty = 2)
        if (plot.adds == TRUE) {
            abline(a = 1, b = -1, col = "darkgrey", lty = 2)
            points(1 - SE, SE, col = "red", pch = 16)
        }
    }
    # Uniform:
    bins <- seq(0, 1, by)
    intervals <- cut(mat[, 1], bins, include.lowest = T, right = F)
    probs <- table(intervals)[intervals]
    if (dim(mat)[1] < 30) 
        warning("Your sample size is low, results must be interpreted with caution.")
    if (sum(table(intervals) < 15) > 0) 
        warning("At least one sutability interval has n < 15, results must be interpreted with caution.")
    if (sum(table(intervals) == 0) > 0) 
        warning(paste("There are", sum(table(intervals) == 0), "interval(s) with zero data, results must be interpreted with caution."))
    if (is.null(deleteBins) == FALSE) {
        toDelete <- levels(intervals)[deleteBins]
        mat <- mat[(intervals %in% toDelete) == FALSE, ]
        probs <- probs[(intervals %in% toDelete) == FALSE]
    }
    uAUC <- c()
    uSE <- c()
    HS <- matrix(nrow = nrow(mat), ncol = rep)
    SP <- matrix(nrow = nrow(mat), ncol = rep)
    for (i in 1:rep) {
        newdata <- mat[sample(nrow(mat), nrow(mat), replace = T, prob = 1/probs), ]
        HS[, i] <- newdata[, 1]
        SP[, i] <- newdata[, 2]
        Wilcox <- wilcox.test(newdata[, 1] ~ newdata[, 2], exact = FALSE)
        uAUC <- c(uAUC, 1 - (Wilcox["statistic"]$statistic/((dim(newdata)[1] - (sum(newdata[, 
            2]))) * sum(newdata[, 2])))[[1]])
        cosa <- ROCR::prediction(newdata[, 1], newdata[, 2])
        temporal <- ROCR::performance(cosa, "sens", "spec")
        diferencia <- abs(temporal@x.values[[1]] - temporal@y.values[[1]])
        uSE <- c(uSE, (temporal@y.values[[1]][which.min(diferencia)] + temporal@x.values[[1]][which.min(diferencia)])/2)
    }
    
    return(list(AUC = AUC[[1]], Se = SE, bins = table(intervals), suit.sim = HS, sp.sim = SP, 
        uAUC = uAUC, uAUC.95CI = quantile(uAUC, c(0.025, 0.5, 0.975)), uSe = uSE, uSe.95CI = quantile(uSE, 
            c(0.025, 0.5, 0.975))))
}
