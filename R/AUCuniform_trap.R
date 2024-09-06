#' @title Calculates the uniform \emph{AUC} and uniform \emph{Se*} by the
#' trapezoidal method.
#'
#' @description This function computes the uniform \emph{AUC} (\emph{uAUC}) and
#' uniform \emph{Se*} (\emph{uSe*}) using the weighted trapezoidal method
#' instead of the weighted bootstrapping method used in \code{AUCuniform} and
#' originally proposed in Jiménez-Valverde (2022). This procedure reduces bias 
#' and improves the coverage of confidence intervals (Jiménez-Valverde 2024). 
#' Additionally, the weights vector associated to each case can be customized. 
#' See Jiménez-Valverde (2024) for details.
#' @param mat A matrix with two columns. The first column must contain the
#' suitability values (i.e., the classification rule); the second column must
#' contain the presences and absences.
#' @param by Size of the suitability intervals (i.e., bins). By default,
#' \code{by} = 0.1.
#' @param deleteBins A vector (e.g., from 1 to 10 if \code{by} = 0.1) with the
#' bins that have to be excluded (1 for [0,0.1), 10 for [0.9,1]) from the
#' calculation of the uniform statistics; \code{NULL} by default.
#' @param w A vector with the weights associated to each case. If \code{NULL}
#' (default), then the uniform \emph{AUC} (\emph{uAUC}) and uniform \emph{Se*}
#' (\emph{uSe*}) are calculated.
#' @param plot Logical. Indicates whether or not the observed ROC curve is
#' plotted (gray dots).
#' @param plot.compare Logical. Indicates whether or not the weighed ROC curve
#' is plotted (black line).
#' @param plot.adds Logical. Indicates whether or not the negative diagonal and
#' the points of equivalence (weighted and unweighted) are added to the ROC
#' plot.
#' @details This function calculates the uniform \emph{AUC} (\emph{uAUC}) and
#' uniform \emph{Se*} (\emph{uSe*}) using the weighted trapezoidal method as
#' suggested in Jiménez-Valverde (2024). A warning message will be shown if
#' the sample size of any bin is zero. Another warning message will be shown if
#' the sample size of any bin is lower than 15. In such case, trimming should be
#' considered using \code{deleteBins} (Jiménez-Valverde 2022). Alternatively,
#' the weights associated to each case can be fully customized with the \code{w}
#' parameter (Jiménez-Valverde 2024). In this case, no warnings regarding
#' sample size issues are raised, and \code{deleteBins} is not used. The
#' \emph{AUC} (non-uniform, unweighted) is estimated non-parametrically by the
#' trapezoidal rule, which is equivalent to the Wilcoxon-based estimation 
#' (Hanley & McNeil 1982) used in \code{AUCuniform}. \emph{Se*} is calculated as 
#' in \code{AUCuniform}.
#' @return A list with the following elements:
#' @return \code{AUC}: the \emph{AUC} value (non-uniform, unweighted), a numeric
#' value between 0 and 1.
#' @return \code{Se}: the \emph{Se*} value (non-uniform, unweighted), a numeric
#' value between 0 and 1.
#' @return \code{bins}: a table with the sample size of each bin (only if
#' \code{w = NULL}).
#' @return \code{uAUC}: the uniform \emph{AUC} value (only if \code{w = NULL}).
#' @return \code{uSe}: the uniform \emph{Se*} value (only if \code{w = NULL}).
#' @return \code{wAUC}: the weighted \emph{AUC} estimated with the vector
#' \code{w}
#' @return \code{wSe}: the weighted \emph{Se*} estimated with the vector
#' \code{w}
#' @examples
#' suit<-rbeta(100, 2, 2) #Generate suitability values
#' random<-runif(100)
#' sp<-ifelse(random < suit, 1, 0) #Generate presence-absence data
#' result<-AUCuniform_trap(cbind(suit, sp), plot = TRUE, plot.compare = TRUE)
#' result$AUC #Get the AUC
#' result$uAUC #Get the uAUC. Note how it is closer to the reference value of
#'             #0.83 since the suitability values are simulated to be
#'             #well-calibrated (see Jimenez-Valverde 2022).
#' @encoding UTF-8
#' @references Hanley, J. A. & McNeil, B. J. (1982). The Meaning and Use of the
#'     Area under a Receiver Operating Characteristic (ROC) Curve.
#'     \emph{Radiology}., 143, 29-36.
#'
#'     Jiménez-Valverde, A. (2022). The uniform AUC: dealing with the
#'     representativeness effect in presence-absence models. \emph{Methods Ecol.
#'     Evol.}, 13, 1224-1236.
#'
#'     Jiménez-Valverde, A. (2024). Improving the uniform AUC (uAUC): towards a 
#'     case-by-case weighting evaluation in species distribution models. 
#'     \emph{In preparation}. 
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics points
#' @export

AUCuniform_trap <- function(mat, by = 0.1, deleteBins = NULL, w = NULL,
    plot = FALSE, plot.compare = FALSE, plot.adds = FALSE) {
    # Non-uniform:
    mat.ord <- mat[order(mat[, 1], decreasing = TRUE), ]
    true.pos <- cumsum(mat.ord[, 2]) / sum(mat.ord[, 2])
    false.pos <- cumsum(1 - mat.ord[, 2]) / sum((1 - mat.ord[, 2]))
    auc.trap <- sum(diff(c(0, false.pos)) * true.pos)
    diferencia <- abs(true.pos - (1 - false.pos))
    se.trap <- (true.pos[which.min(diferencia)] + (1 - false.pos)[which.min(diferencia)]) / 2
    if (plot == TRUE) {
        plot(c(0,false.pos), c(0,true.pos), pch = 16, xlab = "false positive rate",
            ylab = "sensitivity", main = "ROC curve", yaxt = "n", cex.lab = 1.3,
            cex.axis = 1, col="gray")
        axis(side = 2, las = 2, mgp = c(3, 0.75, 0))
        abline(a = 0, b = 1, lty = 2)
    }
    # Uniform:
    if (is.null(w) == TRUE) {
        bins <- seq(0, 1, by)
        intervals <- cut(mat[, 1], bins, include.lowest = TRUE, right = FALSE)
        probs <- as.vector(1 / table(intervals)[intervals])
        probs.ord <- probs[order(mat[, 1], decreasing = TRUE)]
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
            mat.ord <- mat[order(mat[, 1], decreasing = TRUE), ]
            probs.ord <- probs[order(mat[, 1], decreasing = TRUE)]
        }
        wauc <- sum(diff(c(0, cumsum((1 - mat.ord[, 2]) * probs.ord) / sum((1 - mat.ord[, 2]) *
            probs.ord))) * (cumsum(mat.ord[, 2] * probs.ord) / sum(mat.ord[, 2] * probs.ord)))
        diferencia <- abs((cumsum(mat.ord[, 2] * probs.ord) / sum(mat.ord[, 2] * probs.ord)) -
            (1 - (cumsum((1 - mat.ord[, 2]) * probs.ord) / sum((1 - mat.ord[, 2]) * probs.ord))))
        wse <- ((cumsum(mat.ord[, 2] * probs.ord) / sum(mat.ord[, 2] * probs.ord))[which.min(diferencia)] +
            (1 - (cumsum((1 - mat.ord[, 2]) * probs.ord) / sum((1 - mat.ord[, 2]) * probs.ord))))[which.min(diferencia)] / 2
    }
    # Other weights for cases
    if (is.null(w) == FALSE) {
        if (length(w) != dim(mat)[1])
            warning("The number of cases does not match the length of w. The weighted statistics are not calculated and the weighted ROC curve is not plotted.")
        probs.ord <- w[order(mat[, 1], decreasing = TRUE)]
        wauc <- sum(diff(c(0, cumsum((1 - mat.ord[, 2]) * probs.ord) / sum((1 - mat.ord[, 2]) *
            probs.ord))) * (cumsum(mat.ord[, 2] * probs.ord) / sum(mat.ord[, 2] * probs.ord)))
        wse <- ((cumsum(mat.ord[, 2] * probs.ord) / sum(mat.ord[, 2] * probs.ord))[which.min(diferencia)] +
            (1 - (cumsum((1 - mat.ord[, 2]) * probs.ord) / sum((1 - mat.ord[, 2]) * probs.ord))))[which.min(diferencia)] / 2
    }
    # Plots
    if (plot.compare == TRUE) {
        lines(c(0,(cumsum((1 - mat.ord[, 2]) * probs.ord) / sum((1 - mat.ord[, 2]) * probs.ord))),
            c(0,(cumsum(mat.ord[, 2] * probs.ord) / sum(mat.ord[, 2] * probs.ord))),
            pch = 16, col = "black")
    }
    if (plot.adds == TRUE) {
        abline(a = 1, b = -1, col = "darkgrey", lty = 2)
        points(1 - se.trap, se.trap, col = "red", pch = 16)
        points(1 - wse, wse, col = "red", pch = 16)
    }
    # Get results
    if (is.null(w) == TRUE) {
        return(list(AUC = auc.trap, Se = se.trap, bins = table(intervals),
            uAUC = wauc, uSe = wse))
    }
    if (is.null(w) == FALSE) {
        return(list(AUC = auc.trap, Se = se.trap, wAUC = wauc, wSe = wse))
    }
}
