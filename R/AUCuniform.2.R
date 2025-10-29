#' @title Calculates the uniform \emph{AUC} and uniform \emph{Se*} by the
#' direct weighted trapezoidal estimation method.
#'
#' @description This function computes the uniform \emph{AUC} (\emph{uAUC}) and
#' uniform \emph{Se*} (\emph{uSe*}) using the direct weighted trapezoidal
#' estimation method (Jiménez-Valverde 2025), instead of the stratified 
#' bootstrapping with inverse probability weighting method implemented in 
#' \code{AUCuniform} and originally proposed by Jiménez-Valverde (2022). Uniform 
#' statistics are design to account for the representativeness effect 
#' (Jiménez-Valverde 2022). This new method reduces bias and improves the 
#' coverage of confidence intervals relative to the original proposal. 
#' Additionally, the weight vector associated to each case can be customized. 
#' @param mat A matrix with two columns. The first column must contain the
#' the classification rule (e.g., the suitability values); the second column 
#' must contain the presences and absences.
#' @param by The size of the intervals used to divide the classification rule 
#' (i.e., bins width). By default, \code{by} = 0.1. This argument is only used 
#' when \code{w = NULL}.
#' @param deleteBins A vector (e.g., from 1 to 10 if \code{by} = 0.1) with the
#' bins that have to be excluded (1 for [0,0.1), 10 for [0.9,1]) from the
#' calculation of the uniform statistics. The default is \code{NULL}. This 
#' argument is only used when \code{w = NULL}. 
#' @param w A vector with the weights associated with each case. If \code{NULL}
#' (default), each case is weighted by the inverse of the sample size of 
#' its corresponding bin, and the uniform \emph{AUC} (\emph{uAUC}) and uniform 
#' \emph{Se*} (\emph{uSe*}) are calculated (see Details).
#' @param plot Logical. If \code{TRUE}, the standard (unweighted) ROC curve is 
#' plotted (gray dots).
#' @param plot.compare Logical. If \code{TRUE}, the weighed ROC curve is plotted
#' (black line).
#' @param plot.adds Logical. If \code{TRUE}, adds the negative diagonal and the
#' points of equivalence (weighted and unweighted) to the ROC plot.
#' @details This function calculates the uniform \emph{AUC} (\emph{uAUC}) and
#' uniform \emph{Se*} (\emph{uSe*}) using the direct weighted trapezoidal 
#' estimation method proposed by Jiménez-Valverde (2025). To compute the uniform
#' statistics, the \emph{w} parameter must be set to \emph{NULL} (default). The 
#' data set is divided into bins (defined by the parameter \code{by}) based on 
#' the values of the first vector in the input matrix \code{mat} (the 
#' classification rule). Each observation is assigned a weight equal to one 
#' divided by the number of observations in the corresponding bin. Then, the 
#' uniform discrimination statistics are calculated via the direct weighted 
#' trapezoidal estimation method such that, for each threshold, the weighted 
#' true positive and false positive rates are cumulatively updated by summing 
#' the weights of the presences and absences, respectively, with that score 
#' (Jiménez-Valverde 2025). The calculation of the uniform statistics requires
#' the classification rule (\code{mat[,1]}) to range between 0 and 1, and the 
#' value of \code{by} to divide 1 exactly. If any of this conditions are not 
#' match, the function stops. A warning message is displayed if (1) the sample 
#' size is lower than 30, (2) any bin has a sample size of zero, or (3) any bin 
#' has a sample size between 1 and 15. In the latter case, trimming should be 
#' considered using \code{deleteBins}, in which case the uniform statistics are 
#' computed excluding the selected bins. See Jiménez-Valverde (2022) for further 
#' details. 
#' 
#' Alternatively, users may wish to downweight the importance of certain 
#' observations relative to others for reasons unrelated to the 
#' representativeness effect (Jiménez-Valverde 2025). For this purpose, the 
#' weights associated to each case can be fully customized with the \code{w} 
#' parameter (see Examples). The length of the weight vector has to be equal to 
#' \code{dim(mat)[1]}.       
#' 
#' The standard \emph{AUC} (non-uniform, unweighted) is estimated 
#' non-parametrically by the trapezoidal rule, which is equivalent to the 
#' Wilcoxon-based estimation (Hanley & McNeil 1982) used in \code{AUCuniform}. 
#' \emph{Se*} is calculated as in \code{AUCuniform}.
#' @return A list with the following elements:
#' @return \code{AUC}: the standard \emph{AUC} value (unweighted), a numeric 
#' value between 0 and 1.
#' @return \code{Se}: the standard \emph{Se*} value (unweighted), a numeric 
#' value between 0 and 1.
#' @return \code{bins}: a table with the sample size of each bin (returned only 
#' if \code{w = NULL}).
#' @return \code{uAUC}: the uniform \emph{AUC} value (returned only if 
#' \code{w = NULL}).
#' @return \code{uSe}: the uniform \emph{Se*} value (returned only if 
#' \code{w = NULL}).
#' @return \code{wAUC}: the weighted \emph{AUC} estimated with the vector
#' \code{w} (returned only if \code{w} is not \code{NULL}).
#' @return \code{wSe}: the weighted \emph{Se*} estimated with the vector
#' \code{w} (returned only if \code{w} is not \code{NULL}).
#' @return \code{TP}: a vector with the true positive rate for every threshold.
#' @return \code{FP}: a vector with the false positive rate for every threshold.
#' @return \code{TP.W}: a vector with the weighted true positive rate for every
#' threshold.
#' @return \code{FP.W}: a vector with the weighted false positive rate for every
#' threshold.
#' @examples
#' # In this first example, a data set is simulated in such a way that the
#' # classification rule is well-calibrated, i.e., the observed proportion of 
#' # positive cases equates to the simulated probabilities of presence. Since
#' # the objective is to calculate the uAUC to account for the environmental 
#' # representativeness effect (see Jiménez-Valverde 2022), weights are 
#' # automatically calculated and no w vector is needed. 
#' 
#' n <- 1000 # Set the sample size
#' hs <- rbeta(n, 2, 2) # Simulated probabilities (the classification rule)
#' random <- runif(n)
#' sp <- ifelse(random < hs, 1, 0) # Observed presence–absence data
#' 
#' result <- AUCuniform.2(cbind(hs, sp), plot = TRUE, plot.compare = TRUE)
#' 
#' result$AUC  # Get the standard AUC
#' result$uAUC # Get the uniform AUC. Note how it is close to the reference value
#'             # of 0.83 since the probability values  (the classification rule)
#'             # are simulated to be well-calibrated (see Jiménez-Valverde 2022)
#' 
#' # In this second set of examples, the objective is not to calculate the  
#' # uniform AUC, but to assign specific weights to certain observations. These 
#' # examples corresponds to some of those provided in Table 1 of 
#' # Jiménez-Valverde (2025).
#' 
#' hs <- seq(1, 0.05, by = -0.05) # Generate the classification rule 
#' sp <- c(0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0) # Observed presence–absence data
#' 
#' wa <- ifelse(sp == 0, 0.2, 1) # The vector of weights for each case
#' result.a <- AUCuniform.2(cbind(hs, sp), w = wa, plot = TRUE, plot.compare = TRUE)
#' 
#' result.a$AUC  # Get the standard AUC
#' result.a$wAUC # Get the weighted AUC. Since every case within each category of
#'               # sp received the same weight, the weighted AUC value equals the
#'               # standard AUC value
#'
#' wb <- c(rep(1, 19), 0.2) # The vector of weights for each case
#' result.b <- AUCuniform.2(cbind(hs, sp), w = wb, plot = TRUE, plot.compare = TRUE)
#' 
#' result.b$wAUC # Get the weighted AUC. Since a low weight is assigned to an 
#'               # instance of absence associated with a low probability value, 
#'               # the weighted AUC is lower than the standard AUC value.
#'             
#' wc <- c(0.2, rep(1, 19)) # The vector of weights for each case
#' result.c <- AUCuniform.2(cbind(hs, sp), w = wc, plot = TRUE, plot.compare = TRUE)
#' 
#' result.c$wAUC # Get the weighted AUC. Since a low weight is assigned to an 
#'               # instance of absence associated with a high probability value, 
#'               # the weighted AUC is higher than the standard AUC value
#' @encoding UTF-8
#' @references Hanley, J. A. & McNeil, B. J. (1982). The Meaning and Use of the
#'     Area under a Receiver Operating Characteristic (ROC) Curve.
#'     \emph{Radiology}., 143, 29-36.
#'
#'     Jiménez-Valverde, A. (2022). The uniform AUC: dealing with the
#'     representativeness effect in presence-absence models. \emph{Methods Ecol.
#'     Evol.}, 13, 1224-1236.
#'
#'     Jiménez-Valverde, A. (2025). Refining uniform discrimination metrics: 
#'     towards a case-by-case weighting evaluation in species distribution 
#'     models with presence-absence data. \emph{Under review}. 
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics points
#' @importFrom utils head tail
#' @export

AUCuniform.2 <- function(mat, by = 0.1, deleteBins = NULL, w = NULL, plot = FALSE, plot.compare = FALSE, plot.adds = FALSE) {

    if (!plot) {
        plot.compare <- FALSE
        plot.adds <- FALSE
    } else if (!plot.compare) {
        plot.adds <- FALSE
    }

    mat.ord <- mat[order(mat[, 1], decreasing = TRUE), ]

    unique.th <- sort(unique(mat.ord[, 1]), decreasing = TRUE)
    true.pos <- c()
    false.pos <- c()
    n.pos <- sum(mat.ord[, 2])
    n.neg <- length(mat.ord[, 2]) - n.pos
    if (n.pos == 0 || n.neg == 0) {
        stop("The dataset must contain both positive and negative cases.")
    }
    tp <- 0
    fp <- 0
    for (th in unique.th) {
        group <- mat.ord[mat.ord[, 1] == th, , drop = FALSE]
        tp <- tp + sum(group[, 2] == 1)
        fp <- fp + sum(group[, 2] == 0)
        true.pos <- c(true.pos, tp/n.pos)
        false.pos <- c(false.pos, fp/n.neg)
    }
    auc.trap <- sum(diff(false.pos) * (head(true.pos, -1) + tail(true.pos, -1))/2)
    diferencia <- abs(true.pos - (1 - false.pos))
    se.trap <- (true.pos[which.min(diferencia)] + (1 - false.pos)[which.min(diferencia)])/2

    if (is.null(w) == FALSE) {
        if (length(w) != dim(mat)[1])
            stop("The number of cases does not match the length of w.")
        probs.ord <- w[order(mat[, 1], decreasing = TRUE)]
        mat.ord <- cbind(mat.ord, probs.ord)
        unique.th <- sort(unique(mat.ord[, 1]), decreasing = TRUE)
        true.pos.w <- c()
        false.pos.w <- c()
        n.pos.w <- sum(mat.ord[mat.ord[, 2] == 1, 3])
        n.neg.w <- sum(mat.ord[mat.ord[, 2] == 0, 3])
        tp.w <- 0
        fp.w <- 0
        for (th in unique.th) {
            group <- mat.ord[mat.ord[, 1] == th, , drop = FALSE]
            tp.w <- tp.w + sum(group[group[, 2] == 1, 3])
            fp.w <- fp.w + sum(group[group[, 2] == 0, 3])
            true.pos.w <- c(true.pos.w, tp.w/n.pos.w)
            false.pos.w <- c(false.pos.w, fp.w/n.neg.w)
        }
        w.auc <- sum(diff(false.pos.w) * (head(true.pos.w, -1) + tail(true.pos.w, -1))/2)
        diferencia.w <- abs(true.pos.w - (1 - false.pos.w))
        w.se <- (true.pos.w[which.min(diferencia.w)] + (1 - false.pos.w)[which.min(diferencia.w)])/2
    }

    if (is.null(w) == TRUE) {
        if (any(mat[, 1] < 0 | mat[, 1] > 1)) {
            stop("The classification rule is out of the range [0,1]; consider rescaling (min-max normalization) before proceeding.")
        }
        if (abs((1/by) - round(1/by)) > .Machine$double.eps^0.5) {
            stop("'by' must divide 1 exactly (e.g., 0.2, 0.1, 0.05).")
        }
        if (dim(mat)[1] < 30)
            warning("Your sample size is low, results must be interpreted with caution.")
        bins <- seq(0, 1, by)
        intervals <- cut(mat[, 1], bins, include.lowest = TRUE, right = FALSE)
        if (sum(table(intervals) == 0) > 0)
            warning(paste("There are (is)", sum(table(intervals) == 0), "interval(s) with zero data, results must be interpreted with caution."))
        if (sum(table(intervals) < 15) > 0)
            warning(paste("There are (is)", sum(table(intervals) < 15) - sum(table(intervals) == 0), "interval(s) with 0<n<15 data, results must be interpreted with caution."))
        probs <- as.vector(1/table(intervals)[intervals])
        probs.ord <- probs[order(mat[, 1], decreasing = TRUE)]
        if (is.null(deleteBins) == FALSE) {
            toDelete <- levels(intervals)[deleteBins]
            mat <- mat[(intervals %in% toDelete) == FALSE, ]
            probs <- probs[(intervals %in% toDelete) == FALSE]
            mat.ord <- mat[order(mat[, 1], decreasing = TRUE), ]
            probs.ord <- probs[order(mat[, 1], decreasing = TRUE)]
        }
        mat.ord <- cbind(mat.ord, probs.ord)
        unique.th <- sort(unique(mat.ord[, 1]), decreasing = TRUE)
        true.pos.w <- c()
        false.pos.w <- c()
        n.pos.w <- sum(mat.ord[mat.ord[, 2] == 1, 3])
        n.neg.w <- sum(mat.ord[mat.ord[, 2] == 0, 3])
        tp.w <- 0
        fp.w <- 0
        for (th in unique.th) {
            group <- mat.ord[mat.ord[, 1] == th, , drop = FALSE]
            tp.w <- tp.w + sum(group[group[, 2] == 1, 3])
            fp.w <- fp.w + sum(group[group[, 2] == 0, 3])
            true.pos.w <- c(true.pos.w, tp.w/n.pos.w)
            false.pos.w <- c(false.pos.w, fp.w/n.neg.w)
        }
        w.auc <- sum(diff(false.pos.w) * (head(true.pos.w, -1) + tail(true.pos.w, -1))/2)
        diferencia.w <- abs(true.pos.w - (1 - false.pos.w))
        w.se <- (true.pos.w[which.min(diferencia.w)] + (1 - false.pos.w)[which.min(diferencia.w)])/2
    }

    if (plot == TRUE) {
        plot(c(0, false.pos), c(0, true.pos), pch = 16, xlab = "false positive rate", ylab = "sensitivity", main = "ROC curve",
            yaxt = "n", cex.lab = 1.3, cex.axis = 1, col = "gray")
        axis(side = 2, las = 2, mgp = c(3, 0.75, 0))
        abline(a = 0, b = 1, lty = 2)
    }
    if (plot.compare == TRUE) {
        lines(c(0, false.pos.w, 1), c(0, true.pos.w, 1), pch = 16, col = "black")
    }
    if (plot.adds == TRUE) {
        abline(a = 1, b = -1, col = "darkgrey", lty = 2)
        points(1 - se.trap, se.trap, col = "blue", pch = 16)
        points(1 - w.se, w.se, col = "red", pch = 16)
    }
    if (is.null(w) == TRUE) {
        return(list(AUC = auc.trap, Se = se.trap, uAUC = w.auc, uSe = w.se, bins = table(intervals), TP = true.pos,
            FP = false.pos, TP.W = true.pos.w, FP.W = false.pos.w))
    }
    if (is.null(w) == FALSE) {
        return(list(AUC = auc.trap, Se = se.trap, wAUC = w.auc, wSe = w.se, TP = true.pos, FP = false.pos, TP.W = true.pos.w,
            FP.W = false.pos.w))
    }
}
