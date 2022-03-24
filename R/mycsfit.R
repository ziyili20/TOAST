#' csfit: Deconvolution from Known Cell Proportions
#'
#' Deconvolves cell-specific expression using least-squares fit. Input is the
#' heterogeneous sample gene expression of a group of samples and the matching
#' cell-frequencies of the sample. The lower limit for the number of samples
#' needed to deconvolving the cell-specific expression of N cell-types is N+1.
#' For a single color array - the result could be interpreted as the average
#' expression level of a given gene in a cell-type of that group. Multiplied by
#' the frequency of a given cell-type in an individual in the group, it is the
#' amount contributed by that cell type to the overall measured expression on
#' the array.
#'
#'
#' @param G Matrix of gene expression, columns ordered in the same order at the
#' cell-frequency matrix (n by g, n samples, g genes)
#' @param cc Matrix of cell-frequency. (n by k, n samples, k cell-types)
#' @param logRm Exponentiate data for deconvolution stage. Default is FALSE
#' @param logBase Base of logarithm used to determine exponentiation factor.
#' Default is 2
#' @return A list with three attributes:
#' \item{ghat}{A matrix of cell-specific expression for each gene as
#' derived from the coefficients of the fit. (Size: k by g, k cell types, gp
#' genes)}
#' \item{se}{Standard error of the fit coefficients}
#' \item{residuals}{The individual sample residuals.}
#' \item{n}{The number of samples used in the fit.}
#'
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
mycsfit <- function(cc,G,logRm=FALSE,logBase=2) {
     if(logRm == TRUE) {
          G <- logBase^G
     }

     fit1 <- lsfit(cc, G, intercept = FALSE)
     se1 <- ls.diag(fit1)$std.err
     if(logRm == TRUE) {
          ghat <- log(coef(fit1),logBase)
          ghat[is.nan(ghat)] <- 0
          se <- log(se1,logBase)
          res <- list(ghat = ghat, residuals = residuals(fit1),se = se1, n = nrow(G))
     }
     else {
          res <- list(ghat = coef(fit1), residuals = residuals(fit1),se= se1, n = nrow(G))
     }

     # fix dimnames
     dimnames(res$residuals) <- dimnames(G)
     features <- colnames(G)
     colnames(res$ghat) <- features
     colnames(res$se) <- features

     # return
     class(res) <- 'csfit'
     res
}

#' @rawNamespace S3method(coef, csfit)
coef.csfit <- function(object, ...){
     object$ghat
}


