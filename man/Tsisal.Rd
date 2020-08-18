\name{Tsisal}
\alias{Tsisal}
\title{Complete Deconvolution of DNA methylation data based on TOAST and SISAL}
\description{A function to conduct complete reference-free deconvolution on
DNA methylation data. If a full reference or a partial reference panel is
provided, this function also automatically annotate the solved proportions to
known cell types.
}
\usage{
Tsisal(Y_raw, K = NULL, knowRef = NULL, possibleCellNumber = 3:15)
}
\arguments{
   \item{Y_raw}{The DNA methylation 450K array data from complex tissues, rows for CpG sites and columns for samples.}
    \item{K}{The number of pure cell types, we allow users to pre-specify or use our method to estimate.}
    \item{knowRef}{The external reference panel for cell type label assignment.}
    \item{possibleCellNumber}{Range of possible number of cell types. Default is 3:15.}
}
\value{
  \item{estProp}{Estimated proportions.}
  \item{selMarker}{Selected cell type-specific markers.}
  \item{K}{Optional number of cell types.}
}
\references{
Complete deconvolution of DNA methylation signals from complex tissues: a geometric approach.
Weiwei Zhang, Hao Wu and Ziyi Li.
}
\author{
Weiwei Zhang <wwzhangly@163.com>
}
\examples{
### generate a simulation data
knowRef <- matrix(runif(5000*5), 5000, 5)
colnames(knowRef) <- paste0("CellType", 1:5)
Y_raw <- matrix(runif(5000*20), 5000, 20)
rownames(Y_raw) <- paste0("CpG", 1:5000)
colnames(Y_raw) <- paste0("Sample", 1:20)

Tsisal(Y_raw = Y_raw, K = 5, knowRef = knowRef)

## if cell type number is unknown
# Tsisal(Y.raw = Y_raw, K = NULL, knowRef = knowRef, possibleCellNumber = 4:10)
}
\keyword{ compelte deconvolution }% use one of  RShowDoc("KEYWORDS")
