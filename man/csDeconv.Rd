\name{csDeconv}
\alias{csDeconv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Improve reference-free deconvolution using
cross-cell type differential analysis
}
\description{
This function improve the feature selection in
reference-free deconvolution through cross-cell
type differential analysis
}
\usage{
csDeconv(Y_raw, K, FUN, nMarker = 1000,
InitMarker = NULL, TotalIter = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_raw}{
   A G*N matrix, G is the number of features, N is the number of subjects.
}
  \item{K}{
   The number of cell types. Need to be specified a priori.
}
  \item{FUN}{
   The reference-free deconvolution function,
   this function should take Y_raw and K, and
   the return values should be a N by K proportion matrix.
   N is the number of samples and K is the number of cell types.
   Default function is a wrapper of the RefFreeCellMix()
   function from CRAN package RefFreeEWAS.
}
  \item{nMarker}{
   The number of markers used in the deconvolution.
   Default is 1000.
}
  \item{InitMarker}{
   A vector of length L to represent the selection of inital markers.
   L should be equal or smaller than G.  If G is large,
   it is recommended that L is much smaller than G.
   If not specified, the most variable nMarker features will be used.}
  \item{TotalIter}{
   The total number of iterations of applying cross-cell
   type differential analysis. Default is 30.
}
}
\value{
  \item{allProp}{A list of estimated proportions
  from all iterations.}
  \item{allRMSE}{A vector of root mean squared errors
  (RMSE) from all iteratoins.}
  \item{estProp}{A N*K matrix representing the
  mixture proportions
  of K cell types in N subjects,
  chosen from allProp with the smallest RMSE.}
}
\references{
Ziyi Li and Hao Wu.
"Improving reference-free cell composition
estimation by cross-cell type differential analysis".
}
\author{
Ziyi Li <ziyi.li@emory.edu>
}
\examples{
Y_raw <- abs(matrix(runif(10000*20, 0,1),10000,20))
K <- 3

## wrap your reference-free
## deconvolution method into a function
## this function should take Y and K as input
## and output a N by K proprotion matrix
## here we use RefFreeCellMix() as an example
outT <- csDeconv(Y_raw, K)

RefFreeCellMix_wrapper <- function(Y, K){
     outY = RefFreeEWAS::RefFreeCellMix(Y,
     mu0=RefFreeEWAS::RefFreeCellMixInitialize(Y,
     K = K))
     Prop0 = outY$Omega
     return(Prop0)
}

outT <- csDeconv(Y_raw, K,
            FUN = RefFreeCellMix_wrapper)
}
