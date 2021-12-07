\name{MDeconv}
\alias{MDeconv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate cell compoisitons via partial reference-free deconvolution.
}
\description{
This function is the wrapper for TOAST/-P (partial reference-free deconvolution without prior) and TOAST/+P (with prior). It guides cell compoisition estimation through extra biological information, including cell type specific markers and prior knowledge of compoisitons.
}
\usage{
MDeconv(Ymat, SelMarker,
alpha = NULL, sigma = NULL,
epsilon = 0.001, maxIter = 1000,
verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ymat}{
   A gene expression data matrix from complex tissues. Rows represent genes and columns represent samples. Row names (e.g. gene symbol or ID) are needed. There is no need to filter data matrix by marker genes.
}
  \item{SelMarker}{
   A list variable with each element includes cell type-specific markers for this cell type. The marker list can be selected from pure cell type profiles or single cell data using ChooseMarker(). It can also be easily created manually. Please see details section below and example for more information.
}
  \item{alpha}{
   A vector including the prior mean for all cell types.
}
  \item{sigma}{
   A vector including the prior standard deviation for all cell types.
}
  \item{epsilon}{
   A numeric variable to control the level of convergence. With a large epsilon, it converges quickly and the algorithm may not converge well. With a small epsilon, it converges slower and the algorithm may converge "too much". The default value is 1e-3, which we find is a reasonable threshold.
}
  \item{maxIter}{
   Number of maximum iterations.
}
  \item{verbose}{
   A boolean variable of whether to output messages.
}
}
\details{
More about SelMarker: in addition to selecting markers using
ChooseMarker(), we can manually specific marker list. For example,
if we know there are two cell type-specific markers for cell type A
"Gene1" and "Gene2", and two cell type-specific markers for cell type B
"Gene3" and "Gene4", we can create marker list by
SelMarker = list(CellA = c("Gene1","Gene2"), CellB = c("Gene3","Gene4")).

One thing to note is that, the genes in marker list should have matches in
row names of Ymat. If not, the unmatched markers will be removed duing analysis.
}
\value{
A list including
\item{H}{Estimated proportion matrix, rows
for cell types and columns for samples.}
\item{W}{Estimated pure cell type profile matrix,
rows for genes and columns for cell types}
}
\references{
Ziyi Li, Zhenxing Guo, Ying Cheng, Peng Jin, Hao Wu.
"Robust partial reference-free cell compoisiton estimation
from tissue expression profiles."
}
\author{
Ziyi Li <zli16@mdanderson.org>
}
\examples{
# simulate mixed data from complex tissue
# without prior
Wmat <- matrix(abs(rnorm(60*3, 4, 4)), 60, 3)
SelMarker <- list(CellA = 1:20,
                  CellB = 21:40,
                  CellC = 41:60)
for(i in 1:3) {
     Wmat[SelMarker[[i]], i] <- abs(rnorm(20, 50, 10))
}
Hmat <- matrix(runif(3*25), 3, 25)
Hmat <- sweep(Hmat, 2, colSums(Hmat), "/")
Ymat <- Wmat \%*\% Hmat + abs(rnorm(60*25))
rownames(Ymat) <- 1:60

# deconvolution with TOAST/-P (TOAST without prior)
res <- MDeconv(Ymat, SelMarker, verbose = FALSE)
print(dim(Ymat))
cor(t(res$H), t(Hmat))

# supporse we observe the proportions
# for the same tissue from another study
alpha_prior <- rep(0.33, 3)
sigma_prior <- rep(1, 3)

# deconvolution with TOAST/+P (TOAST with prior)
res2 <- MDeconv(Ymat, SelMarker,
                alpha = alpha_prior, sigma = sigma_prior,
                verbose = FALSE)
cor(t(res2$H), t(Hmat))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
