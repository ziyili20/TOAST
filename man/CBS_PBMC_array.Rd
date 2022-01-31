\name{CBS_PBMC_array}
\alias{CBS_PBMC_array}
\docType{data}
\title{
An example dataset for partial reference-free
cell composition estimation from tissue gene expression
}
\description{
The dataset contains 511 microarray gene expressions
for 20 PBMC samples (mixed_all) and PBMC microarray
reference for the matched 511 genes from 5immune
cell types (LM_5). It also contains the true cell
compositions from cell sorting experiment (trueProp)
and prior knowledge of cell compositions for 5 cell types
in PBMC (prior_alpha and prior_sigma).
}
\usage{data("CBS_PBMC_array")}
\references{
Newman, Aaron M., et al. "Robust enumeration
of cell subsets from tissue expression profiles."
Nature methods 12.5 (2015): 453.

Rahmani, Elior, et al. "BayesCCE: a Bayesian
framework for estimating cell-type composition
from DNA methylation without the need for methylation
reference." Genome biology 19.1 (2018): 141.
}
\examples{
data("CBS_PBMC_array")
CBS_PBMC_array$mixed_all[1:5,1:5]
head(CBS_PBMC_array$LM_5,3)
head(CBS_PBMC_array$trueProp,3)
CBS_PBMC_array$prior_alpha
CBS_PBMC_array$prior_sigma
}
\keyword{datasets}
