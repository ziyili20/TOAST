\name{DEVarSelect}
\alias{DEVarSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feature selection for reference-free
deconvolution using cross-cell type differential analysis
}
\description{
This function selects cross-cell type
differential features for reference-free deconvolution.
}
\usage{
DEVarSelect(Y_raw, Prop0, nMarker, bound_negative)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_raw}{
A data matrix containing P features and N samples;
or a SummarizedExperiment object.
}
  \item{Prop0}{
A N by K proportion matrix with K as number of cell types.
}
  \item{nMarker}{
Number of markers selected.
}
\item{bound_negative}{
Whether to bound all negative parameter estimators to zero.
}
}
\value{
Selected markers using cross-cell type differential analysis.
}
\references{
Ziyi Li, Zhijin Wu, Peng Jin, Hao Wu.
"Dissecting differential signals in
high-throughput data from
complex tissues."
}
\author{
Ziyi Li <zli16@mdanderson.org>
}
\examples{
Y_raw <- matrix(runif(5000*20, 0, 1), 5000, 20)
tmp <- matrix(runif(20*4), 20, 4)
Prop0 <- sweep(tmp, 1, rowSums(tmp), "/")
varlist <- DEVarSelect(Y_raw, Prop0,
                       nMarker=1000,
                       bound_negative=FALSE)
}
