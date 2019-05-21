\name{findRefinx_var}
\alias{findRefinx_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
findRefinx_var
}
\description{
Find index for marker genes with largest variance based on raw data.
}
\usage{
findRefinx_var(rawdata, nmarker=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawdata}{
     A data matrix with rows representing features
     and columns represeting samples.
  }
  \item{nmarker}{
     Desired number of markers after selection.
     Default is 1000.
  }
}
\value{
A vector of index for the selected markers.
}
\references{
Ziyi Li, Zhijin Wu, Peng Jin, Hao Wu.
"Dissecting differential signals in
high-throughput data from
complex tissues."
}
\author{
Ziyi Li <ziyi.li@emory.edu>
}
\examples{
Y_raw <- matrix(runif(5000*20, 0, 1), 5000, 20)
idx <- findRefinx_var(Y_raw, nmarker=500)
}
