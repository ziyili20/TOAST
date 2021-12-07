\name{findRefinx}
\alias{findRefinx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
findRefinx
}
\description{
Find index for marker genes with largest
coefficient of variation based on raw data.
}
\usage{
findRefinx(rawdata, nmarker=1000, sortBy = "var")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawdata}{
     A data matrix with rows representing
     features and columns represeting samples;
     or a SummarizedExperiment object.
  }
  \item{nmarker}{
     Desired number of markers after selection.
     Default is 1000.
  }
  \item{sortBy}{
     Desired method to select features.
     "var" represents selecting by largest variance.
     "cv" represents selecting by largest coefficients of variation.
     Default is "var".
  }
}
\value{
A vector of index for the selected markers.
}
\references{
Ziyi Li, Zhijin Wu, Peng Jin, Hao Wu.
"Dissecting differential signals in high-throughput data from
complex tissues."
}
\author{
Ziyi Li <zli16@mdanderson.org>
}
\examples{
Y_raw <- matrix(runif(5000*20, 0, 1), 5000, 20)
idx <- findRefinx(Y_raw, nmarker=500)
idx2 <- findRefinx(Y_raw, nmarker=500, sortBy = "cv")

}
