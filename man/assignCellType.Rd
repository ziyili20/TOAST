\name{assignCellType}
\alias{assignCellType}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Align cell types when reference proportions are known
}
\description{
Align target proportions with reference
proportions by pearson correlation coefficients.
}
\usage{
assignCellType(input,reference)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
     Input proportiona matrix of dimension N by K.
  }
  \item{reference}{
     Reference proportion matrix of dimension N by K.
  }
}
\value{
The aligned proportion matrix, following the cell
type ordering of reference proportion matrix.
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

## generate estimated proportion matrix
estProp <- matrix(abs(runif(50*4,0,1)), 50, 4)
estProp <- sweep(estProp, 1, rowSums(estProp), "/")

## generate reference proportion matrix
refProp <- matrix(abs(runif(50*4,0,1)), 50, 4)
refProp <- sweep(refProp, 1, rowSums(refProp), "/")

estProp_aligned = assignCellType(input = estProp,
reference = refProp)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
