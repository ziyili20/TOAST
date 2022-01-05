\name{myprojectMix}
\alias{myprojectMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Replicate the function myprojectMix() from RefFreeEWAS package
}
\description{
Replicate the function myprojectMix() from RefFreeEWAS
package (https://cran.r-project.org/web/packages/RefFreeEWAS/index.html)
as that package is not in CRAN anymore
}
\usage{myprojectMix(Y, Xmat, nonnegative=TRUE, sumLessThanOne=TRUE, lessThanOne=!sumLessThanOne)}
\arguments{
  \item{Y}{Matrix (m CpGs x n Subjects) of DNA methylation beta values}
  \item{Xmat}{Matrix (m CpGs x K cell types) of cell-type specific methylomes}
  \item{nonnegative}{All coefficients >=0?}
  \item{sumLessThanOne}{Coefficient rows should sum to less than one?}
  \item{lessThanOne}{Every value should be less than one (but possibly sum to value greater than one)?}
}
\details{Function for projecting methylation values (Y) onto space of methyomes (Xmat), with various constraints.
This is the reference-based
method described in Houseman et al. (2012) and also appearing in the minfi package.}
\value{
  Projection coefficients resulting from constrained projection
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
