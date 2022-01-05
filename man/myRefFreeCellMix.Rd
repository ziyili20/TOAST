\name{myRefFreeCellMix}
\alias{myRefFreeCellMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Replicate the function RefFreeCellMix() from RefFreeEWAS package
}
\description{
Replicate the function RefFreeCellMix() from RefFreeEWAS
package (https://cran.r-project.org/web/packages/RefFreeEWAS/index.html)
as that package is not in CRAN anymore}
\usage{myRefFreeCellMix(Y,mu0=NULL,K=NULL,iters=10,Yfinal=NULL,verbose=TRUE)}
\arguments{
  \item{Y}{Matrix (m CpGs x n Subjects) of DNA methylation beta values}
  \item{mu0}{Matrix (m CpGs x K cell types) of *initial* cell-type specific methylomes}
  \item{K}{Number of cell types (ignored if mu0 is provided)}
  \item{iters}{Number of iterations to execute}
  \item{Yfinal}{Matrix (m* CpGs x n Subjects) of DNA methylation beta values on which to base final methylomes}
  \item{verbose}{Report summary of errors after each iteration?}
}
\details{Reference-free decomposition of DNA methylation matrix into cell-type distributions and cell-type
methylomes, Y = Mu Omega^T.  Either an initial estimate of Mu must be provided, or else the number of cell types K,
in which case RefFreeCellMixInitialize will be used to initialize.  Note that the decomposition will be based on Y,
but Yfinal (=Y by default) will be used to determine the final value of Mu based on the last iterated value of Omega.}
\value{
  Object of S3 class RefFreeCellMix, containing the last iteration of Mu and Omega.
}
\references{Houseman, E. Andres, Kile, Molly L., Christiani, David C., et al. Reference-free deconvolution of DNA methylation data and mediation by cell composition effects. BMC bioinformatics, 2016, vol. 17, no 1, p. 259.
}
\author{E. Andres Houseman}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
