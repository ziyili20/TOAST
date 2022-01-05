\name{myRefFreeCellMixInitialize}
\alias{myRefFreeCellMixInitialize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Replicate the function RefFreeCellMixInitialize() from RefFreeEWAS package
}
\description{
Replicate the function RefFreeCellMixInitialize() from RefFreeEWAS
package (https://cran.r-project.org/web/packages/RefFreeEWAS/index.html)
as that package is not in CRAN anymore
}
\usage{myRefFreeCellMixInitialize(Y,K=2,Y.Distance=NULL, Y.Cluster=NULL,
    largeOK=FALSE, dist.method = "euclidean", ...)}
\arguments{
  \item{Y}{Matrix (m CpGs x n Subjects) of DNA methylation beta values}
  \item{K}{Number of cell types}
  \item{Y.Distance}{Distance matrix (object of class "dist") to use for clustering.}
  \item{Y.Cluster}{Hiearchical clustering object (from hclust function)}
  \item{largeOK}{OK to calculate distance matrix for large number of subjects? (See details.)}
  \item{dist.method}{Method for calculating distance matrix}
  \item{\dots}{Additional parameters for hclust function}
}
\details{Initializes the methylome matrix "Mu" for RefFreeCellMix by computing the mean methylation (from Y)
over K clusters of Y, determined by the Y.Cluster object.  If Y.Cluster object does not exist, it will be
created from Y.Distance (using additional clustering parameters if supplied).  If Y.Distance does not exist,
it will be created from t(Y).  As a protection against attempting to fit a very large distance matrix, the
program will stop if the number of columns of Y is > 2500, unless largeOK is explicitly set to TRUE.}
\value{
  An m x K matrix of mean methylation values.
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
