\name{GetPrior}
\alias{GetPrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get prior knowledge for supported tissue types
}
\description{
Users can use this function to get priors provided
in this package. Users can also directly specify tisuse
type in MDeconv function.
}
\usage{
GetPrior(alpha = NULL, sigma = NULL)
}
\arguments{
  \item{alpha}{
  A string chosen from "human pbmc","human
liver", "human brain", "human pancreas", "human skin", or
a nuemric vector for manually specified alpha.
}
  \item{sigma}{
Keep it as NULL for supported tissues. Otherwise a numeric
vector for manuallly specified sigma.
}
}
\value{
\item{alpha_prior}{Prior knowledge for the mean of proportions.}
\item{sigma_prior}{Prior knowledge for the sigma of proprotions.}
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
GetPrior("human pbmc")
GetPrior("human liver")
GetPrior("human brain")
GetPrior("human skin")
GetPrior("human pancreas")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dataset }% use one of  RShowDoc("KEYWORDS")
