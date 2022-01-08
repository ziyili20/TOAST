\name{plotCorr}
\alias{plotCorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Show DE/DM state correlation between cell types

}
\description{
This function generates -log10 transformed p-values for each pair of cell types and calculate
corresponding Pearson correlation and odds ratio - a modification of function 'ggpairs' in package 'GGally'.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotCorr(pval, de.state=NULL, pval.thres=NULL, fdr.thres=NULL, p.size = 0.2,
         p.color = grDevices::adjustcolor( "black", alpha.f = 0.2),
         fig.margin = c(1,1,1,1),
         fig.margin.unit = 'in', line.type = 'dashed', line.color = 'blue')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval}{ matrix of p-values, with rows representing features and columns
  representing cell types.
}
  \item{de.state}{(optional) matrix of DE/DM states (1: DE/DM, 0:non-DE/DM), with rows representing
  features and columns representing cell types.}
  \item{pval.thres}{threshold of p-value to define DE/DM, required if de.state not provided. }
  \item{fdr.thres}{threshold of FDR to define DE/DM, required if de.state not provided and FDR is prefered than p-value.}
  \item{p.size}{point size for scatter plot}
  \item{p.color}{point color for scatter plot}
  \item{fig.margin}{figure margin}
  \item{fig.margin.unit}{unit of figure margin}
  \item{line.type}{line type in scatter plot}
  \item{line.color}{line color in scater plot}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A figure contains scatter plots, Pearson correlaiton and odds ration of -log10 transformed
p-values for each pair of cell types.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Luxiao Chen <luxiao.chen@emory.edu>%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
pval.1 <- runif(1000,0,1)
pval.2 <- pval.1 + rnorm(1000,0,0.01)
pval.2[pval.2 < 0] =0
pval.3 <- runif(1000,0,1)

pval.input <- data.frame('cell.1'=pval.1,
                         'cell.2'=pval.2,
                         'cell.3'=pval.3)

plotCorr(pval = pval.input, pval.thres = 0.05)
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
