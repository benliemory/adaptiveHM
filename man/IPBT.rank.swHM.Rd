\name{IPBT.rank.swHM}
\alias{IPBT.rank.swHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IPBT.rank.swHM(control, treat, hist_var, winSize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{treat}{
%%     ~~Describe \code{treat} here~~
}
  \item{hist_var}{
%%     ~~Describe \code{hist_var} here~~
}
  \item{winSize}{
%%     ~~Describe \code{winSize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (control, treat, hist_var, winSize) 
{
    mean_treat = rowMeans(treat)
    mean_control = rowMeans(control)
    mu = mean_treat - mean_control
    n = dim(control)[2]
    var = numeric(length = nrow(control))
    geneNames = row.names(control)
    row.names(control) = 1:nrow(control)
    s2est = apply(control, 1, var)
    gene_num = length(hist_var)
    order_var = order(hist_var)
    ids = order_var[1:winSize]
    ix = order_var[1:(2 * winSize + 1)]
    var[ids] = bayesHierVar.swHM(control[ix, ], s2est[ix])[1:winSize]
    var[order_sd[(winSize + 1):(gene_num - winSize)]] = sapply((winSize + 
        1):(gene_num - winSize), function(x) {
        ix = order_sd[(x - winSize):(x + winSize)]
        bayesHierVar.swHM(control[ix, ], s2est[ix])[(winSize + 
            1)]
    })
    ids = order_sd[(gene_num - winSize + 1):gene_num]
    ix = order_sd[(gene_num - 2 * winSize):gene_num]
    var[ids] = bayesHierVar.swHM(control[ix, ], s2est[ix])[(winSize + 
        2):(2 * winSize + 1)]
    adjust_t = mu/sqrt(var) * sqrt(2/n)
    Pvalue = 2 * pt(abs(adjust_t), df = 2 * n - 2, lower.tail = FALSE)
    Res = data.frame(adjust_t = adjust_t, Pvalue = Pvalue)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
