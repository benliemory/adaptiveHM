\name{IPBT.rank.stHM}
\alias{IPBT.rank.stHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IPBT.rank.stHM(control, treat, hist_var, groupNumber)
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
  \item{groupNumber}{
%%     ~~Describe \code{groupNumber} here~~
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
function (control, treat, hist_var, groupNumber) 
{
    mean_treat = rowMeans(treat)
    mean_control = rowMeans(control)
    mu = mean_treat - mean_control
    n = dim(control)[2]
    var = numeric(length = nrow(control))
    geneNames = row.names(control)
    row.names(control) = 1:nrow(control)
    breakPoint = quantile(hist_var, (1:groupNumber)/groupNumber)
    groupID = findInterval(hist_var, breakPoint, rightmost.closed = TRUE)
    groupRes = unname(by(control, groupID, bayesHierVar.stHM))
    vecRes = unlist(groupRes)
    var = vecRes[order(as.numeric(names(vecRes)))]
    if (is.null(geneNames) == 0) {
        names(var) = geneNames
    }
    adjust_t = mu/sqrt(var) * sqrt(2/n)
    Pvalue = 2 * pt(abs(adjust_t), df = 2 * n - 2, lower.tail = FALSE)
    Res = data.frame(adjust_t = adjust_t, Pvalue = Pvalue)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
