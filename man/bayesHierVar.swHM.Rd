\name{bayesHierVar.swHM}
\alias{bayesHierVar.swHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bayesHierVar.swHM(x, s2est)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{s2est}{
%%     ~~Describe \code{s2est} here~~
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
function (x, s2est) 
{
    K = ncol(x)
    v = K - 1
    I = nrow(x)
    s2bar = mean(s2est)
    S = var(s2est) * (I - 1)
    Best = (2/v)/(1 + 2/v) * (I - 1)/I + 1/(1 + 2/v) * (2/v) * 
        s2bar^2 * (I - 1)/S
    if (Best < 1) {
        var_hier_est = (1 - Best) * s2est + Best * s2bar
        result = var_hier_est
    }
    else {
        result = rep(s2bar, I)
        names(result) = names(s2est)
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
