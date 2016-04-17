\name{IPBT_rank.main.stHM}
\alias{IPBT_rank.main.stHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IPBT_rank.main.stHM(Control, Treatment, IPBT.prior = FALSE, groupNumber = 20, history = NA, IPBT.id = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Control}{
%%     ~~Describe \code{Control} here~~
}
  \item{Treatment}{
%%     ~~Describe \code{Treatment} here~~
}
  \item{IPBT.prior}{
%%     ~~Describe \code{IPBT.prior} here~~
}
  \item{groupNumber}{
%%     ~~Describe \code{groupNumber} here~~
}
  \item{history}{
%%     ~~Describe \code{history} here~~
}
  \item{IPBT.id}{
%%     ~~Describe \code{IPBT.id} here~~
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
function (Control, Treatment, IPBT.prior = FALSE, groupNumber = 20, 
    history = NA, IPBT.id = NA) 
{
    if ((IPBT.prior == FALSE & all(is.na(history) == 1)) | (IPBT.prior == 
        TRUE & all(is.na(IPBT.id) == 1))) 
        stop("Historical information is missing!\nPlease provide historical data or use IPBT prior!")
    if (IPBT.prior == FALSE) {
        hist_var = apply(history, 1, var)
    }
    if (IPBT.prior == TRUE) {
        data(IPBT3digits)
        data(SampleSize)
        hist_var = IPBT3digits[, IPBT.id]^2
    }
    Bayes.post.prob.res = IPBT.rank.stHM(Control, Treatment, 
        hist_var, groupNumber)
    Adjust_t = Bayes.post.prob.res$adjust_t
    Pvalue_appro = Bayes.post.prob.res$Pvalue
    FDR = p.adjust(Pvalue_appro, method = "BH")
    Output = data.frame(Probe_id = row.names(Bayes.post.prob.res), 
        Adjust_t = Adjust_t, fold_change = rowMeans(Control) - 
            rowMeans(Treatment), P_value_appro = Pvalue_appro, 
        False_Discovery_Rate = FDR)
    Output = Output[order(abs(Output$Adjust_t), decreasing = T), 
        ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
