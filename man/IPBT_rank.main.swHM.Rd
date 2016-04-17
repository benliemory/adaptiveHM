\name{IPBT_rank.main.swHM}
\alias{IPBT_rank.main.swHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detect DE genes by swHM
}
\description{
swHM applies HM in the "neighbouring" genes. Only rank information from historical data is used.
}
\usage{
IPBT_rank.main.swHM(Control, Treatment, IPBT.prior = FALSE, winSize = 50, history = NA, IPBT.id = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Control}{
Data matrix for control samples.
}
  \item{Treatment}{
Data matrix for treatment samples.
}
  \item{IPBT.prior}{
Logical value to specify using external historical data or IPBT prior.
}
  \item{winSize}{
Window size used in swHM. It is better to use GDM to determine optimal number of groups first. 
}
  \item{history}{
External historical data is required to provide if using external historical data.
}
  \item{IPBT.id}{
IPBT prior id is required if using IPBT prior.
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
