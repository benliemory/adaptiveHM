\name{adaptiveHM.main.stHM}
\alias{adaptiveHM.main.stHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detect DE genes by stHM
}
\description{
stHM divides all genes into a few groups and HM is applied in each group while swHM applies HM in the "neighbouring" genes. Only rank information from historical data is used.
}
\usage{
adaptiveHM.main.stHM(Control, Treatment, IPBT.prior = FALSE, groupNumber = 20, history = NA, IPBT.id = NA)
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
  \item{groupNumber}{
Number of groups used in stHM. It is better to use GDM to determine optimal number of groups first. 
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
Ben Li
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
