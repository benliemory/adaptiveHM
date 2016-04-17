\name{GDM_swHM}
\alias{GDM_swHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate GDM for given window size of swHM
}
\description{
Given control data, historical variance (or historical variance rank) and window size of swHM, this function calculate correspoding GDM. 
}
\usage{
GDM_swHM(hist_var,size,control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hist_var}{
historical variance for each gene.
}
  \item{size}{
Given window size for swHM.
}
  \item{control}{
Data matrix for control samples.
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
