\name{bayesHierVar.swHM}
\alias{bayesHierVar.swHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance estimate by empirical Bayes Hierarchical model (swHM)
}
\description{
This function estimates variance for each gene using empirical hierarchical model with the gene's "neighbouring" information for swHM. After deciding the window size, this function will be applied to each gene. To make it more efficient, we calculate each gene's sample variance outside this function. 
}
\usage{
bayesHierVar.swHM(subdata, s2est)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subdata}{
Data matrix for genes included in this hierarchical model.
}
  \item{s2est}{
Sample variance for each gene.
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
