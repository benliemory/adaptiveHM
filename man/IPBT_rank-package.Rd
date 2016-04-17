\name{IPBT_rank-package}
\alias{IPBT_rank-package}
\alias{IPBT_rank}
\docType{package}
\title{
Informative Priors Bayesian Test based on historical ranks (IPBT_rank)
}
\description{
IPBT_rank, using the rank information from historical data, investigates the feasibility and effectiveness of improving hierarchical models (HM) in the problem of detecting differentially expressed genes by applying hierarchical models to similar genes. IPBT_rank includes two different approaches: stratified hierarchical model (stHM) and sliding window hierarchical model (swHM). stHM divides all genes into a few groups and HM is applied in each group while swHM applies HM in the "neighbouring" genes. We also provde a "Group Dividing Metric" (GDM) to determine an optimal group number or window size. IPBT_rank can also be applied to Methylation 450K array. What is more, only using historical rank makes it possible for IPBT_rank to borrow information across different platforms. Similar to IPBT, this package also provides users the possibility of using their own historical data or IPBT priors.96 pre-calculated IPBT priors are built in the package.  
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{IPBT_rank}
\packageIndices{IPBT_rank}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{IPBT_rank}

Maintainer: \packageMaintainer{IPBT_rank}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
