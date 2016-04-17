\name{rankIPBT-package}
\alias{rankIPBT-package}
\alias{rankIPBT}
\docType{package}
\title{
Historical Rank-based Informative Priors Bayesian Test (rankIPBT)
}
\description{
IPBT_rank, using the rank information from historical data, investigates the feasibility and effectiveness of improving hierarchical models (HM) in the problem of detecting differentially expressed genes by applying hierarchical models to similar genes.
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{rankIPBT}
\packageIndices{rankIPBT}
IPBT_rank includes two different approaches: stratified hierarchical model (stHM) and sliding window hierarchical model (swHM). stHM divides all genes into a few groups and HM is applied in each group while swHM applies HM in the "neighbouring" genes. We also provde a "Group Dividing Metric" (GDM) to determine an optimal group number or window size. IPBT_rank can also be applied to Methylation 450K array. What is more, only using historical rank makes it possible for IPBT_rank to borrow information across different platforms. Similar to IPBT, this package also provides users the possibility of using their own historical data or IPBT priors.96 pre-calculated IPBT priors are built in the package.
}
\author{
Ben Li

Maintainer: Ben Li<ben.li@emory.edu>
}
\references{
Lukk, M., et al. (2010) A global map of human gene expression, Nat Biotechnol, 28, 322-324.

Li, B., et al. (2015) Bayesian inference with historical data-based informative priors improves detection of differentially expressed genes, Bioinformatics. 
}

\examples{

}
