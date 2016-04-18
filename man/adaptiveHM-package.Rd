\name{adaptiveHM-package}
\alias{adaptiveHM-package}
\alias{adaptiveHM}
\docType{package}
\title{
Adaptive Hierarchical Model (adaptiveHM)
}
\description{
adaptiveHM, using the rank information from historical data, investigates the feasibility and effectiveness of improving hierarchical models (HM) in the problem of detecting differentially expressed genes by applying hierarchical models to similar genes.
}
\details{

The DESCRIPTION file:
\tabular{ll}{
Package: \tab adaptiveHM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-04-17\cr
License: \tab GPL-2\cr
}
adaptiveHM includes two different approaches: stratified hierarchical model (stHM) and sliding window hierarchical model (swHM). stHM divides all genes into a few groups and HM is applied in each group while swHM applies HM in the "neighbouring" genes. We also provde a "Group Dividing Metric" (GDM) to determine an optimal group number or window size. adaptiveHM can also be applied to Methylation 450K array. What is more, only using historical rank makes it possible for adaptiveHM to borrow information across different platforms. Similar to IPBT, this package also provides users the possibility of using their own historical data or IPBT priors.96 pre-calculated IPBT priors are built in the package.
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

library(adaptiveHM)

##### Load Example Data

data(ExampleData)

##### Using stHM with historical information frin IPBT prior 
##### IPBT id 35 indicates normal heart data

##### Step I Use GDM to determine optimal group number (Step I can be omitted if group number has been decided)

GDMs = GDM_stHM_figure(ExampleData$Control, IPBT.prior=TRUE,IPBT.id=35,groupRanges = 1:200)

##### From GDM figure, we use 50 groups in stHM

##### Step II Obtain DE gene list from Stratified Hierarchial model
DE_gene_lists_stHM = adaptiveHM.main.stHM(ExampleData$Control,ExampleData$Treatment, 
                                    IPBT.prior=TRUE,IPBT.id=35,
                                    groupNumber = 50 )

##### Load 10 samples of normal heart data as external historical information

data(history)

##### Using swHM with external data as historical data

##### Step I Use GDM to determine optimal window size (Step I can be omitted if window size has been decided)

GDMs_swHM = GDM_swHM_figure(ExampleData$Control, IPBT.prior = FALSE,history = history,
                   winSizeRanges = seq(from = 5,to = 200,by = 5) )

##### From GDM figure, we use 10 as window size for swHM

##### Step II Obtain DE gene list from Stratified Hierarchial model
##### Notice: Window size here refers to one side window size. The whole window is 2*window_size + 1 

DE_gene_lists_swHM = adaptiveHM.main.swHM(ExampleData$Control,ExampleData$Treatment,
                                        IPBT.prior = FALSE,history = history,
                                        winSize = 10)


}
