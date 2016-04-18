# adaptiveHM
Historical Rank-based Adaptive Hierarchical model 


## Quick Installation Guide:
If R package "devtools" has been loaded, adaptiveHM can be easily installed by following: 

library(devtools)

install_github(repo = "benliemory/adaptiveHM")

## A Simple Example:

library(adaptiveHM)

\#####  Load Example Data

data(ExampleData)

\##### Using stHM with historical information frin IPBT prior 
\##### IPBT id 35 indicates normal heart data

\##### Step I Use GDM to determine optimal group number (Step I can be omitted if group number has been decided)

GDMs = GDM_stHM_figure(ExampleData$Control, IPBT.prior=TRUE,IPBT.id=35,groupRanges = 1:200)

\##### From GDM figure, we use 50 groups in stHM

\##### Step II Obtain DE gene list from Stratified Hierarchial model
DE_gene_lists_stHM = adaptiveHM.main.stHM(ExampleData$Control,ExampleData$Treatment, 
                                    IPBT.prior=TRUE,IPBT.id=35,
                                    groupNumber = 50 )

\##### Load 10 samples of normal heart data as external historical information

data(history)

\##### Using swHM with external data as historical data

\##### Step I Use GDM to determine optimal window size (Step I can be omitted if window size has been decided)

GDMs_swHM = GDM_swHM_figure(ExampleData$Control, IPBT.prior = FALSE,history = history,
                   winSizeRanges = seq(from = 5,to = 200,by = 5) )

\##### From GDM figure, we use 10 as window size for swHM

\##### Step II Obtain DE gene list from Stratified Hierarchial model
\##### Notice: Window size here refers to one side window size. The whole window is 2*window_size + 1 

DE_gene_lists_swHM = adaptiveHM.main.swHM(ExampleData$Control,ExampleData$Treatment,
                                        IPBT.prior = FALSE,history = history,
                                        winSize = 10)

## Reference
Li B, Sun Z, He Q, Zhu Y, Qin ZS (2015) Bayesian inference with historical data-based informative priors improves detection of differentially expressed genes Bioinformatics (Oxford, England) doi:10.1093/bioinformatics/btv631
