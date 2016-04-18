adaptiveHM.main.stHM <-
function(Control,Treatment, IPBT.prior=FALSE,groupNumber = 20,
                     history=NA,IPBT.id=NA)
{
    if ( (IPBT.prior==FALSE & all(is.na(history)==1)  ) | 
             (IPBT.prior==TRUE &  all(is.na(IPBT.id)==1) )  ) 
        stop("Historical information is missing!\nPlease provide historical data or use IPBT prior!") 
    
    if(IPBT.prior==FALSE)
    {
        hist_var = apply(history,1,var)
    }
    
    if(IPBT.prior==TRUE)
    {   
        data(IPBT3digits)
        hist_var = IPBT3digits[,IPBT.id]^2
    }
    
    mean_treat=rowMeans(Treatment)
    mean_control=rowMeans(Control)
    mu=mean_treat-mean_control
    n = dim(Control)[2]
    
    var=numeric(length = nrow(Control))
    
    row.names(Control)=1:nrow(Control)
    #geneNames = row.names(Control)
    geneNames = names(hist_var)
    
    breakPoint = quantile(hist_var, (1:groupNumber)/groupNumber )
    
    groupID = findInterval(hist_var, breakPoint,rightmost.closed = TRUE)
    
    groupRes = unname(by(Control,groupID,bayesHierVar.stHM) )
    
    vecRes = unlist(groupRes)
    var = vecRes[order(as.numeric(names(vecRes)))]
    
    if(is.null(geneNames)==0)
    {names(var) = geneNames }
    
    Adjust_t = mu/sqrt(var)*sqrt(2/n)
    
    Pvalue_appro = 2*pt(abs(Adjust_t) ,df = 2*n -2, lower.tail = FALSE)
    
    FDR = p.adjust(Pvalue_appro, method = "BH" )
    
    
    
    
    Output = data.frame(Probe_id = geneNames,
                        Adjust_t = Adjust_t,
                        fold_change = rowMeans(Control) - rowMeans(Treatment), 
                        P_value_appro = Pvalue_appro,
                        False_Discovery_Rate = FDR 
                        )
    
    Output = Output[order(abs(Output$Adjust_t),decreasing = T),]
    
    
   
    
    
}
