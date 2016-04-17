IPBT.rank.stHM <-
function(control,treat,hist_var,groupNumber)
    { 
        
        mean_treat=rowMeans(treat)
        mean_control=rowMeans(control)
        mu=mean_treat-mean_control
        n = dim(control)[2]
        
        var=numeric(length = nrow(control))
        
        geneNames = row.names(control)
        row.names(control)=1:nrow(control)
        
        breakPoint = quantile(hist_var, (1:groupNumber)/groupNumber )
        
        groupID = findInterval(hist_var, breakPoint,rightmost.closed = TRUE)
        
        groupRes = unname(by(control,groupID,bayesHierVar.stHM) )

        vecRes = unlist(groupRes)
        var = vecRes[order(as.numeric(names(vecRes)))]
        
        if(is.null(geneNames)==0)
        {names(var) = geneNames }
        
        adjust_t = mu/sqrt(var)*sqrt(2/n)
        
        Pvalue = 2*pt(abs(adjust_t) ,df = 2*n -2, lower.tail = FALSE)
               
        
        Res = data.frame(adjust_t = adjust_t, Pvalue = Pvalue  )
        
                
        
    }
