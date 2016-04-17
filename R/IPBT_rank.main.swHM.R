IPBT_rank.main.swHM <-
function(Control,Treatment, IPBT.prior=FALSE,winSize = 50,
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
            data(SampleSize)
            hist_var = IPBT3digits[,IPBT.id]^2
        }
        
        
        Bayes.post.prob.res = IPBT.rank.swHM(Control,Treatment,hist_var,winSize)
        Adjust_t = Bayes.post.prob.res$adjust_t
        Pvalue_appro = Bayes.post.prob.res$Pvalue
        FDR = p.adjust(Pvalue_appro, method = "BH" )
        
        
        
        
        Output = data.frame(Probe_id = row.names(Bayes.post.prob.res),
                            Adjust_t = Adjust_t,
                            fold_change = rowMeans(Control) - rowMeans(Treatment), 
                            P_value_appro = Pvalue_appro,
                            False_Discovery_Rate = FDR 
        )
        
        Output = Output[order(abs(Output$Adjust_t),decreasing = T),]
        
        
        
        
        
    }
