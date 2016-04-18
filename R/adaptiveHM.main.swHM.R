adaptiveHM.main.swHM <-
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
        
        
        mean_treat=rowMeans(Treatment)
        mean_control=rowMeans(Control)
        mu=mean_treat-mean_control
        n = dim(Control)[2]
        
        var=numeric(length = nrow(Control))
        
        row.names(Control)=1:nrow(Control)
        #geneNames = row.names(Control)
        geneNames = names(hist_var)
        
        s2est=apply(Control,1,var)
        
        gene_num=length(hist_var)             
        order_var=order(hist_var)
        
        ids=order_var[1:winSize]
        ix = order_var[1:(2*winSize+1)]
        var[ids]=bayesHierVar.swHM(Control[ix,],s2est[ix])[1:winSize]
        
        var[order_var[(winSize+1):(gene_num-winSize)]] = sapply((winSize+1):(gene_num-winSize), 
                                                                function(x){ ix = order_var[(x-winSize):(x+winSize)]
                                                                             bayesHierVar.swHM(Control[ix,],s2est[ix])[(winSize+1)]})      
        
        ids=order_var[(gene_num-winSize+1):gene_num]
        ix = order_var[(gene_num-2*winSize):gene_num]
        var[ids]= bayesHierVar.swHM(Control[ix,],s2est[ix])[(winSize+2):(2*winSize+1)]
        
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
