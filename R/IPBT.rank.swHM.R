IPBT.rank.swHM <-
function(control,treat,hist_var,winSize)
    { 
        
        mean_treat=rowMeans(treat)
        mean_control=rowMeans(control)
        mu=mean_treat-mean_control
        n = dim(control)[2]
        
        var=numeric(length = nrow(control))
        
        geneNames = row.names(control)
        row.names(control)=1:nrow(control)
        
        s2est=apply(control,1,var)
                
        gene_num=length(hist_var)             
        order_var=order(hist_var)
               
        ids=order_var[1:winSize]
        ix = order_var[1:(2*winSize+1)]
        var[ids]=bayesHierVar.swHM(control[ix,],s2est[ix])[1:winSize]
        
        var[order_sd[(winSize+1):(gene_num-winSize)]] = sapply((winSize+1):(gene_num-winSize), 
                     function(x){ ix = order_sd[(x-winSize):(x+winSize)]
                                  bayesHierVar.swHM(control[ix,],s2est[ix])[(winSize+1)]})      
        
        ids=order_sd[(gene_num-winSize+1):gene_num]
        ix = order_sd[(gene_num-2*winSize):gene_num]
        var[ids]= bayesHierVar.swHM(control[ix,],s2est[ix])[(winSize+2):(2*winSize+1)]
        
        adjust_t = mu/sqrt(var)*sqrt(2/n)
        
        Pvalue = 2*pt(abs(adjust_t) ,df = 2*n -2, lower.tail = FALSE)
               
        
        Res = data.frame(adjust_t = adjust_t, Pvalue = Pvalue  )
        
                
        
    }
