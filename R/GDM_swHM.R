

GDM_swHM = function(hist_var,size,control)
{   
    # size is the half size of the window
    
    s2est=apply(control,1,var)
    
    gene_num=length(hist_var)  
    groupSD=numeric(gene_num)
    order_sd=order(hist_var)
    
    
    ix = order_sd[1:(2*size+1)]
    groupSD[1:size]=sd(bayesHierVar.swHM(control[ix,],s2est[ix]) )
    
    groupSD[(size+1):(gene_num-size) ] = sapply((size+1):(gene_num-size), 
                 function(x){ ix = order_sd[(x-size):(x+size)]
                              sd(bayesHierVar.swHM(control[ix,],s2est[ix])) } )      
    
    
    ix = order_sd[(gene_num-2*size):gene_num]
    groupSD[(gene_num-size+1):gene_num]=sd(bayesHierVar.swHM(control[ix,],s2est[ix]) )
    
    
    mean(groupSD)
    
    
}


