GDM_stHM_figure <-
function(Control, IPBT.prior=FALSE,groupRanges = 1:200,
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
    

    
    GDM = sapply(groupRanges,function(m){GDM_stHM(hist_var,m,Control)})
    
    plot(groupRanges,GDM,xlab="Group Numbers",ylab = "GDM",col="blue",lwd = 2)
    
    GDM

}
