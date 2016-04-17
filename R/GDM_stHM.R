GDM_stHM <-
function(hist_var,N,control)
{
    breakPoint = quantile(hist_var, (1:N)/N )
    
    groupID = findInterval(hist_var, breakPoint,rightmost.closed = TRUE)
    
    groupRes = unname(by(control,groupID,bayesHierVar.stHM) )
    
    mean(sapply(groupRes,sd))
}
