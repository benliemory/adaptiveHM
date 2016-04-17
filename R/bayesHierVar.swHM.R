bayesHierVar.swHM <-
function(x,s2est){
    K=ncol(x)    #### number of observations
    v=K-1
    I=nrow(x)     ### number of genes
    
    s2bar=mean(s2est)
    S=var(s2est)*(I-1)
    Best=(2/v)/(1+2/v)*(I-1)/I + 1/(1+2/v)*(2/v)*s2bar^2*(I-1)/S
    if(Best<1){
        var_hier_est=(1-Best)*s2est+Best*s2bar
        result=var_hier_est
    } else{
        result=rep(s2bar,I)
        names(result) = names(s2est)
    }
    return(result)
}
