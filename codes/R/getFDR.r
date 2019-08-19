getFDR <- function(pp=NULL,
				  fdr=0.05)
{
	#FDR
    prob0=pp$px0
    S=length(prob0)
    Pvalue=prob0
    sortPvalue=Pvalue[order(Pvalue)]
    newPvalue=rep(0, S)
    positionp1=order(Pvalue)
    tempp=sortPvalue[1]
    newPvalue[1]=tempp/1
    for (s in 2:S)
    {
      tempp=tempp+sortPvalue[s]
      newPvalue[s]=tempp/s
    }
    newPvalueback=newPvalue[order(positionp1)]
    tempX<-ifelse(newPvalueback<=fdr, 1,0)
    
    cpoint=sum(tempX)
    T_BH=ifelse(cpoint==0, 0, sortPvalue[cpoint]) #Threshold
    tempX<-ifelse(Pvalue<=T_BH, 1, 0) ## Consider the tie
    FDR_t=ifelse(sum(tempX)==0, 0, sum(Pvalue[Pvalue<=T_BH])/sum(tempX))
    #print(FDR_t)
    
    return(list(pval=Pvalue, if_sig=tempX))

}