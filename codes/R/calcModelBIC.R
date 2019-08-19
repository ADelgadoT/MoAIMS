
calcModelBIC <- function( loglik, n,signalModel="2S", reg)
{     
    # calculate number of parameters for each case
    
    
    if (reg=="rlm"){
    	npar=1+2+1+2	
    } else{
    	npar=1+10+1+2
    }
    
    if ( signalModel == "2S" ) { 
    	npar <- npar + 3 }
        
    # calculate BIC
    
    penalty=log(n)

    
    val <- -2 * loglik + penalty * npar
    
    return(val)
}
