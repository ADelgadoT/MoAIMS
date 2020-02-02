#fit model

fitEM <- function(bin_count_obj=NULL,
                  min_n_x=50,
                  trunc_val=0.999,
                  k=3,
                  reg="gam",
                  mu_init=0.1, 
                  var_init=0.15,
                  iter_stop=100)
  {
      
    ip=bin_count_obj@ip
    input=bin_count_obj@input
    input_trunc=quantile(input,trunc_val)
    
    
    #exprlev=bin_count_obj@exprlev
    
    
    ###strata input
    Y=ip
    X=input
    Y_freq <- table(Y)
    
    
    # adaptive griding for x
    X_set <- sort( unique(X), decreasing=TRUE )
    ind_X_set <- rep( 0, length(X_set) )
    
    ind_now <- 1
    N_now <- 0
    
    for ( i in 1:length(X_set) )
    {
      N_i <- length( which( X==X_set[i] ) )
      if ( N_now <= min_n_x )
      {
        ind_X_set[i] <- ind_now
        N_now <- N_now + N_i
      } else
      {
        ind_now <- ind_now + 1
        ind_X_set[i] <- ind_now
        N_now <- N_i
      }
    }
    
    X_set_new <- rep( 0, length(X_set) )    
    for ( i in 1:length(unique(ind_X_set)) )
    {
      X_set_new[ind_X_set==i] <- median( X_set[ind_X_set==i] )
    }
    
    X_new <- rep( 0, length(X) )    
    for ( i in 1:length(X_set) )
    {
      X_new[ X==X_set[i] ] <- X_set_new[i]
    }    
    
    
    #after strata, trunc
    X=X_new
    X[which(X>input_trunc)]=input_trunc
    
    
    #mom for background mean
    df_count=data.frame(Y=Y,X=X)
    df_stat=as.data.frame(df_count %>% group_by(X) %>% summarise(mean0=median(Y),var0=stats::mad(Y)^2,n=n())) 
    
    
    mean0=df_stat$mean0
    var0=df_stat$var0
    n_u=df_stat$n
    
    a_u=rep(NA,length(mean0))
    b_u=rep(NA,length(mean0))
    for (i in 1:length(mean0)){
      if( !is.na(var0[i]) && !is.na(mean0[i]) &&  var0[i] > mean0[i] && n_u[i]>=min_n_x) {
        # use MOM only when we have proper mean and variance
        
        a_u[i] <- mean0[i]^2 / (var0[i]-mean0[i])
        b_u[i] <- mean0[i] / (var0[i]-mean0[i]) 
        
      }
    }
    
    par_est=list(mean0_u=mean0,    #mean
                 var0_u=var0,    #var
                 a_u=a_u,    #size
                 b_u = b_u,    #size/mu
                 n_u=n_u,    
                 #exprlev_u=df_stat$exprlev,
                 X_u=df_stat$X,    #x strata
                 Y_val = as.numeric(names(table(Y))), Y_freq = table(Y))       
    
    
    
    #lm fit
    idNA <- which( is.na(par_est$a_u) )
    
    if ( length(idNA)>0 ) {   
      mean0_u <- par_est$mean0_u[-idNA]
      var0_u <- par_est$var0_u[-idNA]
      a_u <- par_est$a_u[-idNA]
      b_u <- par_est$b_u[-idNA]
      n_u <- par_est$n_u[-idNA]
      X_u <- par_est$X_u[-idNA]
      #exprlev_u <- par_est$exprlev_u[-idNA]
      
    } else {      
      mean0_u <- par_est$mean0_u
      var0_u <- par_est$var0_u 
      a_u <- par_est$a_u
      b_u <- par_est$b_u
      n_u <- par_est$n_u
      X_u <- par_est$X_u
      #exprlev_u <- par_est$exprlev_u[-idNA]
    }
    
    a_weight=sum(n_u*a_u)/sum(n_u)
    
  
    Xtrans <- log( X_u + 1 ) #same with log result
        
    
    
    if (reg=="rlm"){  
    	fit <- rlm(log(mean0_u) ~ Xtrans , weights=n_u/sum(n_u) )    
    } else {
    	fit <- gam(log(mean0_u) ~ s(Xtrans),method="REML")
    }       
    
    
    
       
    #coef_val <- coef(fit)
    
    
    #estimate background mu
    oriX=input
    oriX[which(oriX>input_trunc)]=input_trunc

	mu_est <- as.numeric(exp(predict(fit, data.frame(Xtrans = log(oriX+1)))))
	#mu_est <- as.numeric(exp(predict(fit, data.frame(Xtrans = log(oriX+1),exprlev_u=as.factor(exprlev)))))
        
    
    
    
    fitZ0=list(a = a_weight, muEst=mu_est,
        		Y_val=par_est$Y_val, Y_freq=par_est$Y_freq )
        
    pNfit <- mosaics:::.calcPN( Y=Y, k=k, a=fitZ0$a, mu_est=fitZ0$muEst ) 
    
    #1S
	fitZ1_1S <- em1S( fitZ0, Y=Y, 
        			pNfit=pNfit, k=k, mu_init=mu_init, var_init=var_init, iter_stop=iter_stop )
    
    #2S
    fitZ0_update=list(pi0=fitZ1_1S$pi0,a = a_weight, muEst=mu_est,
        			Y_val=par_est$Y_val, Y_freq=par_est$Y_freq )
    
    Y_bd_all <- calcYbdAll( fitZ0_update, k=k )
	fitZ1_2S <- mosaics:::.mosaicsZ1_2S( fitZ0_update, Y=Y, 
        pNfit=pNfit, Y_bd_all=Y_bd_all, k=k )
    
    
    resEst <- new( "MosaicsFitEst",
        pi0=fitZ1_1S$pi0, a=fitZ0$a,muEst=fitZ0$muEst, pNfit=pNfit,
        b=fitZ1_1S$b, c=fitZ1_1S$c,
        p1=fitZ1_2S$p1, b1=fitZ1_2S$b1, c1=fitZ1_2S$c1, b2=fitZ1_2S$b2, c2=fitZ1_2S$c2,
        inputTrunc=input_trunc, analysisType="IO" )
        		     
    
    fitMD1 <- mosaics:::.margDist_1S( mosaicsEst=resEst, 
		tagCount=Y, pNfit=pNfit, k=k)   
	
    loglik_1S <- mosaics:::.logLik( mosaicsEst=resEst, tagCount=Y, 
        pNfit=pNfit, k=k,signalModel="1S")

	fitBIC_1S <- calcModelBIC( 
    	loglik=loglik_1S, n=length(Y),signalModel="1S",reg)  
    
    
    fitMD2 <- mosaics:::.margDist_2S( mosaicsEst=resEst, 
        tagCount=Y, pNfit=pNfit, k=k )
    
    loglik_2S <- mosaics:::.logLik( mosaicsEst=resEst, tagCount=Y, 
        pNfit=pNfit, k=k,signalModel="2S")

	fitBIC_2S <- calcModelBIC( 
    	loglik=loglik_2S, n=length(Y),signalModel="2S",reg)
      
    
              
    fit_obj=new('fitObj')
    fit_obj@mu_b=mu_est
    fit_obj@size_b=a_weight  
    fit_obj@pi_s=1-resEst@pi0
    fit_obj@mu_s=resEst@b/resEst@c
    fit_obj@size_s=resEst@b
    
    
    fit_obj@pi_1s=resEst@p1
    fit_obj@mu_1s=resEst@b1/resEst@c1
    fit_obj@size_1s=resEst@b1
    fit_obj@mu_2s=resEst@b2/resEst@c2
    fit_obj@size_2s=resEst@b2
    
    
    fit_obj@val_BIC_1s=fitBIC_1S
    fit_obj@val_BIC_2s=fitBIC_2S
    
    fit_obj@md1=fitMD1
    fit_obj@md2=fitMD2

    
    
    
    return(fit_obj)
    
}



