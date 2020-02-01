# call enrich regions
callSigRegion<- function(bin_count_obj=NULL,
				  fit_obj=NULL,
				  fdr=0.05,
				  mode="2S")
{
    
    #data1=bin_count_obj@ip
    
    if (mode=="1S"){
    	
		MDZ0 <- fit_obj@md1$MDZ0
    	MDZ1 <- fit_obj@md1$MDZ1 
    
    	denom <- MDZ0*(1-fit_obj@pi_s) + MDZ1*fit_obj@pi_s
    	pH0 = MDZ0 * (1-fit_obj@pi_s) / denom  
    	pH1 = MDZ1 * fit_obj@pi_s / denom

    	if( length(which(is.na(pH0))) > 0 )
    	{
        	pH0[ is.na(pH0) ] = 1
        	pH1[ is.na(pH0) ] = 0
    	}
    
    	pp = list( px1=pH1, px0=pH0 )
    
    } else{
    
    	MDZ0 <- fit_obj@md2$MDZ0
    	MDZ1 <- fit_obj@md2$MDZ1
    	MDZ2 <- fit_obj@md2$MDZ2
    
    
    	# calculate posterior probabilities
    
    	denom <- MDZ0*(1-fit_obj@pi_s) + fit_obj@pi_s * ( MDZ1*fit_obj@pi_1s + MDZ2*(1-fit_obj@pi_1s) )
    	pH0 = MDZ0 * (1-fit_obj@pi_s) / denom  
    	pH1 = MDZ1 * fit_obj@pi_s * MDZ1*fit_obj@pi_1s / denom
    	pH2 = MDZ2 * fit_obj@pi_s * (1-MDZ1*fit_obj@pi_1s) / denom
    

    	if( length(which(is.na(pH0))) > 0 )
    	{
        	pH0[ is.na(pH0) ] = 1
        	pH1[ is.na(pH0) ] = 0
        	pH2[ is.na(pH0) ] = 0
    	}
    
    	pp = list( px2=pH2, px1=pH1, px0=pH0 )
    
    	
    }
    
        
    res_fdr=getFDR(pp,fdr)
    
    df_out=data.frame(bin_id=bin_count_obj@bin_id,
                      ip=bin_count_obj@ip,
                      input=bin_count_obj@input,
                      mu=fit_obj@mu_b,
                      fc=(bin_count_obj@ip+1)/(fit_obj@mu_b+1),
                      if_sig=res_fdr$if_sig,
                      pval=res_fdr$pval)
    
    
    sig_region_obj=new('sigRegion')
    sig_region_obj@sig_region_list=filter(df_out,if_sig==1)

    return(sig_region_obj)
    
}
  








