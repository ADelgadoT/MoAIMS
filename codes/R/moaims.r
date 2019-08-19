#' moaims main function
#'
#' This is the main function of moaims, used to detect enriched regions of MeRIP-Seq. Primary outputs include enriched regions(in BED12 format), GOF(Goodness of Fitting) plots and a summary table of models.
#' @param sample_info_file Required. Sample sheet.
#' @param gtf_file Required. Genome annotation in GTF format.
#' @param strand_specific Required. Sequencing strand protocol. 0 for unstranded, 1 for fr-first, 2 for fr-second.
#' @param is_paired Required. Paired or not.
#' @param work_dir Working directory. Default is current directory.
#' @param proj_name Project name. Default is 'NA'.
#' @param bin_size Bin size. Default is 200.
#' @param expr_thred Threshold for expression value(TPM). Default is 0.5.
#' @param per_thred At least proportion of the samples listed in the sample sheet express the genes of which bins are generated from. This parameter works when setting sep_bin_info FALSE. Default is 0.6.
#' @param min_n_x Minimum number of bins when grouping IP bins with the same Input bin count. Default is 50.
#' @param trunc_val Proportion of Input bins not taken as outliers. Default is 0.999.
#' @param k_down Minimum k. Default is 2.
#' @param k_up Maximum k. Default is 6.
#' @param iter_stop Iteration times. Default is 100.
#' @param fdr False discovery rate. Default is 0.05.
#' @param output_intmd If save intermediate output. Default is TRUE.
#' @param sep_bin_info Get bins of genes accoding genes' expression level in the corresponding sample (set TRUE) or in all the samples listed in sample sheet (set FALSE). This parameter is related to per_thred. Default is TRUE.

#' @export
#' @author Yiqian Zhang
#' @examples
#' moaims(sample_info_file = /absolute/path/to/sample_info_file, 
#'		  gtf_file =/absolute/path/to/gtf_file,
#'		  strand_specific = 1, is_paired = F)
#' @references
#' \itemize{
#' \item Kuan, P.F. et al. A Statistical Framework for the Analysis of ChIP-Seq Data. J Am Stat Assoc 106(495), 891–903 (2011)
#' \item Bao, Y. et al. Accounting for immunoprecipitation efficiencies in the statistical analysis of ChIP-seq data. BMC Bioinformatics 14, 169 (2013)
#' }


moaims <- function(
  sample_info_file=NULL,
  gtf_file=NULL,
  strand_specific=NULL,
  is_paired=NULL, 
  work_dir='.',
  proj_name='NA',
  bin_size=200,
  expr_thred=0.5,
  per_thred=0.6,
  min_n_x=50,
  trunc_val=0.999,
  k_down=2,
  k_up=6,
  iter_stop=100,
  fdr=0.05,
  output_intmd=TRUE,
  sep_bin_info=TRUE  
  ) {
  
  #check working directory
  if(dir.exists(work_dir)){
  	setwd(work_dir)
  } else{
    stop( "Stop! Working directory not found." ) 
  }
  
  #check required parameters
  if ( is.null(sample_info_file) ) { stop( "Stop! Sample sheet is required." ) }
  if ( is.null(gtf_file) ) { stop( "Stop! Annotation GTF is required." ) }
  if ( is.null(strand_specific) ) { stop( "Stop! Strand setting is required. 0 for unstranded, 1 for fr-first, 2 for fr-second." ) }
  if ( is.null(is_paired) ) { stop( "Stop! Pair-or-not setting is required." ) }
  
  
  dir.create(proj_name) 
  setwd(proj_name)
  
  if(output_intmd==TRUE){
      dir.create("intermediate") 
  }


  
  sample_info=read.delim(sample_info_file,sep="\t",header=T,stringsAsFactors = FALSE)
  
  #read sample info
  sample_id=sample_info$SampleID
  ip_bam=sample_info$BamIP
  input_bam=sample_info$BamInput
  

  run_para=list(sample_info_file=paste(getwd(),sample_info_file,sep="/"),
  			  gtf_file=gtf_file,
  			  strand_specific=strand_specific,
  			  is_paired=is_paired,
			  proj_name=proj_name,
			  bin_size=bin_size,
  			  expr_thred=expr_thred,
  			  per_thred=per_thred,
  			  min_n_x=min_n_x,
  			  trunc_val=trunc_val,
  			  k_down=k_down,
  			  k_up=k_up,
  			  fdr=fdr,
  			  output_intmd=output_intmd,
  			  sep_bin_info=sep_bin_info)

  saveRunPara(run_para)

  cat("#Get bins info.\n")
  
  bin_info_obj_all=getBinInfo(
				input_bam_files=input_bam,
  				sample_id=sample_id,
  				gtf_file=gtf_file, 
  				strand_specific=strand_specific, # 0 for unstranded, 1 for fr-first, 2 for fr-second
  				is_paired=is_paired,
  				proj_name=proj_name,
  				bin_size=bin_size,    
  				expr_thred=expr_thred,
  				per_thred=per_thred,
  				output_intmd=output_intmd,
  				sep_bin_info=sep_bin_info  
  )
  
  
  
  #get bin count
  
  fit_obj_all=NULL
  optim_k_all=NULL
  optim_reg_all=NULL
    
  for (i in 1:length(sample_id)){
    print(sample_id[i])
    cat("#Get bin counts.\n")
    
    if(sep_bin_info==TRUE){
      bin_info_obj=bin_info_obj_all[[i]]
    
    } else{
      bin_info_obj=bin_info_obj_all[[1]]
    }
    
  	bin_count_obj=countBin(bin_info_obj=bin_info_obj,
                       ip_bam=ip_bam[i],
                       input_bam=input_bam[i],
                       strand_specific=strand_specific,                       
                       is_paired=is_paired)
    
    if(output_intmd==TRUE){
    	saveBinCount(bin_count_obj,sample_id[i])
    }
    
    
    cat("#Fit models.\n")	
    
    #gam fit
    fit_obj_gam=NULL
    min_BIC_gam=Inf
    optim_k_gam=k_down
    for (k in k_down:k_up){
    
    	fit_obj_tmp=fitEM(bin_count_obj,
                  		  min_n_x=min_n_x,
                  		  trunc_val=trunc_val,
                  		  k=k,
                  		  reg="gam",
                  		  iter_stop=iter_stop)
        if (fit_obj_tmp@val_BIC_2s<min_BIC_gam){
        	fit_obj_gam=fit_obj_tmp
        	min_BIC_gam=fit_obj_tmp@val_BIC_2s
        	optim_k_gam=k
        }
    }
    
    
    #rlm fit
    fit_obj_rlm=NULL
    min_BIC_rlm=Inf
    optim_k_rlm=k_down
    for (k in k_down:k_up){
    	fit_obj_tmp=fitEM(bin_count_obj,
                  min_n_x=min_n_x,
                  trunc_val=trunc_val,
                  k=k,
                  reg="rlm",
                  iter_stop=iter_stop)
        if (fit_obj_tmp@val_BIC_2s<min_BIC_rlm){
        	fit_obj_rlm=fit_obj_tmp
        	min_BIC_rlm=fit_obj_tmp@val_BIC_2s
        	optim_k_rlm=k
        }
    }
    
    if(fit_obj_rlm@val_BIC_2s<fit_obj_gam@val_BIC_2s){
    	fit_obj=fit_obj_rlm
    	optim_k=optim_k_rlm
    	optim_reg="rlm"
    
    } else{ 
    	fit_obj=fit_obj_gam
    	optim_k=optim_k_gam	
    	optim_reg="gam"
    }
    
    fit_obj_all=c(fit_obj_all,fit_obj)
    optim_k_all=c(optim_k_all,optim_k)
	optim_reg_all=c(optim_reg_all,optim_reg)
    
    
    cat("#Save results.\n")
    
    #plotGOF
    plotGOF(fit_obj=fit_obj,
        	bin_count_obj=bin_count_obj,
        	sample_id=sample_id[i],
        	k=optim_k)
    
    
    if(fit_obj@val_BIC_1s<fit_obj@val_BIC_2s){
        #sig region 1S
	    sig_region_obj_1s=callSigRegion(bin_count_obj,fit_obj,fdr=fdr,mode="1S")
	    #output
		saveSigRegion(sig_region_obj_1s,bin_info_obj,sample_id[i])
		
		if(output_intmd==TRUE){
			saveSigBin(sig_region_obj_1s,bin_info_obj,sample_id[i])
		}
  
    } else{
        #sig region 2S
	    sig_region_obj_2s=callSigRegion(bin_count_obj,fit_obj,fdr=fdr,mode="2S")
	    #output 
	    saveSigRegion(sig_region_obj_2s,bin_info_obj,sample_id[i])
	    	    
	    if(output_intmd==TRUE){
			saveSigBin(sig_region_obj_2s,bin_info_obj,sample_id[i])
		}
    }
	
  }
  
  saveFit(fit_obj_all=fit_obj_all,
		  optim_k_all=optim_k_all,
		  optim_reg_all=optim_reg_all,
		  proj_name=proj_name)

  cat("#Done.\n")	  

}
