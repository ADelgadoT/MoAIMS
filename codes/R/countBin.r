# bin count

countBin <- function(bin_info_obj=NULL,
                     ip_bam=NULL,
                     input_bam=NULL,
                     strand_specific=NULL,                       
                     is_paired=NULL)
{

  feat_input=data.frame(GeneID=bin_info_obj@bin_id,
                       Chr=bin_info_obj@chrom,
                       Start=bin_info_obj@start,
                       End=bin_info_obj@end,
                       Strand=bin_info_obj@strand)

  
  bam_files=c(ip_bam,input_bam)
  
  if (strand_specific==0){
  
  	if (is_paired==TRUE){
  		res=featureCounts(files=bam_files,isPaired=TRUE,requireBothEndsMapped=TRUE,
                      annot.ext=feat_input,useMetaFeatures = T,read2pos=5,
                      ignoreDup=T,allowMultiOverlap=T)
  	} else {
  		res=featureCounts(files=bam_files,
                      annot.ext=feat_input,useMetaFeatures = T,read2pos=5,
                      ignoreDup=T,allowMultiOverlap=T)
  	} 
    
  } else if (strand_specific==1) {
  
  	if (is_paired==TRUE){
  		res=featureCounts(files=bam_files,strandSpecific=1, isPaired=TRUE,requireBothEndsMapped=TRUE,
                      annot.ext=feat_input,useMetaFeatures = T,read2pos=5,
                      ignoreDup=T,allowMultiOverlap=T)
  		
  	} else {
  		res=featureCounts(files=bam_files,strandSpecific=1,
                      annot.ext=feat_input,useMetaFeatures = T,read2pos=5,
                      ignoreDup=T,allowMultiOverlap=T)
  	}
    
    
  } else {  
  
  	if (is_paired==TRUE){
  		
  		res=featureCounts(files=bam_files,strandSpecific=2, isPaired=TRUE,requireBothEndsMapped=TRUE,
                      annot.ext=feat_input,useMetaFeatures = T,read2pos=5,
                      ignoreDup=T,allowMultiOverlap=T)
  	} else {
  		res=featureCounts(files=bam_files,strandSpecific=2,
                      annot.ext=feat_input,useMetaFeatures = T,read2pos=5,
                      ignoreDup=T,allowMultiOverlap=T)
  	} 
       
  }
  
  
  df_res=data.frame(res$counts,stringsAsFactors=FALSE)
  df_res$bin_id=rownames(df_res)
  
  
  
  bin_count_obj=new('binCount')
  bin_count_obj@bin_id=df_res[,"bin_id"]
  bin_count_obj@ip=df_res[,1]
  bin_count_obj@input=df_res[,2]
  
  return(bin_count_obj)

}



  
  
  


