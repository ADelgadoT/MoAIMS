
saveRunPara <- function(run_para=NULL){

    out=do.call(rbind.data.frame, run_para)
	rownames(out)=names(run_para)
	write.table(out,"run_para.tsv",sep="\t",col.names = F,quote=F)

}


saveBinInfo <- function(bin_info_obj=NULL,
                        proj_name='NA'){
  
  options(scipen=999)
  feat_input=data.frame(GeneID=bin_info_obj@bin_id,
                        Chr=bin_info_obj@chrom,
                        Start=bin_info_obj@start,
                        End=bin_info_obj@end,
                        Strand=bin_info_obj@strand)
  
  dir_path=paste(getwd(),"intermediate/",sep="/")
  write.table(feat_input,paste(dir_path,"bin_info_",proj_name,".tsv",sep=""),sep="\t",row.names=F,col.names = F,quote=F)
  
}



saveBinCount <- function(bin_count_obj=NULL,
                        sample_id=NULL){
  
  options(scipen=999)
  df_bin_count=data.frame(bin_id=bin_count_obj@bin_id,
                        ip=bin_count_obj@ip,
                        input=bin_count_obj@input)
  
  dir_path=paste(getwd(),"intermediate/",sep="/")
  write.table(df_bin_count,paste(dir_path,"bin_count_",sample_id,".tsv",sep=""),sep="\t",row.names=F,col.names = F,quote=F)
  
}



saveSigRegion <- function(sig_region_obj=NULL,
                     bin_info_obj=NULL,
                     sample_id=NULL){
  
  options(scipen=999)
  sig_region_list=sig_region_obj@sig_region_list
  sig_region_list$bin_id=as.character(sig_region_list$bin_id)
  sig_region_list$gene_id=unlist(lapply(strsplit(sig_region_list$bin_id,"_"), function(x) x[1]))
  sig_region_list$id_num=as.numeric(unlist(lapply(strsplit(sig_region_list$bin_id,"_"), function(x) x[3])))
  
  
  sig_region_list_ext=as.data.frame(sig_region_list %>% 
  group_by(gene_id) %>% 
  arrange(id_num) %>% 
  mutate(id_num_diff = c(1, diff(id_num))) %>% 
  mutate( merge_id = as.character(cumsum(id_num_diff != 1))) %>% 
  ungroup())
  sig_region_list_ext$mbin_id=paste(sig_region_list_ext$gene_id, "mbin", sig_region_list_ext$merge_id,sep = "_")
  
  
  
  bin_info=data.frame(bin_id=bin_info_obj@bin_id,
                      chrom=bin_info_obj@chrom,
                      start=bin_info_obj@start,
                      end=bin_info_obj@end,
                      strand=bin_info_obj@strand)
                      
  bin_info$bin_id=as.character(bin_info$bin_id)
  
  sig_region_list_ext=left_join(sig_region_list_ext,bin_info,by="bin_id")
  
  tmp_out=as.data.frame(sig_region_list_ext %>% group_by(mbin_id) %>% mutate(chromStart=min(start)-1,  #convert to 0-based
                                                       chromEnd=max(end),
                                                       blockSize=paste(c(end-start+1,''),collapse=','),
                                                       merged_size=sum(end-start+1),
                                                       score1=max(ip),
                                                       score2=max(fc),
                                                       score3=-log10(min(pval)),
                                                       blockStarts=ifelse(strand=="+",paste(start-start[1],collapse=','),paste(rev(start-start[n()]),collapse=',')),
                                                       blockCount=n()))
  
  tmp_out$thickStart=tmp_out$chromStart
  tmp_out$thickEnd=tmp_out$chromEnd
  tmp_out$itemRgb=0
  

  #sig bed12
  df_out=distinct(tmp_out[,c("chrom",
                    "chromStart",
                    "chromEnd",
                    "mbin_id",
                    "score1",
                    "strand",
                    "thickStart",
                    "thickEnd",
                    "itemRgb",
                    "blockCount",
                    "blockSize",
                    "blockStarts",
                    "score2","score3")])
    
  write.table(df_out,paste("sig_",sample_id,".bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  
}


saveSigBin <- function(sig_region_obj=NULL,
                       bin_info_obj=NULL,
                       sample_id=NULL){
  
  options(scipen=999)
  bin_info=data.frame(bin_id=bin_info_obj@bin_id,
                      chrom=bin_info_obj@chrom,
                      start=bin_info_obj@start,
                      end=bin_info_obj@end,
                      strand=bin_info_obj@strand)
                      
  bin_info$bin_id=as.character(bin_info$bin_id)
  
  sig_region_list=sig_region_obj@sig_region_list
  sig_region_list$bin_id=as.character(sig_region_list$bin_id)
  sig_region_list_ext=left_join(sig_region_list,bin_info,by="bin_id")
  
  
  tmp_out=as.data.frame(sig_region_list_ext %>% group_by(bin_id) %>% mutate(chromStart=min(start),
                                                       chromEnd=max(end),
                                                       blockSize=paste(c(end-start+1,''),collapse=','),
                                                       score=ip,
                                                       blockStarts=ifelse(strand=="+",paste(start-start[1],collapse=','),paste(rev(start-start[n()]),collapse=',')),
                                                       blockCount=n()))
  
  tmp_out$thickStart=tmp_out$chromStart
  tmp_out$thickEnd=tmp_out$chromEnd
  tmp_out$itemRgb=0
  
  
  df_out=distinct(tmp_out[,c("chrom",
                    "chromStart",
                    "chromEnd",
                    "bin_id",
                    "score",
                    "strand",
                    "thickStart",
                    "thickEnd",
                    "itemRgb",
                    "blockCount",
                    "blockSize",
                    "blockStarts")])
                    
  dir_path=paste(getwd(),"intermediate/",sep="/")
  write.table(df_out,paste(dir_path,"sig_bin_",sample_id,".bed",sep=""),sep="\t",col.names = F,quote=F,row.names=F)
  
  
}


saveFit <- function(fit_obj_all=NULL,
		  			optim_k_all=NULL,
		  			optim_reg_all=NULL,
		  			proj_name='NA',
		  			sample_id_all=NULL){

  pi_s_all=NULL
  val_BIC_1s_all=NULL
  val_BIC_2s_all=NULL
    
  for (i in 1:length(fit_obj_all)){
  	pi_s_all=c(pi_s_all,round(fit_obj_all[[i]]@pi_s,3))
  	val_BIC_1s_all=c(val_BIC_1s_all,round(fit_obj_all[[i]]@val_BIC_1s))
  	val_BIC_2s_all=c(val_BIC_2s_all,round(fit_obj_all[[i]]@val_BIC_2s))
  
  }
    			
  df_out=data.frame(sample_id=sample_id_all,
  					pi_s=pi_s_all,
  					val_BIC_1S=val_BIC_1s_all,val_BIC_2S=val_BIC_2s_all,
  					optim_k=optim_k_all,
  					optim_reg=optim_reg_all)
  
  write.table(df_out,paste("fit_res_",proj_name,".tsv",sep=""),sep="\t",row.names=F,quote=F)
  
}



