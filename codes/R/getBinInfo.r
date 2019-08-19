getBinInfo <- function(
  input_bam_files=NULL,
  sample_id=NULL,
  gtf_file=NULL, 
  strand_specific=NULL, # 0 for unstranded, 1 for fr-first, 2 for fr-second
  is_paired=NULL,
  proj_name='NA',
  bin_size=200,  
  batch_size=1000,    
  expr_thred=0.5,
  per_thred=0.6,
  output_intmd=TRUE,
  sep_bin_info=TRUE 
  ) {
  
  
  #read gtf
  gtf <- import(gtf_file,format="gtf")
  txdb <- makeTxDbFromGRanges(gtf)
  
  #keep standard chromosomes
  txdb=keepStandardChromosomes(txdb)
  
  #filter out chrM
  seqlevels(txdb) <- seqlevels(txdb)[!grepl("M",seqlevels(txdb))]
  
  #return all chromosomes
  #seqlevels(txdb) <- seqlevels0(txdb)
  
  
  #merge exon by gene
  exon_bygene <-  exonsBy(x=txdb,by='gene') 
  merged_gene <- reduce(exon_bygene) 
  gene_len=sum(width(merged_gene))


  #filter genes not expressed  
  if (strand_specific==0){
  
  	if (is_paired==TRUE){ 
  		gene_count <- summarizeOverlaps(exon_bygene, input_bam_files, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, inter.feature=FALSE)  	
  	} else{
  		gene_count <- summarizeOverlaps(exon_bygene, input_bam_files, mode="Union", ignore.strand=TRUE, inter.feature=FALSE)
  	}  
    
  }else if (strand_specific==1) {
    
    if (is_paired==TRUE){ 
    	gene_count <- summarizeOverlaps(exon_bygene, input_bam_files, mode="Union", singleEnd=FALSE, ignore.strand=FALSE, inter.feature=FALSE)  	
    } else {
    	gene_count <- summarizeOverlaps(exon_bygene, input_bam_files, mode="Union",ignore.strand=FALSE, inter.feature=FALSE)
    }
            
  } else {
  
  	if (is_paired==TRUE){ 
  		gene_count <- summarizeOverlaps(exon_bygene, input_bam_files, mode="Union", singleEnd=FALSE, ignore.strand=FALSE, inter.feature=FALSE,preprocess.reads=invertStrand)  		
  	} else{
  		gene_count <- summarizeOverlaps(exon_bygene, input_bam_files, mode="Union", ignore.strand=FALSE, inter.feature=FALSE,preprocess.reads=invertStrand) 	
  	}
        
  }

  
  raw_count <- as.data.frame(assays(gene_count)$counts)


  #get TPM
  tpmx=raw_count / gene_len
  tpm <- as.data.frame(t( t(tpmx) * 1e6 / colSums(tpmx) ))
  colnames(tpm)=sample_id
  tpm$gene_id=rownames(tpm)
  tpm$gene_len=gene_len  
  
  #get gene info. including chrom, gid, gene_type
  df_gtf=as.data.frame(gtf)
  df_ginfo=filter(df_gtf,type=="gene")[,c("seqnames","gene_id","gene_type")]

  #expr_val, gene_id, gene_len, chromsome, gene_type
  df_expr=left_join(tpm,df_ginfo,by="gene_id")
  
  
  bin_info_obj_all=NULL
  
  #require at least 60% input samples tpm>=0.5
  if(sep_bin_info==TRUE){
    for (i in 1:length(sample_id)){
        df_expr_sep=df_expr[,c(i,(dim(df_expr)[2]-3):dim(df_expr)[2])]
        sample_num=1
        gnames=(df_expr_sep %>% mutate(per=(rowSums(.[1:sample_num] >=expr_thred))/sample_num) %>% filter(per>=per_thred))$gene_id
        bin_info_obj=.getBinInfoObj(gnames,merged_gene,batch_size,bin_size)
        
        if(output_intmd==TRUE){
        	saveBinInfo(bin_info_obj,paste0(proj_name,"_",sample_id[i]))
        }
        
        bin_info_obj_all=c(bin_info_obj_all,bin_info_obj)
    }
    
  }else{    
    sample_num=length(input_bam_files)
  	gnames=(df_expr %>% mutate(per=(rowSums(.[1:sample_num] >=expr_thred))/sample_num) %>% filter(per>=per_thred))$gene_id
    bin_info_obj=.getBinInfoObj(gnames,merged_gene,batch_size,bin_size)
    
    if(output_intmd==TRUE){
    	saveBinInfo(bin_info_obj,proj_name)
    }
    
    bin_info_obj_all=c(bin_info_obj_all,bin_info_obj)
    
  }
  
  
  return(bin_info_obj_all)
  

}


.getBinInfoObj <- function(gnames,merged_gene,batch_size=1000,bin_size=200){

  #bin info
  bin_info_all <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(bin_info_all)=c("GeneID","Chr","Start","End","Strand")
  
  batches=floor(length(gnames)/batch_size)    
  
  #start_time <- Sys.time()
  
  for (i in 1:batches){
    #print(i)
    tmpres=lapply(gnames[((i-1)*batch_size+1):(i*batch_size)],function(x) .binByGene(x,merged_gene,bin_size))
    bin_info_all=rbind(bin_info_all,bind_rows(tmpres))
  }
  
  #rest
  tmpres=lapply(gnames[(batches*batch_size+1):length(gnames)],function(x) .binByGene(x,merged_gene,bin_size))
  bin_info_all=rbind(bin_info_all,bind_rows(tmpres))
  
  #end_time <- Sys.time()
  #print(end_time - start_time)
  
  
  out_bin <- new( "binInfo", 
                 bin_id=bin_info_all$GeneID,
                 chrom=bin_info_all$Chr,
                 start=bin_info_all$Start,
                 end=bin_info_all$End,
                 strand=bin_info_all$Strand) 
    
  return(out_bin)

}



.binByGene <- function(gid,gr,bin_size){
  
  bin_info <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(bin_info)=c("GeneID","Chr","Start","End","Strand")
  
  ndf=as.data.frame(gr[[gid]])
  glen=sum(ndf$width)
  strand=as.character(ndf$strand[1])
  chrname=as.character(ndf$seqnames[1])
  exonic_num=nrow(ndf)
  
  #rna coordinate
  rna_coord=c()
  mark_pos=1
  
  if (strand=="+"){
    idx=1:exonic_num
  }else{
    idx=rev(1:exonic_num)
  }
  
  for (j in idx){
    rna_start=mark_pos
    rna_end=mark_pos+ndf[j,"width"]-1
    mark_pos=rna_end+1
    rna_coord=c(rna_coord,rna_start)
    
  }
  
  if (strand=="+"){
    ndf$rna_coord=rna_coord
  }else{
    ndf$rna_coord=rev(rna_coord)
  }
  
  
  #get rna check points
  rna_checkpoint=seq(1,glen,by=bin_size)
  
  #initiate dna check points
  if (strand=="+"){
    dna_checkpoint=c(ndf[1,"start"])
    ndf_idx=c(1)
  }else{
    dna_checkpoint=c(ndf[exonic_num,"end"])
    ndf_idx=c(exonic_num)
  }
  
  #get dna check points
  if (length(rna_checkpoint)>1){    #if glen>bin_size
    
    for (j in 2:length(rna_checkpoint)){
      if (strand=="+"){
        
        if (length(which(ndf$rna_coord==rna_checkpoint[j]))==0){    
          idx=which(order(c(rna_checkpoint[j],ndf$rna_coord))==1)    #insert rna check points to rna_coord
          tmp_idx=idx-1
        } else{
          idx=which(ndf$rna_coord==rna_checkpoint[j])[1]    #rna checkpoints overlap with rna_coord
          tmp_idx=idx
        }
        
        ndf_idx=c(ndf_idx,tmp_idx)
        dna_checkpoint=c(dna_checkpoint,ndf$start[tmp_idx]+(rna_checkpoint[j]-ndf$rna_coord[tmp_idx]))    #start
        
      }else{
        if (length(which(rev(ndf$rna_coord)==rna_checkpoint[j]))==0){
          idx=which(order(c(rna_checkpoint[j],rev(ndf$rna_coord)))==1)
          tmp_idx=exonic_num+1-(idx-1)
        } else{
          idx=which(rev(ndf$rna_coord)==rna_checkpoint[j])[1]    #rna checkpoints overlap with rna_coord
          tmp_idx=exonic_num+1-idx
        }
        
        ndf_idx=c(ndf_idx,tmp_idx)
        dna_checkpoint=c(dna_checkpoint,ndf$end[tmp_idx]-(rna_checkpoint[j]-ndf$rna_coord[tmp_idx]))
      }
    }
    
  }
  
  
  # get bins
  if (strand=="+"){
    # forward strand
    if (length(rna_checkpoint)>1){
      
      for (j in 1:(length(dna_checkpoint)-1)){
        bin_name=paste(gid,"bin",as.character(j),sep="_")
        if (dna_checkpoint[j+1]-dna_checkpoint[j]==bin_size){
          part_start=dna_checkpoint[j]
          part_end=dna_checkpoint[j+1]-1
          bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
        }else{
          #front split part
          part_start=dna_checkpoint[j]
          part_end=ndf[ndf_idx[j],"end"]
          bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
          
          #possible middle part
          if(ndf_idx[j+1]-ndf_idx[j]>1){
            for (k in (ndf_idx[j]+1):(ndf_idx[j+1]-1)){
              part_start=ndf[k,"start"]
              part_end=ndf[k,"end"]
              bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
            }
          }
          
          #possible back split part
          if (ndf[ndf_idx[j+1],"start"]!=dna_checkpoint[j+1]){    #in case back split overlap with next checkpoint
            part_start=ndf[ndf_idx[j+1],"start"]
            part_end=dna_checkpoint[j+1]-1
            bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
            
          }
          
        }
        
      }
      
      
      #fwd last checkpoint
      bin_name=paste(gid,"bin",as.character(length(dna_checkpoint)),sep="_")
      
      if(rev(ndf_idx)[1]==exonic_num){     #checkpoint in the last exonic
        part_start=tail(dna_checkpoint,n=1)
        part_end=ndf[exonic_num,"end"]
        bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
      }else{
        #front split part
        part_start=tail(dna_checkpoint,n=1)
        part_end=ndf[rev(ndf_idx)[1],"end"]
        bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
        
        #rest exonic parts
        for (k in (rev(ndf_idx)[1]+1):exonic_num){
          part_start=ndf[k,"start"]
          part_end=ndf[k,"end"]
          bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
        }
        
      }
      
    }else{
      bin_name=paste(gid,"bin","1",sep="_")    #glen<=bin_size
      idx=1:exonic_num
      for (j in idx){
        part_start=ndf[j,"start"]
        part_end=ndf[j,"end"]
        bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
      }
    }
    
  }else{
    # reverse strand
    if (length(rna_checkpoint)>1){
      
      for (j in 1:(length(dna_checkpoint)-1)){
        bin_name=paste(gid,"bin",as.character(j),sep="_")
        if (dna_checkpoint[j]-dna_checkpoint[j+1]==bin_size){
          part_start=dna_checkpoint[j+1]+1
          part_end=dna_checkpoint[j]
          bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
        }else{
          #front split part
          part_start=ndf[ndf_idx[j],"start"]
          part_end=dna_checkpoint[j]
          bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
          
          #possible middle part
          if(ndf_idx[j]-ndf_idx[j+1]>1){
            for (k in rev((ndf_idx[j+1]+1):(ndf_idx[j]-1))){
              part_start=ndf[k,"start"]
              part_end=ndf[k,"end"]
              bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
            }
          }
          
          #possible back split part
          if (dna_checkpoint[j+1]!=ndf[ndf_idx[j+1],"end"]){
            part_start=dna_checkpoint[j+1]+1
            part_end=ndf[ndf_idx[j+1],"end"]
            bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
            
          }
          
        }
      }
      
      #rev last bin
      bin_name=paste(gid,"bin",as.character(length(dna_checkpoint)),sep="_")
      
      if(rev(ndf_idx)[1]==1){
        part_start=ndf[1,"start"]
        part_end=tail(dna_checkpoint,n=1)
        bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
      }else{
        #front split part
        part_start=ndf[rev(ndf_idx)[1],"start"]
        part_end=tail(dna_checkpoint,n=1)
        bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
        
        #rest exonic parts
        for (k in (rev(ndf_idx)[1]-1):1){
          part_start=ndf[k,"start"]
          part_end=ndf[k,"end"]
          bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
        }
        
      }
      
    }else{
      bin_name=paste(gid,"bin","1",sep="_")  #glen=bin_size
      idx=rev(1:exonic_num)
      for (j in idx){
        part_start=ndf[j,"start"]
        part_end=ndf[j,"end"]
        bin_info[nrow(bin_info) + 1,]=list(bin_name,chrname,part_start,part_end,strand)
      }
    }
    
  }
  
  return(bin_info)
  
}














