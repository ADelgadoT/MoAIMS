plotGOF <- function(fit_obj=NULL,
                    bin_count_obj=NULL,
                    sample_id=NULL,
                    k=3 )
{    
    
    YFreq <- table(bin_count_obj@ip)
    YVal <- as.numeric(names(YFreq))
    
    
    XFreq <- table( bin_count_obj@input)
    XVal <- as.numeric(names(XFreq))
    
        
    N <- sum(YFreq)
    
    
    a <- fit_obj@size_b
    pi0 <- 1-fit_obj@pi_s
    muEst <- fit_obj@mu_b
    bEst <- a / muEst
    
    b <- fit_obj@size_s
    c <- fit_obj@size_s/fit_obj@mu_s
    
    p1 <- fit_obj@pi_1s
    b1 <- fit_obj@size_1s
    c1 <- fit_obj@size_1s/fit_obj@mu_1s
    b2 <- fit_obj@size_2s
    c2 <- fit_obj@size_2s/fit_obj@mu_2s
    
    # simulate one sample signal
    
    probZ_1S <- c( pi0, (1-pi0) )
    set.seed(12345)
    Zsim_1S <- sample( c(0,1), size=N, prob=probZ_1S, replace=TRUE )
    set.seed(12345)
    Ysim_1S <- rnbinom( N, a, bEst/(bEst+1) )
    nS_1S <- length(which( Zsim_1S==1 ))
    Ysim_1S[ Zsim_1S==1 ] <- Ysim_1S[ Zsim_1S==1 ] + rnbinom( nS_1S, b, c/(c+1) ) + k
    
    #YsimZ0Freq <- table( Ysim_1S[ Zsim_1S==0 ] )
    #YsimZ0Val <- as.numeric(names(YsimZ0Freq))
    
    YsimFreq_1S <- table(Ysim_1S)
    YsimVal_1S <- as.numeric(names(YsimFreq_1S))
    
    # simulation two sample signal
    
    probZ_2S <- c( pi0, (1-pi0)*p1, (1-pi0)*(1-p1) )
    set.seed(12345)
    Zsim_2S <- sample( c(0,1,2), size=N, prob=probZ_2S, replace=TRUE )
    set.seed(12345)
    Ysim_2S <- rnbinom( N, a, bEst/(bEst+1) )   
    nS_2S_1 <- length(which(Zsim_2S==1))
    nS_2S_2 <- length(which(Zsim_2S==2))
    Ysim_2S[ Zsim_2S==1 ] <- Ysim_2S[ Zsim_2S==1 ] + rnbinom( nS_2S_1, b1, c1/(c1+1) ) + k
    Ysim_2S[ Zsim_2S==2 ] <- Ysim_2S[ Zsim_2S==2 ] + rnbinom( nS_2S_2, b2, c2/(c2+1) ) + k 
    
    YsimZ0Freq <- table( Ysim_2S[ Zsim_2S==0 ] )
    YsimZ0Val <- as.numeric(names(YsimZ0Freq))    

        
    YsimFreq_2S <- table(Ysim_2S)
    YsimVal_2S <- as.numeric(names(YsimFreq_2S))
    
    # draw GOF plot
    
    x_max=max(c(YVal,YsimVal_1S,YsimVal_2S))
    
    
    df_Y=data.frame(YFreq,Type="Real data")
  	df_Y[,1]=as.numeric(df_Y[,1])
  	colnames(df_Y)[1:2]=c("val","freq")
  
  	df_sim_1s=data.frame(YsimFreq_1S,Type="Simulation_1S")
  	df_sim_1s[,1]=as.numeric(df_sim_1s[,1])
  	colnames(df_sim_1s)[1:2]=c("val","freq")
  	
  	df_sim_2s=data.frame(YsimFreq_2S,Type="Simulation_2S")
  	df_sim_2s[,1]=as.numeric(df_sim_2s[,1])
  	colnames(df_sim_2s)[1:2]=c("val","freq")
  
  	df_plot=rbind(df_Y,df_sim_1s,df_sim_2s)
  	
  	png(paste("GOF_",sample_id,sep=""),width = 500, height = 300) 
  	fig=ggplot(df_plot,aes(x=val+1,y=freq,color=Type))+geom_line(alpha=0.8)+
  		scale_x_continuous(trans = 'log10')+
  		scale_y_continuous(trans = 'log10')+
  		scale_color_manual(values=c("black", "orangered","royalblue"))+
  		xlab("IP bin count")+
  		ylab("Frequency")+
  		theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=12))
 		
 		
 	print(fig)
 	
  	dev.off()
    

}








