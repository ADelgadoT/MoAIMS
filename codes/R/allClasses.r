setClass( Class="binInfo",
          representation=representation(
            bin_id="character",
            chrom="character",
            start="numeric",
            end="numeric",
            strand="character"
          )
)

setClass( Class="binCount",
          representation=representation(
            bin_id="character",
            ip="numeric",
            input="numeric"
          )
)


setClass( Class="fitObj",
          representation=representation(
            mu_b="numeric", 
            size_b="numeric",  
                   
            pi_s="numeric",        
            mu_s="numeric",
            size_s="numeric",
			md1="list",

            pi_1s="numeric",
            mu_1s="numeric",
            size_1s="numeric",
            mu_2s="numeric",  
            size_2s="numeric",  
            md2="list",
            
            val_BIC_1s="numeric",
            val_BIC_2s="numeric"
            
            
          )
)


setClass( Class="sigRegion",
          representation=representation(
            sig_region_list="data.frame"
            
          )
)


# MOSAiCS fitObj

setClass( Class="MosaicsFitEst",
    representation=representation(
        pi0="numeric",
        a="numeric",
        betaEst="numeric",
        muEst="numeric",    
        pNfit="list",
        b="numeric",
        c="numeric",
        p1="numeric",
        b1="numeric",
        c1="numeric",
        b2="numeric",
        c2="numeric",
        inputTrunc="numeric",
        analysisType="character"
    )
)





