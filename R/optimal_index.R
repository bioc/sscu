optimal_index <- function(high_cds_file=NULL,genomic_cds_file=NULL){
  
  optimal_df <- sscu::optimal_codons(high_cds_file = high_cds_file,
                               ref_cds_file = genomic_cds_file)
  
  genomic_cds <- Biostrings::readDNAStringSet(genomic_cds_file)
  
  sscu_index <- sscu::s_index(high_cds_file=high_cds_file,
                        genomic_cds_file=genomic_cds_file)
  
  genomic_cat <- seqinr::s2c(as.character(BiocGenerics::unlist(genomic_cds)))
  genomic_gc3 <- seqinr::GC3(genomic_cat)
  
  optimal_only <- optimal_df[optimal_df$symbol == '+',]
  optimal_four <- data.frame()
  optimal_six <- data.frame()
  
  four_codons_gc <- c( "GTC", 
                       "CCC", 
                       "ACC", 
                       "GCC", 
                       "GGC")
  six_codons_gc <- c( "CGC",  
                      "CTC", 
                      "TCC")
  
  #if(genomic_gc3 >= 0.5){
  #    stop("Genomic GC3 must be lower than 50%!")
  #}else{
  if((TRUE %in% (optimal_only$codon %in% four_codons_gc)) && 
     (TRUE %in% (optimal_only$codon %in% six_codons_gc))){
    optimal_four <- optimal_only[optimal_only$codon %in% four_codons_gc,]
    optimal_six <- optimal_only[optimal_only$codon %in% six_codons_gc,]
    
    corres_four_3rd <- chartr("GCAT","ATGC",substr(optimal_four$codon,3,3))
    corres_four_codon <- paste(
      substr(optimal_four$codon,1,2),corres_four_3rd, sep=''
      ) 
    corres_four <- optimal_df[optimal_df$codon %in% corres_four_codon,]
    
    optimal_four_codon_nr <- optimal_four$high_No_codon
    corres_four_codon_nr <- corres_four$high_No_codon
    
    corres_six_3rd <- chartr("GCAT","ATGC",substr(optimal_six$codon,3,3))
    corres_six_codon <- paste(
      substr(optimal_six$codon,1,2),corres_six_3rd, sep=''
      ) 
    corres_six <- optimal_df[optimal_df$codon %in% corres_six_codon,]
    
    optimal_six_codon_nr <- optimal_six$high_No_codon
    corres_six_codon_nr <- corres_six$high_No_codon
    
    optimal_46_codon_nr <- c(optimal_four_codon_nr, optimal_six_codon_nr)
    corres_46_codon_nr <- c(corres_four_codon_nr, corres_six_codon_nr)
    total_codon46 <- sum(optimal_46_codon_nr,corres_46_codon_nr)
    total_codon4 <- sum(optimal_four_codon_nr,corres_four_codon_nr)
    total_codon6 <- sum(optimal_six_codon_nr,corres_six_codon_nr)
    
    p46 <- (optimal_46_codon_nr/(optimal_46_codon_nr+corres_46_codon_nr))
    p4 <- (optimal_four_codon_nr/(optimal_four_codon_nr+corres_four_codon_nr))
    p6 <- (optimal_six_codon_nr/(optimal_six_codon_nr+corres_six_codon_nr))
    k <- (1-genomic_gc3)/genomic_gc3
    sscu46_values <- vector()
    sscu4_values <- vector()
    sscu6_values <- vector()
    
    for(i in 1:length(p46)){
      sscu46_values[i] <- log(p46[i] * k / (1-p46[i]))
    }
    for(i in 1:length(p4)){
      sscu4_values[i] <- log(p4[i] * k / (1-p4[i]))
    }
    for(i in 1:length(p6)){
      sscu6_values[i] <- log(p6[i] * k / (1-p6[i]))
    }
    
    sscu46 <- sum((sscu46_values*(optimal_46_codon_nr+corres_46_codon_nr)))/
      total_codon46
    sscu4 <- sum((sscu4_values*(optimal_four_codon_nr+corres_four_codon_nr)))/
      total_codon4
    sscu6 <- sum((sscu6_values*(optimal_six_codon_nr+corres_six_codon_nr)))/
      total_codon6
    optimal_four$sscu <- sscu4_values
    optimal_six$sscu <- sscu6_values
    sscu46_sscu_ratio <- sscu46/sscu_index
    #result <- list(sscu46=sscu46,sscu4=sscu4,sscu6=sscu6,optimal_four=optimal_four,corres_four=corres_four,optimal_six=optimal_six,corres_six=corres_six)
    result <- list(
      optimal_mutation_index=sscu46,
      s_index=sscu_index,
      ratio=sscu46_sscu_ratio
      )
    return(result)
  }
  else{
    return(NA)
    #stop("at least one c-ending optimal codon must in the four codon box
    #and six codon box")
  }
  #}
}

#optimal_index(high_cds_file="/home/yu/Data/codon_usage/bee_endosymbionts/sharp_40_highly_dataset/Bbif.ffn",genomic_cds_file="/home/yu/Data/codon_usage/bee_endosymbionts/cds_filtered/Bbif.ffn")
