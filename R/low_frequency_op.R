low_frequency_op <- function(high_cds_file=NULL,genomic_cds_file=NULL,p_cutoff=0.01){
  
  # high_cds_file <- system.file("sequences/Gvag_highly.ffn",package="sscu")
  # genomic_cds_file <- system.file("sequences/Gvag_genome_cds.ffn",package="sscu")
  # p_cutoff <- 0.01
  
  optimal_df <- sscu::op_highly_stats(high_cds_file = high_cds_file,
                               ref_cds_file = genomic_cds_file,
                               p_cutoff = p_cutoff)
  
  # genomic_cds <- Biostrings::readDNAStringSet(genomic_cds_file)
  # 
  # sscu_index <- sscu::s_index(high_cds_file=high_cds_file,
  #                       genomic_cds_file=genomic_cds_file)
  # 
  # genomic_cat <- seqinr::s2c(as.character(BiocGenerics::unlist(genomic_cds)))
  # genomic_gc3 <- seqinr::GC3(genomic_cat)
  
  optimal_only <- optimal_df[optimal_df$symbol == '+',]
  
  low_frequency_optimal <- optimal_only[optimal_only$rscu_high < 0.7,]
  
  corres_3rd <- chartr("GCAT","ATGC",substr(low_frequency_optimal$codon,3,3))
  corres_codon <- paste(
    substr(low_frequency_optimal$codon,1,2),corres_3rd, sep=''
  ) 
  corres <- optimal_df[optimal_df$codon %in% corres_codon,]
  
  # filter1, the RSCU value for low frequency optimal codon is lower than the corresponding codon
  filter1 <- low_frequency_optimal$rscu_high < corres$rscu_high
  # filter2, the corresponding codon is indeen an non-optimal codon
  filter2 <- !(corres_codon %in% optimal_only$codon)
  low_frequency_optimal_filter <- low_frequency_optimal[(filter1 & filter2), ]
  corres_filter <- corres[(filter1 & filter2), ]
  
  result <- list(low_frequency_optimal_codons = low_frequency_optimal_filter,
                 corresponding_codons = corres_filter)
  return(result)
}
