op_highly_stats <- function(high_cds_file=NULL,ref_cds_file=NULL,p_cutoff=0.01){
  
  calc_frequency <- function(cds_file=NULL){
    
    cds <- Biostrings::readDNAStringSet(cds_file)
    cat <- unlist(cds)
    nuc_freq <- Biostrings::trinucleotideFrequency(cat,step=3)
    aa_freq <- Biostrings::alphabetFrequency(
      Biostrings::translate(cat,if.fuzzy.codon="X")
    )
    nuc_exp <- c()
    
    for(i in 1:64){
      codon <- attr(nuc_freq[i],"names")
      aa <- Biostrings::GENETIC_CODE[attr(GENETIC_CODE,"names")==codon]
      aa_num <- aa_freq[attr(aa_freq,"names")==aa]
      table_gc <- table(GENETIC_CODE)
      codon_expect <- aa_num/table_gc[attr(table_gc,"dimnames")$GENETIC_CODE==aa]
      attr(codon_expect,"names") <- codon
      nuc_exp <- c(nuc_exp,codon_expect)
    }
    return(list(nuc_freq=nuc_freq,nuc_exp=nuc_exp))
  }
  
  reference_result <- calc_frequency(cds_file=ref_cds_file)
  high_result <- calc_frequency(cds_file=high_cds_file)
  
  p <- numeric()
  symbol <- character()
  rscu_high <- numeric()
  rscu_ref <- numeric()
  codon <- character()
  aa_for_codon <- character()
  
  for(i in 1:64){
    codon[i] <- attr(high_result[[1]][i],"names")
    aa_for_codon[i] <- GENETIC_CODE[attr(GENETIC_CODE,"names")==codon[i]]
    x <- matrix(c(high_result[[1]][i],high_result[[2]][i],
                  reference_result[[1]][i],reference_result[[2]][i]),ncol=2)
    p[i] <- stats::chisq.test(x,simulate.p.value = TRUE)$p.value  
    
    rscu_high[i] <- high_result[[1]][i]/high_result[[2]][i]
    rscu_ref[i] <- reference_result[[1]][i]/reference_result[[2]][i]
    
    if(is.na(rscu_high[i])){
      symbol[i] <- 'NA'
    }
    else{
      if(rscu_high[i]>rscu_ref[i]){
        if(p[i]<p_cutoff){
          symbol[i] <- '+'
        }
        else{
          symbol[i] <- 'NA'
        }
      }
      else{
        if(p[i]<p_cutoff){
          symbol[i] <- '-'
        }
        else{
          symbol[i] <- 'NA'
        }
      }
    }
  }
  
  out_list <- data.frame(codon=codon,
                         aa=aa_for_codon,
                         rscu_high=round(rscu_high,digits=2),
                         rscu_ref=round(rscu_ref,digits=2),
                         high_No_codon=high_result[[1]],
                         high_expect_No_codon=round(high_result[[2]],digits=0),
                         ref_No_codon=reference_result[[1]],
                         ref_expect_No_codon=round(reference_result[[2]],digits=0),
                         p_value=round(p,digits=3),
                         symbol=symbol)
  codon_level <- c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG",
                   "TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG",
                   "CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
                   "CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG",
                   "ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG",
                   "AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG",
                   "GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG",
                   "GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG")
  out_list$codon  <- factor(out_list$codon, levels = codon_level)
  sorted_out <- out_list[order(out_list$codon),]
  return(sorted_out)
}

#df <- optimal_codons(high_cds_file = "/home/yu/Data/codon_usage/bee_endosymbionts/sharp_40_highly_dataset/Bin2.ffn",ref_cds_file = "/home/yu/Data/codon_usage/bee_endosymbionts/cds_filtered/Bin2.ffn",p_cutoff = 0.05)
