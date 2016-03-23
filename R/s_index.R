s_index <- function(high_cds_file=NULL,genomic_cds_file=NULL,gc3=NULL){
  
  if(is.null(high_cds_file))
    stop("You must give filepath for the highly expressed CDS!")
  else if(!file.exists(high_cds_file))
    stop("Check your filepath for highly expressed CDS. 
             The file is not exist!")
  else if(is.null(genomic_cds_file)){
    if(is.numeric(gc3))
      genomic_gc3 <- gc3
    else  
      stop("You must either give filepath for the genomic CDS file, 
                 or give the genomic gc3 content!")
  }
  else {
    genomic_cds <- Biostrings::readDNAStringSet(genomic_cds_file)
    genomic_cat <- seqinr::s2c(as.character(BiocGenerics::unlist(genomic_cds)))
    genomic_gc3 <- seqinr::GC3(genomic_cat)
  }
  
  
  
  highly_cds <- Biostrings::readDNAStringSet(high_cds_file)
  highly_cat <- BiocGenerics::unlist(highly_cds)
  nuc_freq <- Biostrings::trinucleotideFrequency(highly_cat,step=3)
  
  tac_num <- nuc_freq[attr(nuc_freq,"names")=="TAC"]
  tat_num <- nuc_freq[attr(nuc_freq,"names")=="TAT"]
  aac_num <- nuc_freq[attr(nuc_freq,"names")=="AAC"]
  aat_num <- nuc_freq[attr(nuc_freq,"names")=="AAT"]
  atc_num <- nuc_freq[attr(nuc_freq,"names")=="ATC"]
  att_num <- nuc_freq[attr(nuc_freq,"names")=="ATT"]
  ttc_num <- nuc_freq[attr(nuc_freq,"names")=="TTC"]
  ttt_num <- nuc_freq[attr(nuc_freq,"names")=="TTT"]
  total_num <- tac_num+tat_num+aac_num+aat_num+atc_num+att_num+ttc_num+ttt_num
  
  s_tyr <- log((tac_num/(tac_num+tat_num))*((1-genomic_gc3)/genomic_gc3)*1/
                 (1-(tac_num/(tac_num+tat_num))))
  s_asn <- log((aac_num/(aac_num+aat_num))*((1-genomic_gc3)/genomic_gc3)*1/
                 (1-(aac_num/(aac_num+aat_num))))
  s_ile <- log((atc_num/(atc_num+att_num))*((1-genomic_gc3)/genomic_gc3)*1/
                 (1-(atc_num/(atc_num+att_num))))
  s_phe <- log((ttc_num/(ttc_num+ttt_num))*((1-genomic_gc3)/genomic_gc3)*1/
                 (1-(ttc_num/(ttc_num+ttt_num))))
  
  s <- (s_tyr*(tac_num+tat_num)+s_asn*(aac_num+aat_num)+s_ile*
          (atc_num+att_num)+s_phe*(ttc_num+ttt_num))/total_num
  attr(s,"names") <- NULL
  
  return(s)
  
}
