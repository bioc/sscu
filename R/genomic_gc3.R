genomic_gc3 <- function(inputfile){
  
  genomic_cds <- Biostrings::readDNAStringSet(inputfile)
  genomic_cat <- seqinr::s2c(as.character(BiocGenerics::unlist(genomic_cds)))
  genomic_gc3 <- seqinr::GC3(genomic_cat)    
  
  return(genomic_gc3)
}