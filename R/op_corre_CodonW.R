op_corre_CodonW <- function(genomic_cds_file=NULL, correspondence_file=NULL){

genomic_cds <- read.fasta(file = genomic_cds_file)
#genomic_cds <- read.fasta(file = genomic_cds_file)
num_seq <- length(genomic_cds)
rscu_df <- data.frame()

aa_codons <- list(
    c(5:8), #thr
    c(21:24), #pro
    c(37:40), #ala
    c(41:44), #gly
    c(45:48), #val
    c(10,12,53:56), #ser
    c(9,11,25:28), #arg
    c(29:32,61,63), #leu
    c(13,14,16), #ile
    c(1,3), #lys
    c(2,4), #asn
    c(17,19), #gln
    c(18,20), #his
    c(33,35), #glu
    c(34,36), #asp
    c(50,52), #tyr
    c(58,60), #cys
    c(62,64) #phe
)

logic_vector <- c()
for(i in 1:num_seq){
    rscu <- uco(genomic_cds[[i]], index="rscu")
    codons <- uco(genomic_cds[[i]], index="eff")
    codon_for_aa <- c()
    for(h in 1:18){
        codon_for_aa <- c(codon_for_aa, sum(codons[aa_codons[[h]]]))
    }
    codon_family_num <- sum(codon_for_aa > 0)
    codon_count <- sum(codons)
    logic_vector <- c(logic_vector, codon_count > 50 && codon_family_num > 9)
    rscu_df <- rbind(rscu_df,as.vector(rscu))
}

nc_df <- read.table(file = correspondence_file, header=TRUE)
#nc_df <- read.table(file = '/home/yu/Data/codon_usage/bee_endosymbionts/correspondence_analysis/Gvag/Gvag.out', header=TRUE)

correlations_nc <- c()
correlations_nc_p <- c()
for(j in 1:64){
    cortest_nc <- cor.test(rscu_df[,j][logic_vector],as.numeric(as.character(nc_df$Nc[logic_vector])),
                          method="spearman", exact=FALSE)
    correlation_nc <- cortest_nc$estimate
    correlation_nc_p <- cortest_nc$p.value
    correlations_nc <- c(correlations_nc, correlation_nc)
    correlations_nc_p <- c(correlations_nc_p, correlation_nc_p)
}

names(correlations_nc) <- c("aaa","aac","aag","aat","aca","acc","acg","act",
                            "aga","agc","agg","agt","ata","atc","atg","att",
                            "caa","cac","cag","cat","cca","ccc","ccg","cct",
                            "cga","cgc","cgg","cgt","cta","ctc","ctg","ctt",
                            "gaa","gac","gag","gat","gca","gcc","gcg","gct",
                            "gga","ggc","ggg","ggt","gta","gtc","gtg","gtt",
                            "taa","tac","tag","tat","tca","tcc","tcg","tct",
                            "tga","tgc","tgg","tgt","tta","ttc","ttg","ttt")


optimal_codons_nc <- c()
for(i in 1:18){
    logic_vec <- correlations_nc_p[aa_codons[[i]]] < 0.05/sum(correlations_nc_p[aa_codons[[i]]])
    min_corres <- min(correlations_nc[aa_codons[[i]]][logic_vec])
    optimal_codon_nc <- correlations_nc[aa_codons[[i]]][match(min_corres, correlations_nc[aa_codons[[i]]])]
    optimal_codons_nc <- c(optimal_codons_nc, optimal_codon_nc)
}
optimal_codon_chr_nc <- attributes(na.omit(optimal_codons_nc))$name

return (optimal_codon_chr_nc);
}

# correlative_test_codonw(genomic_cds_file = '/home/yu/Data/codon_usage/bee_endosymbionts/cds_filtered/Gvag.ffn', correspondence_file = '/home/yu/Data/codon_usage/bee_endosymbionts/correspondence_analysis/Gvag/Gvag.out')
