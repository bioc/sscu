\name{low_frequency_op}
\alias{low_frequency_op}
\title{
  the function identify low frequency optimal codons
}
\description{
  The function low_frequency_op identify the low frequency optimal codons. Based on the previous study, in some species, the optimal codons identified by the op_highly function has bit strange patterns: the optimal codons do have lower prequency than the non-optimal codons. It occurs most in the mutation-shifting species such as G. vaginalis. The function can identify these low frequency optimal codons.   
}
\usage{
low_frequency_op(high_cds_file = NULL, genomic_cds_file = NULL, p_cutoff=0.01)
}
\arguments{
  \item{high_cds_file}{
	a character vector for the filepath of the highly expressed genes
  }
  \item{genomic_cds_file}{
	a character vector for the filepath of the whole genome cds file
  }
  \item{p_cutoff}{
	a numeric vector to set the cutoff of p value for the chi.square test, default is set to 0.01
  }
}
\details{

  The function low_frequency_op identify the low frequency optimal codons. Based on the previous study (Sun Y et al. Switches in genomic GC content drives shifts of optimal codons under sustained selection on synonymous sites. Genome Biol Evol. 2016 Aug 18.), in some species, the optimal codons identified by the op_highly function has bit strange patterns: the optimal codons do have lower prequency than the non-optimal codons. It occurs most in the mutation-shifting species such as G. vaginalis. The function can identify these low frequency optimal codons. 
  
  The function first calculated all the optimal codons statistics by the function op_highly_stats in the package, then tried to find the low frequency optimal codons with the following settings: 1) it has to be an optimal codon 2) The RSCU value for the optimal codon is lower than 0.7 (this is the quite arbitrary setting) 3) the RSCU value for the optimal codon is lower than the corresponding non-optimal codons, which has the same first two nucleotide as the optimal codon but the third position experience the point mutation transition. With this detailed filters, we can identify the low frequency optimal codons in the given genome.

  The argument high_cds_file should specific the path for the highly expressed gene dataset. It is up to the users how to define which dataset of highly expressed genes. Some studies use the expression data, or Nc value to divide genes into highly/lowly sets. Other studies use a specific dataset, such as only including the very highly expressed genes (ribosomal genes). In the example, I used the ribosomal genes as the representative for the highly expressed genes. 

  The arguments, genomic_cds_file, is used to calculate the optimal codons and statistics for highly and all genes.

  Same as the op_highly and op_highly_stats, the p value cutoff was set as 0.01.
}
\value{
  a list is returned
  \item{low_frequency_optimal_codons}{a dataframe with statistics for the low frequency optimal codons}
  \item{corresponding_codons}{a dataframe with statistics for the corresponding codons}
}
\references{
  Sun Y et al. Switches in genomic GC content drives shifts of optimal codons under sustained selection on synonymous sites. Genome Biol Evol. 2016 Aug 18.
  
  unpublished paper from Yu Sun
}
\author{
Yu Sun
}
\seealso{
  the op_highly and op_highly_stats function in the same package
}
\examples{
# ----------------------------------------------- #
#     Lactobacillus kunkeei example               #
# ----------------------------------------------- #

  # Here is an example to load the data included in the sscu package
  # input the two multifasta files to calculate sscu 
   low_frequency_op(high_cds_file=system.file("sequences/Gvag_highly.ffn",package="sscu"),genomic_cds_file=system.file("sequences/Gvag_genome_cds.ffn",package="sscu"))

  # if you want to load your own data, you just specify the file path for your input as these examples
  # low_frequency_op(high_cds_file="/home/yu/Data/codon_usage/bee_endosymbionts/sharp_40_highly_dataset/Bin2.ffn",genomic_cds_file="/home/yu/Data/codon_usage/bee_endosymbionts/cds_filtered/Bin2.ffn")

}

