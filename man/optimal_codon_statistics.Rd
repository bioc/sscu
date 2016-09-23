\name{optimal_codon_statistics}
\alias{optimal_codon_statistics}
\alias{optimal_codons_table}
\title{
  statistical table for the optimal codons
}
\description{
  Optimal codons are defined as codons significantly enriched in the highly expressed genes compared to the lowly expressed genes, or other set of appropriate reference genes. In another word, these codons were favored by translational selection. This function calculate the optimal codon list with p-values, thus user could have a general idea of which codons were preferred by selection in the genome.
}
\usage{
  optimal_codon_statistics(high_cds_file = NULL,ref_cds_file = NULL,p_cutoff = 0.05)
}
\arguments{
  \item{high_cds_file}{
	a character vector for the filepath of the highly expressed genes
}
  \item{ref_cds_file}{
	a character vector for the filepath of the reference cds file
}
  \item{p_cutoff}{
	a numeric vector to set the cutoff of p value for the chi.square test, default is set to 0.05
}
}
\details{
  Optimal codons are defined as codons significantly enriched in the highly expressed genes compared to the lowly expressed genes, or other set of appropriate reference genes. In another word, these codons were favored by translational selection, which was strongest among highly expressed genes. This function calculate the optimal codon list with p-values, thus user could have a general idea of which codons were preferred by selection in the genome.
  
  The argument high_cds_file should specific the path for the highly expressed gene dataset. It is up to the users how to define which dataset of highly expressed genes. Some studies use the expression data, or Nc value to divide genes into highly/lowly sets. Other studies use a specific dataset, such as only including the very highly expressed genes (ribosomal genes).

  The argument ref_cds_file should specific the path for the lowly expressed gene dataset, or any appropriate dataset. In Sharp PM paper (Forces that influence the evolution of codon bias), he used the all gene data set as neutral reference and also get a list of optimal codons. 

  The argument p_cutoff set the cutoff for p values in the chi.square test. Only codons are significantly enriched in the highly expressed genes are marked with + symbol in the ouotput tables. The codons are significantly lower presented in the highly expressed genes are marked with - symbol. The codons are not significantly differently presented compared to the reference dataset are marked with NA symbol.

  The function also output the rscu value for the high expressed dataset and reference dataset.
}
\value{
  a dataframe is returned
  \item{rscu_high}{rscu value for the highly expressed dataset}
  \item{rscu_ref}{rscu value for the reference dataset}
  \item{high_No_codon}{number of codons found in the highly expressed dataset}
  \item{high_expect_No_codon}{number of expected codons in the highly expressed dataset}
  \item{ref_No_codon}{number of codons found in the reference dataset}
  \item{ref_expect_No_codon}{number of expected codons in the reference dataset}
  \item{p_value}{p value for the chi.square test}
  \item{symbol}{codons are significantly enriched in the highly expressed genes are marked with +; codons are significantly lower presented in the highly expressed genes are marked with -; codons are not significantly differently presented compared to the reference dataset are marked with NA }
}
\author{
  Yu Sun
}
\seealso{
  \code{\link{uco}} in seqinr library for rscu calculation.  
}
\examples{
# ----------------------------------------------- #
#     Lactobacillus kunkeei example               #
# ----------------------------------------------- #

  # Here is an example to load the data included in the sscu package
  optimal_codon_statistics(high_cds_file=system.file("sequences/L_kunkeei_highly.ffn",package="sscu"),ref_cds_file=system.file("sequences/L_kunkeei_genome_cds.ffn",package="sscu"))

  # if you want to set the p value cutoff as 0.01
  optimal_codon_statistics(high_cds_file=system.file("sequences/L_kunkeei_highly.ffn",package="sscu"),ref_cds_file=system.file("sequences/L_kunkeei_genome_cds.ffn",package="sscu"),p_cutoff=0.01)

  # if you want to load your own data, you just specify the file path for your input as these examples
  # optimal_codon_statistics(high_cds_file = "/home/yu/Data/codon_usage/bee_endosymbionts/sharp_40_highly_dataset/Bin2.ffn",ref_cds_file = "/home/yu/Data/codon_usage/bee_endosymbionts/cds_filtered/Bin2.ffn")

}

