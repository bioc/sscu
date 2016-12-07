\name{op_highly}
\alias{op_highly}
\alias{optimal_codons}
\alias{selected_codons}
\title{
  Identify optimal codons by using the highly expressed genes method
}
\description{
  Optimal codons can be defined as codons significantly enriched in the highly expressed genes compared to the lowly expressed genes, or other set of appropriate reference genes. In another word, these codons were favored by translational selection. This function calculate the optimal codon list, thus user could have a general idea of which codons were preferred by selection in the genome.
}
\usage{
  op_highly(high_cds_file = NULL,ref_cds_file = NULL,p_cutoff = 0.01)
}
\arguments{
  \item{high_cds_file}{
	a character vector for the filepath of the highly expressed genes
}
  \item{ref_cds_file}{
	a character vector for the filepath of the reference cds file
}
  \item{p_cutoff}{
	a numeric vector to set the cutoff of p value for the chi.square test, default is set to 0.01
}
}
\details{
  Optimal codons can be defined as codons significantly enriched in the highly expressed genes compared to the lowly expressed genes, or other set of appropriate reference genes. In another word, these codons were favored by translational selection, which was strongest among highly expressed genes. This function calculate the optimal codon list with p-values, thus user could have a general idea of which codons were preferred by selection in the genome.
  
  The argument high_cds_file should specific the path for the highly expressed gene dataset. It is up to the users how to define which dataset of highly expressed genes. Some studies use the expression data, or Nc value to divide genes into highly/lowly sets. Other studies use a specific dataset, such as only including the very highly expressed genes (ribosomal genes).

  The argument ref_cds_file should specific the path for the lowly expressed gene dataset, or any appropriate dataset. In Sharp PM paper (Forces that influence the evolution of codon bias), he used the all gene data set as neutral reference and also get a list of optimal codons. 

  The argument p_cutoff set the cutoff for p values in the chi.square test. Only codons are significantly enriched in the highly expressed genes are marked with + symbol in the ouotput tables. The codons are significantly lower presented in the highly expressed genes are marked with - symbol. The codons are not significantly differently presented compared to the reference dataset are marked with NA symbol.
}
\value{
  a character vector for all the optimal codons is returned
}
\references{
  Sharp PM, Emery LR, Zeng K. 2010. Forces that influence the evolution of codon bias. Philos Trans R Soc Lond B Sci. 365:1203-1212.
}
\author{
  Yu Sun
}
\examples{
# ----------------------------------------------- #
#     Lactobacillus kunkeei example               #
# ----------------------------------------------- #

  # Here is an example to load the data included in the sscu package
  op_highly(high_cds_file=system.file("sequences/L_kunkeei_highly.ffn",package="sscu"),ref_cds_file=system.file("sequences/L_kunkeei_genome_cds.ffn",package="sscu"))

  # if you want to set the p value cutoff as 0.05
  op_highly(high_cds_file=system.file("sequences/L_kunkeei_highly.ffn",package="sscu"),ref_cds_file=system.file("sequences/L_kunkeei_genome_cds.ffn",package="sscu"),p_cutoff=0.05)

  # if you want to load your own data, you just specify the file path for your input as these examples
  # op_highly(high_cds_file = "/home/yu/Data/codon_usage/bee_endosymbionts/sharp_40_highly_dataset/Bin2.ffn",ref_cds_file = "/home/yu/Data/codon_usage/bee_endosymbionts/cds_filtered/Bin2.ffn")

}

