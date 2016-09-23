\name{op_corre_CodonW}
\alias{op_corre_CodonW}
\title{
  Identify optimal codons by using the correlative method from Hershberg & Petrov, the input file is from CodonW
}
\description{
  The function identify the optimal codons based on the correlative method from Hershberg & Petrov.
  This method take the whole genome into consideration, and predict the optimal codons by making the correlation between 
  the frequency of each codon within each gene and the overall codon bias (Nc or Nc').
  The input file include the correspondence analysis output file from the program CodonW (to get the Nc value), 
  and the genomic cds file (to get the codon usage information for each gene). 
}
\usage{
  op_corre_CodonW(genomic_cds_file=NULL, correspondence_file=NULL)
}
\arguments{
  \item{genomic_cds_file}{
	a character vector for the filepath of the whole genome cds file
  }
  \item{correspondence_file}{
	a character vector for the filepath of the correspondence file from the CodonW program
  }
}
\details{
  The function identify the optimal codons based on the correlative method from Hershberg & Petrov.
  This method take the whole genome into consideration, and predict the optimal codons by making the correlation between 
  the frequency of each codon within each gene and the overall codon bias (Nc or Nc').
  The input file include the correspondence analysis output file from the program CodonW (to get the Nc value), 
  and the genomic cds file (to get the codon usage information for each gene). 
  
  For further details regard how to use CodonW, you can refer to the site http://codonw.sourceforge.net/.
  Note, you must input the same genomic cds file to CodonW and to the op_corre_CodonW funtion, so that the order and number of
  genes are consistent in the files. 
  
}
\value{
  a character vector for all the optimal codons is returned
}
\references{
  Hershberg R, Petrov DA. 2009. General rules for optimal codon choice. Plos Genet. 5:e1001115.
  http://codonw.sourceforge.net/
}
\author{
Yu Sun
}
\examples{
# ----------------------------------------------- #
#     Lactobacillus kunkeei example               #
# ----------------------------------------------- #

  # Here is an example to calculate the genomic gc3
  # input the one multifasta files to calculate genomic gc3 
   op_corre_CodonW(genomic_cds_file=system.file("sequences/Gvag_genome_cds.ffn",package="sscu"),correspondence_file=system.file("correlative_test/Gvag.codonw",package="sscu"))

}

