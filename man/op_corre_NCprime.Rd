\name{op_corre_NCprime}
\alias{op_corre_NCprime}
\title{
  Identify optimal codons by using the correlative method from Hershberg & Petrov, the input file is from NCprime
}
\description{
  The function identify the optimal codons based on the correlative method from Hershberg & Petrov.
  This method take the whole genome into consideration, and predict the optimal codons by making the correlation between 
  the frequency of each codon within each gene and the overall codon bias (Nc or Nc').
  The input file include the output file from the program ENCprime (to get the Nc and Nc' value), 
  and the genomic cds file (to get the codon usage information for each gene). 
}
\usage{
  op_corre_NCprime(genomic_cds_file=NULL, nc_file=NULL)
}
\arguments{
  \item{genomic_cds_file}{
	a character vector for the filepath of the whole genome cds file
  }
  \item{nc_file}{
	a character vector for the filepath of the correspondence file from the CodonW program
  }
}
\details{
  The function identify the optimal codons based on the correlative method from Hershberg & Petrov.
  This method take the whole genome into consideration, and predict the optimal codons by making the correlation between 
  the frequency of each codon within each gene and the overall codon bias (Nc or Nc').
  The input file include the output file from the program ENCprime (to get the Nc and Nc' value), 
  and the genomic cds file (to get the codon usage information for each gene). 
  
  For further details regard how to use ENCprime, you can refer to the site https://github.com/jnovembre/ENCprime.
  Note, you must input the same genomic cds file to ENCprime and to the op_corre_NCprime funtion, so that the order and number of
  genes are consistent in the two files. 
  
}
\value{
  a character vector for all the optimal codons is returned
}
\references{
  Hershberg R, Petrov DA. 2009. General rules for optimal codon choice. Plos Genet. 5:e1001115.
  Novembre JA. 2002. Accounting for background necleotide composition when measuring codon usage bias. Mol Biol Evol. 19: 1390-1394.
  https://github.com/jnovembre/ENCprime
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
   op_corre_NCprime(genomic_cds_file=system.file("sequences/LbDelBA1_genome_cds.ffn",package="sscu"),nc_file=system.file("correlative_test/LbDelBA1.NCprime",package="sscu"))

}

