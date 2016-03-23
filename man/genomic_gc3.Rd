\name{genomic_gc3}
\alias{genomic_gc3}
\title{
  genomic gc3 for an multifasta genomic file
}
\description{
  The function calculates the genomic gc3 for an multifasta genomic CDS file. The function first concatenated all 
  the CDS sequences in the file into one long CDS string, than calculated the gc3 from the GC3 function in seqinr package.
  You can also use the function to calculate the gc3 for a single gene, or a set of genes, depends what content you put in 
  the input file.
}
\usage{
  genomic_gc3(inputfile)
}
\arguments{
  \item{inputfile}{
	a character vector for the filepath of the whole genome cds file
}
}
\details{
  The function calculates the genomic gc3 for an multifasta genomic CDS file. The function first concatenated all 
  the CDS sequences in the file into one long CDS string, than calculated the gc3 from the GC3 function in seqinr package.
  You can also use the function to calculate the gc3 for a single gene, or a set of genes, depends what content you put in 
  the input file. The result can be used as input for the s_index calculation.
}
\value{
  a numeric vector genomic_gc3 is returned
}
\author{
Yu Sun
}
\seealso{
  \code{\link{GC3}} in seqinr library
}
\examples{
# ----------------------------------------------- #
#     Lactobacillus kunkeei example               #
# ----------------------------------------------- #

  # Here is an example to calculate the genomic gc3
  # input the one multifasta files to calculate genomic gc3 
  genomic_gc3(system.file("sequences/L_kunkeei_genome_cds.ffn",package="sscu"))

}

