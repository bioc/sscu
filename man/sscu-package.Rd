\name{sscu-package}
\alias{sscu-package}
\docType{package}
\title{Strength of Selected Codon Usage}
\description{
  The package sscu (Strength of Selected Codon Usage) calculates the selective profile in codon usage in bacteria species. First of all, the package calculate the s index based on Paul Sharp's method . The method take into account of background mutation rate, and focus only on codons with universal translational advantages in all bacterial species. Thus the s index is comparable among different species. In addition, detailed optimal codons (selected codons) information can be calculated by optimal_codons function, so the users will have a more accurate selective scheme for each codons. Furthermore, we added one more function optimal_index in the package. The function has similar mathematical formula as s index, but focus on the estimates the amount of GC-ending optimal codon for the highly expressed genes in the four and six codon boxes. The function takes into account of background mutation rate, and it is comparable with the s index. However, since the set of GC-ending optimal codons are likely to be different among different species, the index can not be compared among different species.
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{sscu}
\packageIndices{sscu}
}

\author{
Yu Sun

Maintainer: Yu Sun <sunyu1357@gmail.com>
}

\references{
  Sharp PM, Bailes E, Grocock RJ, Peden JF, Sockett RE (2005). Variation in the strength of selected codon usage bias among bacteria. Nucleic Acids Research.
}

