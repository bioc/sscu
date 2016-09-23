\name{sscu-package}
\alias{sscu-package}
\docType{package}
\title{Strength of Selected Codon Usage}
\description{
  The package can calculate the indexes for selective stength in codon usage in bacteria species. (1) The package can calculate the strength of selected codon usage bias (sscu, also named as s_index) based on Paul Sharp's method. The method take into account of background mutation rate, and focus only on four pairs of codons with universal translational advantages in all bacterial species. Thus the sscu index is comparable among different species. (2) Translational accuracy selection can be inferred from Akashi's test. The test tabulating all codons into four categories with the feature as conserved/variable amino acids and optimal/non-optimal codons. (3) Optimal codon lists (selected codons) can be calculated by either op_highly function (by using the highly expressed genes compared with all genes to identify optimal codons biased used in the highly expressed genes), or op_corre_CodonW/op_corre_NCprime function (by correlative method developed by Hershberg & Petrov). Users will have a list of optimal codons for further analysis, such as input to the Akashi's test. (4) The detailed codon usage information, such as RSCU value, number of optimal codons in the highly/all gene set, as well as the genomic gc3 value, can be calculate by the optimal_codon_statistics and genomic_gc3 function. (5) Furthermore, we added one test function proportion_index in the package. The function focus on the proportion of optimal codon against its corresponding non-optimal codons for the the four and six codon boxes. 
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
  Sharp PM, Emery LR, Zeng K. 2010. Forces that influence the evolution of codon bias. Philos Trans R Soc Lond B Sci. 365:1203-1212.
  Hershberg R, Petrov DA. 2009. General rules for optimal codon choice. Plos Genet. 5:e1001115.
  Akashi H. Synonymous codon usage in Drosophila melanogaster: natural selection and translational accuracy. Genetics 1994 Mar;136(3):927-35.
  http://drummond.openwetware.org/Akashi's_Test.html
  Novembre JA. 2002. Accounting for background necleotide composition when measuring codon usage bias. Mol Biol Evol. 19: 1390-1394.
  https://github.com/jnovembre/ENCprime
  http://codonw.sourceforge.net/
}

