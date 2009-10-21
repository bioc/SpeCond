\name{expSetSpeCondExample}
\Rdversion{1.1}
\alias{expSetSpeCondExample}
\docType{data}
\title{
An ExpressionSet example object used in the SpeCond package}
\description{
\code{expSetSpeCondExample} is an ExpressionSet example object used as an example for the SpeCond package. This ExpressionSet only contains an expression matrix and the phenoData. This object has only the purpose of illustrating how SpeCond can be used with an ExpressionSet input object. 
}
\usage{data(expSetSpeCondExample)}
\format{
  The format is:
Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	2 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:2] "Tissue names" "Experience number"
  .. .. ..@ data             :'data.frame':	64 obs. of  2 variables:
  .. .. .. ..$ Tissue: Factor w/ 32 levels "Adrenal_cortex",..: 23 23 5 5 1 1 16 16 32 32 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:64] "S_1" "S_2" "S_3" "S_4" ...
  .. .. .. ..$ Exp   : Factor w/ 2 levels "Exp1","Exp2": 1 2 1 2 1 2 1 2 1 2 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:64] "S_1" "S_2" "S_3" "S_4" ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{getMatrixFromExpressionSet}}
\examples{
data(expSetSpeCondExample)
}
\keyword{datasets}
