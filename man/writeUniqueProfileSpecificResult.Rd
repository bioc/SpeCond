\name{writeUniqueProfileSpecificResult}
\alias{writeUniqueProfileSpecificResult}
\title{Write the specific profiles from the SpeCond analysis}
\description{
 Produces a text file with the unique specific profiles among the conditions detected by the SpeCond analysis.
}
\usage{
writeUniqueProfileSpecificResult(L.specific.result, file.name.unique.profile =
 "specific.unique_profile.txt", full.list.gene = FALSE)
}
\arguments{
 \item{L.specific.result}{ the \code{L.specific.result} list of the included in the result of the main \code{SpeCond} function: generalResult$specificResult$L.specific.result}
  \item{file.name.unique.profile}{ the name of the produced file containing the gene's profiles}
  \item{full.list.gene}{ If TRUE, the last column correspond to the gene's names which have the profile described in the row}
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{SpeCond}, \code{getProfile}, \code{writeSpeCondResult}, \code{writeGeneResult}}
\examples{
library(SpeCond)
data(expressionSpeCondExample)
##Perform the condition specific detection analysis with SpeCond()
generalResult=SpeCond(expressionSpeCondExample, 
  param.detection=NULL, multitest.correction.method="BY", prefix.file="E", 
  print.hist.pv=TRUE, fit1=NULL, fit2=NULL, specificOutlierStep1=NULL)
specificResult=generalResult$specificResult

##write the result file
writeUniqueProfileSpecificResult(L.specific.result=specificResult$L.specific.result,
  file.name.unique.profile="Example_specific_unique_profile.txt", full.list.gene=FALSE)
}
