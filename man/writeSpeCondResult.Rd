\name{writeSpeCondResult}
\alias{writeSpeCondResult}
\title{Write in text files the main result of the SpeCond function}
\description{\code{writeSpeCondResult} produces three text files: 
- The table of the gene detected as specific and in which condition they are specific (0: no specific, 1: specific up-regulated, -1:specific down-regulated). The default name is file.name.profile="specific_profile.txt".
- The list of the specific genes. The default name is: "list_specific_probeset.txt".
- The table of the unique specific profiles detected. The default name is: "specific_unique_profile.txt".
}
\usage{
writeSpeCondResult(L.specific.result, file.name.profile =
 "specific_profile.txt", file.specific.gene = "list_specific_gene.txt", 
 file.name.unique.profile = "specific_unique_profile.txt")
}
\arguments{
  \item{L.specific.result}{The \code{L.specific.result} list of the included in the result of the main \code{SpeCond} function: generalResult$specificResult$L.specific.result}
  \item{file.name.profile}{The name of the produced file containing the gene's profiles}
  \item{file.specific.gene}{The name of the produced file containing the list of the specific genes}
  \item{file.name.unique.profile}{The name of the produced file containing the unique gene's profiles}
}
\author{Florence Cavalli, florence@ebi.ac.uk
}
\seealso{ \code{SpeCond}, \code{getProfile}, \code{writeUniqueProfileSpecifcResult}, 
\code{writeGeneResult}}
\examples{
library(SpeCond)
data(expressionSpeCondExample)
##Perform the condition specific detection analysis with SpeCond()
generalResult=SpeCond(expressionSpeCondExample, param.detection=NULL, 
 multitest.correction.method="BY", prefix.file="E", print.hist.pv=TRUE, fit1=NULL, 
 fit2=NULL, specificOutlierStep1=NULL)
 specificResult=generalResult$specificResult

##write the SpeCond results files
 writeSpeCondResult(specificResult$L.specific.result,file.name.profile=
 "Example_specific_profile.txt", file.specific.gene="Example_list_specific_gene.txt",
  file.name.unique.profile="Example_specific_unique_profile.txt")
}
