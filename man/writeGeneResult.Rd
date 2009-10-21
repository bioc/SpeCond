\name{writeGeneResult}
\alias{writeGeneResult}
\title{Write a condition-specific analysis result text file}
\description{
 \code{writeGeneResult} produces a text file containing the list of gene, if they have been detected as tissue-specific or not (S/N), for how many tissues in total, how many tissue as up-regulated, how many tissue as down-regulated, in which tissues for up-regulated and down-regulated.
}
\usage{
writeGeneResult(L.specific.result, file.name.result.gene = "gene_summary_result.txt",
gene.names = NULL)
}
\arguments{
 \item{L.specific.result}{the \code{L.specific.result} list of the included in the result of the main \code{SpeCond} function: generalResult$specificResult$L.specific.result}
  \item{file.name.result.gene}{the name of the produced file containing the list of specific genes an thier specific detection}
  \item{gene.names}{vector of gene's names to select a suset of genes. The default is NULL, all genes from the input matrix in \code{SpeCond} function are used}
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{SpeCond},\code{getProfile},\code{writeSpeCondResult},\code{writeUniqueProfileSpecificResult}}
\examples{
library(SpeCond)
data(expressionSpeCondExample)
##Perform the condition specific detection analysis with SpeCond()
generalResult=SpeCond(expressionSpeCondExample, param.detection=NULL, multitest.correction.method="BY",
  prefix.file="E", print.hist.pv=TRUE, fit1=NULL, fit2=NULL, specificOutlierStep1=NULL)
specificResult=generalResult$specificResult
##write the result file
writeGeneResult(specificResult$L.specific.result, file.name.result.gene=
  "Example_gene_summary_result.txt", gene.names=rownames(expressionSpeCondExample)[1:10])
}
