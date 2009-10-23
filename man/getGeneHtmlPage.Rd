\name{getGeneHtmlPage}
\alias{getGeneHtmlPage}
\title{ Visualise for each gene the condition-specific detection result from SpeCond}
\description{
\code{getGeneHtmlPage} generates html results pages for a set of genes as well as an index page. The index allows to navigate between the gene result pages.
}
\usage{
getGeneHtmlPage(expressionMatrix, specificResult, name.index.html = "index.html",
prefix.file = NULL, outdir="Single_result_pages", force=TRUE, gene.html = NULL, 
gene.html.ids = c(1:10))
}
\arguments{
  \item{expressionMatrix}{the matrix of expression values initially used}
  \item{specificResult}{the \code{sp_list} class object result of the \code{getSpecificProbeset} function}
  \item{name.index.html}{ the name of the html index, by default is index.html}
  \item{prefix.file}{ a prefix added to the generated file(s) and \code{outdir} directory name to linked to the index file. The default is NULL, the \code{prefix.file}  attribute of \code{specificResult}  is used}  
  \item{outdir}{ the name of the directory in which the generated files will be created. The default is "Single_result_pages"}
  \item{force}{ if TRUE the \code{outdir} will be cleaned then written if it already exists. The default is TRUE}
  \item{gene.html}{a vector of gene names for which you want to create html pages, same as the row names of the expressionMatrix object. The default is NULL (the values of the \emph{gene.html.ids} argument will be used)}
  \item{gene.html.ids}{a vector of  integer corresponding to the row numbers in the expressionMatrix object of the genes for which you want to create html pages. The default is the 10 first rows (or the number of row of the expressionMatrix if inferior to 10)}
}
\details{
The main file \code{name.index.html} is created in the current directory. The result page(s) to which it points are created in the \code{outdir} directory.
 If both \code{gene.html} and \code{gene.html.ids} are set to NULL, the gene html pages for every gene in the expressionMatrix object will be generated 
It is useful to change the prefix when you create a new index as well as changing the \code{name.index.html} value. As you may want to get index with the same genes but different parameters set and plots so using a different specificResult object. It is possible to use \code{gene.html} or \code{gene.html.ids} to select a list of gene.
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{getFullHtmlSpeCondResult}}
\examples{
library(SpeCond)
data(expressionSpeCondExample)
##Perform the condition specific detection analysis with SpeCond()
generalResult=SpeCond(expressionSpeCondExample, param.detection=NULL, multitest.correction.method="BY",
  prefix.file="E", print.hist.pv=TRUE, fit1=NULL, fit2=NULL, specificOutlierStep1=NULL)
specificResult=generalResult$specificResult
##Produce the Gene html page results for the first 20 genes using the specificResult object
genePageInfo=getGeneHtmlPage(expressionSpeCondExample, specificResult, name.index.html=
  "index_example_SpeCond_Results.html", outdir="Single_result_pages_dir", gene.html.ids=c(1:20))
}
