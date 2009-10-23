\name{getFullHtmlSpeCondResult}
\alias{getFullHtmlSpeCondResult}
\title{ Visualisation function of the SpeCond analysis results}
\description{
 \code{getFullHtmlSpeCondResult} generates a full result html page.
}
\usage{
getFullHtmlSpeCondResult(SpeCondResult=NULL, L.specific.result = NULL, 
param.detection = NULL, page.name = "SpeCond_result",
page.title = "Condition-specific analysis results", prefix.file = NULL, 
outdir="General_Result", force=TRUE, sort.condition = "all", 
gene.page.info=NULL, heatmap.profile = TRUE, heatmap.expression = FALSE, 
heatmap.unique.profile = FALSE, expressionMatrix = NULL)
}
\arguments{
  \item{SpeCondResult}{the \code{sp_list} class object result of the \code{SpeCond} functions}
  \item{L.specific.result}{List of results present in the \code{sp_list} class specificResult object, see \code{SpeCond} or \code{getSpecificResult} functions}
  \item{param.detection}{The parameter matrix used by the SpeCond detection procedure}
  \item{page.name}{The name of the result html page. The default is "SpeCond\_result"}
  \item{page.title}{The title of the result html page. The default is "Condition-specific analysis results"}
  \item{prefix.file}{a prefix added to the generated file(s) and the \code{outdir} directory name to linked them to the full result html page, by default is NULL, the \code{prefix.file}  attribute of the \code{SpeCondResult}  is used. It is useful to change the prefix when you create a new result page. As you may want to get results with different parameter sets and plots so using a different \code{SpeCondResult} or \code{L.specific.result} objects}
  \item{outdir}{ the name of the directory in which the generated files will be created. The default is "General_result"}
  \item{force}{ if TRUE the \code{outdir} will be cleaned then written if it already exists. The default is TRUE}
  \item{sort.condition}{If the condition must sorted in the barplot presented the number of specific genes by condition. Can table the values: positive", "negative", "all": the conditions are sorted respectively by the number of specific genes detected as up-regulated, down-regulated or both}
\item{gene.page.info}{ the result of the \code{getGeneHtmlPage} function. Enables the creation of links between this full result page and the single result pages created by the previous function. The default is "NULL"; no links are created}
  \item{heatmap.profile}{TRUE/FALSE, to print or not a heatmap showing the specific profile of the genes. The default is FALSE}
  \item{heatmap.expression}{TRUE/FALSE, to print or not a heatmap showing the expression of the genes. The default is FALSE}
  \item{heatmap.unique.profile}{TRUE/FALSE, to print or not a heatmap showing the unique specific profile. The default is FALSE}
  \item{expressionMatrix}{Must not be NULL if heatmap.expression=TRUE, must be the same as the input expression matrix. The default is NULL}
}
\details{
Either \code{SpeCondResult} or \code{L.specific.result} can be specified to use this function. If you use \code{L.specific.result} you ahve tp define \code{prefix.file}.
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{\code{getGeneHtmlPage}}
\examples{
library(SpeCond)
data(expressionSpeCondExample)
##Perform the condition specific detection analysis with SpeCond()
generalResult=SpeCond(expressionSpeCondExample, param.detection=NULL, multitest.correction.method="BY",
  prefix.file="E", print.hist.pv=TRUE, fit1=NULL, fit2=NULL, specificOutlierStep1=NULL)
specificResult=generalResult$specificResult
##Produce the general html page results
getFullHtmlSpeCondResult(SpeCondResult=generalResult, param.detection=	
  specificResult$param.detection, page.name="Example_SpeCond_results",
  page.title="Tissue specific results", sort.condition="all", heatmap.profile=TRUE,
  heatmap.expression=FALSE, heatmap.unique.profile=FALSE, expressionMatrix=expressionSpeCondExample)

##Produce the Gene html page results for the first 20 genes using the specificResult object to be able to link 
## these pages to the table result in the general html page
specificResult=generalResult$specificResult
genePageInfo=getGeneHtmlPage(expressionSpeCondExample, specificResult, name.index.html=
  "index_example_SpeCond_Results.html",outdir="Single_result_pages_example", gene.html.ids=c(1:20))

##Produce the general html page results
getFullHtmlSpeCondResult(L.specific.result=specificResult$L.specific.result,
  param.detection=specificResult$param.detection, page.name="Example_SpeCond_results2",
  page.title="Tissue specific results", prefix.file="S2", sort.condition="all", heatmap.profile=TRUE,
  heatmap.expression=FALSE, heatmap.unique.profile=FALSE, expressionMatrix=Mexp,  
  gene.page.info=genePageInfo)
}
