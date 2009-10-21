\name{getSpecificResult}
\alias{getSpecificResult}
\title{Detect the condition-specific genes for the second step on the SpeCond procedure}
\description{
Perform the second detection step of the SpeCond procedure. Use the second fitting (without prior and ignoring the outliers detected in the first step) of the gene expression value with a mixture of normal distribution results and a set of rules to detect the outliers. It returns the outliers detected as specifically expressed for each gene.
}
\usage{
getSpecificResult(expressionMatrix, fit2 = NULL, param.detection = NULL, 
specificOutlierStep1 = NULL, multitest.correction.method = "BY", 
prefix.file = NULL, print.hist.pv = FALSE)
}
\arguments{
  \item{expressionMatrix}{the gene expression matrix (genes * conditions)}
  \item{fit2}{The result of \code{fitNoPriorWithExclusion} containing the parameter of the mixture normal model of the expression data ignoring the outliers detected in the first step of the procedure}
  \item{param.detection}{the parameter for the detection, a vector with the names ("per","md","mlk","rsd","pv") or the second row of the matrix obtained by \code{getDefaultParameter} or \code{createParameterMatrix}}
  \item{specificOutlierStep1}{the list of outliers detected by the first step procedure, result of the \code{getSpecificOtuliersStep1} function}
  \item{multitest.correction.method}{the multitest correction method. The default is "BY", for the possible values see \code{p.adjust}}
  \item{prefix.file}{a prefix added to the generated file. The default is NULL but as to be set. It is useful to change the prefix when you perform a new analysis. As you may want to compare the results with different parameters set}
  \item{print.hist.pv}{ a logical (TRUE/FALSE) whether to print in a pdf file the (non-adjusted) p-value histogram; the default is FALSE}
}
\value{
An object of class \code{sp_list}
  \item{prefix.file}{the prefix used for this analysis. It will be used by default in the function \code{getGeneHtmlPage}}
  \item{fit}{the fitting parameters used by the detection i.e. the argument fit2}
  \item{param.detection}{the parameters used for the detection i.e. the argument parm.detection}
\item{L.specific.result}{Full detection results (It will be used by the \code{getFullHtmlSpeCondResult}). This list contains 7 attributes:}
\item{M.specific.all}{matrix of 0: not selective, 1: selective up-regulated, -1: selective down-regulated; same dimensions as the input expression values matrix}
\item{M.specific}{same as M.specific.all but reduced to the specific genes. NULL if no gene has been detected as specific}
\item{M.specific.sum.row}{Number of conditions in which the gene is specific}
\item{M.specific.sum.column}{Number of specific genes by conditions}
\item{L.pv}{list of all genes with a matrix of conditions and the corresponding p-values (if the gene is specific)}
\item{specific}{vector of size the number of genes with "Not specific" or "Specific" according to the specificity of the gene}
\item{L.condition.specific.id}{list of the specific genes with a vector of column numbers (condition ids), for which the gene is specific}
  \item{L.null}{a list of vectors of 1 and 0 representing the null distribution. The length of the vector for each gene corresponds to the number of normal distributions fitting the gene expression value. The list is sorted as the gene order in the input expression matrix}
  \item{L.mlk}{a list of vectors containing the min log-likelihood computed between normal distribution components. NULL if the mixture model of the gene is composed of only one component or if the proportion of all components is superior to the per.2 parameter}
  \item{L.rsd}{a list of vectors containing the standard deviation ratio computed between normal distribution components. NULL if the mixture model of the gene is composed of only one component}
  \item{identic.row.ids}{row number(s) from the initial input matrix which contain identical values for all conditions. These rows are not considered in the analysis}
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{fitNoPriorwithExclusion},\code{SpeCond},\code{getSpecificResult}}
\examples{
library(SpeCond)
data(expressionSpeCondExample)
##Perform the SpeCond analysis step by step
param.detection=getDefaultParameter()
param.detection

fit1=fitPrior(expressionSpeCondExample, param.detection=param.detection)

specificOutlierStep1=getSpecificOutliersStep1(expressionSpeCondExample, fit=fit1$fit1,
 param.detection, multitest.correction.method="BY", prefix.file="run1_Step1",
 print.hist.pv=FALSE)

fit2=fitNoPriorWithExclusion(expressionSpeCondExample, specificOutlierStep1=specificOutlierStep1,
  param.detection=param.detection)

specificResult=getSpecificResult(expressionSpeCondExample, fit=fit2,
  specificOutlierStep1=specificOutlierStep1, param.detection, 
  multitest.correction.method="BY", prefix.file="run1_Step2",
  print.hist.pv=FALSE)

}
