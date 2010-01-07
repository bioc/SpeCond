\name{getSpecificOutliersStep1}
\alias{getSpecificOutliersStep1}
\title{Detect the condition-specific as outliers in for the first step on the SpeCond procedure}
\description{
Perform the first detection step of the SpeCond procedure. Use the fitting of the gene expression value with a mixture of normal distribution results and a set of rules to detect the outliers. It returns the outliers detected as specifically expressed for each gene.
}
\usage{
getSpecificOutliersStep1(expressionMatrix, fit1 = NULL, 
param.detection = NULL, multitest.correction.method = "BY", 
prefix.file = NULL, print.hist.pv = FALSE)
}
\arguments{
  \item{expressionMatrix}{the gene expression matrix (genes * conditions)}
  \item{fit1}{the result of \code{fitPrior} containing the parameter of the mixture normal model of the expression data}
  \item{param.detection}{the parameter for the detection, a vector with the names ("per","md","mlk","rsd","pv") or the first row of the matrix obtained by \code{getDefaultParameter} or 
  \code{createParameterMatrix}}
  \item{multitest.correction.method}{the multitest correction method. The default is "BY", for the possible values see \code{p.adjust}}
  \item{prefix.file}{a prefix added to the generated file. The default is NULL but has to be set. It is useful to change the prefix when you perform a new analysis. As you may want to compare the results with different parameters set.}
  \item{print.hist.pv}{to print in a pdf file the (non-adjusted) p-value histogram}
}
\value{
A list of size the number of rows (genes) in the expressionMatrix. If the gene has outlier expression, the column number of this outlier is stored, NULL if not.
}
\details{
Frist essential method to obtain the matrix of expression value from your ExpressionSet to apply the SpeCond procedure step by step using the following function \code{fitPrior}, \code{fitNoPriorwithExclusion}, \code{getSpecificOutliersStep1}, \code{getSpecificResult}. The returned matrix will be the expressionMatrix argument of the above function
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{fitPrior},\code{SpeCond},\code{getSpecificResult}}

\examples{
library(SpeCond)
data(expressionSpeCondExample)
##Perform the SpeCond analysis step by step
param.detection=getDefaultParameter()
param.detection

fit1=fitPrior(expressionSpeCondExample, param.detection=param.detection)

specificOutlierStep1=getSpecificOutliersStep1(expressionSpeCondExample, 
 fit=fit1$fit1, param.detection, multitest.correction.method="BY",
 prefix.file="run1_Step1", print.hist.pv=FALSE)

##then use fitNoPriorWithExclusion() and getSpecificResult()
}
