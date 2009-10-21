\name{getMatrixFromExpressionSet}
\alias{getMatrixFromExpressionSet}
\title{Obtain the expression matrix from an ExpressionSet object}
\description{
\code{getMatrixFromExpressionSet} method returns an matrix of expression values from an ExpressionSet object. It takes into consideration the need of summarizing the samples values by conditions to perform the SpeCond analysis 
}
\usage{
getMatrixFromExpressionSet(expSet, condition.factor = NULL, 
                       condition.method = c("mean", "median","max"))
}
\arguments{
  \item{expSet}{an ExpressionSet object}
  \item{condition.factor}{a factor object of length equal to the number of columns (samples) of the ExpressionSet object specifying which sample(s) belong to which condition (condition.factor levels); can be extracted from the phenoData}
  \item{condition.method}{the method (mean, median or max) to summarise the samples by conditions (defined by the condition.factor vector)}
}
\details{
For each level of the condition.factor, the expression values of the ExpressionSet object are computed using the condition.method method. 
If there is only one sample for a condition the expression value is not changed
if condition.factor is NULL, the expression matrix of the ExpressionSet object will simply be extracted using exprs()}
\value{
A matrix of expression values of size (number of row in the ExpressionSet * number of level of the condition.factor) 
}
\references{
Biobase
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{SpeCond}}
\examples{
library(SpeCond)
data(expSetSpeCondExample)
expSetSpeCondExample
f_Tissues=factor(paste("Tissue_",rep(1:32,each=2),sep=""))
f_Tissues
Mexp=getMatrixFromExpressionSet(expSetSpeCondExample,
  condition.factor=f_Tissues,condition.method="mean")
## or
Mexp=getMatrixFromExpressionSet(expSetSpeCondExample,
  condition.factor=expSetSpeCondExample$Tissue,condition.method="mean")
}
