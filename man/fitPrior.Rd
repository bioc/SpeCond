\name{fitPrior}
\alias{fitPrior}
\title{ Fit the expression value profiles with a mixture of normal components }
\description{
\code{firPrior} performs a clustering of expression values for each gene profile using the mclust function. This results to a mixture of normal distribution components (from 1 to 3 components) fitting the expression values.
}
\usage{
fitPrior(expressionMatrix, param.detection = NULL, lambda = 1, 
beta = 6, evaluation.lambda.beta = FALSE)
}
\arguments{
  \item{expressionMatrix}{the expression value matrix, genes*conditions}
  \item{param.detection}{the matrix of parameters as obtained by \code{getDefaultParameter} or \code{createParamterMatrix}. It must contain positive values for "lambda" and "beta". If NULL, the function \code{getDefaultParameter} will be used}
  \item{lambda}{positive value, it influences the choice of models by affecting the selection of one, two or three normal distributions, thus introducing some weight on the effect of number of parameters to be defined. The default is 1, the model uses the BIC value taking into account the log-likelihood value}
  \item{beta}{positive value, it  influences the prior applied during the determination of the variance of the normal distributions. It is important for this fitting since it allows the model to capture isolated outliers. The default value is 6}
  \item{evaluation.lambda.beta}{if TRUE, an extra attribute will be return indicating for how many gene the lambda and beta parameters change the number of normal component chosen to fit the expression values}
}
\value{
 \item{fit1}{list of the gene as first attributes, for each gene a list of three attributes:}
 \item{G}{number of normal components fitting the data}
 \item{NorMixParam}{the parameters of each normal component: proportion, mean and standard deviation for the gene}
 \item{classification}{the normal component id in  which the expression values of the gene are attributed}
}
\details{
 If \code{evaluation.lambda.beta} is TRUE an additional attribute \code{G.lambda.beta.effect} is returned. It is a matrix presenting the number of time the values of G (number of normal components for a particular gene) has changed between \code{lambda}=0 and the \code{lambda.1} value and between \code{beta}=0 and the \code{beta.1} value.
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{fitNoPriorWithExclusion}, \code{SpeCond}}
\examples{
library(SpeCond)
data(expressionSpeCondExample)
##Perform the SpeCond analysis step by step
param.detection=getDefaultParameter()
param.detection

fit1=fitPrior(expressionSpeCondExample, param.detection=param.detection)

##then use getSpecificOutliersStep1(), fitNoPriorWithExclusion() and getSpecificResult()
}
