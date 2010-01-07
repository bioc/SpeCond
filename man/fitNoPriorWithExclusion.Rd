\name{fitNoPriorWithExclusion}
\alias{fitNoPriorWithExclusion}
\title{ Fit the expression values profile with a mixture of normal components ignoring outliers}
\description{
\code{firPrior} performs a clustering of expression values for each gene profile using the mclust function ignoring the outliers (detected by the first step of the SpeCond prcedure) present in the SpecificOutlierStep1 argument . This results to a mixture of normal distribution components (from 1 to 3 components) fitting the expression values.
}

\usage{
fitNoPriorWithExclusion(expressionMatrix, specificOutlierStep1 = FALSE, 
param.detection = NULL, lambda = 1, beta = 0)
}
\arguments{
  \item{expressionMatrix}{ the expression value matrix, genes*conditions}
  \item{specificOutlierStep1}{the list of outliers detected by the first step procedure, result of the
 \code{getSpecificOutliersStep1} function or an attritube of the \code{SpeCond} result object. These outliers won't be taken into account for the mixture normal modelling performed by this function}
  \item{param.detection}{ the matrix of parameters as obtained by \code{getDefaultParameter} or
 \code{createParamterMatrix}. It must contain positive values for "lambda" and "beta". If NULL, the function \code{getDefaultParameter} will be used}
  \item{lambda}{ positive value, it influences the choice of models by affecting the selection of one, two or three normal distributions, thus introducing some weight on the effect of number of parameters to be defined. The default is 1, the model uses the BIC value taking into account the log-likelihood value}
  \item{beta}{Should be equal to 0; prior is put on the variance determination of the normal distribution}
}
\value{
 \item{fit2}{list of the gene as first attributes, for each gene a list of three attributes:}
 \item{G}{number of normal components fitting the data}
 \item{NorMixParam}{the parameters of each normal component: proportion, mean and standard deviation for the gene}
 \item{classification}{the normal component id in  which the expression values of the gene are attributed}
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{ \code{fitPrior}, \code{SpeCond}}
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

fit2=fitNoPriorWithExclusion(expressionSpeCondExample, 
specificOutlierStep1=specificOutlierStep1,
param.detection=param.detection)

##then use getSpecificResult()
}
