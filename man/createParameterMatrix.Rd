\name{createParameterMatrix}
\alias{createParameterMatrix}
\title{ Create or Modify the SpeCond argument parameters}
\description{\code{createParameterMatrix} creates and/or modifies the param.detection matrix used as argument in the SpeCond function. If parm.detection is NULL the param.detection matrix used is the one containing the default parameter values, as obtained by \code{getDefaultParameter}. The remaining arguments enable to change the values of the param.detection matrix.
}
\usage{
createParameterMatrix(param.detection = NULL, beta.1 = NULL, beta.2 = NULL,
lambda.1 = NULL, lambda.2 = NULL, per.1 = NULL, per.2 = NULL, md.1 = NULL, 
md.2 = NULL, mlk.1 = NULL, mlk.2 = NULL, rsd.1 = NULL, rsd.2 = NULL, 
pv.1 = NULL, pv.2 = NULL)
}
\arguments{
  \item{param.detection}{ a matrix of 2 rows and 7 columns as the result of \code{getDefaultParameter}}
  \item{beta.1}{ Influences the prior applied during the determination of the variance of the normal distributions. It is necessary in the first fitting step to allow the model to capture isolated outliers.}
  \item{beta.2}{The normal use of SpeCond does not prior on Step2: must be set to 0}
  \item{lambda.1}{ Influences the choice of models by affecting the selection of one, two or three normal distributions, thus introducing some weight on the effect of number of parameters to be defined. The default is 1, the model uses the BIC value taking into account the log-likelihood value. }
  \item{lambda.2}{ Same as \code{lambda.1} for the second step of the SpeCond function}
  \item{per.1}{ percentage threshold: this is the percentage of conditions that can be detected as specific. As per increases a larger number of expression values per genes can be identified as specific. The default is 0.1}
  \item{per.2}{ percentage threshold: This is the final percentage of condition that can be detected as specific. As per increases a larger number of expression values per genes can be identified as specific. The default is 0.3}
  \item{md.1}{  median difference: this is the minimum value between the median values of two mixture components that is allowed to identify one of them as representing outliers, i.e. possibly not part of the null distribution. This corresponds to a biological fact; specific expression that corresponds to noise should not be detected as specific}
  \item{md.2}{ Same as \code{md.1} for the second step of the SpeCond function. For consistency should be equal to \code{md.1}}
  \item{mlk.1}{minimum log-likelihood: enables the identification of clusters of conditions that are well separated from the others in the model. If the gene mlk value>mlk, the mixture component can be detected as outlier (i.e. not part of the null distribution)}
  \item{mlk.2}{same as \code{mlk.2} for the second step of the SpeCond function}
  \item{rsd.1}{minimum of standard deviation ratio: enables the identification of clusters of conditions that are extremely spread out compared to the distribution clustering of most expression values. If the gene rsd values< rsd the mixture component can be detected as outlier (i.e. not part of the null distribution)}
  \item{rsd.2}{same as \code{rsd.1} for the second step of the SpeCond function}
  \item{pv.1}{p-value threshold to detect a condition as specific for a given gene}
  \item{pv.2}{same as \code{pv.1} for the second step of the SpeCond function, for consistency should be equal to pv.1}
}

\value{
  param.detection: a matrix of 2 row and 7 columns. The rows "Step1 "and "Step2" correspond respectively to the first and second set of parameters for the SpeCond function. The parameters (columns) are: lambda, beta, per, md, mlk, rsd. See the \code{createParameterMatrix} documentation for more details about the parameters.
} 
\author{Florence Cavalli, florence@ebi.ac.uk}
\section{Warning}{
The SpeCond code is based on:
beta.2=0
md.1=md.2
per.1<=per.2
pv.1=pv.2
}
\seealso{getDefaultParameter}
\examples{
##Get the default parameters and changing the mlk.1 value to 10:
param.detection2=createParameterMatrix(mlk.1=10)
param.detection2
## Modify param.detection2 with mlk.1 value to 15 and rsd.2 value to 0.2
param.detection2B=createParameterMatrix(param.detection=param.detection2,
mlk.1=10, rsd.2=0.2) 
param.detection2B
}
