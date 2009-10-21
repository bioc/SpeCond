\name{getDefaultParameter}
\alias{getDefaultParameter}
\title{ Get the default parameter to use SpeCond function}
\description{
This function returns the matrix param.detection default argument for the SpeCond function
}
\usage{
getDefaultParameter()
}
\value{
  param.detection: a matrix of 2 row and 7 columns. The rows "Step1 "and "Step2" correspond respectively to the first and second set of parameters for the SpeCond function. The parameters (columns) are: lambda, beta, per, md, mlk, rsd. See the createParameterMatrix documentation for more details about the parameters.
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{\code{createParameterMatrix}}
\examples{
  param.detection=getDefaultParameter()
  param.detection
  }
