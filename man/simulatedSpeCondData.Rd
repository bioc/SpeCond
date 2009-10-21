\name{simulatedSpeCondData}
\Rdversion{1.1}
\alias{simulatedSpeCondData}
\docType{data}
\title{
An example of  simulated expression matrix used in the SpeCond package
}
\description{
\code{simulatedSpeCondData} is a expression value matrix used as an example for the SpeCond package. The expression values were randomly generated from three different normal distributions.
}
\usage{data(simulatedSpeCondData)}
\format{ A matrix of 600 rows and 30 columns}
\details{
The default expression values for each probeset is randomly generated from a normal distribution of mean=7 and sd=0.6. The probesets 1 to 100 have specific expression values for the conditions 10, 20 and 30 coming from a normal distribution of mean=11, sd=0.5. The probesets 200 to 300 have specific expression values for the conditions 9, 18 and 27 coming from a normal distribution of mean=13, sd=0.4.
This data set is used to show the ipmportance and the effect of the paramters in the SpeCond detection. See the SpCond vignette for more detailsy
}
\examples{
data(simulatedSpeCondData)
}
\keyword{datasets}
