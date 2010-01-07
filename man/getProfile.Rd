\name{getProfile}
\alias{getProfile}
\title{ Create the condition-specific profile of specific matrix result from SpeCond}
\description{
\code{getProfile} converts a matrix of 0,1,-1 values in a matrix of one columns. Each row is transformed to a character chain of the values separated by comma.
}
\usage{
getProfile(M.specific)
}
\arguments{
 \item{M.specific}{Is a matrix result present in the SpeCond object result: 
 		      generalResult\$specificResult\$L.specific.result\$M.specific}
}
\value{
 \item{M.specific.profile }{a matrix of number of row as the M.specific matrix x 2 columns. The first column "profile" is the profile: character chain of the values in M.specific separated by commas. The second column of the 2 columns: "sum.row" is the number of condition in which the genes is specific (up or down regulated)}
  \item{M.specific.profile.unique}{a matrix of number of unique profile * number of conditions. The columns order is the same as M.specific}
  \item{M.specific.profile.table}{a matrix of number of unique profile *2. The columns are: profile, nb.gene. The first column is the profile: character chain of the unique rows in M.specific separated by commas. The second column is the number of genes (rows) from M.specific which have this profile}
}
\author{Florence Cavalli, florence@ebi.ac.uk}
\seealso{\code{SpeCond}, \code{writeSpeCondResult}, \code{writeUniqueProfileSpecifcResult},
\code{writeGeneResult}}

\examples{
library(SpeCond)
data(expressionSpeCondExample)
dim(expressionSpeCondExample)
##Perform the condition specific detection analysis with SpeCond()
generalResult=SpeCond(expressionSpeCondExample, param.detection=NULL, 
  multitest.correction.method="BY", prefix.file="E", print.hist.pv=TRUE, fit1=NULL,
  fit2=NULL, specificOutlierStep1=NULL)

##get the profiles for each gene
L.specific.result.profile=getProfile(generalResult$specificResult$L.specific.result
                                                      $M.specific)
##or
specificResult=generalResult$specificResult
L.specific.result.profile=getProfile(specificResult$L.specific.result$M.specific)
}
