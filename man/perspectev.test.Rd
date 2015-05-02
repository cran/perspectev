\name{perspectev.test}
\alias{perspectev.test}
\title{
Test for irreducibility of relationship between higher level traits and survivorship
}
\description{
Performs permutation tests by permuting upper level labels between lower levels, recalculating upper trait value, and taking the correlation between upper level trait and survivorship. This process is repeated until a null distribution is generated. This is then compared against observed covariance to give a p value for the null hypothesis that a relationship between trait and survivorship is explainable by random aggregations of lower level traits.}
\usage{
perspectev.test(data,iterations=1000,cores=1,traitfun=mcpRange,vlist=NULL)
}
\arguments{
  \item{data}{
Dataframe in perspectev format (see below).}
  \item{iterations}{
Number of iterations to perform. At least 1000 is recommended, though can be slow.
}
  \item{cores}{
Number of cores over which to parallelize the test.
}
\item{traitfun}{
Function for calculating trait values at each level.
}
\item{vlist}{
Optional variable list for trait function.
}

}
\value{
    \item{correlation_permuted}{Correlations between trait and survivorship obtained from permuted upper levels (Si)}
    \item{correlation_observed}{Observed correlation between upper level trait and survivorship (Ri)}
    \item{pvalue}{Portion of permuted genus correlations (S) larger than observed value (R)}
    \item{permuted_quantiles}{Matrix of interquartile trait values obtained from each upper level permutation}
}
\examples{
	data(testData)

  	data = perspectev.read(testData,extinctionAge=5,occurrenceAge="Age",
  	upper="Genus",lower="Species",t1="Lat",t2="Long",traitfun=mcpRange,projection=FALSE)

  	#4 iterations chosen out of convenience - use more!
	mcpTest  = perspectev.test(data,4,1,traitfun=mcpRange)
	mcpSim  = perspectev.simulate(data,4,1,traitfun=mcpRange)
	perspectev.plot(mcpTest,list(mcpSim),c("S1"),"Test")
}

\author{
    Kenneth B. Hoehn <perspectev@gmail.com>
}

