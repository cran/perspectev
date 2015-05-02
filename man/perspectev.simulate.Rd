\name{perspectev.simulate}
\alias{perspectev.simulate}
\title{
Simulate species/genus level selection at different intensities and observe resulting p values.
}
\description{
Perform simulations of natural selection on a data set, and observe what p values result. Perform linear transformation of trait vs. survivorship (using provided or fitted slope and intercept) and then transform into survviorship probabilites using inv.logit function.
}
\usage{
perspectev.simulate(data,simulations,cores,traitfun=mcpRange,vlist=NULL,
	binary=NA,intercept=NA,slope=NA,level=NA,noise=0,fit=FALSE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{data}{
Dataframe in perspectev format. See ?perspectev.read for more information.}
  \item{simulations}{
Number of simulations to perform. At least 1000 is recommended but can be slow.}
  \item{cores}{
Number of cores over which to parallelize the script.}
\item{traitfun}{
Function for calculating trait value at each level.
}
\item{vlist}{
Optional variable list for trait function.
}
  \item{intercept}{
Intercept of the regression between range and survivorship. Can be manually provided or, if not provided, is fit to the data set. Will attempt to set automatically if unspecified.
}
  \item{slope}{
Slope of the regression between range and survivorship. Can be manually provided or, if not provided, is fit to the data set. Will attempt to set automatically if unspecified.
}
  \item{level}{
Simulate either upper or lower level selection. Accepts either "lower" or "upper" respectively. Will set to 'upper' if unspecified.
}
  \item{noise}{
Extra stochasticity to the simulation. 
}
  \item{fit}{
Fit simulation to observed parameters? Will do automatically if intercept or slope are not specified.
}
\item{binary}{
Is survivorship binary or continuous? Will attempt to set automatically if unspecified.
}
}

\value{

\item{correlation_simulated}{Correlation between trait and surivovorship from permuted groups at each simulation.}
\item{correlation_observed}{Observed correlation between trait and surivovorship from simulated groups.}
\item{pvalue}{p value observed under simulation of selection.}
\item{level}{Level of simulated selection.}
\item{intercept}{Intercept of simulation.}
\item{slope}{Slope of the simulation.}
\item{noise}{Noise level selected for simulation.}
\item{fitted_model}{Was the model fitted to the data set?}
}


\author{
    Kenneth B. Hoehn <perspectev@gmail.com>
}
\note{
See original publication and tutorial for more information on simulation.
}
