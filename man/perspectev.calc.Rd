\name{perspectev.calc}
\alias{perspectev.calc}
\title{
Calculate upper and lower level trait values and survivorship
}
\description{
Calculates upper and lower level trait values and survivorship given a data set and trait functions}
\usage{
perspectev.calc(data,traitfun,vlist=list())
}
\arguments{
 \item{data}{
Dataframe in perspectev format (see below).}
\item{traitfun}{
Function for calculating trait value at each level.
}
\item{vlist}{
Optional variable list for trait function.
}
}

\value{
    \item{upper}{Table of survivorship and trait value at upper level (scaled to mean 0 and variance 1)}
    \item{lower}{Table of survivorship and trait value at lower level (scaled to mean 0 and variance 1)}
    \item{stats}{[1]: correlation between upper level trait and survivorship, [2-4]: quantiles of upper level trait values (unscaled)}
    
}

\author{
    Kenneth B. Hoehn <perspectev@gmail.com>
}

