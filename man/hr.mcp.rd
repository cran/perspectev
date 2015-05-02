\name{hr.mcp}
\alias{hr.mcp}
\title{
Home range estimation
}
\description{
Compute and plot a minimum convex polygon estimate of "home range." Source code borrowed from wild1 package.
}
\usage{
hr.mcp(x,y=NULL,n.min=50,plot=TRUE,add=FALSE,ID=NULL,...)
}

\arguments{
  \item{x}{
        A vector, matrix, or data frame of coordinates.
}
  \item{y}{
        A numeric vector (if \code{x} is a vector) or \code{NULL} 
}
  \item{n.min}{
        If \code{x} includes fewer than \code{n.min} points, the function will return \code{NA} and a warning.
}
  \item{plot}{
        Parameter nonfunctional for perspectev implementation of this function.
}
  \item{add}{
        Logical. Add to existing plot on current output device?  Default is \code{FALSE}.
}
  \item{ID}{
        Required argument from Polygons
}
  \item{\dots}{
        Optional arguments for \code{\link{plot}} or \code{\link{points}} (see also \code{\link{par}}).
}
}
\details{
See plot.Polygons for another way of plotting polygon objects.
}
\value{
Returns an object of class Polygons (sp).
}

\author{
Glen A. Sargeant\cr
U.S. Geological Survey\cr
Northern Prairie Wildlife Research Center\cr
\email{glen_sargeant@usgs.gov}

Edited by Kenneth B. Hoehn\cr
\email{perspectev@gmail.com}

}
\keyword{ aplot }
\keyword{ hplot }

