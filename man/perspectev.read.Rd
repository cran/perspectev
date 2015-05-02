\name{perspectev.read}
\alias{perspectev.read}
\title{
Convert from PBDB to perspectev input formats
}
\description{
Convert occurrence data (each row a separate occurrence) into standardized perspectev format.}
\usage{
perspectev.read(data,extinctionAge,occurrenceAge,upper,lower,t1,t2,
  traitfun=mcpRange,vlist=NULL,trim=TRUE,projection=FALSE)
}
\arguments{
  \item{data}{
Occurrence data.
}
  \item{occurrenceAge}{
Column desigation for the age of the occurrence.
}
  \item{extinctionAge}{
Age of the extinction. E.g. 65.5 for K-Pg mass extinction if occurrenceAge is in mya.
}
  \item{upper}{
Column designation for upper level names.
}
  \item{lower}{
Column designation for lower level names.
}
  \item{t1}{
Column designation for trait-data 1.
}
  \item{t2}{
Column designation for trait-data 2.
}
  \item{trim}{
Remove lower levels that give NA's from trait function (recommended)?
}
  \item{projection}{
Convert coordinates to a Mollweide projection?
} 
\item{traitfun}{
Function for calculating trait values at each level. Not used if trim==FALSE.
}
\item{vlist}{
Optional variable list for trait function.
}
}
\details{
See tutorial for details. Note that this will assign any lower level entity found both before/during and after the extinctionAge parameter will be listed as a survivor, while any lower level entity only found in the beginning will be listed as a victim. It is important to ensure your data set does not span too far back, or many of the species that went extinct before the extinction will be listed as victims. This is not good.

If survivorship is to be given probabilistically rather than as a binary, run this function with extinctionAge=0, and input the survivorship probabilities into the Survivorship column of the resulting data frame, which lists each unique species on a separate row.
}
\examples{
  data(testData)

  data = perspectev.read(testData,extinctionAge=5,occurrenceAge="Age",
  upper="Genus",lower="Species",t1="Lat",t2="Long",traitfun=mcpRange,projection=FALSE)
}

\author{
    Kenneth B. Hoehn <perspectev@gmail.com>
}
