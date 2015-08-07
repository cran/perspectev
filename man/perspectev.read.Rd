\name{perspectev.read}
\alias{perspectev.read}
\title{
Convert occurrence data to perspectev input format
}
\description{
Convert occurrence data (each row a separate occurrence) into standardized perspectev format.}
\usage{
perspectev.read(data,upper,lower,traits,survProb,traitfun,extinctionAge=NULL,
  occurrenceAge=NULL,vlist=NULL,trim=TRUE,
  projection=FALSE,quiet=FALSE)
}
\arguments{
  \item{data}{
Occurrence data. This must list each occurrence on a separate row, with upper/lower level classification specified in separate columns, as well as one or more columns specifying the trait-data for that occurrence (e.g. latitude and longitude). For example see data(occurrence_data).
}
  \item{upper}{
Column name designation for upper level names.
}
  \item{lower}{
Column name designation for lower level names.
}
  \item{traits}{
  A column name (string) or vector of column names which correspond to the trait data under investigation. The order in which trait names are specified is important and corresponds to the t1/t2 designation of trait functions.
}
  \item{survProb}{
  Designates lower level survivorship probability, may be either:
  1) A dataframe with three columns: Upper name, lower name, survivorship probability.
      Preferred option - compatible with continuous survivorship probability.

  2) A dataframe with two columns: Upper name, lower name. 
      Specifies lower units named as survivorship probability = 1, and all others as 0.

  3) A vector of names of surviving lower level units. 
      Upper and lower names must be separated by a space.

  If unspecified will assign survivorship based on occurrenceAge and extinctionAge.
}
  \item{occurrenceAge}{
Column name desigation for the age of the occurrence. Only necessary if survProb unspecified, and will be ignored if survProb is specified.
}
  \item{extinctionAge}{
Age of the extinction in backward time. E.g. 65.5 for K-Pg mass extinction if occurrenceAge is in mya. Occurrences with ages greater than or equal to extinctionAge will be used to calculate trait values, and lower level units with occurrences before and after extinctionAge are marked as survivors. Lower level units with occurrences only after extinctionAge are excluded. Only necessary if survProb unspecified, and will be ignored if survProb is specified.
}
\item{traitfun}{
Function for calculating trait values at each level. Not necessary if trim=FALSE.
}
\item{vlist}{
Optional variable list for trait function.
}
  \item{trim}{
Remove lower levels that give NA's from trait function? (recommended)
}
  \item{projection}{
Convert coordinates to a Mollweide projection?
} 
  \item{quiet}{
Suppress information messages?
} 
}

\value{

\item{Survivorship}{Dataframe containing upper/lower level classification and survivorship probability for each lower level unit.}
\item{Traits}{List containing trait data vectors for each lower level unit.}
\item{Key}{Shows indexes of trait values in Traits list, corresponds to order used in "traits" parameter.}

}

\details{

Function specifies traits under analysis and survivorship probabilities of lower level units. See tutorial (http://evolve.zoo.ox.ac.uk/Evolve/Perspectev.html) for examples.

Traits are specified by supplying the column names of trait data values (e.g. Latitude and Longitude) as a vector. The order in which the traits are specified is important, as this will be the basis of the index in trait functions (first in vector is first in trait list, second in vector is second in trait list, etc..). You can check this after the fact using the Key value. If trim=TRUE, supply a trait function using the traitfun parameter, and this function will remove all lower units which give an NA value in the trait function (e.g. having < 3 distinct coordinates for MCP range).

Survivorship probability may be specified in four ways:

1) The best way is to supply a dataframe to survProb. This data frame should be of column format: <Upper Name> <Lower Name> <Survivorship Probability>, and should specify a survivorship probability for all lower level units in the occurrence data.

2) Similar to 1), but a dataframe of column format: <Upper Name> <Lower Name>. This will assign all lower units specified as having a survivorship probability of 1, and all others as 0.

3) Similar to 2), but a vector of names formed by pasting Upper and Lower names together, separated by a space (e.g. c("Genus1 Species1", "Genus2 Species2")) 

4) Specifying occurrenceAge and extinctionAge. occurrenceAge is the column name for the age of the occurrences, and extinctionAge is the age of the extinction (both in backward time). All units with occurrences that only occur after the extinction are discarded. All units that have some occurrences with ages < extinctionAge, and some occurrences >= extinctionAge (found before and after extinction) are assigned survivorship = 1. All those which only have occurrences >= extinctionAge are assigned survivorship = 0. Note that this will assign any lower level entity found both before/during and after the extinctionAge parameter as a survivor, while any lower level entity only found in the beginning will be listed as a victim. Specifying these parameters is only necessary if survProb is unspecified, and they will be ignored if it is.

In all cases it is important to ensure your data set does not span too far back, because all occurrences provided, unless otherwise stated (see 4, above), will be used to calculate trait values. This is may or may not be good.
}

\examples{
  data(testData)
  survivorship = data.frame(cbind(testData$Genus,testData$Species,runif(nrow(testData),0,1)))
  data = perspectev.read(testData,survProb=survivorship,
  upper="Genus",lower="Species",traits=c("Lat","Long"),traitfun=mcpRange,
  projection=FALSE)
}

\author{
    Kenneth B. Hoehn <perspectev@gmail.com>
}
