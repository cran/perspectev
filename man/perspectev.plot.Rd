\name{perspectev.plot}
\alias{perspectev.plot}
\title{
Plotting function for permutation analysis
}
\description{
Make plots to summarize permutation results and power analysis results
}
\usage{
perspectev.plot(observed, simulated, names, title)
}
\arguments{
  \item{observed}{
Object returned from perspectev.test().
}
  \item{simulated}{
List of up to eight objects obtained from perspectev.simulate().
}
  \item{names}{
Names of each object in the list provided to "simulated" in order.
}
  \item{title}{
Title of the plot. Optional.
}
}

\author{
    Kenneth B. Hoehn <perspectev@gmail.com>
}
