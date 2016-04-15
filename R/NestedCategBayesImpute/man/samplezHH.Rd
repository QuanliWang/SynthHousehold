\name{samplezHH}
\alias{samplezHH}
\title{
Update household level component indexes.
}
\description{
Update household level component indexes.
}
\usage{
samplezHH(phi, data, w, pi, S, HHdata, lambda)
}
\arguments{
  \item{phi}{
  phi parameters.
}
  \item{data}{
  Input data.
}
  \item{w}{
  w parameter.
}
  \item{pi}{
  pi parameter.
}
  \item{S}{
  household size variable.
}
  \item{HHdata}{
  Household level variables. Currently hard-coded to two variables in the paper.
}
  \item{lambda}{
  a list of lambdas for household level variables
}
}
\details{
Need to add more details by Monica or Michael.
}
\value{
A list with two variables.
  \item{zHH }{component indexes for input households}
  \item{z_HH_Individuals }{component indexes for input individuals}
}

\author{
Quanli Wang
}

\keyword{sampler}
\keyword{household level}
