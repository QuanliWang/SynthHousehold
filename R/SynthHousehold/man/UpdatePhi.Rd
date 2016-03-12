\name{UpdatePhi}
\alias{UpdatePhi}
\title{
Update phi parameter
}
\description{
Update phi parameter.
}
\usage{
UpdatePhi(IndividualData_all, z_Individual_all, K, L, p, d, maxd)
}
\arguments{
  \item{IndividualData_all}{
  data matrix for all individuals, including invididuals from impossible households.
}
  \item{z_Individual_all}{
 component indicator matrix for all indiciduals, including individuals from impossible households.
}
  \item{K}{
  Maximum number of components at the household level.
}
  \item{L}{
  Maximum number of components at the individual level.
}
  \item{p}{
  Number of individual level variables in the model.
}
  \item{d}{
  A vector for number of levels of each individual level variable.
}
  \item{maxd}{
  Maximum number of levels from all individual level variables.
}
}

\value{
Updated phi parameter.
}

\author{
Quanli Wang
}

\keyword{ sampler }
