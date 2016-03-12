\name{UpdateLambda}
\alias{UpdateLambda}
\title{
Update lambda parameters
}
\description{
Update lambda parameters.
}
\usage{
UpdateLambda(dHH, K, z_HH_all, HHdata_all)
}
\arguments{
  \item{dHH}{
  A vector for number of levels for each household level variable.
}
  \item{K}{
  Maximum number of components at the household level.
}
  \item{z_HH_all}{
  Component indicators at the household level for all households, including all impossible households.
}
  \item{HHdata_all}{
  Data matrix for household level variables for all households, including all impossible households.
}
}

\value{
Updated phi parameter.
}

\author{
Quanli Wang
}

\keyword{ sampler }

