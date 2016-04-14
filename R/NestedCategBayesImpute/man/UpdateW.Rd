\name{UpdateW}
\alias{UpdateW}
\title{
Update W and v parameters
}
\description{
Update W and v parameters.
}
\usage{
UpdateW(beta, z_Individual_all, K, L)
}
\arguments{
  \item{beta}{
  Beta parameter.
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
}

\value{
A list containing updated W and v parameters.
}

\author{
Quanli Wang
}

\keyword{sampler}

