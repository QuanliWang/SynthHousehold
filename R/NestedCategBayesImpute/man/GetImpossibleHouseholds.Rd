\name{GetImpossibleHouseholds}
\alias{GetImpossibleHouseholds}
\title{
Generate the desired number of impossible households for given household size and model parameters.
}
\description{
Generate the desired number of impossible households for given household size and model parameters. Also generate possible household data (synthetic data) when required.
}
\usage{
GetImpossibleHouseholds(d, ACS_count, lambda, w, phi, pi, howmany, n, synindex, format2)
}
\arguments{
  \item{d}{
  number of levels for each variable.
}
  \item{ACS_count}{
  number of households for different household size. Need to change this name later.
}
  \item{lambda}{lambda parameters for the household level variables.
}
  \item{w}{
  w parameter.
}
  \item{phi}{
  phi parameter.
}
  \item{pi}{
  pi parameter.
}
  \item{howmany}{
  number households to be generated at each batch. It is used as a tuning parameter for performance.
}
  \item{n}{
  number of household in the original input data.
}
  \item{synindex}{
  indicates if a possible household data(synthetic data) is needed.
}
  \item{format2}{
  indicator for whether or not the household head is moved to the household level (TRUE) or included as an individual within the household (FALSE) . 
}
}

\value{
A list of variales  relating to impossible household/individuals and synthetic data. Really confusing here and need more clarifiction.
  \item{Individuals_extra }{A data matrix for individual members from the impossible households.}
  \item{HHdata_extra }{A data matrix for impossible household level variables.}
  \item{z_HHdata_individual_extra }{A data matrix for impossible household level component indicators. The household level indicators are expanded to individual members within this data matrix.}
  \item{z_HH_extra }{A vector for the household level compoment indicators.}
  \item{hh_size_new }{A vector for the number of impossible households for different household size.}
  \item{synIndividuals_all }{A data matrix for individual members from the possible(synthetic) households when synindex is non-zero. NULL otherwise.}
}

\author{
Quanli Wang
}

\keyword{sampler}
\keyword{impossible household}
\keyword{synthetic data}
