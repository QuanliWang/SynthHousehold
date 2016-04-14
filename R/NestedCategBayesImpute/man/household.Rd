\name{household}
\alias{household}
\docType{data}
\title{
example dataset to be used for library development.
}
\description{
This is the example data set that is used in Monica's paper (not exactly the same, but a new random sample of same size).
}
\usage{data("household")}

\format{
  A data frame with 10000 observations on the following 9 variables.\describe{
    \item{\code{Hhindex}}{Household index}
    \item{\code{pernum}}{Ask Monica to fill it in}
    \item{\code{sex}}{Gender. 1 = male, 2 = female}
    \item{\code{race}}{Recoded general race code. 1 = white alone, 2 = black alone,3 = American Indian/Alaska Native alone,4 = Asian or Pacific Islander alone,5 = other, 6 = two or more races.}
    \item{\code{sthn}}{Gender. 1 = male, 2 = female}
    \item{\code{age}}{Age. 1=0(lessthenoneyearold),2=1,...,94 = 93}
    \item{\code{relate}}{1 = head/householder, 2 = spouse, 3 = child, 4 = child-in-law, 5 = parent, 6 = parent-in- law, 7 = sibling, 8 = sibling-in-law,9 = grandchild, 10 = other relatives,11 = partner, friend, visitor,12 = other non-relatives}
    \item{\code{ownership}}{Ownership of dwelling. 1 = owned or being bought (loan), 2 = rented.}
    \item{\code{householdsize}}{Household size. 2=2people,3=3people,4=4people.}
  }
}
\details{
Monica please provide the details here about the dataset.
}
\source{
From Monica. Need real sources.
}
\references{
Jingchen Hu, Jerome P. Reiter and Quanli Wang. Dirichlet process mixture models for nested categorical data. Vol 00(000). DOI: 0000.
}
\examples{
data(household)
}
\keyword{datasets}
\keyword{household data}
\keyword{example}
