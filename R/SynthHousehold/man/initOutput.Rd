\name{initOutput}
\alias{initOutput}

\title{
Initilize output data structure for saved parameters.
}
\description{
Initilize output data structure for saved parameters.
}
\usage{
initOutput(data, hyper, mc)
}

\arguments{
  \item{data}{
  Input data, including household info, and level of each variable.
}
  \item{hyper}{
  Hyper parameters for priors.
}
  \item{mc}{
  MCMC parameters.
}
}
\details{
 Currently, the program is harded code to the specific data format used in the the paper. Needs to be generalized once a standard foramt is defined and the general format to define structural zeros is finalized.
}
\value{
A list of output parameters to be saved.
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
}

