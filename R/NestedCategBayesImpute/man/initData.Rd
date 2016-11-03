\name{initData}
\alias{initData}
\title{
Initialize the input data structure.
}
\description{
Initialize the input data structure.
}
\usage{
initData(household, individual_variable_index, household_variable_index)
}
\arguments{
  \item{household}{
  A data matrix for the input household data.
}
  \item{individual_variable_index}{
  Vector of column indexes for the individual-level variables.
}
  \item{household_variable_index}{
  Vector of column indexes for the group-level variables.
}
}
%\details{
%Curently the variables and the level of each variable are harded coded. This needs to be changed later on once it is generalized.
%}
\value{
A list object including all the necessary data variables needed by the  model.
}

\author{
Quanli Wang
}

