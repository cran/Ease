% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UserFunctions.R
\name{setMetapopulation}
\alias{setMetapopulation}
\title{Setting a metapopulation}
\usage{
setMetapopulation(populations, migMat = matrix(1))
}
\arguments{
\item{populations}{a list of \code{Population} objects}

\item{migMat}{a migration matrix}
}
\value{
a \code{Metapopulation} object
}
\description{
A metapopulation is a set of population(s) (from 1) that are simulated
with potential migration between them. Only genotypes can migrate, i.e.
adult individuals.
}
\details{
The construction of a \code{Metapopulation} object requires only two
arguments (one optional). The first is a population(s) list, defined
from the population class. The second is a migration matrix, which
connects the populations together. This matrix is a probability matrix
(square with the sum of the rows equal to 1, whose size is equal to the
number of populations) where each value corresponds to the proportion
of individuals (genotypes) that disperse from their source population
(row) to their target population (column).
}
\examples{
# Definition of a population in its simplest form:
DL <- list(dl = c("A", "a"))
HL <- list(hl = c("B", "b"))
mutations <- list(
  mutation(from = "A", to = "a", rate = 1e-3),
  mutation(from = "B", to = "b", rate = 1e-3)
)
genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
pop <- setPopulation(
  name = "A",
  size = 1000,
  dioecy = TRUE,
  genomeObj = genomeObj,
  selectionObj = setSelectNeutral(genomeObj),
  mutMatrixObj = setMutationMatrix(genomeObj, mutations = mutations)
)
metapop <- setMetapopulation(populations = list(pop))
metapop <- simulate(metapop, nsim = 10, seed = 123)
# Other examples available in the documentation of the package

}
\author{
Ehouarn Le Faou
}
