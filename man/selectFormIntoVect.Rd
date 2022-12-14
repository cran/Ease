% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelFunctions.R
\name{selectFormIntoVect}
\alias{selectFormIntoVect}
\title{Conversion of selection formulas}
\usage{
selectFormIntoVect(selectFormula, genomeObj, haplo = FALSE)
}
\arguments{
\item{selectFormula}{a list of selection formulas}

\item{genomeObj}{a \code{Genome} object}

\item{haplo}{logical indicating whether the selection should apply to
haplotypes (in the case of gametic selection for example)}
}
\value{
a vector of fitness values
}
\description{
Conversion of a list of selection formulas into a genotypic (or haplotypic)
fitness vector associated with a \code{Genome} object.
}
\author{
Ehouarn Le Faou
}
