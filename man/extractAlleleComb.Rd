% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelFunctions.R
\name{extractAlleleComb}
\alias{extractAlleleComb}
\title{Extract the allele combination}
\usage{
extractAlleleComb(xVect)
}
\arguments{
\item{xVect}{allelic combination extracted from a selection formula.}
}
\value{
the list of alleles that must be present in the genotype to match
the input allelic combination
}
\description{
Conversion of an allelic combination defined in a selection formula into the
vector listing the alleles present (alleles that must be in the homozygous
state appear 2 times, 1 time for heterozygous).
}
\author{
Ehouarn Le Faou
}
