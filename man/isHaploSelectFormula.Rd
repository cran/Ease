% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelFunctions.R
\name{isHaploSelectFormula}
\alias{isHaploSelectFormula}
\title{Are there any allelic combinations including homozygosity}
\usage{
isHaploSelectFormula(selectFormula)
}
\arguments{
\item{selectFormula}{a list of selection formula}
}
\value{
logical indicating if there are homozygotes
}
\description{
Test if there are homozygotes in the specified allelic combinations of
a list of selection formulas
}
