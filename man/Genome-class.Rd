% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenomeClass.R
\docType{class}
\name{Genome-class}
\alias{Genome-class}
\title{\code{Genome} class}
\description{
The \code{Genome} class allows to define all the characteristics of the
genome which will be used as a basis for the construction of transition
matrices from one generation to another in simulations of the model.
}
\details{
A genome includes the list of all possible haplotypes and genotypes
resulting from the combination of the alleles defined in input.
As the \code{Ease} package was originally built for population genetics
simulations including both diploid and haploid loci, it is necessary
that both types of loci are defined. Despite this, the user can define
only diploid or only haploid loci if they wish. If no diploid locus is
defined, one is automatically generated with only one allele, thus not
influencing the simulation. The same applies if no haploid locus is defined.

Each locus is described by a vector of factors which are the names of
the possible alleles at that locus. All diploid (resp. haploid) loci
thus defined are grouped in a list, called \code{listDipLoci} (resp.
\code{listHapLoci}). Therefore, a \code{Genome} class object has two lists
of loci defined in this way, one for diploid loci, one for haploid loci.
The alleles and loci (diploid and haploid) must all have different
names so that no ambiguity can persist.

If several are defined, the order of the diploid loci in the list is not
trivial. The rates of two-to-one combinations between them must indeed be
defined by the vector \code{recRate}. For example, if three diploid loci
are defined, \code{recRate} must be of length 2, the first of its values
defining the recombination rate between the first and second loci, the
second of its values the recombination rate between the second and third
loci. For example, if we want to define two groups of two loci that are
linked to each other but are on two different chromosomes, we can define
a \code{recRate = c(0.1, 0.5, 0.1)}. The first two loci are thus relatively
linked (recombination rate of 0.1), as are the last two loci. On the other
hand, the recombination rate of 0.5 between the second and third loci
ensures that the two groups are independent.
}
\section{Slots}{

\describe{
\item{\code{listHapLoci}}{a list of haploid loci}

\item{\code{listDipLoci}}{a list of diploid loci}

\item{\code{recRate}}{a two-by-two recombination rate vector}

\item{\code{nbHL}}{the number of haploid loci}

\item{\code{nbDL}}{the number of diploid loci}

\item{\code{listLoci}}{the list of all loci}

\item{\code{haplotypesHL}}{haplotypes of haploid loci only}

\item{\code{haplotypesDL}}{haplotypes of diploid loci only}

\item{\code{haplotypes}}{haplotypes of all loci}

\item{\code{alleles}}{the vector of all the alleles}

\item{\code{nbAlleles}}{the number of alleles}

\item{\code{nbHaplo}}{the number of haplotypes}

\item{\code{IDhaplotypes}}{IDs of haplotypes}

\item{\code{genotypes}}{the list of genotypes}

\item{\code{nbGeno}}{the number of genotypes}

\item{\code{IDgenotypes}}{IDs of genotypes}

\item{\code{IDgenome}}{ID of the genome}
}}

\author{
Ehouarn Le Faou
}
