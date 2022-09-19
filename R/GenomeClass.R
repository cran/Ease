# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                 Genome CLASS                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Validity check ----

#' The validity check for the \code{Genome} class
#'
#' @param object a \code{Genome} object
#'
#' @return A logical corresponding to whether the object is a correct
#' \code{Genome} object.
#'
#' @author Ehouarn Le Faou
#'
check.genome <- function(object) {

  # Checking the allele enumeration vector
  if (any(unlist(lapply(object@listHapLoci, class)) != "factor")) {
    stop(paste(
      "The enumeration list of the alleles of the loci must be",
      "a list of loci where each of them is associated with a",
      "vector containing as factors the names of each of its",
      "alleles."
    ))
  }
  if (any(unlist(lapply(object@listDipLoci, class)) != "factor")) {
    stop(paste(
      "The enumeration list of the alleles of the loci must be",
      "a list of loci where each of them is associated with a",
      "vector containing as factors the names of each of its",
      "alleles."
    ))
  }

  # Checking the recombination vector
  if (!identical(object@recRate, numeric())) {
    if (any(object@recRate < 0 | object@recRate > 1)) {
      stop("The recombination rate(s) must be a value(s) between 0 and 1.")
    }
    if (length(object@recRate) != (length(object@listDipLoci) - 1)) {
      stop(paste(
        "The number of recombination rates should be equal to the",
        "number of diploid loci minus one."
      ))
    }
  } else {
    if (length(object@listDipLoci) > 1) {
      stop(paste(
        "The number of recombination rates should be equal to the",
        "number of diploid loci minus one."
      ))
    }
  }

  # Checking the names of loci/alleles
  loci <- names(object@listLoci)
  if ("" %in% loci) {
    stop("Each loci must be named.")
  }
  alleles <- as.character(unlist(object@listLoci))
  if (length(union(loci, loci)) != length(loci) ||
    length(union(alleles, alleles)) != length(alleles)) {
    stop("Each loci and allele must have a unique name.")
  }

  return(TRUE)
}

## Class definition ----

#' \code{Genome} class
#'
#' The \code{Genome} class allows to define all the characteristics of the
#' genome which will be used as a basis for the construction of transition
#' matrices from one generation to another in simulations of the model.
#'
#' A genome includes the list of all possible haplotypes and genotypes
#' resulting from the combination of the alleles defined in input.
#' As the \code{Ease} package was originally built for population genetics
#' simulations including both diploid and haploid loci, it is necessary
#' that both types of loci are defined. Despite this, the user can define
#' only diploid or only haploid loci if they wish. If no diploid locus is
#' defined, one is automatically generated with only one allele, thus not
#' influencing the simulation. The same applies if no haploid locus is defined.
#'
#' Each locus is described by a vector of factors which are the names of
#' the possible alleles at that locus. All diploid (resp. haploid) loci
#' thus defined are grouped in a list, called \code{listDipLoci} (resp.
#' \code{listHapLoci}). Therefore, a \code{Genome} class object has two lists
#' of loci defined in this way, one for diploid loci, one for haploid loci.
#' The alleles and loci (diploid and haploid) must all have different
#' names so that no ambiguity can persist.
#'
#' If several are defined, the order of the diploid loci in the list is not
#' trivial. The rates of two-to-one combinations between them must indeed be
#' defined by the vector \code{recRate}. For example, if three diploid loci
#' are defined, \code{recRate} must be of length 2, the first of its values
#' defining the recombination rate between the first and second loci, the
#' second of its values the recombination rate between the second and third
#' loci. For example, if we want to define two groups of two loci that are
#' linked to each other but are on two different chromosomes, we can define
#' a \code{recRate = c(0.1, 0.5, 0.1)}. The first two loci are thus relatively
#' linked (recombination rate of 0.1), as are the last two loci. On the other
#' hand, the recombination rate of 0.5 between the second and third loci
#' ensures that the two groups are independent.
#'
#' @slot listHapLoci a list of haploid loci
#' @slot listDipLoci a list of diploid loci
#' @slot recRate a two-by-two recombination rate vector
#' @slot nbHL the number of haploid loci
#' @slot nbDL the number of diploid loci
#' @slot listLoci the list of all loci
#' @slot haplotypesHL haplotypes of haploid loci only
#' @slot haplotypesDL haplotypes of diploid loci only
#' @slot haplotypes haplotypes of all loci
#' @slot alleles the vector of all the alleles
#' @slot nbAlleles the number of alleles
#' @slot nbHaplo the number of haplotypes
#' @slot IDhaplotypes IDs of haplotypes
#' @slot genotypes the list of genotypes
#' @slot nbGeno the number of genotypes
#' @slot IDgenotypes IDs of genotypes
#'
#' @examples
#' # Definition of loci:
#' LD <- list(dl = as.factor(c("A", "a")))
#' HL <- list(hl = as.factor(c("B", "b")))
#'
#' # Creation of the object:
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#'
#' # Showing:
#' genomeObj
#'
#' # Printing:
#' print(genomeObj)
#'
#' @author Ehouarn Le Faou
#'
#' @export
setClass("Genome",
  representation(
    recRate = "numeric",
    listHapLoci = "list",
    listDipLoci = "list",
    # Attributes not specified by the user
    nbHL = "numeric",
    nbDL = "numeric",
    listLoci = "list",
    haplotypesHL = "data.frame",
    haplotypesDL = "data.frame",
    haplotypes = "list",
    alleles = "character",
    nbAlleles = "numeric",
    nbHaplo = "numeric",
    IDhaplotypes = "character",
    genotypes = "list",
    nbGeno = "numeric",
    IDgenotypes = "character"
  ),
  validity = check.genome
)

## Initialize method ----

#' Initialize method for the \code{Genome} class
#'
#' @param .Object a \code{Genome} object
#' @param listHapLoci a list of haploid loci
#' @param listDipLoci a list of diploid loci
#' @param recRate a two-by-two recombination rate vector
#'
#' @return A \code{Genome} object
#'
#' @author Ehouarn Le Faou
#'
setMethod("initialize", "Genome", function(.Object, listHapLoci, listDipLoci,
                                           recRate) {

  # Definition of attributes
  .Object@recRate <- recRate
  .Object@listDipLoci <- listDipLoci
  .Object@listHapLoci <- listHapLoci

  # Warnings
  if (length(.Object@listHapLoci) == 0) {
    warning(paste(
      "No haploid locus has been set. By construction it is ",
      "necessary that there is at least 1 haploid locus, so it ",
      "has been defined with only one allele (this will not ",
      "affect the simulations)."
    ))
  }
  if (length(.Object@listDipLoci) == 0) {
    warning(paste(
      "No diploid locus has been set. By construction it is ",
      "necessary that there is at least 1 haploid locus, so it ",
      "has been defined with only one allele (this will not ",
      "affect the simulations)."
    ))
  }

  # Completion if no haploid/diploid/both locus is defined
  if (length(.Object@listDipLoci) == 0) {
    .Object@listDipLoci <- list(NoneDL = as.factor(c("NoneD")))
  }
  if (length(.Object@listHapLoci) == 0) {
    .Object@listHapLoci <- list(NoneHL = as.factor(c("NoneH")))
  }

  # Calculation of other attributes
  .Object@listLoci <- c(.Object@listHapLoci, .Object@listDipLoci)
  .Object@nbHL <- length(.Object@listHapLoci)
  .Object@nbDL <- length(.Object@listDipLoci)

  # Validity of the object
  validObject(.Object)

  # Allele enumeration
  HLalleles <- as.character(unlist(.Object@listHapLoci))
  DLalleles <- as.character(unlist(.Object@listDipLoci))
  .Object@alleles <- c(HLalleles, DLalleles)
  .Object@nbAlleles <- length(.Object@alleles)

  # Haplotyping
  haplo <- haplotyping(.Object)
  .Object@haplotypesHL <- haplo$haplotypesHL
  .Object@haplotypesDL <- haplo$haplotypesDL
  .Object@haplotypes <- haplo$haplotypes
  .Object@nbHaplo <- nrow(.Object@haplotypes$DL)
  .Object@IDhaplotypes <- sapply(
    1:.Object@nbHaplo,
    function(x) {
      IDhaplotypeGeneration(
        .Object@haplotypes$DL[x, ],
        .Object@haplotypes$HL[x, ]
      )
    }
  )

  # Genotyping
  .Object@genotypes <- genotyping(.Object)
  .Object@nbGeno <- nrow(.Object@genotypes$DL1)
  .Object@IDgenotypes <- sapply(
    1:.Object@nbGeno,
    function(x) {
      IDgenotypeGeneration(
        .Object@genotypes$DL1[x, ],
        .Object@genotypes$DL2[x, ],
        .Object@genotypes$HL[x, ]
      )
    }
  )

  return(.Object)
})

## Show method ----

#' Show method for the \code{Genome} class
#'
#' @param object a \code{Genome} object
#'
#' @return No return value, only a display.
#'
#' @author Ehouarn Le Faou
#'
setMethod("show", "Genome", function(object) {
  nbAlleleHL <- sapply(object@listHapLoci, length)
  nbAlleleDL <- sapply(object@listDipLoci, length)
  catn("-=-=-=-=-=-= GENOME OBJECT =-=-=-=-=-=-")
  cat(" # ", object@nbHL)
  if (object@nbHL == 1) {
    catn(" haploid locus, with ", nbAlleleHL, " allele(s)", sep = "")
  } else {
    catn(" haploid loci, with respectively ",
      Reduce(function(x, y) {
        paste0(x, ", ", y)
      }, nbAlleleHL[-object@nbHL]),
      " and ", nbAlleleHL[object@nbHL], " allele(s)",
      sep = ""
    )
  }
  cat(" # ", object@nbDL)
  if (object@nbDL == 1) {
    catn(" diploid locus, with ", nbAlleleDL, " allele(s)", sep = "")
  } else {
    catn(" diploid loci, with respectively ",
      Reduce(function(x, y) {
        paste0(x, ", ", y)
      }, nbAlleleDL[-object@nbDL]),
      " and ", nbAlleleDL[object@nbDL], " allele(s)",
      sep = ""
    )
  }
  catn(" # ", object@nbHaplo, "haplotypes")
  catn(" # ", object@nbGeno, "genotypes")

  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  catn("(use print for a list of haplotypes and genotypes)")
})


## Print method ----

#' Print method for the \code{Genome} class
#'
#' @param x a \code{Genome} object
#' @param ... Ignored.
#'
#' @return No return value, only a display.
#'
#' @author Ehouarn Le Faou
#'
setMethod("print", "Genome", function(x, ...) {
  catn("-=-=-=-=-=-= GENOME OBJECT =-=-=-=-=-=-")
  catn("              in details")
  catn()
  catn(" #  ", x@nbHL, " haploid loci:")
  for (i in 1:x@nbHL) {
    cat("      - '", names(x@listHapLoci)[i], "' with ", sep = "")
    cat(listing(x@listHapLoci[[i]]))
    catn(" alleles")
  }
  catn()
  catn(" #  ", x@nbDL, " diploid loci:")
  for (i in 1:x@nbDL) {
    cat("      - '", names(x@listDipLoci)[i], "' with ", sep = "")
    cat(listing(x@listDipLoci[[i]]))
    catn(" alleles")
  }
  catn()
  catn(" #  ", x@nbHaplo, " haplotypes:")
  print(paste("  ", paste(1:x@nbHaplo, ")", sep = ""),
    x@IDhaplotypes,
    sep = " "
  ),
  quote = FALSE, row.names = FALSE
  )
  catn()
  catn(" #  ", x@nbGeno, " genotypes:")
  print(paste("  ", paste(1:x@nbGeno, ")", sep = ""),
    x@IDgenotypes,
    sep = " "
  ),
  quote = FALSE, row.names = FALSE
  )
  catn()
  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
})
