# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                             MutationMatrix CLASS                             #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# The MutationMatrix class is a class that converts allelic mutation matrices
# into a haplotypic mutation matrix and checks whether the input allelic
# mutation matrices are compatible with an object of class 'genome'.
# The creation of a separate class for this mutation matrix is justified by the
# need to have these allelic mutation matrices as input, which must therefore
# be verified.

## Validity check ----

#' The validity check associated with the \code{MutationMatrix} class
#'
#' @param object an object of class \code{MutationMatrix}
#'
#' @return A logical corresponding to whether \code{x} is a correct
#' \code{MutationMatrix} object.
#'
#' @author Ehouarn Le Faou
#'
check.mutationMatrix <- function(object) {
  for (i in 1:length(object@mutHapLoci)) {
    is.correct.transition.matrix(object@mutHapLoci[[i]], "mutation", "hapLoci")
  }
  for (i in 1:length(object@mutDipLoci)) {
    is.correct.transition.matrix(object@mutDipLoci[[i]], "mutation", "dipLoci")
  }
  if (length(object@mutHapLoci) != object@nbHL) {
    stop(paste(
      "The number of mutation matrices for haploid loci is not",
      "equal to their number."
    ))
  }
  if (length(object@mutDipLoci) != object@nbDL) {
    stop(paste(
      "The number of mutation matrices for diploid loci is not",
      "equal to their number."
    ))
  }
  sizeMutMatHL <- unlist(lapply(object@mutHapLoci, nrow))
  sizeMutMatDL <- unlist(lapply(object@mutDipLoci, nrow))
  if (any(sizeMutMatHL != object@nbAlHL) | any(sizeMutMatDL != object@nbAlDL)) {
    stop(paste(
      "The size of the mutation matrices for each loci should be",
      "equal to the number of alleles"
    ))
  }
  return(TRUE)
}

## Class definition ----

#' Mutation matrix
#'
#' A mutation matrix is used to simulate mutations that affect loci. An object
#' of the class \code{MutationMatrix} does not only contain a (haplotypic)
#' mutation matrix. It also contains the attributes necessary for the
#' construction and easy-to-read display of this matrix.
#'
#' The mutation matrix itself is a square matrix of size equal to the number of
#' haplotypes. It is a probability matrix in that the sum of the values in
#' each row is equal to 1. For a given haplotype, the row associated with it
#' describes the probabilistic proportions that lead by mutation of this
#' haplotype to the production of the other haplotypes (and of itself if there
#' are no mutations).
#'
#' @slot mutHapLoci a list of haploid locus by locus allelic mulation matrices.
#' @slot mutDipLoci a list of diploid locus by locus allelic mulation matrices.
#' @slot mutLoci a list concatenating \code{mutHapLoci} and \code{mutDipLoci}
#' @slot nbAlDL a vector of the number(s) of alleles at each haploid locus
#' @slot nbAlHL a vector of the number(s) of alleles at each diploid locus
#' @slot mutationMatrix the haplotypic mutation matrix
#' @slot nbHaplo the number of haplotypes
#' @slot nbDL the number of diploid loci
#' @slot nbHL the number of haploid loci
#' @slot haplotypes the enumeration of haplotypes
#' @slot IDgenome ID of the associated genome
#'
#' @author Ehouarn Le Faou
#'
#' @export
setClass("MutationMatrix",
  representation(
    mutHapLoci = "list",
    mutDipLoci = "list",
    mutLoci = "list",
    nbAlDL = "numeric",
    nbAlHL = "numeric",
    mutationMatrix = "matrix",
    nbHaplo = "numeric",
    nbDL = "numeric",
    nbHL = "numeric",
    haplotypes = "list",
    IDgenome = "character"
  ),
  validity = check.mutationMatrix
)

## Initialize method ----

#' Initialize method for the \code{MutationMatrix} class
#'
#' @param .Object a \code{MutationMatrix} object
#' @param genomeObj a \code{Genome} object
#' @param mutHapLoci a list of haploid locus by locus allelic mulation matrices.
#' @param mutDipLoci a list of diploid locus by locus allelic mulation matrices.
#'
#' @return A \code{MutationMatrix} object
#'
#' @author Ehouarn Le Faou
#'
setMethod("initialize", "MutationMatrix", function(.Object, genomeObj,
                                                   mutHapLoci, mutDipLoci) {
  # Definition of attributes
  .Object@mutHapLoci <- mutHapLoci
  .Object@mutDipLoci <- mutDipLoci
  names(.Object@mutHapLoci) <- names(genomeObj@listHapLoci)
  names(.Object@mutDipLoci) <- names(genomeObj@listDipLoci)
  .Object@nbAlDL <- unlist(lapply(genomeObj@listDipLoci, length))
  .Object@nbAlHL <- unlist(lapply(genomeObj@listHapLoci, length))
  .Object@haplotypes <- genomeObj@haplotypes
  .Object@nbHaplo <- genomeObj@nbHaplo
  .Object@nbDL <- genomeObj@nbDL
  .Object@nbHL <- genomeObj@nbHL
  .Object@IDgenome <- genomeObj@IDgenome

  # Validity of the object
  validObject(.Object)

  # Column and row naming for allelic mutation matrices
  for (i in 1:.Object@nbHL) {
    colnames(.Object@mutHapLoci[[i]]) <- genomeObj@listHapLoci[[i]]
    rownames(.Object@mutHapLoci[[i]]) <- genomeObj@listHapLoci[[i]]
  }
  for (i in 1:.Object@nbDL) {
    colnames(.Object@mutDipLoci[[i]]) <- genomeObj@listDipLoci[[i]]
    rownames(.Object@mutDipLoci[[i]]) <- genomeObj@listDipLoci[[i]]
  }
  # Compilation of both haploid and diploid loci mutation matrices
  .Object@mutLoci <- list(HL = .Object@mutHapLoci, DL = .Object@mutDipLoci)

  # Haplotypic mutation matrix
  .Object@mutationMatrix <- matrix(0, .Object@nbHaplo, .Object@nbHaplo)
  for (i in 1:.Object@nbHaplo) {
    for (j in 1:.Object@nbHaplo) {
      mut <- 1
      for (k in 1:.Object@nbDL) {
        mut <- mut * .Object@mutLoci$DL[[k]][
          which(.Object@haplotypes$DL[i, k] == genomeObj@listDipLoci[[k]]),
          which(.Object@haplotypes$DL[j, k] == genomeObj@listDipLoci[[k]])
        ]
      }
      for (k in 1:.Object@nbHL) {
        mut <- mut * .Object@mutLoci$HL[[k]][
          which(.Object@haplotypes$HL[i, k] == genomeObj@listHapLoci[[k]]),
          which(.Object@haplotypes$HL[j, k] == genomeObj@listHapLoci[[k]])
        ]
      }
      .Object@mutationMatrix[i, j] <- mut
    }
  }
  colnames(.Object@mutationMatrix) <- genomeObj@IDhaplotypes
  rownames(.Object@mutationMatrix) <- genomeObj@IDhaplotypes

  return(.Object)
})

## Show method ----

#' Show method for the \code{MutationMatrix} class
#'
#' @param object a \code{MutationMatrix} object
#'
#' @return No return value, only a display.
#'
#' @author Ehouarn Le Faou
#'
#' @export
setMethod("show", "MutationMatrix", function(object) {
  catn("-=-=-=- MUTATION MATRIX OBJECT -=-=-=-")

  catn(" #  Haplotypic mutation matrix:")
  print(object@mutationMatrix)

  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  catn("(use print to access the allelic mutation matrices)")
})


## Print method ----

#' Print method for the \code{MutationMatrix} class
#'
#' @param x a \code{MutationMatrix} object
#' @param ... there are no more parameters.
#'
#' @author Ehouarn Le Faou
#'
#' @return No return value, only a display.
#'
#' @export
setMethod("print", "MutationMatrix", function(x, ...) {
  catn("-=-=-=- MUTATION MATRIX OBJECT -=-=-=-")
  catn("              in details")
  catn()
  if (x@nbHL == 1) {
    catn(" #  ", x@nbHL, " haploid locus allelic matrix:")
  } else {
    catn(" #  ", x@nbHL, " haploid loci allelic matrices:")
  }
  for (i in seq_len(x@nbHL)) {
    catn()
    catn(names(x@mutHapLoci)[i], ":")
    print(x@mutHapLoci[[i]])
  }
  catn()
  if (x@nbDL == 1) {
    catn(" #  ", x@nbDL, " diploid locus allelic matrix:")
  } else {
    catn(" #  ", x@nbDL, " diploid loci allelic matrices:")
  }
  for (i in seq_len(x@nbDL)) {
    catn()
    catn(names(x@mutDipLoci)[i], ":")
    print(x@mutDipLoci[[i]])
  }
  catn()
  catn(" #  Haplotypic mutation matrix:")
  print(x@mutationMatrix)

  catn()
  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
})
