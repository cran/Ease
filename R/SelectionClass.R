# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                               Selection CLASS                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Validity check ----

#' The validity check for the \code{Selection} class
#'
#' @param object \code{Selection} object
#'
#' @return A logical corresponding to whether the object is a correct
#' \code{Selection} object.
#'
#' @author Ehouarn Le Faou
#'
check.selection <- function(object) {
  lenVectsGeno <- c(
    object@indFit,
    object@gamProdFit
  )
  lenVectsHaplo <- object@gamFit
  if (any(sapply(lenVectsGeno, length) != object@nbGeno)) {
    stop(paste(
      "All selection vectors that apply to individuals (directly or",
      "via gamete production) must be equal in length to the number",
      "of genotypes in the associated genome."
    ))
  }
  if (any(sapply(lenVectsHaplo, length) != object@nbHaplo)) {
    stop(paste(
      "All selection vectors that apply to gametes must be equal in",
      "length to the number of haplotypes in the associated genome."
    ))
  }
  return(TRUE)
}


## Class definition ----

#' \code{Selection} class
#'
#' Class used to generate objects that manage the selection in the simulations.
#'
#' An object of type \code{Selection} is an object which describes the set of
#' fitnesses which will be taken into account in the simulations. The
#' selection according to these fitnesses can be applied at three levels:
#' at the level of the individual, at the level of the production of
#' gametes and at the level of the gametes themselves.
#' Selection is therefore genotypic in the first two cases (each genotype
#' is associated with a fitness value) and haplotypic in the third (each
#' haplotype is associated with a fitness value).
#'
#' @slot genome a \code{Genome} object
#' @slot IDhaplotypes IDs of haplotypes
#' @slot IDgenotypes IDs of genotypes
#' @slot IDgenome ID of the associated genome
#' @slot nbHaplo the number of haplotypes
#' @slot nbGeno the number of genotypes
#' @slot gamFit the list of gametes' fitness
#' @slot indFit the list of individuals' fitness
#' @slot gamProdFit the list of gamete production fitness
#' @slot sOnInds a logical indicating whether a selection on individuals
#' has been configured by the user
#' @slot sOnGams a logical indicating whether a selection on gametes
#' has been configured by the user
#' @slot sOnGamsProd a logical indicating whether a selection on gamete
#' production has been configured by the user
#'
#' @author Ehouarn Le Faou
#'
#' @export
setClass("Selection",
  representation(
    genome = "Genome",
    IDhaplotypes = "character",
    IDgenotypes = "character",
    IDgenome = "character",
    nbHaplo = "numeric",
    nbGeno = "numeric",
    gamFit = "list",
    indFit = "list",
    gamProdFit = "list",
    sOnInds = "logical",
    sOnGams = "logical",
    sOnGamsProd = "logical"
  ),
  validity = check.selection
)

## Initialize method ----

#' Initialize method for the \code{Selection} class
#'
#' @param .Object a \code{Selection} object
#' @param genomeObj a \code{Genome} object
#'
#' @return A \code{Selection} object
#'
#' @author Ehouarn Le Faou
#'
setMethod("initialize", "Selection", function(.Object, genomeObj) {
  .Object@genome <- genomeObj
  .Object@IDhaplotypes <- genomeObj@IDhaplotypes
  .Object@IDgenotypes <- genomeObj@IDgenotypes
  .Object@nbHaplo <- genomeObj@nbHaplo
  .Object@nbGeno <- genomeObj@nbGeno
  .Object@IDgenome <- genomeObj@IDgenome

  .Object@gamFit <- list(
    female = rep(1, genomeObj@nbHaplo),
    male = rep(1, genomeObj@nbHaplo)
  )
  names(.Object@gamFit$female) <- genomeObj@IDhaplotypes
  names(.Object@gamFit$male) <- genomeObj@IDhaplotypes
  .Object@indFit <- list(
    ind = rep(1, genomeObj@nbGeno),
    female = rep(1, genomeObj@nbGeno),
    male = rep(1, genomeObj@nbGeno)
  )
  names(.Object@indFit$ind) <- genomeObj@IDgenotypes
  names(.Object@indFit$female) <- genomeObj@IDgenotypes
  names(.Object@indFit$male) <- genomeObj@IDgenotypes
  .Object@gamProdFit <- list(
    female = rep(1, genomeObj@nbGeno),
    male = rep(1, genomeObj@nbGeno)
  )
  names(.Object@gamProdFit$female) <- genomeObj@IDgenotypes
  names(.Object@gamProdFit$male) <- genomeObj@IDgenotypes

  .Object@sOnInds <- FALSE
  .Object@sOnGams <- FALSE
  .Object@sOnGamsProd <- FALSE
  validObject(.Object)

  return(.Object)
})


## Show method ----

#' Show method for the \code{Selection} class
#'
#' @param object a \code{Selection} object
#'
#' @return No return value, only a display.
#'
#' @author Ehouarn Le Faou
#'
#' @export
setMethod("show", "Selection", function(object) {
  catn("-=-=-=-=-=- SELECTION OJBECT =-=-=-=-=-")
  catn(" #  On individuals: ", c("NO", "YES")[as.integer(object@sOnInds) + 1])
  catn(" #  On gametes: ", c("NO", "YES")[as.integer(object@sOnGams) + 1])
  catn(
    " #  On gamete production: ",
    c("NO", "YES")[as.integer(object@sOnGamsProd) + 1]
  )
  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  catn("(use print to access the fitness values)")
})

## Print method ----

#' Print method for the \code{Selection} class
#'
#' @param x a \code{Selection} object
#' @param ... there are no more parameters.
#'
#' @return No return value, only a display.
#'
#' @author Ehouarn Le Faou
#'
#' @export
setMethod("print", "Selection", function(x, ...) {
  catn("-=-=-=-=-=- SELECTION OJBECT =-=-=-=-=-")
  catn("              in details")
  catn()
  if (!x@sOnInds & !x@sOnGams & !x@sOnGamsProd) {
    catn("No selection defined.")
    catn()
  }
  if (x@sOnInds) {
    tbp <- cbind(t(t(x@indFit[["ind"]])), t(t(x@indFit[["female"]])), t(t(x@indFit[["male"]])))
    colnames(tbp) <- c("Individuals", "Female", "Male")
    print(tbp)
    catn()
  }
  if (x@sOnGams) {
    catn(" #  On gametes: ")
    tbp <- cbind(t(t(x@gamFit[["female"]])), t(t(x@gamFit[["male"]])))
    colnames(tbp) <- c("Female gamete", "Male gamete")
    print(tbp)
    catn()
  }
  if (x@sOnGamsProd) {
    catn(" #  On gamete production: ")
    tbp <- cbind(t(t(x@gamProdFit[["female"]])), t(t(x@gamProdFit[["male"]])))
    colnames(tbp) <- c("Female gamete", "Male gamete")
    print(tbp)
    catn()
  }
  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
})
