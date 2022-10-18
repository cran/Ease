# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                              Population CLASS                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Validity check ----

#' The validity check for the \code{Population} class
#'
#' @param object a \code{Population} object
#'
#' @return a boolean corresponding to whether the object is a correct
#' \code{Population} object.
#'
#' @author Ehouarn Le Faou
#'
check.population <- function(object) {
  if (object@size < 0 | object@size != as.integer(object@size)) {
    stop("The population size must be a positive integer.")
  }
  if (object@size == 0) {
    stop("A necessarily empty population (size = 0) is not valid.")
  }
  if (object@initPopSize < 0 | object@initPopSize != as.integer(object@initPopSize)) {
    stop("The initial populaton size must be a positive integer.")
  }
  if (object@selfRate > 1 | object@selfRate < 0) {
    stop("The selfing rate should be between 0 and 1.")
  }
  if (object@growthRate < 0) {
    stop(paste(
      "The growth rate of a population is a positive real."
    ))
  }
  if (object@selection@IDgenome != object@genome@IDgenome) {
    stop(paste(
      "The genome used to build the Selection object and the",
      "genome given as input do not match."
    ))
  }
  if (object@mutMat@IDgenome != object@genome@IDgenome) {
    stop(paste(
      "The genome used to build the MutationMatrix object and the",
      "genome given as input do not match."
    ))
  }
}


## Class definition ----

#' Population
#'
#' The \code{Population} class allows for the collection of the parameters necessary
#' to characterise a biological population. It is an essentially useful class
#' in that no method associated with the population class can simulate its
#' dynamics. To do this, it is necessary to use the Metapopulation class,
#' which takes as input a list of populations (from one). The Population
#' class is also used to check that each of these parameters is compatible
#' with each other.
#'
#' Thus to build an object of class \code{Ease}, it is necessary to have
#' defined an object \code{Genome}, as well as an object \code{MutationMatrix}
#' and an object \code{Selection} (even if it is neutral, see
#' \link[Ease]{setSelectNeutral}).
#'
#' @slot name the name of the population.
#' @slot size the size of the population.
#' @slot dioecy logical indicating whether the population is dioecious or not
#' (hermaphrodite).
#' @slot selfRate the selfing rate of the population
#' @slot demography logical indicating whether the population has stochastic
#' demography (this does not include migration), i.e. non-constant size and
#' potentially population growth or decay, depending on the situation it is in.
#' @slot growthRate growth rate of the population.
#' @slot initGenoFreq A row matrix of the size of the genotype number
#' describing the initial allele frequencies common to all simulations
#' @slot genome a \code{Genome} object
#' @slot initPopSize initial population size, knowing that if the demography
#' is extinct, the initial population size will automatically be set equal to
#' the population size.
#' @slot selection a \code{Selection} object
#' @slot mutMat a \code{MutationMatrix} object
#'
#' @author Ehouarn Le Faou
#'
#' @export
setClass("Population",
  representation(
    name = "character",
    size = "numeric",
    dioecy = "logical",
    selfRate = "numeric",
    demography = "logical",
    growthRate = "numeric",
    initGenoFreq = "matrix",
    genome = "Genome",
    initPopSize = "numeric",
    selection = "Selection",
    mutMat = "MutationMatrix"
  ),
  validity = check.population
)


## Initialize method ----

#' Initialize method for the \code{Population} class
#'
#' @param .Object a \code{Population} object
#' @param name the name of the population.
#' @param size the size of the population.
#' @param dioecy logical indicating whether the population is dioecious or not
#' (hermaphrodite).
#' @param selfRate the selfing rate of the population
#' @param demography logical indicating whether the population has stochastic
#' demography (this does not include migration), i.e. non-constant size and
#' potentially population growth or decay, depending on the situation it is in.
#' @param growthRate growth rate of the population.
#' @param initGenoFreq A row matrix of the size of the genotype number
#' describing the initial allele frequencies common to all simulations
#' @param genomeObj a \code{Genome} object
#' @param initPopSize initial population size, knowing that if the demography
#' is extinct, the initial population size will automatically be set equal to
#' the population size.
#' @param selectionObj a \code{Selection} object
#' @param mutMatrixObj a \code{MutationMatrix} object
#'
#' @return a \code{Population} object
#'
#' @author Ehouarn Le Faou
#'
setMethod("initialize", "Population", function(.Object, name, size, dioecy,
                                               selfRate, demography, growthRate,
                                               initGenoFreq, genomeObj,
                                               initPopSize, selectionObj,
                                               mutMatrixObj) {

  # Definition of attributes
  .Object@name <- name
  .Object@size <- size
  .Object@selfRate <- selfRate
  .Object@dioecy <- dioecy
  .Object@demography <- demography
  .Object@growthRate <- growthRate
  .Object@genome <- genomeObj
  .Object@initPopSize <- initPopSize
  .Object@selection <- selectionObj
  .Object@mutMat <- mutMatrixObj

  if (is.null(initPopSize)) .Object@initPopSize <- size

  # Initial genotype frequencies
  if (is.null(initGenoFreq)) {
    initGenoFreq <- matrix(c(1, rep(0, .Object@genome@nbGeno - 1)), 1)
  } else {
    if (length(initGenoFreq) != .Object@genome@nbGeno) {
      stop(paste0(
        "The vector of initial genotype frequencies must be the ",
        "same length as the number of genotypes in the input ",
        "genome (", .Object@genome@nbGeno, ")."
      ))
    }
    if (round(sum(initGenoFreq), 10) != 1) {
      stop(paste(
        "The sum of the initial genotype frequencies must be",
        "equal to 1."
      ))
    }
    initGenoFreq <- t(as.matrix(initGenoFreq))
  }
  .Object@initGenoFreq <- initGenoFreq
  colnames(.Object@initGenoFreq) <- .Object@genome@IDgenotypes

  # Warnings
  if (!.Object@demography) {
    if (.Object@growthRate != 0) {
      warning(paste(
        "The population growth rate is not taken into account in",
        "a population where demography is not enabled."
      ))
      .Object@growthRate <- 0
    }
    if (.Object@initPopSize != .Object@size) {
      warning(paste(
        "The initial population size is necessarily equal to",
        "the population size when demography is not enabled."
      ))
      .Object@initPopSize <- .Object@size
    }
  }

  # Validity of the object
  validObject(.Object)

  return(.Object)
})


## Show method ----

#' Show method for the \code{Population} class
#'
#' @param object a \code{Population} object
#'
#' @author Ehouarn Le Faou
#'
setMethod("show", "Population", function(object) {
  catn("-=-=-=-=-= Population OBJECT =-=-=-=-=-")
  cat(paste0(
    "Population '", object@name, "' of ", object@size, " ",
    c("hermaphroditic", "dioecious")[as.integer(object@dioecy) + 1],
    " individuals"
  ))
  if (!object@dioecy) {
    cat("\n   with a ")
    tbp <- round(object@selfRate * 100, 2)
    if (tbp == object@selfRate * 100) {
      catn(paste0(tbp, "% selfing rate."))
    } else {
      catn(paste0("~", tbp, "% selfing rate."))
    }
  } else {
    catn()
  }
  if (!object@demography) {
    catn("There is no demography.")
  } else {
    catn(paste0(
      "The demography is active, with a growth rate of ",
      object@growthRate, "."
    ))
  }
  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  catn("(use print for details on the genome)")
})


## Print method ----

#' Print method for the \code{Population} class
#'
#' @param x a \code{Population} object
#' @param ... the other parameter is \code{frame}, which is a logic indicating
#' whether the frame surrounding the display of the population characteristics
#' should be displayed or not.
#'
#' @author Ehouarn Le Faou
#'
setMethod("print", "Population", function(x, ...) {
  args <- list(...)
  if ("frame" %in% names(args)) {
    if (args$frame) {
      catn("-=-=-=-=-= Population OBJECT =-=-=-=-=-")
      catn("              in details")
    }
  } else {
    catn("-=-=-=-=-= Population OBJECT =-=-=-=-=-")
    catn("              in details")
  }
  cat(paste0(
    "Population '", x@name, "' of ", x@size, " ",
    c("hermaphroditic", "dioecious")[as.integer(x@dioecy) + 1],
    " individuals"
  ))
  if (!x@dioecy) {
    cat("\n   with a ")
    tbp <- round(x@selfRate * 100, 2)
    if (tbp == x@selfRate * 100) {
      catn(paste0(tbp, "% selfing rate."))
    } else {
      catn(paste0("~", tbp, "% selfing rate."))
    }
  } else {
    catn()
  }
  if (!x@demography) {
    catn("There is no demography.")
  } else {
    catn(paste0(
      "The demography is active, with a growth rate of ",
      x@growthRate, " and an initial population size of ", x@initPopSize, "."
    ))
  }
  catn("The initial genotypes frequency are: ")
  tbp <- as.vector(x@initGenoFreq)
  names(tbp) <- x@genome@IDgenotypes
  print(tbp)

  catn("Selection: ")
  if (!x@selection@sOnInds & !x@selection@sOnGams & !x@selection@sOnGamsProd) {
    catn("No selection defined.")
  }
  if (x@selection@sOnInds) {
    if (length(x@selection@indFit[["ind"]]) > 0) {
      catn("   - On individuals: ")
      tbp <- t(t(x@selection@indFit[["ind"]]))
      colnames(tbp) <- "Fitness"
      print(tbp)
    } else {
      tbp <- cbind(t(t(x@selection@indFit[["female"]])), t(t(x@selection@indFit[["male"]])))
      colnames(tbp) <- c("Female", "Male")
      print(tbp)
    }
  }
  if (x@selection@sOnGams) {
    catn("   - On gametes: ")
    tbp <- cbind(t(t(x@selection@gamFit[["female"]])), t(t(x@selection@gamFit[["male"]])))
    colnames(tbp) <- c("Female gamete", "Male gamete")
    print(tbp)
  }
  if (x@selection@sOnGamsProd) {
    catn("   - On gamete production: ")
    tbp <- cbind(t(t(x@selection@gamProdFit[["female"]])), t(t(x@selection@gamProdFit[["male"]])))
    colnames(tbp) <- c("Female gamete", "Male gamete")
    print(tbp)
  }


  if ("frame" %in% names(args)) {
    if (args$frame) {
      catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
    }
  } else {
    catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  }
})
