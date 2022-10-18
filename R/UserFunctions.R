# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                               User Functions                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Mutation matrix ----


#' Definition of a mutation
#'
#' Utility function to easily generate a mutation matrix
#' (see \link[Ease]{setMutationMatrix}).
#'
#' Mutation occurs from one allele to another at a specific rate. Please
#' take care to define alleles as traits, that these alleles are present
#' in the genome you are using and that the alleles are associated with
#' the same locus.
#'
#' @param from name of the original allele
#' @param to name of the mutant allele
#' @param rate rate at which the mutation occurs
#'
#' @return A standardised list of input parameters that will be used by the
#' function \link[Ease]{setMutationMatrix} to generate the mutation matrix.
#'
#' @examples
#' ### Example with two loci, each with two alleles ###
#'
#' # Definition of the genome
#' DL <- list(dl = c("A", "a"))
#' HL <- list(hl = c("B", "b"))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#'
#' # The mutation function allows each transition from one allele to
#' # another to be defined individually, to produce the mutation matrix
#' # as follows:
#' mutMatrixObj <- setMutationMatrix(genomeObj,
#'   mutations = list(
#'     mutation(from = "A", to = "a", rate = 0.1),
#'     mutation(from = "B", to = "b", rate = 0.1)
#'   )
#' )
#'
#' @export
mutation <- function(from, to, rate) {
  return(list(from = from, to = to, rate = rate))
}


#' Setting the mutation matrix
#'
#' Generation of the mutation matrix associated with the genome given as input.
#' A mutation matrix is used to simulate mutations that affect loci. An object
#' of the class \code{MutationMatrix} does not only contain a (genotypic)
#' mutation matrix. It also contains the attributes necessary for the
#' construction and easy-to-read display of this matrix.
#' The mutation matrix itself is a square matrix of size equal to the number of
#' genotypes. It is a probability matrix in that the sum of the values in
#' each row is equal to 1. For a given genotype, the row associated with it
#' describes the probabilistic proportions that lead by mutation of this
#' genotype to the production of the other genotypes (and of itself if there
#' are no mutations).
#'
#' There are three ways to define the mutation matrix associated with a
#' \code{Genome} class object.
#'
#' 1) By giving two lists of allelic mutation matrices \code{mutHapLoci} and
#' \code{mutDipLoci}, for haploid and diploid loci respectively. Each of
#' these lists contains as many matrices as there are loci. These matrices
#' are transition matrices (squares, with the sum of the rows equal to 1)
#' of size equal to the number of alleles at the locus concerned.
#'
#' 2) By giving a forward and a backward allelic mutation rate
#' (\code{forwardMut} and \code{backwardMut} respectively). The generated
#' mutation matrices will thus be defined with the same rates for all loci. A
#' forward mutation rate means that the transition from one allele to another
#' is done in the order in which they were defined when the Genome class object
#' was created, and in the other direction for the backward rate.
#'
#' 3) By giving a list of \code{mutations} generated through the
#' \link[Ease]{mutation} function.
#'
#' @param genomeObj a \code{Genome} object
#' @param ... see details.
#'
#' @return a \code{MutationMatrix} object
#'
#' @examples
#' ### Example with two loci, each with two alleles ###
#'
#' # Definition of the genome
#' DL <- list(dl = c("A", "a"))
#' HL <- list(hl = c("B", "b"))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#'
#' # Three ways to define the same mutation matrix associated with the
#' # genome defined above:
#'
#' # 1) Mutation matrix from matrices
#' mutHapLoci <- list(matrix(c(0.99, 0.01, 0.01, 0.99), 2))
#' mutDipLoci <- list(matrix(c(0.99, 0.01, 0.01, 0.99), 2))
#' # One can then define the MutationMatrix class object:
#' setMutationMatrix(genomeObj,
#'   mutHapLoci = mutHapLoci,
#'   mutDipLoci = mutDipLoci
#' )
#'
#' # 2) Mutation matrix from mutation rates
#' mutMatrixObj <- setMutationMatrix(genomeObj, forwardMut = 0.1)
#' # or by adding a backward mutation rate:
#' mutMatrixObj <- setMutationMatrix(genomeObj,
#'   forwardMut = 1e-3,
#'   backwardMut = 1e-4
#' )
#'
#' # 3) Mutation matrix from single mutation definition
#' mutMatrixObj <- setMutationMatrix(genomeObj,
#'   mutations = list(
#'     mutation(from = "A", to = "a", rate = 0.1),
#'     mutation(from = "B", to = "b", rate = 0.1)
#'   )
#' )
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setMutationMatrix <- function(genomeObj, ...) {
  args <- list(...)
  argsNames <- names(args)
  if ("mutHapLoci" %in% argsNames | "mutDipLoci" %in% argsNames) {
    if (is.null(args[["mutHapLoci"]])) {
      mutHapLoci <- list()
    } else {
      mutHapLoci <- args[["mutHapLoci"]]
    }
    if (is.null(args[["mutDipLoci"]])) {
      mutDipLoci <- list()
    } else {
      mutDipLoci <- args[["mutDipLoci"]]
    }
    return(new("MutationMatrix", genomeObj, mutHapLoci, mutDipLoci))
  } else if ("forwardMut" %in% argsNames | "backwardMut" %in% argsNames) {
    if (is.null(args[["forwardMut"]])) {
      forwardMut <- 0
    } else {
      forwardMut <- args[["forwardMut"]]
    }
    if (is.null(args[["backwardMut"]])) {
      backwardMut <- 0
    } else {
      backwardMut <- args[["backwardMut"]]
    }

    if (forwardMut < 0 | forwardMut > 1 | backwardMut < 0 | backwardMut > 1) {
      stop("Forward and backward mutation rates should be between 0 and 1.")
    }
    mutHapLoci <- lapply(genomeObj@listHapLoci, function(alleles) {
      mutMatRates(alleles, forwardMut, backwardMut)
    })
    mutDipLoci <- lapply(genomeObj@listDipLoci, function(alleles) {
      mutMatRates(alleles, forwardMut, backwardMut)
    })
    return(new("MutationMatrix", genomeObj, mutHapLoci, mutDipLoci))
  } else if ("mutations" %in% argsNames) {
    rates <- lapply(args[["mutations"]], function(x) x$rate)

    if (any(rates < 0) | any(rates > 1)) {
      stop("Mutation rates should be between 0 and 1.")
    }

    mutMats <- mutMatFriendly(genomeObj, args[["mutations"]])
    return(new(
      "MutationMatrix", genomeObj, mutMats[["mutHapLoci"]],
      mutMats[["mutDipLoci"]]
    ))
  }
  stop("Invalid input. See ?setMutationMatrix")
}

## Genome ----

#' Setting the genome
#'
#' Generation of a genome class object from the list of haploid loci
#' and diploid loci. Each loci is defined by a factor vector that enumerates
#' its alleles.
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
#' @param listHapLoci a list of haploid loci
#' @param listDipLoci a list of diploid loci
#' @param recRate a two-by-two recombination rate vector
#'
#' @return a \code{Genome} object
#'
#' @examples
#' DL <- list(dl = c("A", "a"))
#' HL <- list(hl = c("B", "b"))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setGenome <- function(listHapLoci = list(), listDipLoci = list(),
                      recRate = numeric()) {
  if (all(sapply(listHapLoci, inherits, what = "character"))) {
    listHapLoci <- lapply(listHapLoci, as.factor)
  }
  if (all(sapply(listDipLoci, inherits, what = "character"))) {
    listDipLoci <- lapply(listDipLoci, as.factor)
  }
  return(new("Genome", listHapLoci, listDipLoci, recRate))
}


## Selection ----

#' Setting the selection
#'
#' Generation of a neutral class \code{Selection} object. It can be used as a
#' basis for adding selection layers with the \code{setSelectOnInds},
#' \code{setSelectOnGametes} or \code{setSelectOnGametesProd} functions, or
#' if the model is neutral.
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
#' @param genomeObj a \code{Genome} object
#'
#' @return a \code{Selection} object
#'
#' @examples
#' ### Example with two loci, each with two alleles ###
#' # Definition of the diploid locus
#' DL <- list(dl = c("A", "a"))
#' # Definition of the haploid locus
#' HL <- list(hl = c("B", "b"))
#' # Definition of the object of Genome class
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#' genomeObj
#'
#' ### Exemple with more diploid loci ###
#' # Definition of the diploid loci
#' DL <- list(
#'   dl1 = c("A", "a"),
#'   dl2 = c("B", "b"),
#'   dl3 = c("C", "c")
#' )
#' # Definition of the haploid locus
#' HL <- list(hl = c("D", "d"))
#' # Definition of the object of Genome class, with in addition the necessary
#' # definition of recombination rates between loci:
#' genomeObj <- setGenome(
#'   listHapLoci = HL, listDipLoci = DL,
#'   recRate = c(0.1, 0.5)
#' )
#' # Here we have a 0.1 recombination rate between dl1 and dl2 and a 0.5
#' # recombination rate between dl2 and dl3. It is as if dl1 and dl2 were linked,
#' # for example on the same chromosome, and that dl2 (and dl1 by consequence)
#' # and dl3 were independent, for example on different chromosomes.
#'
#' genomeObj
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setSelectNeutral <- function(genomeObj) {
  return(new("Selection", genomeObj))
}


#' Setting the selection on individuals
#'
#' Generation of an object of the \code{Selection} class which defines a
#' selection among the individuals either by adding this type of selection
#' to an already existing \code{SelectionObj} object (parameter
#' \code{selectionObj}) or by creating one.
#'
#' @param genomeObj a \code{Genome} object
#' @param indFit a genotypic fitness vector for all individuals (whether or not
#' they are hermaphordite)
#' @param femaleFit a genotypic fitness vector for females only (only if the
#' population is dioecious)
#' @param maleFit a genotypic fitness vector for males only (only if the
#' population is dioecious)
#' @param selectionObj a \code{Selection} object (in the case where the
#' selection on individuals is overlaid on an existing \code{Selection} object)
#'
#' @return a \code{Selection} object
#'
#' @examples
#' DL <- list(dl = c("A", "a"))
#' HL <- list(hl = c("B", "b"))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#' selectionObj <- setSelectOnInds(
#'   genomeObj = genomeObj,
#'   indFit = c(1, 1, 1, 1, 0.5, 0)
#' )
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#' @importFrom methods validObject
#'
#' @export
setSelectOnInds <- function(genomeObj = NULL, indFit = c(), femaleFit = c(),
                            maleFit = c(), selectionObj = NULL) {
  if (is.null(selectionObj)) {
    if (is.null(genomeObj)) {
      stop(paste(
        "If no Selection class object is specified here, it is",
        "necessary to enter a Genome class object for the definition",
        "of a new Selection class object."
      ))
    }
    selectionObj <- new("Selection", genomeObj)
  } else {
    if (!is.null(genomeObj)) {
      warning(paste(
        "The Genome class object given as input is not used",
        "because a selection object is specified."
      ))
    }
    genomeObj <- selectionObj@genome
  }
  indFit <- selectInputTreatment(indFit, genomeObj)
  femaleFit <- selectInputTreatment(femaleFit, genomeObj)
  maleFit <- selectInputTreatment(maleFit, genomeObj)

  if (all(sapply(list(indFit, femaleFit, maleFit), length) == 0)) {
    return(selectionObj)
  }
  if (length(indFit) > 0 && length(maleFit) == 0 && length(femaleFit) == 0) {
    selectionObj@indFit[["female"]] <- indFit
    selectionObj@indFit[["male"]] <- indFit
    selectionObj@indFit[["ind"]] <- indFit
  } else {
    if (length(indFit) > 0) {
      selectionObj@indFit[["ind"]] <- indFit
    }
    if (length(femaleFit) > 0) {
      selectionObj@indFit[["female"]] <- femaleFit
    }
    if (length(maleFit) > 0) {
      selectionObj@indFit[["male"]] <- maleFit
    }
  }
  names(selectionObj@indFit[["female"]]) <- selectionObj@IDgenotypes
  names(selectionObj@indFit[["male"]]) <- selectionObj@IDgenotypes
  names(selectionObj@indFit[["ind"]]) <- selectionObj@IDgenotypes
  selectionObj@sOnInds <- TRUE

  validObject(selectionObj)

  return(selectionObj)
}

#' Setting the selection on gametes
#'
#' Generation of an object of the \code{Selection} class which defines a
#' selection among the individuals either by adding this type of selection
#' to an already existing \code{SelectionObj} object (parameter
#' \code{selectionObj}) or by creating one.
#'
#' @param genomeObj a \code{Genome} object
#' @param gamFit an haplotypic fitness vector for all individuals
#' @param femaleFit an haplotypic fitness vector for females only
#' @param maleFit an haplotypic fitness vector for males only
#' @param selectionObj a \code{Selection} object (in the case where the
#' selection on individuals is overlaid on an existing \code{Selection} object)
#'
#' @return a \code{Selection} object
#'
#' @examples
#' DL <- list(dl = c("A", "a"))
#' HL <- list(hl = c("B", "b"))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#' selectionObj <- setSelectOnGametes(
#'   genomeObj = genomeObj,
#'   gamFit = c(1, 1, 0.5, 0)
#' )
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#' @importFrom methods validObject
#'
#' @export
setSelectOnGametes <- function(genomeObj = NULL, gamFit = c(), femaleFit = c(),
                               maleFit = c(), selectionObj = NULL) {
  if (is.null(selectionObj)) {
    if (is.null(genomeObj)) {
      stop(paste(
        "If no Selection class object is specified here, it is",
        "necessary to enter a Genome class object for the definition",
        "of a new Selection class object."
      ))
    }
    selectionObj <- new("Selection", genomeObj)
  } else {
    if (!is.null(genomeObj)) {
      warning(paste(
        "The Genome class object given as input is not used",
        "because a selection object is specified."
      ))
    }
    genomeObj <- selectionObj@genome
  }

  gamFit <- selectInputTreatment(gamFit, genomeObj, haplo = TRUE)
  femaleFit <- selectInputTreatment(femaleFit, genomeObj, haplo = TRUE)
  maleFit <- selectInputTreatment(maleFit, genomeObj, haplo = TRUE)

  if (all(sapply(list(gamFit, femaleFit, maleFit), length) == 0)) {
    return(selectionObj)
  }
  if (((length(gamFit) != 0) & (length(femaleFit) != 0)) |
    ((length(gamFit) != 0) & (length(maleFit) != 0))) {
    warning(paste(
      "A definition of both selection for all gametes",
      "('indFit') and differentiating for male and female gametes",
      "is ambiguous. Only the definition of selection for all",
      "gametes is retained."
    ))
  }
  if (length(gamFit) > 0) {
    selectionObj@gamFit[["female"]] <- gamFit
    selectionObj@gamFit[["male"]] <- gamFit
  } else {
    if (length(femaleFit) > 0) {
      selectionObj@gamFit[["female"]] <- femaleFit
    }
    if (length(maleFit) > 0) {
      selectionObj@gamFit[["male"]] <- maleFit
    }
  }
  names(selectionObj@gamFit[["female"]]) <- selectionObj@IDhaplotypes
  names(selectionObj@gamFit[["male"]]) <- selectionObj@IDhaplotypes
  selectionObj@sOnGams <- TRUE

  validObject(selectionObj)

  return(selectionObj)
}

#' Setting the selection on gamete production
#'
#' Generation of an object of the \code{Selection} class which defines a
#' selection on the gamete production either by adding this type of selection
#' to an already existing \code{SelectionObj} object (parameter
#' \code{selectionObj}) or by creating one.
#'
#' @param genomeObj a \code{Genome} object
#' @param indProdFit a genotypic fitness vector for all individuals
#' @param femProdFit a genotypic fitness vector for females only
#' @param maleProdFit a genotypic fitness vector for males only
#' @param selectionObj a \code{Selection} object (in the case where the
#' selection on individuals is overlaid on an existing \code{Selection} object)
#'
#' @return a \code{Selection} object
#'
#' @examples
#' DL <- list(dl = c("A", "a"))
#' HL <- list(hl = c("B", "b"))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#' selectionObj <- setSelectOnGametesProd(
#'   genomeObj = genomeObj,
#'   indProdFit = c(1, 1, 1, 1, 0.5, 0)
#' )
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#' @importFrom methods validObject
#'
#' @export
setSelectOnGametesProd <- function(genomeObj = NULL, indProdFit = c(),
                                   femProdFit = c(), maleProdFit = c(),
                                   selectionObj = NULL) {
  if (is.null(selectionObj)) {
    if (is.null(genomeObj)) {
      stop(paste(
        "If no Selection class object is specified here, it is",
        "necessary to enter a Genome class object for the definition",
        "of a new Selection class object."
      ))
    }
    selectionObj <- new("Selection", genomeObj)
  } else {
    if (!is.null(genomeObj)) {
      warning(paste(
        "The Genome class object given as input is not used",
        "because a selection object is specified."
      ))
    }
    genomeObj <- selectionObj@genome
  }

  indProdFit <- selectInputTreatment(indProdFit, genomeObj)
  femProdFit <- selectInputTreatment(femProdFit, genomeObj)
  maleProdFit <- selectInputTreatment(maleProdFit, genomeObj)

  if (all(sapply(list(indProdFit, femProdFit, maleProdFit), length) == 0)) {
    return(selectionObj)
  }
  if (((length(indProdFit) != 0) & (length(femProdFit) != 0)) |
    ((length(indProdFit) != 0) & (length(maleProdFit) != 0))) {
    stop(paste(
      "Prohibited definition of both selection at the individual",
      "level and for at least one of the two sexes."
    ))
  }
  if (length(indProdFit) > 0) {
    selectionObj@gamProdFit[["female"]] <- indProdFit
    selectionObj@gamProdFit[["male"]] <- indProdFit
  } else {
    if (length(femProdFit) > 0) {
      selectionObj@gamProdFit[["female"]] <- femProdFit
    }
    if (length(maleProdFit) > 0) {
      selectionObj@gamProdFit[["male"]] <- maleProdFit
    }
  }
  names(selectionObj@gamProdFit[["female"]]) <- selectionObj@IDgenotypes
  names(selectionObj@gamProdFit[["male"]]) <- selectionObj@IDgenotypes
  selectionObj@sOnGamsProd <- TRUE
  validObject(selectionObj)
  return(selectionObj)
}


## Population ----

#' Setting a population
#'
#' Generation of a population by providing all the necessary ingredients for
#' its definition, including a genome, a mutation matrix and a selection regime.
#'
#' A population is defined strictly by a name, a size, a sexual system
#' (dioecy or hermaphodite), and the three objects defined previously:
#' genome, mutation matrix and selection. In addition to that, it is
#' possible to define
#'  - a selfing rate (by default equal to 0)
#'  - a vector of initial genotypic frequencies
#'  - a demography
#'
#' Two demographic regimes are possible: no demography, i.e. a fixed population
#' size, or demography, i.e. a population where the size fluctuates
#' stochastically. The boolean argument `demography` is used to define whether
#' there should be stochasticity. For a fixed population size, it is therefore
#' sufficient to define that `demography = FALSE` (default) and to set the
#' desired population size with the `popSize` parameter.
#'
#' For a fluctuating demography, `demography` must be `TRUE` and three other
#' parameters are then needed: the initial population size (`initPopSize`),
#' the population growth rate (`growthRate`) and the carrying capacity of the
#' population (the population size, `popSize`).
#'
#' It is also possible to avoid defining a population size altogether, by
#' setting off the genetic drift (`drift` parameter). This will allow the
#' model to be simulated deterministically.
#'
#' @param name the name of the population
#' @param size the population size
#' @param dioecy logical indicating whether the simulated population is
#' dioecious or hermaphroditic
#' @param genomeObj a \code{Genome} object
#' @param mutMatrixObj a \code{MutationMatrix} object
#' @param selectionObj a \code{Selection} object
#' @param selfRate the selfing rate
#' @param demography a logic indicating whether the population should have
#' a demography (stochasticity in the number of individuals present in the
#' population + logistic growth with carrying capacity equal to the \code{size}
#' parameter)
#' @param growthRate a \code{Genome} object
#' @param initPopSize the initial size of the population. It is necessarily
#' equal to \code{size} if the population has no \code{demography}.
#' @param initGenoFreq a vector of the size of the genotype number
#' describing the initial allele frequencies common to all simulations
#'
#' @return a \code{Population} object
#'
#' @examples
#' # Definition of a population in its simplest form:
#' DL <- list(dl = c("A", "a"))
#' HL <- list(hl = c("B", "b"))
#' mutations <- list(
#'   mutation(from = "A", to = "a", rate = 1e-3),
#'   mutation(from = "B", to = "b", rate = 1e-3)
#' )
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#' pop <- setPopulation(
#'   name = "A",
#'   size = 1000,
#'   dioecy = TRUE,
#'   genomeObj = genomeObj,
#'   selectionObj = setSelectNeutral(genomeObj),
#'   mutMatrixObj = setMutationMatrix(genomeObj, mutations = mutations)
#' )
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setPopulation <- function(name, size, dioecy, genomeObj, mutMatrixObj, selectionObj,
                          selfRate = 0, demography = F, growthRate = 0,
                          initPopSize = NULL, initGenoFreq = NULL) {
  if (is.null(initPopSize)) {
    if (demography) {
      warning(paste(
        "As demography is enabled but the population size is not",
        "set, it is defined by default as equal to the regular",
        "population size (size parameter)."
      ))
    }
    initPopSize <- size
  } else {
    if (!demography & (initPopSize != size)) {
      warning(paste(
        "Without demography, the initial population size",
        "(initPopSize parameter) is necessarily",
        "equal to the regular population size (size",
        "parameter)."
      ))
      initPopSize <- size
    }
  }

  return(new(
    "Population", name, size, dioecy, selfRate, demography, growthRate,
    initGenoFreq, genomeObj, initPopSize, selectionObj, mutMatrixObj
  ))
}

## Metapopulation ----

#' Setting a metapopulation
#'
#' A metapopulation is a set of population(s) (from 1) that are simulated
#' with potential migration between them. Only genotypes can migrate, i.e.
#' adult individuals.
#'
#' The construction of a \code{Metapopulation} object requires only two
#' arguments (one optional). The first is a population(s) list, defined
#' from the population class. The second is a migration matrix, which
#' connects the populations together. This matrix is a probability matrix
#' (square with the sum of the rows equal to 1, whose size is equal to the
#' number of populations) where each value corresponds to the proportion
#' of individuals (genotypes) that disperse from their source population
#' (row) to their target population (column).
#'
#' @param populations a list of \code{Population} objects
#' @param migMat a migration matrix
#'
#' @return a \code{Metapopulation} object
#'
#' @examples
#' # Definition of a population in its simplest form:
#' DL <- list(dl = c("A", "a"))
#' HL <- list(hl = c("B", "b"))
#' mutations <- list(
#'   mutation(from = "A", to = "a", rate = 1e-3),
#'   mutation(from = "B", to = "b", rate = 1e-3)
#' )
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
#' pop <- setPopulation(
#'   name = "A",
#'   size = 1000,
#'   dioecy = TRUE,
#'   genomeObj = genomeObj,
#'   selectionObj = setSelectNeutral(genomeObj),
#'   mutMatrixObj = setMutationMatrix(genomeObj, mutations = mutations)
#' )
#' metapop <- setMetapopulation(populations = list(pop))
#' metapop <- simulate(metapop, nsim = 10, seed = 123)
#' # Other examples available in the documentation of the package
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setMetapopulation <- function(populations, migMat = matrix(1)) {
  return(new("Metapopulation", populations, migMat))
}

#' Getting the simulation results
#'
#' @param metapop a \code{Metapopulation} objects
#'
#' @return A data.frame where each line corresponds to a simulation. The
#' results include :
#' - the last generation reached (the threshold or the generation that first
#' verified at least one of the stopping conditions)
#' - the final population size
#' - the genotype frequencies
#' - allelic frequencies
#' - the reason(s) for the stop (either the threshold was reached, i.e.
#' \code{unstopped} or the stop condition(s) that was (were) reached, in the
#' form of boolean values
#' - Average fitness (individual, gamete production and gametic)
#'
#' @author Ehouarn Le Faou
#'
#' @export
getResults <- function(metapop) {
  if (identical(dim(metapop@results), c(0L, 0L))) {
    warning(paste(
      "No results available for this metapopulation. Please use",
      "the simulate method to generate them."
    ))
  } else {
    return(metapop@results)
  }
}

#' Getting the simulation results
#'
#' @param metapop a \code{Metapopulation} objects
#'
#' @return A list where each item is associated with a simulation. Each of
#' these elements consists of a list of data.frames, one per population.
#' These data.frames consist of the same columns as the results
#' (see \link[Ease]{getResults} documentation), except that they do not
#' include the stop conditions.
#'
#' @author Ehouarn Le Faou
#'
#' @export
getRecords <- function(metapop) {
  if (identical(metapop@records, list())) {
    warning(paste(
      "No records available for this metapopulation. Please use",
      "the simulate method and activate the recording parameter",
      "to generate them."
    ))
  } else {
    return(metapop@records)
  }
}


#' Getting the custom output
#'
#' @param metapop a \code{Metapopulation} objects
#'
#' @return The list generated through the custom result function, if at least
#' it was specified during the simulation of the \code{Metapopulation}.
#'
#' @author Ehouarn Le Faou
#'
#' @export
getCustomOutput <- function(metapop) {
  if (is.null(unlist(metapop@customOutput))) {
    warning(paste(
      "No custom output available for this metapopulation. Please use",
      "the simulate method and use the nameOutFunct parameter",
      "to generate them."
    ))
  } else {
    return(metapop@customOutput)
  }
}
