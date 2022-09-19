# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                               User Functions                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#' Setting the mutation matrix
#'
#' Generation of the mutation matrix associated with the genome given as input
#' and specifying the mutation matrices locus by locus by the arguments
#' \code{mutHapLoci} and \code{mutDipLoci}.
#'
#' A mutation matrix is used to simulate mutations that affect loci. An object
#' of the class \code{MutationMatrix} does not only contain a (genotypic)
#' mutation matrix. It also contains the attributes necessary for the
#' construction and easy-to-read display of this matrix.
#'
#' The mutation matrix itself is a square matrix of size equal to the number of
#' genotypes. It is a probability matrix in that the sum of the values in
#' each row is equal to 1. For a given genotype, the row associated with it
#' describes the probabilistic proportions that lead by mutation of this
#' genotype to the production of the other genotypes (and of itself if there
#' are no mutations).
#'
#' @param genomeObj a \code{Genome} object
#' @param mutHapLoci a list of haploid locus by locus allelic mulation matrices.
#' @param mutDipLoci a list of diploid locus by locus allelic mulation matrices.
#'
#' @return a \code{MutationMatrix} object
#'
#' @examples
#' ### Example with two loci, each with two alleles ###
#' # Definition of the diploid locus
#' LD <- list(dl = as.factor(c("A", "a")))
#' # Definition of the haploid locus
#' HL <- list(hl = as.factor(c("B", "b")))
#' # Definition of the object of Genome class
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#' # The mutation matrices can be defined as follows:
#' mutHapLoci <- list(matrix(c(0.9, 0.1, 0.1, 0.9), 2))
#' mutDipLoci <- list(matrix(c(0.99, 0.01, 0.01, 0.99), 2))
#' # One can then define the MutationMatrix class object:
#' setMutationMatrix(genomeObj, mutHapLoci, mutDipLoci)
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setMutationMatrix <- function(genomeObj, mutHapLoci, mutDipLoci) {
  return(new("MutationMatrix", genomeObj, mutHapLoci, mutDipLoci))
}

#' Setting the mutation matrix by rates
#'
#' Generation of the mutation matrix associated with the genome given as input
#' and specifying the forward and backward mutation rates.
#'
#' A mutation matrix is used to simulate mutations that affect loci. An object
#' of the class \code{MutationMatrix} does not only contain a (genotypic)
#' mutation matrix. It also contains the attributes necessary for the
#' construction and easy-to-read display of this matrix.
#'
#' The mutation matrix itself is a square matrix of size equal to the number of
#' genotypes. It is a probability matrix in that the sum of the values in
#' each row is equal to 1. For a given genotype, the row associated with it
#' describes the probabilistic proportions that lead by mutation of this
#' genotype to the production of the other genotypes (and of itself if there
#' are no mutations).
#'
#' @param genomeObj a \code{Genome} object
#' @param forwardMut the forward mutation rate
#' @param backwardMut the backward mutation rate
#'
#' @return a \code{MutationMatrix} object
#'
#' @examples
#' ### Example with two loci, each with two alleles ###
#' # Definition of the diploid locus
#' LD <- list(dl = as.factor(c("A", "a")))
#' # Definition of the haploid locus
#' HL <- list(hl = as.factor(c("B", "b")))
#' # Definition of the object of Genome class
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#' # The mutation matrices can be defined as follows:
#' mutMatrixObj <- setMutationMatrixByRates(genomeObj, forwardMut = 1e-3)
#' mutMatrixObj
#'
#' # One can also add a backward mutation rate:
#' mutMatrixObj <- setMutationMatrixByRates(genomeObj,
#'   forwardMut = 1e-3,
#'   backwardMut = 1e-4
#' )
#' mutMatrixObj
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setMutationMatrixByRates <- function(genomeObj, forwardMut = 0,
                                     backwardMut = 0) {
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
}


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
#' LD <- list(dl = as.factor(c("A", "a")))
#' HL <- list(hl = as.factor(c("B", "b")))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setGenome <- function(listHapLoci = list(), listDipLoci = list(),
                      recRate = numeric()) {
  return(new("Genome", listHapLoci, listDipLoci, recRate))
}

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
#' LD <- list(dl = as.factor(c("A", "a")))
#' # Definition of the haploid locus
#' HL <- list(hl = as.factor(c("B", "b")))
#' # Definition of the object of Genome class
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#' genomeObj
#'
#' ### Exemple with more diploid loci ###
#' # Definition of the diploid loci
#' LD <- list(
#'   dl1 = as.factor(c("A", "a")),
#'   dl2 = as.factor(c("B", "b")),
#'   dl3 = as.factor(c("C", "c"))
#' )
#' # Definition of the haploid locus
#' HL <- list(hl = as.factor(c("D", "d")))
#' # Definition of the object of Genome class, with in addition the necessary
#' # definition of recombination rates between loci:
#' genomeObj <- setGenome(
#'   listHapLoci = HL, listDipLoci = LD,
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
#' LD <- list(dl = as.factor(c("A", "a")))
#' HL <- list(hl = as.factor(c("B", "b")))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
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
  }
  if (all(sapply(list(indFit, femaleFit, maleFit), length) == 0)) {
    return(selectionObj)
  }
  if (length(indFit) > 0 && length(maleFit) == 0 && length(femaleFit) == 0) {
    selectionObj@femindFit <- indFit
    selectionObj@maleindFit <- indFit
    selectionObj@indFit <- indFit
  } else {
    if (length(indFit) > 0) {
      selectionObj@indFit <- indFit
    }
    if (length(femaleFit) > 0) {
      selectionObj@femindFit <- femaleFit
    }
    if (length(maleFit) > 0) {
      selectionObj@maleindFit <- maleFit
    }
  }
  names(selectionObj@femindFit) <- selectionObj@IDgenotypes
  names(selectionObj@maleindFit) <- selectionObj@IDgenotypes
  names(selectionObj@indFit) <- selectionObj@IDgenotypes
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
#' LD <- list(dl = as.factor(c("A", "a")))
#' HL <- list(hl = as.factor(c("B", "b")))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
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
  }

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
    selectionObj@femgamFit <- gamFit
    selectionObj@malegamFit <- gamFit
  } else {
    if (length(femaleFit) > 0) {
      selectionObj@femgamFit <- femaleFit
    }
    if (length(maleFit) > 0) {
      selectionObj@malegamFit <- maleFit
    }
  }
  names(selectionObj@femgamFit) <- selectionObj@IDhaplotypes
  names(selectionObj@malegamFit) <- selectionObj@IDhaplotypes
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
#' LD <- list(dl = as.factor(c("A", "a")))
#' HL <- list(hl = as.factor(c("B", "b")))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
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
  }
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
    selectionObj@femProdFit <- indProdFit
    selectionObj@maleProdFit <- indProdFit
  } else {
    if (length(femProdFit) > 0) {
      selectionObj@femProdFit <- femProdFit
    }
    if (length(maleProdFit) > 0) {
      selectionObj@maleProdFit <- maleProdFit
    }
  }
  names(selectionObj@femProdFit) <- selectionObj@IDgenotypes
  names(selectionObj@maleProdFit) <- selectionObj@IDgenotypes
  selectionObj@sOnGamsProd <- TRUE
  validObject(selectionObj)
  return(selectionObj)
}


#' Setting the parameters for simulating a model
#'
#' Last step before the simulation of the model, the creation of a \code{Ease}
#' class object makes it possible to gather all the ingredients for a complete
#' parameterization of a model and then to simulate it.
#'
#' The \code{Ease} class is used to manage the simulations by handling the
#' objects needed to build the model. Thus to build an object of class
#' \code{Ease}, it is necessary to have defined an object \code{Genome},
#' as well as an object \code{MutationMatrix} and an object \code{Selection}
#' (even if it is neutral, see \link[Ease]{setSelectNeutral}).
#'
#' Once simulated with the method \code{simulate}, an Ease object contains
#' the results of the simulations and the records of these last ones if the
#' parameter \code{recording} of \code{simulate} was fixed at \code{TRUE}.
#' To obtain these results and records, it is necessary to use the functions
#' \link[Ease]{getResults} and \link[Ease]{getRecords}.
#'
#' @param N the population size
#' @param threshold the maximum number of generations
#' @param dioecy logical indicating whether the simulated population is
#' dioecious or hermaphroditic
#' @param mutMatrixObj a \code{MutationMatrix} object
#' @param genomeObj a \code{Genome} object
#' @param selectionObj a \code{Selection} object
#' @param stopCondition list of vectors that each describe the alleles that
#' must be fixed to define a stop condition. Each of these stop conditions
#' will therefore be associated with a stop condition
#' @param initGenoFreq A vector of the size of the genotype number
#' describing the initial allele frequencies common to all simulations
#' @param selfRate the selfing rate
#'
#' @return an \code{Ease} object
#'
#' @examples
#'
#' library(tidyr)
#' library(ggplot2)
#'
#' ### Simple example of a succession of allele replacements to each other
#' # in a deterministic way and simply by mutations (selection is neutral) ###
#'
#' # Let's put a single diploid locus with 5 alleles:
#' LD <- list(dl1 = as.factor(c("a1", "a2", "a3", "a4", "a5", "a6")))
#' HL <- list(hl = as.factor("noHL"))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#'
#' # The only possible mutations are mutations from the a1 to a2 allele,
#' # from a2 to a3, etc.:
#' mutMatrixObj <- setMutationMatrixByRates(
#'   genomeObj = genomeObj,
#'   forwardMut = 1e-2
#' )
#' # The matrix thus constructed looks like this:
#' mutMatrixObj
#'
#' # The selection is neutral:
#' selectionObjNeutral <- setSelectNeutral(genomeObj = genomeObj)
#'
#' # We can thus define an Ease object with a population size of 100
#' # sex-separated individuals and a threshold of 700 generations:
#' mod <- setEase(
#'   N = 100, threshold = 700, dioecy = TRUE,
#'   mutMatrixObj = mutMatrixObj,
#'   genomeObj = genomeObj,
#'   selectionObj = selectionObjNeutral
#' )
#'
#' # For the simulation we shut down the drift to have a deterministic
#' # evolution of allelic frequencies:
#' mod <- simulate(mod, recording = TRUE, verbose = TRUE, drift = FALSE)
#'
#' # We recover the \code{data.frame} of the simulation record and we modify
#' # a little the organization of the data:
#' records <- getRecords(mod)[[1]]
#' records <- gather(records, "allele", "freqAllele", 45:50)
#'
#' # Then we display the evolution of the allelic frequencies of each of the
#' # alleles of the locus:
#' ggplot(records, aes(x = gen, y = freqAllele, color = allele)) +
#'   geom_line() +
#'   ylim(0, 1)
#'
#'
#' ### Example of simulation of cyto-nuclear Bateson-Dobzhansky-Muller
#' # incompatibilities (BDMI) ###
#'
#' # Two loci: a haploid locus and a diploid locus. Each has two alleles,
#' # an ancestral allele in upper case and a derived allele in lower case:
#' LD <- list(dl = as.factor(c("A", "a")))
#' HL <- list(hl = as.factor(c("B", "b")))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#'
#' # The mutation rate to derived alleles is set at 1e-2:
#' mutMatrixObj <- setMutationMatrixByRates(genomeObj = genomeObj, forwardMut = 1e-2)
#'
#' # The two derived alleles a and b are incompatible and therefore impose
#' # a fitness cost on their carrier:
#' selectionObj <- setSelectOnInds(
#'   genomeObj = genomeObj,
#'   indFit = c(1, 1, 1, 1, 0.5, 0)
#' )
#'
#' # We can then define the Ease object by specifying the population size (100),
#' # the maximum generation threshold (1e6), that we want the individuals to be
#' # hermaphroditic (dioecy = FALSE) and that they reproduce at 50% by
#' # selfing. The simulation stops if one of the two derived alleles is fixed:
#' mod <- setEase(
#'   N = 100, threshold = 1e6, dioecy = FALSE, selfRate = 0.5,
#'   stopCondition = list("a", "b"),
#'   mutMatrixObj = mutMatrixObj,
#'   genomeObj = genomeObj,
#'   selectionObj = selectionObj
#' )
#'
#' # The model is simulated:
#' mod <- simulate(mod, nsim = 10, verbose = TRUE)
#'
#' # And the results plotted:
#' plot(mod)
#'
#' @author Ehouarn Le Faou
#'
#' @importFrom methods new
#'
#' @export
setEase <- function(N, threshold, dioecy, mutMatrixObj, genomeObj,
                    selectionObj, stopCondition = list(),
                    initGenoFreq = matrix(), selfRate = NA_real_) {
  return(new(
    "Ease", N, threshold, dioecy, mutMatrixObj, genomeObj,
    selectionObj, stopCondition, selfRate, initGenoFreq
  ))
}


#' Retrieving simulation results
#'
#' A simple function to retrieve the results of a simulation of an
#' \code{Ease} object. The results can be given as a list of \code{data.frames}
#' distinguishing: the parameters, the genotypic frequencies (unless
#' they have not been recorded), the generations, and the stop conditions.
#'
#' @param easeObj an \code{Ease} object
#' @param asList logical defining whether the results should be output as a list
#'
#' @return A \code{data.frame} corresponding to the results of the simulations.
#' By default, the results are presented in the form of \code{data.frame}, each
#' row corresponding to a simulation. The columns correspond, in order, to:
#' - the parameters of the simulation: the population size, the generation
#'   threshold, the sexual system and finally the self-fertilisation rate (also
#'   present in dioecy, but always equal to 0);
#' - the genotypic frequencies at the end of the simulation (of females,
#'   then males, or only individuals in the case of hermaphrodites)
#' - the allelic frequencies at the end of the simulation;
#' - the generation at which the simulation stopped;
#' - the stop conditions, one column being dedicated to each of them, they
#'   indicate by logic if the simulation has reached these stop conditions;
#' - the fitness of the individuals at the end of the simulation.
#' The parameters \code{includeParams} and \code{includeFitness} of the
#' \code{simulate} method can delete the input parameters and fitness
#' of the result \code{data.frame} respectively.
#'
#' @examples
#' LD <- list(dl = as.factor(c("A", "a")))
#' HL <- list(hl = as.factor(c("B", "b")))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#' selectionObj <- setSelectOnInds(
#'   genomeObj = genomeObj,
#'   indFit = c(1, 1, 1, 1, 0.5, 0)
#' )
#' mutMatrixObj <- setMutationMatrixByRates(genomeObj = genomeObj, forwardMut = 1e-2)
#' mod <- setEase(
#'   N = 100, threshold = 1e6, dioecy = FALSE, selfRate = 0.5,
#'   stopCondition = list("a", "b"),
#'   mutMatrixObj = mutMatrixObj,
#'   genomeObj = genomeObj,
#'   selectionObj = selectionObj
#' )
#'
#' mod <- simulate(mod, nsim = 10)
#'
#' # As a single \code{data.frame} :
#' getResults(mod)
#'
#' # As a list :
#' getResults(mod, asList = TRUE)
#'
#' @author Ehouarn Le Faou
#'
#' @export
getResults <- function(easeObj, asList = FALSE) {
  if (identical(easeObj@results, list())) {
    stop(paste(
      "No results available for this object (use the simulate",
      "method to generate them)."
    ))
  }
  if (asList) {
    return(easeObj@results)
  } else {
    return(Reduce(cbind, easeObj@results))
  }
}




#' Retrieving simulation records
#'
#' A simple function to retrieve the records of a simulation of an
#' \code{Ease} object.
#'
#' @param easeObj an \code{Ease} object
#'
#' @return A list of \code{data.frame} corresponding to the records of the
#' simulations.
#'
#' @examples
#' LD <- list(dl = as.factor(c("A", "a")))
#' HL <- list(hl = as.factor(c("B", "b")))
#' genomeObj <- setGenome(listHapLoci = HL, listDipLoci = LD)
#' selectionObj <- setSelectOnInds(
#'   genomeObj = genomeObj,
#'   indFit = c(1, 1, 1, 1, 0.5, 0)
#' )
#' mutMatrixObj <- setMutationMatrixByRates(genomeObj = genomeObj, forwardMut = 1e-2)
#' mod <- setEase(
#'   N = 100, threshold = 1e6, dioecy = FALSE, selfRate = 0.5,
#'   stopCondition = list("a", "b"),
#'   mutMatrixObj = mutMatrixObj,
#'   genomeObj = genomeObj,
#'   selectionObj = selectionObj
#' )
#'
#' mod <- simulate(mod, recording = TRUE, nsim = 10)
#'
#' getRecords(mod)
#'
#' @author Ehouarn Le Faou
#'
#' @export
getRecords <- function(easeObj) {
  if (identical(easeObj@records, list())) {
    stop(paste(
      "No records available for this object (use the simulate",
      "method and the parameter recording = TRUE to generate them)."
    ))
  }
  return(easeObj@records)
}
