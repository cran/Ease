# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                              Metapopulation CLASS                            #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#' The validity check for the \code{Metapopulation} class
#'
#' @param object a \code{Metapopulation} object
#'
#' @return a boolean corresponding to whether the object is a correct
#' \code{Metapopulation} object.
#'
#' @author Ehouarn Le Faou
#'
check.metapopulation <- function(object) {
  is.correct.transition.matrix(object@migMat, "migration", "migMat")
  if (nrow(object@migMat) != object@nbPop) {
    stop(paste(
      "The size of the migration matrix does not correspond to the",
      "population number given as input."
    ))
  }
  if (length(unique(object@names)) != object@nbPop) {
    stop(paste(
      "The names of the populations in the metapopulation should",
      "all be different."
    ))
  }

  if (!all(object@dioecies == object@dioecies[[1]])) {
    stop(paste(
      "All populations should have the same sexual system",
      "(dioecy parameters)."
    ))
  }

  if (length(unique(unlist(object@genomeIDs))) != 1) {
    stop(paste(
      "All populations in the metapopulation must have been created",
      "with the same genome."
    ))
  }

  return(TRUE)
}


#' Metapopulation
#'
#' The class \code{Metapopulation} is used to centralise the information
#' relating to the populations that we want to simulate, as well as to define
#' the migration conditions between them if there are several. This class
#' is thus defined by a list of objects \code{Population} and a migration
#' matrix.
#'
#' @slot populations list of objects \code{Population}
#' @slot nbPop number of populations
#' @slot names names of the populations
#' @slot migMat migration matrix between population (if there is more than one)
#' @slot sizes sizes of the populations
#' @slot dioecies sexual systems of the populations (they must all be the same)
#' @slot selfRates selfing rates of the populations
#' @slot demographies demography parameter of the populations
#' @slot growthRates growth rates of the populations
#' @slot initPopSizes initial population sizes of the populations
#' @slot initGenoFreqs initial genotypic frequencies of the populations
#' @slot genome a \code{Genome} object
#' @slot genomeIDs ID of the \code{Genome} object
#' @slot mutMat a \code{MutationMatrix} object
#' @slot selection a \code{Selection} object
#' @slot recMat recombination matrix
#' @slot meiosisMat a meiosis matrix
#' @slot haploCrossMat an haplotype crossing matrix
#' @slot haploCrossMatNamed an haplotype crossing matrix with names of
#' genotypes instead of their indices
#' @slot gametogenesisMat a gametogenesis matrix
#' @slot alleleFreqMat a matrix for calculating allelic frequencies
#' @slot rawOutputSimul raw output of the simulation function, its refinement
#' is done directly afterwards in the \code{simulate} method
#' @slot stopCondition list of stop conditions for the simulation (if required)
#' @slot IDstopCondition names of stop conditions. They are given an arbitrary
#' name if none is given by the user.
#' @slot results data.frame.
#' @slot records list.
#' @slot customOutput list.
#'
#' @author Ehouarn Le Faou
#'
#' @export
setClass("Metapopulation",
  representation(
    populations = "list",
    nbPop = "numeric",
    names = "list",
    migMat = "matrix",
    sizes = "list",
    dioecies = "list",
    selfRates = "list",
    demographies = "list",
    growthRates = "list",
    initPopSizes = "list",
    initGenoFreqs = "list",
    genome = "Genome",
    genomeIDs = "list",
    mutMat = "MutationMatrix",
    selection = "list",
    recMat = "matrix",
    meiosisMat = "matrix",
    haploCrossMat = "matrix",
    haploCrossMatNamed = "data.frame",
    gametogenesisMat = "matrix",
    alleleFreqMat = "matrix",
    rawOutputSimul = "list",
    stopCondition = "list",
    IDstopCondition = "character",
    results = "data.frame",
    records = "list",
    customOutput = "list"
  ),
  validity = check.metapopulation
)


## Initialize method ----

#' Initialize method for the \code{Metapopulation} class
#'
#' @param .Object a \code{Metapopulation} object
#' @param populations list of \code{Population} object(s)
#' @param migMat migration matrix
#'
#' @return a \code{Metapopulation} object
#'
#' @author Ehouarn Le Faou
#'
setMethod("initialize", "Metapopulation", function(.Object, populations, migMat) {

  # Definition of attributes
  .Object@populations <- populations
  .Object@nbPop <- length(populations)

  .Object@names <- lapply(populations, function(x) x@name)
  .Object@sizes <- lapply(populations, function(x) x@size)
  .Object@dioecies <- lapply(populations, function(x) x@dioecy)
  .Object@selfRates <- lapply(populations, function(x) x@selfRate)
  .Object@demographies <- lapply(populations, function(x) x@demography)
  .Object@growthRates <- lapply(populations, function(x) x@growthRate)
  .Object@initPopSizes <- lapply(populations, function(x) x@initPopSize)
  .Object@genomeIDs <- lapply(populations, function(x) x@genome@IDgenome)
  .Object@genome <- populations[[1]]@genome
  .Object@initGenoFreqs <- lapply(populations, function(x) x@initGenoFreq)
  .Object@selection <- lapply(populations, function(x) x@selection)

  # Matrices
  .Object@mutMat <- populations[[1]]@mutMat
  .Object@recMat <- recombinationMatrix(.Object@genome)
  .Object@meiosisMat <- meiosisMatrix(.Object@genome)
  .Object@haploCrossMat <- haploCrossMatrix(.Object@genome)

  # Management of duplicate population names
  if (length(.Object@names) != length(unique(.Object@names))) {
    for (n in .Object@names) {
      wn <- which(.Object@names == n)
      if (length(wn) > 1) {
        .Object@names[wn] <- paste0(n, "(", 1:length(wn), ")")
      }
    }
  }

  haploCrossMatNamed <- .Object@haploCrossMat
  haploCrossMatNamed <-
    apply(.Object@haploCrossMat, 2, function(x) {
      .Object@genome@IDgenotypes[x]
    })
  haploCrossMatNamed <- as.data.frame(haploCrossMatNamed, 2, as.factor)
  rownames(haploCrossMatNamed) <- rownames(.Object@haploCrossMat)
  .Object@haploCrossMatNamed <- haploCrossMatNamed

  .Object@gametogenesisMat <- .Object@recMat %*% .Object@meiosisMat %*%
    .Object@mutMat@mutationMatrix
  .Object@alleleFreqMat <- alleleFreqMatGeneration(.Object@genome)


  if (is.default.matrix(migMat)) {
    .Object@migMat <- matrix(0, .Object@nbPop, .Object@nbPop)
    diag(.Object@migMat) <- 1
  } else {
    .Object@migMat <- migMat
  }

  # Validity of the object
  validObject(.Object)

  rownames(.Object@migMat) <- .Object@names
  colnames(.Object@migMat) <- .Object@names


  return(.Object)
})

## Simulate method ----

#' Simulate method for the \code{Metapopulation} class
#'
#' Performing simulations of an Metapopulation object. The returned object is the same
#' Metapopulation object completed with the results and records if they have been
#' activated.
#'
#' @param object a \code{Metapopulation} object
#' @param nsim the number of simulation to perform
#' @param seed the RNG seed to be fixed (allows exact reproduction of
#' results)
#' @param threshold maximum duration of a simulation (in generations)
#' @param includefreqGeno a logical indicating whether to include genotype
#' frequencies in the results
#' @param recording a logical indicating whether to record all mutations, i.e.
#' to record allelic and genotypic frequencies along the simulations
#' @param recordGenGap the number of generations between two records during
#' simulation, if the record parameter is TRUE. Whatever the value of this
#' parameter, both the first and the last generation will be included in
#' the record
#' @param drift a logical indicating whether genetic drift should be
#' considered (i.e. whether deterministic simulations are performed or not)
#' @param includeParams a logical indicating whether the parameters should be
#' included in the result data.frame (can be useful when compiling multiple
#' result tables)
#' @param includeFitness a logical indicating whether the mean fitness should
#' be included in the result data.frame (can be useful when compiling multiple
#' result tables)
#' @param verbose logical determining if the progress of the simulations should
#' be displayed or not (useful in case of many simulations)
#' @param stopCondition list of vectors that each describe the allele(s) that
#' must be fixed to define a stop condition. Each of these vectors
#' will therefore be associated with a stop condition
#' @param nameOutFunct name of the custom output function. This function is
#' called each generation in each population of a simulation and systematically
#' returns a list with the first element being a logic that indicates whether
#' something should be saved. If so, the second element of this
#' list will be saved.If the customOutFunct parameter is null (default),
#' there will be no custom output.
#'
#' @return An \code{Metapopulation} object from which we can now extract the results
#' (or the records if recording = TRUE) with the getResults and getRecords
#' functions.
#'
#' @author Ehouarn Le Faou
#'
#' @export
#'
setMethod("simulate", "Metapopulation", function(object, nsim = 1,
                                                 seed = NULL,
                                                 threshold = 500,
                                                 includefreqGeno = TRUE,
                                                 recording = FALSE,
                                                 recordGenGap = 1,
                                                 drift = TRUE,
                                                 includeParams = TRUE,
                                                 includeFitness = TRUE,
                                                 verbose = FALSE,
                                                 stopCondition = list(),
                                                 nameOutFunct = "outFunct") {
  # Warnings
  if (!drift & nsim > 1) {
    warning(paste(
      "It is useless to repeat simulations without the drift",
      "(they will all be identical). The number of simulations",
      "has therefore been fixed at 1."
    ))
    nsim <- 1
  }

  # Stop condition -> THE ALLELE NAMES MUST ALL BE DIFFERENT
  if (length(union(object@genome@alleles, unlist(stopCondition))) !=
    length(object@genome@alleles)) {
    stop(paste(
      "The allele names given in the list of stop conditions do not",
      "correspond to the alleles in the metapopulation genome."
    ))
  }

  if (!identical(stopCondition, list())) {
    IDstopCondition <- names(stopCondition)
    if (identical(IDstopCondition, NULL)) {
      IDstopCondition <- paste0("stop", 1:length(stopCondition))
    }
    if ("" %in% IDstopCondition) {
      w <- which(IDstopCondition == "")
      IDstopCondition[w] <- paste0("stop", 1:length(w))
    }
    if (length(union(IDstopCondition, object@genome@alleles)) !=
      length(IDstopCondition) + length(object@genome@alleles)) {
      stop(paste(
        "The names of the stop conditions cannot be similar to the",
        "names of the alleles (to avoid redundancy between the names",
        "of the columns of simulation results)"
      ))
    }
    object@IDstopCondition <- IDstopCondition
    object@stopCondition <- rep(
      list(rep(NA, length(object@genome@alleles))),
      length(stopCondition)
    )
    for (i in 1:length(stopCondition)) {
      object@stopCondition[[i]][sapply(stopCondition[[i]], function(x) {
        which(object@genome@alleles == x)
      })] <- 1
    }
    names(object@stopCondition) <- IDstopCondition
  } else {
    object@stopCondition <- list()
  }


  if (verbose) {
    message("-=-=-=-=-=-=-=-=- METAPOPULATION SIMULATION =-=-=-=-=-=-=-=-=-")
    message("..............................................................")
  }
  set.seed(seed)
  object@rawOutputSimul <- METAPOP_SIMULATION(
    nbPop = object@nbPop,
    ids = object@names,
    migMat = object@migMat,
    nsim = nsim,
    verbose = verbose,
    recording = recording,
    recordGenGap = recordGenGap,
    drift = drift,
    nbHaplo = object@genome@nbHaplo,
    nbGeno = object@genome@nbGeno,
    idGeno = object@genome@IDgenotypes,
    nbAlleles = object@genome@nbAlleles,
    idAlleles = object@genome@alleles,
    nbLoci = (object@genome@nbDL + object@genome@nbHL),
    initGenoFreq = object@initGenoFreqs,
    meiosisMat = object@meiosisMat,
    gametogenesisMat = object@gametogenesisMat,
    popSize = object@sizes,
    threshold = threshold,
    dioecy = object@dioecies[[1]],
    selfRate = object@selfRates,
    stopCondition = object@stopCondition,
    IDstopCondition = object@IDstopCondition,
    haploCrossMat = object@haploCrossMat,
    alleleFreqMat = object@alleleFreqMat,
    gamFit = lapply(object@selection, function(x) x@gamFit),
    indFit = lapply(object@selection, function(x) x@indFit),
    gamProdFit = lapply(object@selection, function(x) x@gamProdFit),
    demography = object@demographies,
    growthRate = object@growthRates,
    initPopSize = object@initPopSizes,
    nameOutFunct = nameOutFunct
  )

  if (verbose) {
    message("''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''")
    message("Done.")
    message("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
    message()
  }

  # Formatting of results
  res <- lapply(object@rawOutputSimul, function(x) {
    lapply(x, function(y) {
      y$results
    })
  })
  resFormated <- list()
  for (s in res) {
    resFormated <- c(resFormated, lapply(s, function(x) {
      rowResultGen(x)
    }))
  }
  resFormated <- Reduce(rbind, resFormated)
  resFormated <- as.data.frame(resFormated)
  resFormated$simulation <- rep(1:nsim, each = object@nbPop)
  resFormated$population <- as.factor(unlist(rep(object@names, nsim)))
  object@results <- resFormated

  # Formatting of records
  if (recording) {
    rec <- lapply(object@rawOutputSimul, function(x) {
      lapply(x, function(y) {
        y$records
      })
    })
    recFormated <- list()
    for (s in rec) {
      recFormated <- c(recFormated, list(lapply(s, function(x) {
        rowResultGen(x)
      })))
    }
    names(recFormated) <- paste0("s", 1:nsim)
    recFormated <- lapply(recFormated, function(x) {
      lapply(x, as.data.frame)
    })
    object@records <- recFormated
  } else {
    object@records <- list()
  }

  # Formatting of custom output
  object@customOutput <- lapply(object@rawOutputSimul, function(x) {
    lapply(x, function(y) {
      y$custom
    })
  })

  return(object)
})


## Show method ----


#' Show method for the \code{Metapopulation} class
#'
#' @param object a \code{Metapopulation} object
#'
#' @author Ehouarn Le Faou
#'
setMethod("show", "Metapopulation", function(object) {
  catn("-=-=-=-= Metapopulation OBJECT =-=-=-=-")
  popSizes <- unlist(object@sizes)
  catn(paste0(object@nbPop, " populations of size ", listing(popSizes), ","))
  catn(paste0("forming a metapopulation of ", sum(popSizes), " individuals."))
  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  catn("(use print to get population details)")
})

## Print method ----

#' Print method for the \code{Metapopulation} class
#'
#' @param x a \code{Metapopulation} object
#' @param ... Ignored.
#'
#' @author Ehouarn Le Faou
#'
setMethod("print", "Metapopulation", function(x, ...) {
  catn("-=-=-=-= Metapopulation OBJECT =-=-=-=-")
  catn("             in details")
  catn()
  catn(" #  Populations: ")
  catn()
  print(x@populations[[1]], frame = F)
  for (pop in x@populations[-1]) {
    catn("~-~-~")
    print(pop, frame = F)
  }

  if (x@nbPop > 1) {
    catn()
    catn(" #  Migration matrix: ")
    catn()
    print(x@migMat)
  }

  catn()
  catn(" #  Haploptypes: ")
  catn()
  haplo <- x@genome@IDhaplotypes
  names(haplo) <- 1:length(haplo)
  print(haplo)

  catn()
  catn(" #  Genotypes: ")
  catn()
  geno <- x@genome@IDgenotypes
  names(geno) <- 1:length(geno)
  print(geno)

  catn()
  catn(" #  Matrices involved in gametogenesis: ")
  catn()
  catn("   - Mutation matrix: ")
  print(x@mutMat@mutationMatrix)
  catn()
  catn("   - Recombination matrix: ")
  print(x@recMat)
  catn()
  catn("   - Meiosis matrix: ")
  print(x@meiosisMat)
  catn()
  catn("   - Final gametogenesis matrix: ")
  print(x@gametogenesisMat)
  catn()

  catn(" #  Haplotypes crossing matrix: ")
  catn()
  print(x@haploCrossMatNamed)
  catn()

  catn(" #  Allele frequencies from genotype frequencies matrix: ")
  catn()
  print(x@alleleFreqMat)
  catn()


  catn("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
})
