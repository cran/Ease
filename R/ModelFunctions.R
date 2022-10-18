# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                               MODEL FUNCTIONS                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Genome, Geno- and Haplotypes ----

### Identifiers ----

#' Genome identifier
#'
#' Generation of the input genome ID, i.e. the concatenation in string form
#' of the names of the loci and alleles constituting this genome.
#'
#' @param listLoci the list of all loci
#' @param alleles the vector of all the alleles
#'
#' @return The genome ID as a character string.
#'
#' @author Ehouarn Le Faou
#'
IDgenomeGeneration <- function(listLoci, alleles) {
  IDlist <- unlist(lapply(listLoci, length))
  IDvect <- paste(rep(names(IDlist), IDlist), alleles, sep = ":")
  return(Reduce(function(x, y) {
    paste(x, y, sep = "/")
  }, IDvect))
}

#' Genotype identifier
#'
#' Generation of the input genotype ID, i.e. the concatenation in string form
#' of the names of the alleles constituting these haplotypes (the two from
#' the diploid genome and the one from the haploid genome).
#'
#' @param dl1 the first diploid haplotype as a character (or factors) vector.
#' @param dl2 the second diploid haplotype as a character (or factors) vector.
#' @param hl the haploid haplotype as a character (or factors) vector.
#'
#' @return The genotype ID as a character string.
#'
#' @author Ehouarn Le Faou
#'
IDgenotypeGeneration <- function(dl1, dl2, hl = NULL) {
  if (is.null(hl)) {
    haploStr <- unlist(lapply(list(dl1, dl2), function(x) {
      as.character(Reduce(paste0, x))
    }))
    return(as.character(paste0(haploStr[[1]], "/", haploStr[[2]])))
  }
  haploStr <- unlist(lapply(list(dl1, dl2, hl), function(x) {
    as.character(Reduce(paste0, x))
  }))
  return(as.character(paste0(
    haploStr[[1]], "/",
    haploStr[[2]], "||",
    haploStr[[3]]
  )))
}

#' Haplotype identifier
#'
#' Generation of the input haplotype ID, i.e. the concatenation in string form
#' of the names of the alleles constituting this haplotype.
#'
#' @param dl the diploid haplotype as a character (or factors) vector.
#' @param hl the haploid haplotype as a character (or factors) vector.
#'
#' @return The haplotype ID as a character string.
#'
#' @author Ehouarn Le Faou
#'
IDhaplotypeGeneration <- function(dl, hl) {
  haploStr <- unlist(lapply(list(dl, hl), function(x) {
    as.character(Reduce(paste0, x))
  }))
  return(as.character(paste0(haploStr[[1]], "||", haploStr[[2]])))
}

### Generation ----

#' Haplotyping
#'
#' Generation of haplotypes associated with a \code{Genome} object.
#'
#' The generated haplotypes are output as a list of three enumeration
#' in the form of matrices of alleles (each row corresponding to an haplotype,
#' each column to a locus). The first enumeration corresponds to haplotypes
#' considering only haploid loci, the second only diploid loci and the third
#' all loci (with two matrices, 1 for haploid loci, 1 for diploid loci).
#'
#' @param genomeObj a \code{Genome} object
#'
#' @return A list of matrices describing haplotypes in rows.
#'
#' @author Ehouarn Le Faou
#'
haplotyping <- function(genomeObj) {
  haplotypesHL <- expand.grid(
    genomeObj@listLoci[1:genomeObj@nbHL]
  )
  haplotypesDL <- expand.grid(
    genomeObj@listLoci[(genomeObj@nbHL + 1):(genomeObj@nbHL + genomeObj@nbDL)]
  )
  nbHaploHL <- nrow(haplotypesHL)
  nbHaploDL <- nrow(haplotypesDL)

  haplotypes <- list()
  haplotypes$DL <- haplotypesDL[rep(seq_len(nbHaploDL), nbHaploHL),
    seq_len(genomeObj@nbDL),
    drop = FALSE
  ]
  haplotypes$HL <- haplotypesHL[rep(seq_len(nbHaploHL), each = nbHaploDL),
    seq_len(genomeObj@nbHL),
    drop = FALSE
  ]

  if (is.null(names(genomeObj@listHapLoci))) {
    colnames(haplotypes$HL) <- paste0("HL", seq_len(genomeObj@nbHL))
  } else {
    colnames(haplotypes$HL) <- names(genomeObj@listHapLoci)
  }
  if (is.null(names(genomeObj@listDipLoci))) {
    colnames(haplotypes$DL) <- paste0("DL", seq_len(genomeObj@nbDL))
  } else {
    colnames(haplotypes$DL) <- names(genomeObj@listDipLoci)
  }
  rownames(haplotypes$HL) <- seq_len(nrow(haplotypes$HL))
  rownames(haplotypes$DL) <- seq_len(nrow(haplotypes$DL))

  return(list(
    haplotypesHL = haplotypesHL, haplotypesDL = haplotypesDL,
    haplotypes = haplotypes
  ))
}

#' Genotyping
#'
#' Generation of genotypes associated with a \code{Genome} object.
#'
#' The output genotypes are described as a list of three matrices. A
#' genotype consists of two diploid haplotypes (first two matrices)
#' and one haploid haplotype (third matrix), which are read at the
#' same row number on all three matrices.
#'
#' @param genomeObj a \code{Genome} object.
#'
#' @return A list of matrices describing genotypes in rows.
#'
#' @author Ehouarn Le Faou
#'
genotyping <- function(genomeObj) {

  # Only diploid loci genotypes
  genotypesDL <- list(
    DL1 = matrix(factor(), 0, genomeObj@nbDL),
    DL2 = matrix(factor(), 0, genomeObj@nbDL)
  )

  colnames(genotypesDL$DL1) <- names(genomeObj@haplotypes$DL)
  colnames(genotypesDL$DL2) <- names(genomeObj@haplotypes$DL)

  k <- 1
  for (i in 1:nrow(genomeObj@haplotypesDL)) {
    for (j in k:nrow(genomeObj@haplotypesDL)) {
      genotypesDL$DL1 <- rbind(
        genotypesDL$DL1,
        as.matrix(genomeObj@haplotypes$DL[i, ])
      )
      genotypesDL$DL2 <- rbind(
        genotypesDL$DL2,
        as.matrix(genomeObj@haplotypes$DL[j, ])
      )
    }
    k <- k + 1
  }
  nbGenoDL <- nrow(genotypesDL$DL1)
  rownames(genotypesDL$DL1) <- 1:nbGenoDL
  rownames(genotypesDL$DL2) <- 1:nbGenoDL

  # Diploid+haploid loci genotypes
  genotypes <- list(
    DL1 = matrix(0, 0, genomeObj@nbDL),
    DL2 = matrix(0, 0, genomeObj@nbDL),
    HL = matrix(0, 0, genomeObj@nbHL)
  )
  for (i in 1:nrow(genomeObj@haplotypesHL)) {
    genotypes$DL1 <- rbind(
      genotypes$DL1,
      genotypesDL$DL1
    )
    genotypes$DL2 <- rbind(
      genotypes$DL2,
      genotypesDL$DL2
    )
    genotypes$HL <- rbind(
      genotypes$HL,
      matrix(rep(as.matrix(genomeObj@haplotypesHL[i, ]), nbGenoDL),
        ncol = genomeObj@nbHL,
        byrow = TRUE
      )
    )
  }
  genomeObj@nbGeno <- nrow(genotypes$DL1)
  colnames(genotypes$DL1) <- names(genomeObj@haplotypes$DL)
  colnames(genotypes$DL2) <- names(genomeObj@haplotypes$DL)
  colnames(genotypes$HL) <- names(genomeObj@haplotypes$HL)
  for (i in 1:3) {
    rownames(genotypes[[i]]) <- 1:genomeObj@nbGeno
  }
  genotypes$DL1 <- as.data.frame(unclass(genotypes$DL1),
    stringsAsFactors = TRUE
  )
  genotypes$DL2 <- as.data.frame(unclass(genotypes$DL2),
    stringsAsFactors = TRUE
  )
  genotypes$HL <- as.data.frame(unclass(genotypes$HL),
    stringsAsFactors = TRUE
  )
  return(genotypes)
}


## Mutation matrix ----

#' Mutation matrix from rates
#'
#' Generation of a mutation matrix from the allele enumeration vector of a loci
#' and the forward and backward mutation rates.
#'
#' See \code{MutationMatrix} for more details on mutation matrices.
#'
#' @param alleles allele enumeration vector of a locus
#' @param forwardMut forward mutation rate
#' @param backwardMut backward mutation rate
#'
#' @return An allelic mutation matrix (probability matrix which associates to
#' each allele in a row the probability of mutating or not to the other alleles
#' of the locus in question).
#'
#' @author Ehouarn Le Faou
#'
mutMatRates <- function(alleles, forwardMut, backwardMut) {
  le <- length(alleles)
  if (le == 1) {
    return(matrix(1, 1, 1))
  } else if (le == 2) {
    return(matrix(c(1 - forwardMut, forwardMut, backwardMut, 1 - backwardMut),
      2, 2,
      byrow = TRUE
    ))
  } else {
    return(matrix(c(
      1 - forwardMut, forwardMut, rep(0, le - 2),
      rep(c(
        backwardMut, 1 - forwardMut - backwardMut, forwardMut,
        rep(0, le - 2)
      ), le - 2),
      backwardMut, 1 - backwardMut
    ), le, le, byrow = TRUE))
  }
}


#' Individual mutation definition to allelic mutation matrices
#'
#' Translation of the list of individually defined mutations into allelic
#' mutation matrices which are then used to generate the genotypic mutation
#' matrix.
#'
#' @param genomeObj a \code{Genome} object
#' @param mutations list of mutations defined individually with the function
#' \link[Ease]{mutation}
#'
#' @return A list of the two list of allelic mutation matrices, for haploid
#' and diploid loci respectively.
#'
#' @author Ehouarn Le Faou
#'
mutMatFriendly <- function(genomeObj, mutations) {

  # Blank matrices definition
  mutHapLoci <- lapply(genomeObj@listHapLoci, function(alleles) {
    mutMatRates(alleles, 0, 0)
  })
  mutDipLoci <- lapply(genomeObj@listDipLoci, function(alleles) {
    mutMatRates(alleles, 0, 0)
  })

  # Naming
  for (i in 1:genomeObj@nbHL) {
    colnames(mutHapLoci[[i]]) <- genomeObj@listHapLoci[[i]]
    rownames(mutHapLoci[[i]]) <- genomeObj@listHapLoci[[i]]
  }
  for (i in 1:genomeObj@nbDL) {
    colnames(mutDipLoci[[i]]) <- genomeObj@listDipLoci[[i]]
    rownames(mutDipLoci[[i]]) <- genomeObj@listDipLoci[[i]]
  }

  # Translation of mutations in matrices
  for (mut in mutations) {
    fromHapLoci <- unlist(lapply(genomeObj@listHapLoci, function(x) mut[["from"]] %in% x))
    fromDipLoci <- unlist(lapply(genomeObj@listDipLoci, function(x) mut[["from"]] %in% x))

    toHapLoci <- unlist(lapply(genomeObj@listHapLoci, function(x) mut[["to"]] %in% x))
    toDipLoci <- unlist(lapply(genomeObj@listDipLoci, function(x) mut[["to"]] %in% x))

    comparison <- which(c(fromHapLoci, fromDipLoci) & c(toHapLoci, toDipLoci))

    if (length(comparison) == 0) {
      stop(paste(
        "In particular, check that the mutations defined are between",
        "alleles of the same loci, and that the names of the alleles",
        "are correct."
      ))
    } else if (length(comparison) > 1) {
      stop(paste("Ambiguous input."))
    }

    if (names(comparison) %in% names(genomeObj@listHapLoci)) {
      indexLoci <- which(fromHapLoci & toHapLoci)
      indexFrom <- which(genomeObj@listHapLoci[[indexLoci]] == mut$from)
      indexTo <- which(genomeObj@listHapLoci[[indexLoci]] == mut$to)
      mutHapLoci[[indexLoci]][indexFrom, indexTo] <- mut$rate
    } else {
      indexLoci <- which(fromDipLoci & toDipLoci)
      indexFrom <- which(genomeObj@listDipLoci[[indexLoci]] == mut$from)
      indexTo <- which(genomeObj@listDipLoci[[indexLoci]] == mut$to)
      mutDipLoci[[indexLoci]][indexFrom, indexTo] <- mut$rate
    }
  }
  # Generation of diagonal values
  for (i in seq_len(length(mutHapLoci))) {
    diag(mutHapLoci[[i]]) <- 0
    diag(mutHapLoci[[i]]) <- 1 - apply(mutHapLoci[[i]], 1, sum)
  }
  for (i in seq_len(length(mutDipLoci))) {
    diag(mutDipLoci[[i]]) <- 0
    diag(mutDipLoci[[i]]) <- 1 - apply(mutDipLoci[[i]], 1, sum)
  }
  return(list(mutHapLoci = mutHapLoci, mutDipLoci = mutDipLoci))
}


## Other matrices ----

#' Recombination matrix generation
#'
#' Generation of the recombination matrix associated to a \code{Genome} object.
#'
#' A recombination matrix is a square matrix of size equal to the number of
#' genotypes. It is a probability matrix in that the sum of the values in each
#' row is equal to 1. For a given genotype, the row associated with it
#' describes the probabilistic proportions that lead by recombination
#' between diploid loci to the production of the other genotypes (and of
#' itself if there are no mutations).
#'
#' @param genomeObj a \code{Genome} object
#'
#' @return A recombination matrix (probability matrix which associates to each
#' genotype in a row the probability of recombining or not and of becoming
#' another genotype or remaining the same).
#'
#' @author Ehouarn Le Faou
#'
recombinationMatrix <- function(genomeObj) {
  if (genomeObj@nbDL == 1) {
    recombinationMat <- matrix(0, genomeObj@nbGeno, genomeObj@nbGeno)
    diag(recombinationMat) <- rep(1, genomeObj@nbGeno)
  } else {
    recPatterns <- expand.grid(
      as.list(as.data.frame(matrix(c(FALSE, TRUE), 2, genomeObj@nbDL - 1)))
    )
    probRecPatterns <- apply(
      recPatterns, 1,
      function(ch) {
        prod(c(
          genomeObj@recRate[ch],
          (1 - genomeObj@recRate)[!ch]
        ))
      }
    )
    recombinationMat <- matrix(0, genomeObj@nbGeno, genomeObj@nbGeno)
    for (k in 1:genomeObj@nbGeno) {
      haplo1Ref <- unlist(genomeObj@genotypes$DL1[k, ])
      haplo2Ref <- unlist(genomeObj@genotypes$DL2[k, ])
      haploHLRef <- unlist(genomeObj@genotypes$HL[k, ])
      for (i in 1:nrow(recPatterns)) {
        haplo1 <- haplo1Ref
        haplo2 <- haplo2Ref
        rp <- unlist(recPatterns[i, ])
        for (j in 1:ncol(recPatterns)) {
          if (rp[j]) {
            haplo1prov <- c(haplo1[1:j], haplo2[(j + 1):genomeObj@nbDL])
            haplo2 <- c(haplo2[1:j], haplo1[(j + 1):genomeObj@nbDL])
            haplo1 <- haplo1prov
          }
        }
        index <- which(genomeObj@IDgenotypes ==
          IDgenotypeGeneration(haplo1, haplo2, haploHLRef) |
          genomeObj@IDgenotypes ==
            IDgenotypeGeneration(haplo2, haplo1, haploHLRef))
        recombinationMat[k, index] <- recombinationMat[k, index] +
          probRecPatterns[i]
      }
    }
  }
  colnames(recombinationMat) <- genomeObj@IDgenotypes
  rownames(recombinationMat) <- genomeObj@IDgenotypes
  return(recombinationMat)
}

#' Meiosis matrix generation
#'
#' Generation of the meiosis matrix associated to a \code{Genome} object.
#'
#' A meiosis matrix is a matrix where the number of rows is equal to the number
#' of genotypes and the number of columns to the number of haplotypes. It is a
#' matrix that allows to pass from parental genotypes to gametic haplotypes
#' by meiosis. It is a probability matrix in that the sum of the values in each
#' row is equal to 1. For a given genotype, the row associated with it
#' describes the probabilistic proportions that lead by meiosis to the
#' production of the other genotypes (and of itself if there are no mutations).
#'
#' @param genomeObj a \code{Genome} object
#'
#' @return A meiosis matrix (probability matrix that associates to each
#' genotype in a row the probability of producing each of the possible
#' haplotypes).
#'
#' @author Ehouarn Le Faou
#'
meiosisMatrix <- function(genomeObj) {
  meiosisMat <- matrix(0, genomeObj@nbGeno, genomeObj@nbHaplo)
  for (i in 1:genomeObj@nbGeno) {
    for (j in 1:genomeObj@nbHaplo) {
      ID1 <- IDgenotypeGeneration(
        genomeObj@genotypes$DL1[i, ],
        genomeObj@genotypes$HL[i, ]
      )
      ID2 <- IDgenotypeGeneration(
        genomeObj@genotypes$DL2[i, ],
        genomeObj@genotypes$HL[i, ]
      )
      ID <- IDgenotypeGeneration(
        genomeObj@haplotypes$DL[j, ],
        genomeObj@haplotypes$HL[j, ]
      )
      if (ID1 == ID) {
        meiosisMat[i, j] <- meiosisMat[i, j] + 0.5
      }
      if (ID2 == ID) {
        meiosisMat[i, j] <- meiosisMat[i, j] + 0.5
      }
    }
  }
  colnames(meiosisMat) <- genomeObj@IDhaplotypes
  rownames(meiosisMat) <- genomeObj@IDgenotypes
  return(meiosisMat)
}



#' Haplotype crossing matrix generation
#'
#' Generation of the haplotype crossing matrix associated to a \code{Genome}
#' object.
#'
#' A crossover matrix is a square matrix of size equal to the number of
#' haplotypes. It describes for each combination of two gametic haplotypes
#' the genotype index resulting from their syngamy. In the general case it
#' is not a symmetrical matrix (it is if a single haploid locus with a single
#' allele is defined), because the transmission of haploid loci is only
#' maternal, therefore non-symmetrical as is the transmission of diploid
#' loci. It is therefore necessary to enter the haplotype frequencies of
#' male gametes in the columns and the haplotype frequencies of female
#' gametes in the rows during the calculations (this is done in the
#' simulations).
#'
#' @param genomeObj a \code{Genome} object
#'
#' @return An haplotype crossing matrix.
#'
#' @author Ehouarn Le Faou
#'
haploCrossMatrix <- function(genomeObj) {
  haploCrossMat <- matrix(0, genomeObj@nbHaplo, genomeObj@nbHaplo)
  for (fem in 1:genomeObj@nbHaplo) {
    for (male in 1:genomeObj@nbHaplo) {
      ID1 <- IDgenotypeGeneration(
        genomeObj@haplotypes$DL[fem, ],
        genomeObj@haplotypes$DL[male, ],
        genomeObj@haplotypes$HL[fem, ]
      )
      ID2 <- IDgenotypeGeneration(
        genomeObj@haplotypes$DL[male, ],
        genomeObj@haplotypes$DL[fem, ],
        genomeObj@haplotypes$HL[fem, ]
      )
      haploCrossMat[fem, male] <- which(genomeObj@IDgenotypes == ID1 |
        genomeObj@IDgenotypes == ID2)
    }
  }
  colnames(haploCrossMat) <- paste0(genomeObj@IDhaplotypes, "(Mal)")
  rownames(haploCrossMat) <- paste0(genomeObj@IDhaplotypes, "(Fem)")
  return(haploCrossMat)
}



#' Generation of the matrix for calculating allelic frequencies
#'
#' Generates a matrix that allows to go from genotypic frequencies
#' to allelic frequencies.
#'
#' An allele frequency matrix is a matrix with rows equal to the number of
#' genotypes and columns equal to the number of alleles. By multiplying a
#' row matrix of genotype frequencies we obtain a row matrix of associated
#' allele frequencies.
#'
#' @param genomeObj a \code{Genome} object
#'
#' @return A matrix for calculating allelic frequencies from genotypes
#' frequencies.
#'
#' @author Ehouarn Le Faou
#'
alleleFreqMatGeneration <- function(genomeObj) {
  HLalleles <- as.character(unlist(genomeObj@listHapLoci))
  DLalleles <- as.character(unlist(genomeObj@listDipLoci))
  ploidy <- c(rep("HL", length(HLalleles)), rep("DL", length(DLalleles)))

  IDlocus <- c(
    rep(1:genomeObj@nbHL, unlist(lapply(genomeObj@listHapLoci, length))),
    rep(1:genomeObj@nbDL, unlist(lapply(genomeObj@listDipLoci, length)))
  )

  alleleFreqMat <- matrix(0, genomeObj@nbGeno, length(genomeObj@alleles))

  for (i in 1:genomeObj@nbGeno) {
    for (j in 1:length(genomeObj@alleles)) {
      if (ploidy[j] == "DL") {
        alleleFreqMat[i, j] <- (
          as.integer(genomeObj@genotypes$DL1[i, IDlocus[j]]
          == genomeObj@alleles[j]) +
            as.integer(genomeObj@genotypes$DL2[i, IDlocus[j]]
            == genomeObj@alleles[j])) / 2
      } else {
        alleleFreqMat[i, j] <-
          as.integer(genomeObj@genotypes$HL[i, IDlocus[j]]
          == genomeObj@alleles[j])
      }
    }
  }
  colnames(alleleFreqMat) <- genomeObj@alleles
  rownames(alleleFreqMat) <- genomeObj@IDgenotypes
  return(alleleFreqMat)
}



## Selection formulas handling ----

#' Extract the allele combination
#'
#' Conversion of an allelic combination defined in a selection formula into the
#' vector listing the alleles present (alleles that must be in the homozygous
#' state appear 2 times, 1 time for heterozygous).
#'
#' @param xVect allelic combination extracted from a selection formula.
#'
#' @return the list of alleles that must be present in the genotype to match
#' the input allelic combination
#'
#' @author Ehouarn Le Faou
#'
extractAlleleComb <- function(xVect) {
  yVect <- c()
  for (x in xVect) {
    y <- gsub("h", "", x)
    y <- gsub("\\(|\\)", "", y)
    if (identical(x, y)) {
      yVect <- c(yVect, x)
    } else {
      yVect <- c(yVect, rep(y, 2))
    }
  }
  return(yVect)
}


#' Treatment of a selection formula
#'
#' Conversion of the factors of a selection formula into the list of
#' corresponding allelic combinations
#'
#' @param factors formula factors (right-hand members)
#' @param genomeObj a \code{Genome} object for the test (see \code{checking}
#' parameter)
#' @param checking logical indicating whether a test verifying the
#' compatibility of input factors with the genome
#'
#' @return A list of vectors enumerating the allelic combinations that
#' correspond to the factors
#'
#' @author Ehouarn Le Faou
#'
selection.form.treatment <- function(factors, genomeObj = NULL, checking = FALSE) {
  # Treatment of the factors
  allAlleleCombsRaw <- colnames(factors)
  alleleCombRaw <- sapply(allAlleleCombsRaw, strsplit, split = ":")
  alleleCombList <- lapply(alleleCombRaw, extractAlleleComb)

  if (checking) {
    # Checking the validity of the expression
    test <- length(union(unlist(alleleCombList), genomeObj@alleles)) ==
      genomeObj@nbAlleles
    if (!test) {
      stop(paste(
        "At least one of the selection definition formulas is not",
        "correctly defined."
      ))
    }
  }
  return(alleleCombList)
}


#' Is this haplo/geno-type affected ?
#'
#' Determination for a given genotype or haplotype whether it
#' contains the allelic combination under selection
#'
#' @param refDNAtype the reference allelic combination (that of a genotype or a
#' haplotype
#' @param selDNAtype the selected allelic combination
#'
#' @return a logic indicating whether the reference genotype or haplotype
#' is affected by the allelic combination under selection
#'
#' @author Ehouarn Le Faou
#'
isAffected <- function(refDNAtype, selDNAtype) {
  test <- TRUE
  for (i in seq_len(length(selDNAtype))) {
    a <- selDNAtype[i]
    id <- which(names(refDNAtype) == names(a))
    if (!identical(id, integer(0))) {
      if (!(refDNAtype[id] == a)) {
        test <- FALSE
      }
    } else {
      test <- FALSE
    }
  }
  return(test)
}

#' Which alleles are homozygous in the input?
#'
#' Determine which alleles are at least once input as homozygous in the formula.
#'
#' @param formula a selection formula
#'
#' @return the enumeration of alleles that appear at least once homozygous
#'
#' @author Ehouarn Le Faou
#'
#' @import stats
whichHomoz <- function(formula) {
  fact <- attr(terms(formula), which = "factors")
  alleleComb <- selection.form.treatment(fact)
  countFact <- lapply(alleleComb, table)
  countFactVect <- lapply(countFact, function(x) {
    vect <- as.vector(x)
    names(vect) <- names(x)
    vect
  })
  names(countFactVect) <- NULL
  allCountFact <- unlist(countFactVect)
  homozAlleleOccur <- unique(names(allCountFact)[allCountFact == 2])
  return(homozAlleleOccur)
}

#' Are there any allelic combinations including homozygosity
#'
#' Test if there are homozygotes in the specified allelic combinations of
#' a selection formula
#'
#' @param formula a selection formula
#'
#' @return logical indicating if there are homozygotes
#'
#' @author Ehouarn Le Faou
#'
areThereHomoz <- function(formula) {
  return(!identical(whichHomoz(formula), character(0)))
}

#' Are there any allelic combinations including homozygosity
#'
#' Test if there are homozygotes in the specified allelic combinations of
#' a list of selection formulas
#'
#' @param selectFormula a list of selection formula
#'
#' @return logical indicating if there are homozygotes
#'
isHaploSelectFormula <- function(selectFormula) {
  return(!any(unlist(lapply(selectFormula, areThereHomoz))))
}



#' Conversion of selection formulas
#'
#' Conversion of a list of selection formulas into a genotypic (or haplotypic)
#' fitness vector associated with a \code{Genome} object.
#'
#' @param selectFormula a list of selection formulas
#' @param genomeObj a \code{Genome} object
#' @param haplo logical indicating whether the selection should apply to
#' haplotypes (in the case of gametic selection for example)
#'
#' @return a vector of fitness values
#'
#' @author Ehouarn Le Faou
#'
#' @import stats
selectFormIntoVect <- function(selectFormula, genomeObj, haplo = FALSE) {
  if (haplo) {
    if (!isHaploSelectFormula(selectFormula)) {
      stop(paste(
        "Selection formulas are not compatible with the definition",
        "of fitness for haplotypes (they include homozygous",
        "allelic combinations)"
      ))
    }
    nbDNAtype <- genomeObj@nbHaplo
    DNAtype <- genomeObj@haplotypes
  } else {
    nbDNAtype <- genomeObj@nbGeno
    DNAtype <- genomeObj@genotypes
  }

  haploidAlleles <- unlist(genomeObj@listHapLoci)
  homozAllelesForm <- unlist(lapply(selectFormula, whichHomoz))

  if (any(haploidAlleles %in% homozAllelesForm)) {
    stop(paste(
      "The alleles of the haploid loci cannot be homozygous in",
      "the selection formulas."
    ))
  }


  # Converting the list of selection formulas into the selection list
  listSelect <- list()
  for (formula in selectFormula) {
    factors <- attr(terms(formula), which = "factors")
    listSelect <- c(listSelect, list(list(
      alleleComb = selection.form.treatment(factors, genomeObj, checking = TRUE),
      fitnessEffect = eval(parse(text = rownames(factors)[1]))
    )))
  }
  names(listSelect) <- NULL


  fitTypes <- c()
  for (i in seq_len(nbDNAtype)) {
    DNAtypeAlleles <- as.character(
      unlist(lapply(DNAtype, function(x) x[i, ]))
    )
    fitType <- 1
    for (select in listSelect) {
      fit <- select$fitnessEffect
      for (sel in select$alleleComb) {
        refDNAtype <- sapply(
          unique(DNAtypeAlleles),
          function(x) length(which(DNAtypeAlleles == x))
        )
        selDNAtype <- sapply(
          unique(sel),
          function(x) length(which(sel == x))
        )

        if (isAffected(refDNAtype, selDNAtype)) {
          fitType <- fitType * fit
        }
      }
    }
    fitTypes <- c(fitTypes, fitType)
  }

  return(fitTypes)
}

#' Treatment of selection formulas
#'
#' Determines whether an entry for the selection is a list of selection
#' formulas or just a vector. If it is a list of formulas, turns them
#' into a vector. If it is a vector, does nothing.
#'
#' @param selectInput a selection input
#' @param genomeObj a \code{Genome} object
#' @param haplo logical indicating whether the selection should apply to
#' haplotypes (in the case of gametic selection for example)
#'
#' @return a vector of fitness values
#'
#' @author Ehouarn Le Faou
#'
selectInputTreatment <- function(selectInput, genomeObj, haplo = FALSE) {
  if (inherits(selectInput, "list")) {
    if (all(sapply(selectInput, inherits, what = "formula"))) {
      return(selectFormIntoVect(selectInput, genomeObj, haplo))
    } else {
      stop(paste(
        "The only possible input for the selection is a vector of",
        "fitness values or a list of selection formulas. "
      ))
    }
  }
  return(selectInput)
}

## Results ----


#' Custom output function
#'
#' Allow to produce a custom output for a simulation.
#'
#' This function is called each generation in each population of a simulation
#' and systematically returns a list with the first element being a logic that
#' indicates whether something should be saved. If so, the second element
#' of this list will be saved.
#'
#' By default the save nothing function, but it can be changed by the user as
#' an argument in the \code{simulate} method of the \code{Metapopulation}
#' class.
#'
#' @param pop list of some characteristics of the population :
#' - customOutput : list of all previous savings
#' - gen : generation
#' - freqGeno : list of genotypic frequency matrices (matrix 1 x # genotypes).
#' The list is constructed as follows: if the population is hermaphroditic it
#' has only one element "ind", if the population is dioecious it has three
#' elements, "female", "male" and "ind" which correspond respectively to
#' the genotypic frequencies of the females, the males and the average of
#' the two (assuming a sex ratio of 50:50).
#' - freqHaplo : list of genotypic frequency matrices (matrix 1 x # haplotypes).
#' The list is constructed in the same way as for genotypic frequencies (see
#' above).
#' - freqAlleles : list of allelic frequency matrices (matrix 1 x # alleles).
#' The list is constructed in the same way as for genotypic frequencies (see
#' above).
#'
#' @author Ehouarn Le Faou
#'
#' @export
outFunct <- function(pop) {
  return(list(FALSE))
}

#' Processing a result (or record) list
#'
#' @param x list of result or record
#'
#' @return Merges the column names of the matrices making up the list with
#' the names of the matrices, then merges the matrices together.
#'
#' @author Ehouarn Le Faou
#'
rowResultGen <- function(x) {
  namesRes <- names(x)
  x <- sapply(1:length(x), function(i) {
    colnames(x[[i]]) <- paste0(namesRes[i], colnames(x[[i]]))
    x[[i]]
  })
  return(Reduce(cbind, x))
}
