## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(Ease)

## -----------------------------------------------------------------------------
DL <- list(dl = c("A", "a"))
HL <- list(hl = c("B", "b"))
genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)

## -----------------------------------------------------------------------------
genomeObj

## -----------------------------------------------------------------------------
print(genomeObj)

## -----------------------------------------------------------------------------
mutMatrixObj <- setMutationMatrix(
  genomeObj = genomeObj,
  mutHapLoci = list(matrix(c(
    0.95, 0.05,
    0.03, 0.97
  ), 2, byrow = TRUE)),
  mutDipLoci = list(matrix(c(
    0.9, 0.1,
    0.09, 0.91
  ), 2, byrow = TRUE))
)
mutMatrixObj

## -----------------------------------------------------------------------------
mutMatrixObj <- setMutationMatrix(genomeObj = genomeObj, forwardMut = 0.01)
mutMatrixObj

## -----------------------------------------------------------------------------
mutMatrixObj <- setMutationMatrix(
  genomeObj = genomeObj,
  mutations = list(
    mutation(from = "A", to = "a", rate = 0.01),
    mutation(from = "B", to = "b", rate = 0.01)
  )
)
mutMatrixObj

## -----------------------------------------------------------------------------
selectionObj <- setSelectNeutral(genomeObj = genomeObj)

## -----------------------------------------------------------------------------
selectionObj

## -----------------------------------------------------------------------------
print(selectionObj)

## -----------------------------------------------------------------------------
indFit <- list(
  1 - h * s ~ a:b,
  1 - s ~ h(a):b
)
s <- 0.8
h <- 0.5
selectionObj <- setSelectOnInds(
  genomeObj = genomeObj,
  indFit = indFit
)

## -----------------------------------------------------------------------------
selectionObj

## -----------------------------------------------------------------------------
print(selectionObj)

## -----------------------------------------------------------------------------
indProdFit <- list(
  1 - h * s ~ a:b,
  1 - s ~ h(a):b
)
s <- 0.8
h <- 0.5
selectionObj <- setSelectOnGametesProd(
  genomeObj = genomeObj,
  indProdFit = indProdFit
)

## -----------------------------------------------------------------------------
gamFit <- list(
  1 - s ~ a:b
)
s <- 0.8
h <- 0.5
selectionObj <- setSelectOnGametes(
  genomeObj = genomeObj,
  gamFit = gamFit
)

## -----------------------------------------------------------------------------
indFit <- list(
  1 - h * s ~ a:b,
  1 - s ~ h(a):b
)
indProdFit <- list(
  1 - h * s ~ a:b,
  1 - s ~ h(a):b
)
gamFit <- list(
  1 - s ~ a:b
)
s <- 0.8
h <- 0.5
selectionObj <- setSelectOnInds(genomeObj = genomeObj, indFit = indFit)
selectionObj <- setSelectOnGametesProd(indProdFit = indProdFit, selectionObj = selectionObj)
selectionObj <- setSelectOnGametes(femaleFit = gamFit, selectionObj = selectionObj)
print(selectionObj)

## -----------------------------------------------------------------------------
DL <- list(dl = c("A", "a"))
HL <- list(hl = c("B", "b"))
mutations <- list(
  mutation(from = "A", to = "a", rate = 1e-3),
  mutation(from = "B", to = "b", rate = 1e-3)
)
genomeObj <- setGenome(listHapLoci = HL, listDipLoci = DL)
pop <- setPopulation(
  name = "A",
  size = 1000,
  dioecy = TRUE,
  genomeObj = genomeObj,
  selectionObj = setSelectNeutral(genomeObj),
  mutMatrixObj = setMutationMatrix(genomeObj, mutations = mutations)
)
pop

## -----------------------------------------------------------------------------
print(pop)

## -----------------------------------------------------------------------------
customOutFunct <- function(pop) {
  if (pop$freqAlleles[4] < 0.1 | pop$freqAlleles[4] > 0.9) {
    return(list(TRUE, a = list(gen = pop$gen, freq = pop$freqAlleles[4])))
  }
  return(list(FALSE))
}

## -----------------------------------------------------------------------------
DL <- list(dl = c("A", "a"))
genomeObj <- setGenome(listDipLoci = DL)
print(genomeObj)

## -----------------------------------------------------------------------------
mutMatrixObj <- setMutationMatrix(
  genomeObj = genomeObj,
  forwardMut = 1e-3,
  backwardMut = 1e-3
)
print(mutMatrixObj)

## -----------------------------------------------------------------------------
# Selection in population 1
indFit <- list(
  1 + h * s ~ a,
  1 + s ~ h(a)
)
h <- 0.5
s <- 0.05
selectionObj1 <- setSelectOnInds(
  genomeObj = genomeObj,
  indFit = indFit
)
print(selectionObj1)

# Selection in population 2
indFit <- list(
  1 + h * s ~ A,
  1 + s ~ h(A)
)
h <- 0.5
s <- 0.05
selectionObj2 <- setSelectOnInds(
  genomeObj = genomeObj,
  indFit = indFit
)
print(selectionObj2)

## -----------------------------------------------------------------------------
migMat <- matrix(c(
  0.995, 0.005,
  0.005, 0.995
), 2, 2)

## -----------------------------------------------------------------------------
metapop <- setMetapopulation(
  populations = list(
    setPopulation(
      name = "pop1", size = 500, dioecy = F, genomeObj = genomeObj,
      selectionObj = selectionObj1, mutMatrixObj = mutMatrixObj, selfRate = 0.5
    ),
    setPopulation(
      name = "pop2", size = 500, dioecy = F, genomeObj = genomeObj,
      selectionObj = selectionObj2, mutMatrixObj = mutMatrixObj, selfRate = 0.5
    )
  ),
  migMat = migMat
)
metapop

## -----------------------------------------------------------------------------
print(metapop)

## -----------------------------------------------------------------------------
nsim <- 100
metapop <- simulate(metapop,
  nsim = nsim, seed = 123,
  recording = TRUE, recordGenGap = 50,
  threshold = 800
)
metapopDeterminist <- simulate(metapop,
  drift = FALSE, seed = 123,
  recording = TRUE, threshold = 800
)

## -----------------------------------------------------------------------------
rec <- getRecords(metapop)

recMean <- Reduce(
  function(x, y) {
    x$pop1 <- x$pop1 + y$pop1
    x$pop2 <- x$pop2 + y$pop2
    x
  },
  rec
)
recMean$pop1 <- recMean$pop1 / nsim
recMean$pop2 <- recMean$pop2 / nsim

recDeterminist <- getRecords(metapopDeterminist)$s1

## ---- dpi = 100, fig.width = 7, fig.asp = 2-----------------------------------
par(mfrow = c(2, 1))

### Allelic frequency
plot(c(),
  xlim = c(0, 800), ylim = c(0, 1), col = "blue",
  xlab = "Generation\n(blue for A, red for B)",
  ylab = "Frequency of a"
)

# Stochastic
for (i in seq_len(nsim)) {
  lines(rec[[i]]$pop1$gen, rec[[i]]$pop1$a, lty = "longdash", col = "lightblue")
  lines(rec[[i]]$pop2$gen, rec[[i]]$pop2$a, lty = "longdash", col = "lightpink")
}
lines(recMean$pop1$gen, recMean$pop1$a, lty = "longdash", col = "blue")
lines(recMean$pop2$gen, recMean$pop2$a, lty = "longdash", col = "red")

# Deterministic
lines(recDeterminist$pop1$gen, recDeterminist$pop1$a, lty = "longdash", col = "blue", lwd = 2)
lines(recDeterminist$pop2$gen, recDeterminist$pop2$a, lty = "longdash", col = "red", lwd = 2)

# Mean stochastic
lines(recMean$pop1$gen, recMean$pop1$a, col = "blue")
lines(recMean$pop2$gen, recMean$pop2$a, col = "red")


### Mean fitness
plot(c(),
  xlim = c(0, 800), ylim = c(1, 1 + s), col = "blue",
  xlab = "Generation\n(blue for A, red for B)",
  ylab = "Frequency of a"
)

# Stochastic
for (i in seq_len(nsim)) {
  lines(rec[[i]]$pop1$gen, rec[[i]]$pop1$indMeanFit, lty = "longdash", col = "lightblue")
  lines(rec[[i]]$pop2$gen, rec[[i]]$pop2$indMeanFit, lty = "longdash", col = "lightpink")
}
lines(recMean$pop1$gen, recMean$pop1$indMeanFit, lty = "longdash", col = "blue")
lines(recMean$pop2$gen, recMean$pop2$indMeanFit, lty = "longdash", col = "red")

# Deterministic
lines(recDeterminist$pop1$gen, recDeterminist$pop1$indMeanFit, lty = "longdash", col = "blue", lwd = 2)
lines(recDeterminist$pop2$gen, recDeterminist$pop2$indMeanFit, lty = "longdash", col = "red", lwd = 2)

# Mean stochastic
lines(recMean$pop1$gen, recMean$pop1$indMeanFit, col = "blue")
lines(recMean$pop2$gen, recMean$pop2$indMeanFit, col = "red")

