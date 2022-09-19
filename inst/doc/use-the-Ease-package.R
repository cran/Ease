## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(Ease)

## -----------------------------------------------------------------------------
LD = list(dl = as.factor(c("A", "a")))
HL = list(hl = as.factor(c("B", "b")))
genomeObj = setGenome(listHapLoci = HL, listDipLoci = LD)

## -----------------------------------------------------------------------------
genomeObj

## -----------------------------------------------------------------------------
print(genomeObj)

## -----------------------------------------------------------------------------
mutMatrixObj = setMutationMatrix(genomeObj = genomeObj,
                                 mutHapLoci = list(matrix(c(0.95, 0.05, 0.03, 0.97), 2, byrow = T)),
                                 mutDipLoci = list(matrix(c(0.9, 0.1, 0.09, 0.91), 2, byrow = T)))
mutMatrixObj

## -----------------------------------------------------------------------------
mutMatrixObj = setMutationMatrixByRates(genomeObj = genomeObj, forwardMut = 1e-2)
mutMatrixObj

## -----------------------------------------------------------------------------
selectionObj = setSelectNeutral(genomeObj = genomeObj)

## -----------------------------------------------------------------------------
selectionObj

## -----------------------------------------------------------------------------
print(selectionObj)

## -----------------------------------------------------------------------------
s = 0.8
h = 0.5
selectionObj = setSelectOnInds(genomeObj = genomeObj, indFit = c(1, 1, 1, 1, 1 - h*s, 1 - s))

## -----------------------------------------------------------------------------
selectionObj

## -----------------------------------------------------------------------------
print(selectionObj)

## -----------------------------------------------------------------------------
selectionObj = setSelectOnGametesProd(genomeObj = genomeObj, indProdFit = c(1, 1, 1, 1, 1 - h*s, 1 - s))

## -----------------------------------------------------------------------------
selectionObj = setSelectOnGametes(genomeObj = genomeObj, femaleFit = c(1, 1, 1 - s, 1 - s))

## -----------------------------------------------------------------------------
s = 0.8
h = 0.5
selectionObj = setSelectOnInds(genomeObj = genomeObj, 
                               indFit = c(1, 1, 1, 1, 1 - h*s, 1 - s))
selectionObj = setSelectOnGametesProd(indProdFit = c(1, 1, 1, 1, 1 - h*s, 1 - s),
                                      selectionObj = selectionObj)
selectionObj = setSelectOnGametes(femaleFit = c(1, 1, 1 - s, 1 - s),
                                  selectionObj = selectionObj)
print(selectionObj)

## ---- echo = FALSE------------------------------------------------------------
s = 0.8
h = 0.5
selectionObj = setSelectOnInds(genomeObj = genomeObj, indFit = c(1, 1, 1, 1, 1 - h*s, 1 - s))

## -----------------------------------------------------------------------------
mod = setEase(N = 100, threshold = 1e6, dioecy = F, selfRate = 0.5,
            stopCondition = list(nucleo = "a", cyto = "b"),
            mutMatrixObj = mutMatrixObj,
            genomeObj = genomeObj,
            selectionObj = selectionObj)

## -----------------------------------------------------------------------------
mod = simulate(mod, nsim = 50, recording = T, seed = 123)

## ---- fig.width = 7, fig.asp = 0.6--------------------------------------------
plot(mod)

