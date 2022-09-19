
library(testthat)
library(Ease)

# NO LOCUS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that(
  desc = "Genome: no locus",
  code =
    {
      DL = list(dl = as.factor(c("A", "a")))
      expect_warning(expect_warning(setGenome()))
    }
)


# ONLY DIPLOID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


test_that(
  desc = "Genome: only diploid, 1 locus",
  code =
    {
      DL = list(dl = as.factor(c("A", "a")))
      expect_warning(setGenome(listDipLoci = DL))
    }
)

test_that(
  desc = "Genome: only diploid, 2 locus",
  code =
    {
      DL = list(dl1 = as.factor(c("A", "a")),
                dl2 = as.factor(c("B", "b")))
      expect_warning(setGenome(listDipLoci = DL,
                               recRate = 0.5))
    }
)

test_that(
  desc = "Genome: only diploid, 3 locus",
  code =
    {
      DL = list(dl1 = as.factor(c("A", "a")),
                dl2 = as.factor(c("B", "b")),
                dl3 = as.factor(c("C", "c")))
      expect_warning(setGenome(listDipLoci = DL,
                               recRate = c(0.5, 0.2)))
    }
)

# ONLY HAPLOID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that(
  desc = "Genome: only haploid, 1 locus",
  code =
    {
      HL = list(hl = as.factor(c("A", "a")))
      expect_warning(setGenome(listHapLoci = HL))
    }
)

test_that(
  desc = "Genome: only haploid, 2 locus",
  code =
    {
      HL = list(hl1 = as.factor(c("A", "a")),
                hl2 = as.factor(c("B", "b")))
      expect_warning(setGenome(listHapLoci = HL))
    }
)

test_that(
  desc = "Genome: only haploid, 3 locus",
  code =
    {
      HL = list(hl1 = as.factor(c("A", "a")),
                hl2 = as.factor(c("B", "b")),
                hl3 = as.factor(c("C", "c")))
      expect_warning(setGenome(listHapLoci = HL))
    }
)


# COMBINATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


test_that(
  desc = "Genome: combination, 5 loci",
  code =
    {
      DL = list(dl1 = as.factor(c("A", "a")),
                dl2 = as.factor(c("C", "c")))
      HL = list(hl1 = as.factor(c("B", "b")),
                hl2 = as.factor(c("D", "d")),
                hl3 = as.factor(c("E", "e")))
      expect_s4_class(object = setGenome(listHapLoci = HL, listDipLoci = DL,
                                         recRate = 0.5),
                      class = "Genome")
    }
)








