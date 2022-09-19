
#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace std;


//' Matrix product
//'
//' Standard matrix product. The first matrix must have as many columns as the
//' second has rows.
//'
//' @param MAT1 matrix
//' @param MAT2 matrix
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix MATRIX_PRODUCT(NumericMatrix MAT1, NumericMatrix MAT2)
{
	int nRow1(MAT1.nrow());
	int nCol1(MAT1.ncol());
	int nCol2(MAT2.ncol());

	NumericMatrix MAT_res(nRow1, nCol2);
	for (int i = 0; i < nRow1; i++) {
		for (int j = 0; j < nCol2; j++) {
			float sum(0);
			for (int k = 0; k < nCol1; k++) {
				sum += MAT1(i, k) * MAT2(k, j);
			}
			MAT_res(i, j) = sum;
		}
	}
	return MAT_res;
}

//' Matrix row bind
//'
//' Binding of two matrices by their row. Both matrices must have the same
//' number of rows.
//'
//' @param MAT1 matrix
//' @param MAT2 matrix
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix ROW_BIND(NumericMatrix MAT1, NumericMatrix MAT2)
{
	int nCol(MAT1.ncol());
	int nRow1(MAT1.nrow());
	int nRow2(MAT2.nrow());
	NumericMatrix MAT_res(nRow1 + nRow2, nCol);
	for (int i = 0; i < nCol; i++) {
		for (int j = 0; j < nRow1; j++) {
			MAT_res(j, i) = MAT1(j, i);
		}
		for (int j = 0; j < nRow2; j++) {
			MAT_res(j + nRow1, i) = MAT2(j, i);
		}
	}
	return MAT_res;

}


//_________________________________________________________________
//These functions are used to generate a multinomial random draw
//using the rmultinom function implemented in R


//' Importing the multinomial draw function from R into C++
//'
//' @param size number of random vectors to draw.
//' @param probs event probabilites
//' @param N number of trials
//'
//' @author Ehouarn Le Faou
//'
IntegerVector rmultinom_1(unsigned int& size, NumericVector& probs, unsigned int& N) {
	IntegerVector outcome(N);
	R::rmultinom(size, probs.begin(), N, outcome.begin());
	return outcome;
}

//' Multinomial draw
//'
//' @param n number of trials
//' @param size number of random vectors to draw.
//' @param probs event probabilites
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix rmultinom_rcpp(unsigned int& n, unsigned int& size, NumericVector& probs) {
	unsigned int N = probs.length();
	NumericMatrix sim(N, n);
	for (unsigned int i = 0; i < n; i++) {
		sim(_, i) = rmultinom_1(size, probs, N);
	}
	return sim;
}

//' Genetic drift
//'
//' Simulates genetic drift by drawing in a multinomial distribution with the
//' population size as the number of trials and the genotype frequencies as probabilities
//'
//' @param freq genotypic frequencies
//' @param N population size
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix DRIFT(NumericMatrix freq, int N)
{
	NumericMatrix freqP(1, freq.ncol());
	for (int j = 0; j < freq.ncol(); j++) { freqP(0, j) = freq(0, j); }

	NumericMatrix ind(freq.ncol(), 1);
	unsigned int n(1);
	unsigned int size(N);

	ind = rmultinom_rcpp(n, size, freqP);

	NumericMatrix freq2(1, freq.ncol());
	for (int i(0); i < freq.ncol(); i++) { freq2(0, i) = ind(i, 0) / N; }

	return freq2;
}

//' Should the simulation stop?
//'
//' Determination from the allele frequencies if the simulation should stop.
//'
//' @param freqAlleles allelic frequencies
//' @param stopCondition list of stop conditions
//'
//' @author Ehouarn Le Faou
//'
bool HAVE_TO_STOP(NumericMatrix freqAlleles, List stopCondition)
{
	if (stopCondition.size() == 0)
	{
		return false;
	}
	bool stopGlobal(false);
	for (int k = 0; k < stopCondition.size(); k++)
	{
		bool stop(true);
		NumericVector SC(stopCondition[k]);
		LogicalVector test = !is_na(SC);
		for (int i = 0; i < SC.length(); i++)
		{
			if (test(i)) { stop = stop && (freqAlleles(0, i) == SC(i)); }
		}
		stopGlobal = stopGlobal || stop;
	}
	return stopGlobal;
}

//' What stop conditions has the simulation reached?
//'
//' Determination of the stop condition index or indices that the simulation
//' has reached.
//'
//' @param freqAlleles allelic frequencies
//' @param stopCondition list of stop conditions
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix WHICH_STOP(NumericMatrix freqAlleles, List stopCondition)
{
	NumericMatrix idStop(1, stopCondition.size());
	for (int k = 0; k < stopCondition.size(); k++)
	{
		bool stop(true);
		NumericVector SC(stopCondition[k]);
		LogicalVector test = !is_na(SC);
		for (int i = 0; i < SC.length(); i++)
		{
			if (test(i)) { stop = stop && (freqAlleles(0, i) == SC(i)); }
		}
		idStop(0, k) = stop;
	}
	return idStop;
}

//' Standardisation of a matrix
//'
//' Divides each cell of a matrix by its sum, so that the sum of all these
//' cells becomes 1.
//'
//' @param x matrix
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix STANDARDISATION(NumericMatrix x)
{
	double sum(0);
	for (int i = 0; i < x.nrow(); i++) {
		for (int j = 0; j < x.ncol(); j++) {
			sum += x(i, j);
		}
	}
	NumericMatrix xStandard(x.nrow(), x.ncol());
	for (int i = 0; i < x.nrow(); i++) {
		for (int j = 0; j < x.ncol(); j++) {
			xStandard(i, j) = x(i, j) / sum;
		}
	}
	return xStandard;
}

//' Specify the seed for the RNG
//'
//' Sets the seed for random number generation to ensure repeatability of
//' simulations.
//'
//' @author Ehouarn Le Faou
//'
void SET_SEED(double seed) {
	Rcpp::Environment base_env("package:base");
	Rcpp::Function set_seed_r = base_env["set.seed"];
	set_seed_r(std::floor(std::fabs(seed)));
}

//' Crossing of male and female gametes
//'
//' Crossing of male and female gametes by knowing their respective genotype
//' frequencies and thanks to the crossing matrix which allows to determine
//' which genotypes result from each crossing between gametic haplotypes.
//'
//' @param nbHaplo number of haplotypes
//' @param nbGeno number of genotypes
//' @param freqHaploFemale female genotype frequencies
//' @param freqHaploMale male genotype frequencies
//' @param haploCrossMat haplotypes crossing matrix
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix CROSSING(int nbHaplo, int nbGeno, NumericMatrix freqHaploFemale, NumericMatrix freqHaploMale, NumericMatrix haploCrossMat)
{
	NumericMatrix freqHF(nbHaplo, 1);
	NumericMatrix freqHM(1, nbHaplo);
	for (int i = 0; i < nbHaplo; i++) {
		freqHF(i, 0) = freqHaploFemale(0, i);
		freqHM(0, i) = freqHaploMale(0, i);
	}
	NumericMatrix freqGenoCrossed(MATRIX_PRODUCT(freqHF, freqHM));

	NumericMatrix freqGeno(1, nbGeno);
	for (int i = 0; i < nbHaplo; i++) {
		for (int j = 0; j < nbHaplo; j++) {
			freqGeno(0, haploCrossMat(i, j) - 1) += freqGenoCrossed(i, j);
		}
	}
	return freqGeno;
}

//' Selfing
//'
//' Reproduction of the population by self-fertilization (only if the
//' population is hermaphroditic).
//'
//' @param nbHaplo number of haplotypes
//' @param nbGeno number of genotypes
//' @param freqGeno genotype frequencies
//' @param gametogenesisMat gametogenesis matrix
//' @param haploCrossMat haplotypes crossing matrix
//' @param femgamFit fitness of female gametes
//' @param malegamFit fitness of male gametes
//' @param femProdFit fitness for female gamete production
//' @param maleProdFit fitness for male gamete production
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix SELFING(int nbHaplo, int nbGeno, NumericMatrix freqGeno, NumericMatrix gametogenesisMat, NumericMatrix haploCrossMat, NumericVector femProdFit, NumericVector maleProdFit, NumericVector femgamFit, NumericVector malegamFit)
{
	NumericMatrix freqHaploMale(1, nbHaplo);
	NumericMatrix freqHaploFemale(1, nbHaplo);
	NumericMatrix freqGenoOffspringProv(1, nbGeno);
	NumericMatrix freqGenoOffspring(1, nbGeno);

	for (int i = 0; i < nbGeno; i++)
	{
		NumericMatrix freqGenoFemale(1, nbGeno);
		NumericMatrix freqGenoMale(1, nbGeno);
		freqGenoFemale(0, i) = femProdFit(i);
		freqGenoMale(0, i) = maleProdFit(i);
		freqHaploFemale = MATRIX_PRODUCT(freqGenoFemale, gametogenesisMat);
		freqHaploMale = MATRIX_PRODUCT(freqGenoMale, gametogenesisMat);
		for (int j = 0; j < nbHaplo; j++)
		{
			freqHaploFemale(0, j) = freqHaploFemale(0, j) * femgamFit(j);
			freqHaploMale(0, j) = freqHaploMale(0, j) * malegamFit(j);
		}
		freqGenoOffspringProv = CROSSING(nbHaplo, nbGeno, freqHaploFemale, freqHaploMale, haploCrossMat);

		for (int j = 0; j < nbGeno; j++)
		{
			freqGenoOffspring(0, j) = freqGenoOffspring(0, j) + freqGenoOffspringProv(0, j) * freqGeno(0, i);
		}
	}
	return STANDARDISATION(freqGenoOffspring);
}

//' Outcrossing
//'
//' Reproduction of the population by outcrossing.
//'
//' @param nbHaplo number of haplotypes
//' @param nbGeno number of genotypes
//' @param freqGenoFemale female genotype frequencies
//' @param freqGenoMale male genotype frequencies
//' @param gametogenesisMat gametogenesis matrix
//' @param haploCrossMat haplotypes crossing matrix
//' @param femgamFit fitness of female gametes
//' @param malegamFit fitness of male gametes
//' @param femProdFit fitness for female gamete production
//' @param maleProdFit fitness for male gamete production
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix OUTCROSSING(int nbHaplo, int nbGeno, NumericMatrix freqGenoFemale, NumericMatrix freqGenoMale, NumericMatrix gametogenesisMat, NumericMatrix haploCrossMat, NumericVector femProdFit, NumericVector maleProdFit, NumericVector femgamFit, NumericVector malegamFit)
{
	NumericMatrix freqGenoFemaleSelected(1, nbGeno);
	NumericMatrix freqGenoMaleSelected(1, nbGeno);
	NumericMatrix freqHaploMale(1, nbHaplo);
	NumericMatrix freqHaploFemale(1, nbHaplo);
	for (int i = 0; i < nbGeno; i++)
	{
		freqGenoFemaleSelected(0, i) = freqGenoFemale(0, i) * femProdFit(i);
		freqGenoMaleSelected(0, i) = freqGenoMale(0, i) * maleProdFit(i);
	}
	freqHaploFemale = MATRIX_PRODUCT(freqGenoFemaleSelected, gametogenesisMat);
	freqHaploMale = MATRIX_PRODUCT(freqGenoMaleSelected, gametogenesisMat);
	for (int i = 0; i < nbHaplo; i++)
	{
		freqHaploFemale(0, i) = freqHaploFemale(0, i) * femgamFit(i);
		freqHaploMale(0, i) = freqHaploMale(0, i) * malegamFit(i);
	}
	NumericMatrix freqGenoOffspring(CROSSING(nbHaplo, nbGeno, freqHaploFemale, freqHaploMale, haploCrossMat));
	return STANDARDISATION(freqGenoOffspring);
}

//' Reproduction of an hermaphroditic population
//'
//' Reproduction of a hermaphroditic population by selfing and outcrossing,
//' knowing their respective frequencies.
//'
//' @param nbHaplo number of haplotypes
//' @param nbGeno number of genotypes
//' @param selfRate selfing rate
//' @param freqGeno genotype frequencies
//' @param gametogenesisMat gametogenesis matrix
//' @param haploCrossMat haplotypes crossing matrix
//' @param femgamFit fitness of female gametes
//' @param malegamFit fitness of male gametes
//' @param femProdFit fitness for female gamete production
//' @param maleProdFit fitness for male gamete production
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix REPRODUCTION_HERMA(int nbHaplo, int nbGeno, const double selfRate, NumericMatrix freqGeno, NumericMatrix gametogenesisMat, NumericMatrix haploCrossMat,
	NumericVector femgamFit, NumericVector malegamFit, NumericVector femProdFit, NumericVector maleProdFit)
{
	if (selfRate == 1)
	{
		NumericMatrix freqGenoOffspring(SELFING(nbHaplo, nbGeno, freqGeno, gametogenesisMat, haploCrossMat, femProdFit, maleProdFit, femgamFit, malegamFit));
		return STANDARDISATION(freqGenoOffspring);
	}
	if (selfRate == 0)
	{
		NumericMatrix freqGenoOffspring(OUTCROSSING(nbHaplo, nbGeno, freqGeno, freqGeno, gametogenesisMat, haploCrossMat, femProdFit, maleProdFit, femgamFit, malegamFit));
		return STANDARDISATION(freqGenoOffspring);
	}
	NumericMatrix freqOut(OUTCROSSING(nbHaplo, nbGeno, freqGeno, freqGeno, gametogenesisMat, haploCrossMat, femProdFit, maleProdFit, femgamFit, malegamFit));
	NumericMatrix freqSelf(SELFING(nbHaplo, nbGeno, freqGeno, gametogenesisMat, haploCrossMat, femProdFit, maleProdFit, femgamFit, malegamFit));
	NumericMatrix freqGenoOffspring(1, nbGeno);
	for (int i = 0; i < nbGeno; i++) { freqGenoOffspring(0, i) = selfRate * freqSelf(0, i) + (1 - selfRate) * freqOut(0, i); }
	return STANDARDISATION(freqGenoOffspring);
}


//' Reproduction of a dioecious population
//'
//' @param nbHaplo number of haplotypes
//' @param nbGeno number of genotypes
//' @param freqGenoFemale female genotype frequencies
//' @param freqGenoMale male genotype frequencies
//' @param gametogenesisMat gametogenesis matrix
//' @param haploCrossMat haplotypes crossing matrix
//' @param femgamFit fitness of female gametes
//' @param malegamFit fitness of male gametes
//' @param femProdFit fitness for female gamete production
//' @param maleProdFit fitness for male gamete production
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix REPRODUCTION_DIOECY(int nbHaplo, int nbGeno, NumericMatrix freqGenoFemale, NumericMatrix freqGenoMale, NumericMatrix gametogenesisMat, NumericMatrix haploCrossMat,
	NumericVector femgamFit, NumericVector malegamFit, NumericVector femProdFit, NumericVector maleProdFit)
{
	NumericMatrix freqGenoOffspring(OUTCROSSING(nbHaplo, nbGeno, freqGenoFemale, freqGenoMale, gametogenesisMat, haploCrossMat, femProdFit, maleProdFit, femgamFit, malegamFit));
	return STANDARDISATION(freqGenoOffspring);
}

//' Selection on individuals
//'
//' Selection on individuals (selection on gametes and gamete production is
//' carried out during reproduction).
//'
//' @param nbGeno number of genotypes
//' @param freqGeno genotype frequencies
//' @param fitness fitness of individuals
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix INDIVIDUAL_SELECTION(int nbGeno, NumericMatrix freqGeno, NumericVector fitness)
{
	NumericMatrix freqGenoSelected(freqGeno);
	for (int i = 0; i < nbGeno; i++) { freqGenoSelected(0, i) = freqGeno(0, i) * fitness(i); }
	return STANDARDISATION(freqGenoSelected);
}


//' Simulate the model
//'
//' @param recording a boolean indicating whether to record all mutations, i.e.
//' to record allelic and genotypic frequencies along the simulations
//' @param recordGenGap the number of generations between two records during
//' simulation, if the record parameter is TRUE. Whatever the value of this
//' parameter, both the first and the last generation will be included in
//' the record
//' @param drift a boolean indicating whether genetic drift should be
//' considered (i.e. whether deterministic simulations are performed or not)
//' @param nbHaplo number of haplotypes
//' @param nbGeno number of genotypes
//' @param nbAlleles number of alleles for each loci
//' @param initGenoFreq initial genotype frequencies
//' @param gametogenesisMat gametogenesis matrix
//' @param N population size
//' @param threshold threshold for simulations
//' @param dioecy whether the population is dioecious or not (hermaphrodism)
//' @param selfRate selfing rate (only for hermaphroditic population)
//' @param stopCondition list of stop conditions
//' @param haploCrossMat haplotypes crossing matrix
//' @param alleleFreqMat matrix for calculating allelic frequencies
//' @param femgamFit fitness of female gametes
//' @param malegamFit fitness of male gametes
//' @param femindFit fitness of female individuals
//' @param maleindFit fitness of male individuals
//' @param indFit fitness of individuals (only for hermaphroditic population)
//' @param femProdFit fitness for female gamete production
//' @param maleProdFit fitness for male gamete production
//'
//' @author Ehouarn Le Faou
//'
List SIMULATION(bool recording, int recordGenGap, bool drift, int nbHaplo, int nbGeno, int nbAlleles, NumericMatrix initGenoFreq, NumericMatrix gametogenesisMat, int N, int threshold, bool dioecy,
	const double selfRate, List stopCondition, NumericMatrix haploCrossMat, NumericMatrix alleleFreqMat,
	NumericVector femgamFit, NumericVector malegamFit, NumericVector femindFit, NumericVector maleindFit, NumericVector indFit, NumericVector femProdFit, NumericVector maleProdFit)
{

	if (!dioecy) {

		NumericMatrix freqGeno(initGenoFreq);
		NumericMatrix freqAlleles(MATRIX_PRODUCT(freqGeno, alleleFreqMat));

		NumericMatrix moniFreqGeno(freqGeno);
		NumericMatrix moniFreqAlleles(freqAlleles);
		NumericVector moniGenerations(0);
		moniGenerations.push_back(0);

		int gen(0);
		while (gen < threshold && !HAVE_TO_STOP(freqAlleles, stopCondition))
		{
			freqGeno = REPRODUCTION_HERMA(nbHaplo, nbGeno, selfRate, freqGeno, gametogenesisMat, haploCrossMat, femgamFit, malegamFit, femProdFit, maleProdFit);
			freqGeno = INDIVIDUAL_SELECTION(nbGeno, freqGeno, indFit);
			if (drift) { freqGeno = DRIFT(freqGeno, N); }
			freqAlleles = MATRIX_PRODUCT(freqGeno, alleleFreqMat);
			if (recording && gen > 0 && (gen % recordGenGap) == 0)
			{
				moniFreqGeno = ROW_BIND(moniFreqGeno, freqGeno);
				moniFreqAlleles = ROW_BIND(moniFreqAlleles, freqAlleles);
				moniGenerations.push_back(gen);
			}
			gen++;
		}
		if (moniGenerations(moniGenerations.size() - 1) != gen)
		{
			moniFreqGeno = ROW_BIND(moniFreqGeno, freqGeno);
			moniFreqAlleles = ROW_BIND(moniFreqAlleles, freqAlleles);
			moniGenerations.push_back(gen);
		}
		List res = List::create(Named("gen") = gen, Named("freqGeno") = freqGeno, Named("freqHaplo") = freqGeno, Named("freqAlleles") = freqAlleles, Named("moniFreqGeno") = moniFreqGeno, Named("moniFreqAlleles") = moniFreqAlleles, Named("moniGenerations") = moniGenerations);
		return res;

	}
	else
	{
		NumericMatrix freqGeno(initGenoFreq);
		NumericMatrix freqGenoFemale(initGenoFreq);
		NumericMatrix freqGenoMale(initGenoFreq);
		NumericMatrix freqAlleles(MATRIX_PRODUCT(freqGeno, alleleFreqMat));

		NumericMatrix moniFreqGenoFemale(freqGeno);
		NumericMatrix moniFreqGenoMale(freqGeno);
		NumericMatrix moniFreqAlleles(freqAlleles);
		NumericVector moniGenerations(0);
		moniGenerations.push_back(0);

		int gen(0);
		while (gen < threshold && !HAVE_TO_STOP(freqAlleles, stopCondition))
		{
			freqGeno = REPRODUCTION_DIOECY(nbHaplo, nbGeno, freqGenoFemale, freqGenoMale, gametogenesisMat, haploCrossMat, femgamFit, malegamFit, femProdFit, maleProdFit);

			freqGenoFemale = INDIVIDUAL_SELECTION(nbGeno, freqGeno, femindFit);
			if (drift) { freqGenoFemale = DRIFT(freqGenoFemale, N); }
			freqGenoMale = INDIVIDUAL_SELECTION(nbGeno, freqGeno, maleindFit);
			if (drift) { freqGenoMale = DRIFT(freqGenoMale, N); }

			for (int i = 0; i < nbGeno; i++)
			{
				freqGeno(0, i) = freqGenoFemale(0, i) + freqGenoMale(0, i);
			}
			freqAlleles = MATRIX_PRODUCT(STANDARDISATION(freqGeno), alleleFreqMat);
			if (recording && gen > 0 && (gen % recordGenGap) == 0)
			{
				moniFreqGenoFemale = ROW_BIND(moniFreqGenoFemale, freqGenoFemale);
				moniFreqGenoMale = ROW_BIND(moniFreqGenoMale, freqGenoMale);
				moniFreqAlleles = ROW_BIND(moniFreqAlleles, freqAlleles);
				moniGenerations.push_back(gen);
			}
			gen++;
		}
		if (moniGenerations(moniGenerations.size() - 1) != gen)
		{
			moniFreqGenoFemale = ROW_BIND(moniFreqGenoFemale, freqGenoFemale);
			moniFreqGenoMale = ROW_BIND(moniFreqGenoMale, freqGenoMale);
			moniFreqAlleles = ROW_BIND(moniFreqAlleles, freqAlleles);
			moniGenerations.push_back(gen);
		}

		List res = List::create(Named("gen") = gen, Named("freqGenoFemale") = freqGenoFemale, Named("freqGenoMale") = freqGenoMale, Named("freqAlleles") = freqAlleles, Named("moniFreqGenoFemale") = moniFreqGenoFemale, Named("moniFreqGenoMale") = moniFreqGenoMale, Named("moniFreqAlleles") = moniFreqAlleles, Named("moniGenerations") = moniGenerations);
		return res;
	}
}

//' Simulate the model multiple times
//'
//' Simulation of the model multiple times with all the necessary parameters
//' set.
//'
//' @param nsim number of simulations
//' @param recording a boolean indicating whether to record all mutations, i.e.
//' to record allelic and genotypic frequencies along the simulations
//' @param recordGenGap the number of generations between two records during
//' simulation, if the record parameter is TRUE. Whatever the value of this
//' parameter, both the first and the last generation will be included in
//' the record
//' @param drift a boolean indicating whether genetic drift should be
//' considered (i.e. whether deterministic simulations are performed or not)
//' @param nbHaplo number of haplotypes
//' @param nbGeno number of genotypes
//' @param nbAlleles number of alleles for each loci
//' @param initGenoFreq initial genotype frequencies
//' @param gametogenesisMat gametogenesis matrix
//' @param N population size
//' @param threshold threshold for simulations
//' @param dioecy whether the population is dioecious or not (hermaphrodism)
//' @param selfRate selfing rate (only for hermaphroditic population)
//' @param stopCondition list of stop conditions
//' @param haploCrossMat haplotypes crossing matrix
//' @param alleleFreqMat matrix for calculating allelic frequencies
//' @param femgamFit fitness of female gametes
//' @param malegamFit fitness of male gametes
//' @param femindFit fitness of female individuals
//' @param maleindFit fitness of male individuals
//' @param indFit fitness of individuals (only for hermaphroditic population)
//' @param femProdFit fitness for female gamete production
//' @param maleProdFit fitness for male gamete production
//' @param verbose boolean determining if the progress of the simulations should be displayed or not (useful in case of many simulations)
//'
//' @return Each simulation is described by data.frame where each line corresponds to a generation (by default,
//' otherwise according to the \code{recordGenGap} parameter of the \code{simulate} method). The first one(s) describe(s)
//' the genotypic frequencies at the end of the simulations of the individuals in the case of a hermaphroditic population
//' or of the females then of the males in the case of a dioecious population. Then there are in order: the allelic
//' frequencies at the end of the simulations, the generations where the simulations stopped, the stop conditions reached.
//' Then there are the lists of records where each simulation is described by a data.frame where each line is a generation.
//' The record lists correspond to the genotypic frequencies, the allelic frequencies and finally the generations.
//' 
//' @author Ehouarn Le Faou
//'
// [[Rcpp::export]]
List SIMULATION_MULTIPLE(int nsim, bool recording, int recordGenGap, bool drift, int nbHaplo, int nbGeno, int nbAlleles, NumericMatrix initGenoFreq, NumericMatrix gametogenesisMat, int N, int threshold, bool dioecy,
	const double selfRate, List stopCondition, NumericMatrix haploCrossMat, NumericMatrix alleleFreqMat,
	NumericVector femgamFit, NumericVector malegamFit, NumericVector femindFit, NumericVector maleindFit, NumericVector indFit, NumericVector femProdFit, NumericVector maleProdFit, bool verbose)
{
	int nbSC = stopCondition.size();

	if (!dioecy) {
		NumericMatrix memGenoFreq(nsim, nbGeno);
		NumericMatrix memAlleles(nsim, nbAlleles);
		NumericMatrix memGen(nsim, 1);
		NumericMatrix memStopCondition(nsim, nbSC);
		List moniFreqGeno;
		List moniFreqAlleles;
		List moniGenerations;

		Progress prog(nsim, verbose);

		for (int i = 0; i < nsim; i++)
		{
			List res;
			res = SIMULATION(recording, recordGenGap, drift, nbHaplo, nbGeno, nbAlleles, initGenoFreq, gametogenesisMat, N, threshold, dioecy, selfRate, stopCondition,
				haploCrossMat, alleleFreqMat, femgamFit, malegamFit, femindFit, maleindFit, indFit, femProdFit, maleProdFit);

			NumericMatrix IDstops = WHICH_STOP(res["freqAlleles"], stopCondition);
			for (int j = 0; j < nbSC; j++) { memStopCondition(i, j) = IDstops(0, j); }
			NumericMatrix freqGeno = res["freqGeno"];
			for (int j = 0; j < nbGeno; j++) { memGenoFreq(i, j) = freqGeno(0, j); }
			NumericMatrix freqAlleles = res["freqAlleles"];
			for (int j = 0; j < nbAlleles; j++) { memAlleles(i, j) = freqAlleles(0, j); }
			memGen(i, 0) = res["gen"];

			if (recording)
			{
				moniFreqGeno.push_back(res["moniFreqGeno"]);
				moniFreqAlleles.push_back(res["moniFreqAlleles"]);
				moniGenerations.push_back(res["moniGenerations"]);
			}

			prog.increment(); 
		}

		return List::create(Named("genoFreq") = memGenoFreq, Named("allelesFreq") = memAlleles, Named("gen") = memGen, Named("stopCondition") = memStopCondition,
			Named("moniFreqGeno") = moniFreqGeno, Named("moniFreqAlleles") = moniFreqAlleles, Named("moniGenerations") = moniGenerations);

	}
	else
	{
		NumericMatrix memGenoFreqFemale(nsim, nbGeno);
		NumericMatrix memGenoFreqMale(nsim, nbGeno);
		NumericMatrix memAlleles(nsim, nbAlleles);
		NumericMatrix memGen(nsim, 1);
		NumericMatrix memStopCondition(nsim, nbSC);
		List moniFreqGenoFemale;
		List moniFreqGenoMale;
		List moniFreqAlleles;
		List moniGenerations;

		Progress prog(nsim, verbose);

		for (int i = 0; i < nsim; i++)
		{
			List res;
			res = SIMULATION(recording, recordGenGap, drift, nbHaplo, nbGeno, nbAlleles, initGenoFreq, gametogenesisMat, N, threshold, dioecy, selfRate, stopCondition,
				haploCrossMat, alleleFreqMat, femgamFit, malegamFit, femindFit, maleindFit, indFit, femProdFit, maleProdFit);

			NumericMatrix IDstops = WHICH_STOP(res["freqAlleles"], stopCondition);
			for (int j = 0; j < nbSC; j++) { memStopCondition(i, j) = IDstops(0, j); }
			NumericMatrix freqGenoFemale = res["freqGenoFemale"];
			for (int j = 0; j < nbGeno; j++) { memGenoFreqFemale(i, j) = freqGenoFemale(0, j); }
			NumericMatrix freqGenoMale = res["freqGenoMale"];
			for (int j = 0; j < nbGeno; j++) { memGenoFreqMale(i, j) = freqGenoMale(0, j); }
			NumericMatrix freqAlleles = res["freqAlleles"];
			for (int j = 0; j < nbAlleles; j++) { memAlleles(i, j) = freqAlleles(0, j); }
			memGen(i, 0) = res["gen"];

			if (recording)
			{
				moniFreqGenoFemale.push_back(res["moniFreqGenoFemale"]);
				moniFreqGenoMale.push_back(res["moniFreqGenoMale"]);
				moniFreqAlleles.push_back(res["moniFreqAlleles"]);
				moniGenerations.push_back(res["moniGenerations"]);
			}

			prog.increment();
		}

		return List::create(Named("genoFreqFemale") = memGenoFreqFemale, Named("genoFreqMale") = memGenoFreqMale, Named("allelesFreq") = memAlleles, Named("gen") = memGen, Named("stopCondition") = memStopCondition,
			Named("moniFreqGenoFemale") = moniFreqGenoFemale, Named("moniFreqGenoMale") = moniFreqGenoMale, Named("moniFreqAlleles") = moniFreqAlleles, Named("moniGenerations") = moniGenerations);
	}
}




