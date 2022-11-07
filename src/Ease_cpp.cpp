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

//' Matrix col bind
//'
//' Binding of two matrices by their columns. Both matrices must have the same
//' number of columns.
//'
//' @param MAT1 matrix
//' @param MAT2 matrix
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix COL_BIND(NumericMatrix MAT1, NumericMatrix MAT2)
{
	int nCol1(MAT1.ncol());
	int nCol2(MAT2.ncol());
	int nRow(MAT1.nrow());
	NumericMatrix MAT_res(nRow, nCol1 + nCol2);
	for (int i = 0; i < nRow; i++) {
		for (int j = 0; j < nCol1; j++) {
			MAT_res(i, j) = MAT1(i, j);
		}
		for (int j = 0; j < nCol2; j++) {
			MAT_res(i, j + nCol1) = MAT2(i, j);
		}
	}
	return MAT_res;
}

//' Logical matrix col bind
//'
//' Binding of two logical matrices by their columns. Both matrices must have the same
//' number of columns.
//'
//' @param MAT1 logical matrix
//' @param MAT2 logical matrix
//'
//' @author Ehouarn Le Faou
//'
LogicalMatrix BOOL_COL_BIND(LogicalMatrix MAT1, LogicalMatrix MAT2)
{
	int nCol1(MAT1.ncol());
	int nCol2(MAT2.ncol());
	int nRow(MAT1.nrow());
	LogicalMatrix MAT_res(nRow, nCol1 + nCol2);
	for (int i = 0; i < nRow; i++) {
		for (int j = 0; j < nCol1; j++) {
			MAT_res(i, j) = MAT1(i, j);
		}
		for (int j = 0; j < nCol2; j++) {
			MAT_res(i, j + nCol1) = MAT2(i, j);
		}
	}
	return MAT_res;
}

//' Sum of a matrix
//'
//' Sum of all the cells of one matrix.
//'
//' @param MAT matrix
//'
//' @author Ehouarn Le Faou
//'
double SUM_MAT(NumericMatrix MAT)
{
	double res = 0;
	for (int i = 0; i < MAT.nrow(); i++) {
		for (int j = 0; j < MAT.ncol(); j++) {
			res += MAT(i, j);
		}
	}
	return res;
}

//' Mean of two matrices
//'
//' Calculates the average of the values of each cell one by one.
//'
//' @param MAT1 matrix
//' @param MAT2 matrix
//'
//' @author Ehouarn Le Faou
//'
NumericMatrix MEAN_MATS(NumericMatrix MAT1, NumericMatrix MAT2)
{
	int nCol(MAT1.ncol());
	int nRow(MAT1.nrow());
	NumericMatrix MAT_res(nRow, nCol);
	for (int i = 0; i < nRow; i++) {
		for (int j = 0; j < nCol; j++) {
			MAT_res(i, j) = (MAT1(i, j) + MAT2(i, j))/2;
		}
	}
	return MAT_res;
}

//' Binding matrices on matrices
//'
//' Links to each of the matrices in the REF_LIST list the
//' matrices in the ADD_LIST list, respectively.
//'
//' @param REF_LIST list of reference matrices
//' @param ADD_LIST list of matrices to ve added
//'
//' @author Ehouarn Le Faou
//'
List APPEND_RECORD_MATRIX_LIST(List REF_LIST, List ADD_LIST)
{
	int k = REF_LIST.size();
	List RES_LIST = List::create();
	for (int i = 0; i < k; i++)
	{
		RES_LIST.push_back(ROW_BIND(REF_LIST(i), ADD_LIST(i)));
	}
	RES_LIST.names() = REF_LIST.names();
	return RES_LIST;
}

//' Binding vectors on vectors
//'
//' Links to each of the vectors in the REF_LIST list the
//' vectors in the ADD_LIST list, respectively.
//'
//' @param REF_LIST list of reference vectors
//' @param ADD_LIST list of vectors to ve added
//'
//' @author Ehouarn Le Faou
//'
List APPEND_RECORD_VECTOR_LIST(List REF_LIST, List ADD_LIST)
{
	int k = REF_LIST.size();
	List RES_LIST = List::create();
	for (int i = 0; i < k; i++)
	{
		NumericVector vect = REF_LIST(i);
		vect.push_back(ADD_LIST(i));
		RES_LIST.push_back(vect);
	}
	RES_LIST.names() = REF_LIST.names();
	return RES_LIST;
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
NumericMatrix rmultinom_rcpp(unsigned int& n, unsigned int& size, NumericVector& probs) {
	unsigned int N = probs.length();
	NumericMatrix sim(N, n);
	for (unsigned int i = 0; i < n; i++) {
		sim(_, i) = rmultinom_1(size, probs, N);
	}
	return sim;
}


IntegerVector rpois_rcpp(unsigned int& n, NumericVector& lambda) {
	unsigned int lambda_i = 0;
	IntegerVector sim(n);
	for (unsigned int i = 0; i < n; i++) {
		sim(i) = R::rpois(lambda[lambda_i]);
		// update lambda_i to match next realized value with correct mean
		lambda_i++;
		// restart lambda_i at 0 if end of lambda reached
		if (lambda_i == lambda.length()) {
			lambda_i = 0;
		}
	}
	return sim;
}

double rpois_simul(double lambda) {
	NumericVector lambda_vect = { lambda };
	unsigned int n = 1;
	IntegerVector draw = rpois_rcpp(n, lambda_vect);
	return (double)draw(0);
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

	return STANDARDISATION(freq2);
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
LogicalMatrix WHICH_STOP(NumericMatrix freqAlleles, List stopCondition)
{
	LogicalMatrix idStop(1, stopCondition.size());
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
	return freqGenoSelected;
}

double COMPUTE_MEAN_FITNESS(NumericMatrix freqs, NumericVector fitness)
{
	int nCol(freqs.ncol());
	double res = 0;
	for (int j = 0; j < nCol; j++)
	{
		res += freqs(0, j) * fitness(j);
	}	
	return res;
}

List COMPUTE_MEAN_FITNESS_LIST(List freqs_list, List fitness_list)
{
	int k = freqs_list.size();
	List res_list;
	for (int i = 0; i < k; i++)
	{
		res_list.push_back(COMPUTE_MEAN_FITNESS(freqs_list(i), fitness_list(i)));
	}
	res_list.names() = freqs_list.names();
	return res_list;
}


class Population
{
	// Access specifier
public:

	// Constant population attributes
	CharacterVector id;
	int nbHaplo;
	int nbGeno;
	CharacterVector idGeno;
	int nbAlleles;
	CharacterVector idAlleles;
	int popSize;
	bool dioecy;
	double selfRate;
	double growthRate;
	bool demography;
	int initPopSize;
	NumericMatrix initGenoFreq;
	// Constant matrices
	NumericMatrix meiosisMat;
	NumericMatrix gametogenesisMat;
	NumericMatrix haploCrossMat;
	NumericMatrix alleleFreqMat;
	// Constant selection vectors
	List gamFit;
	List indFit;
	List gamProdFit;
	// Simulation-related constants
	int threshold;
	bool recording;
	int recordGenGap;
	bool drift;
	List stopCondition;
	CharacterVector IDstopCondition;
	String nameOutFunct;
	// Simulation-related variables
	int gen;
	double popSize_current;
	List freqGeno;
	List freqHaplo;
	List freqAlleles;
	List indMeanFitness;
	List gamProdMeanFitness;
	List gamMeanFitness;
	List moniFreqGeno;
	List moniFreqAlleles;
	List moniIndMeanFitness;
	List moniGamProdMeanFitness;
	List moniGamMeanFitness;
	List customOut;
	NumericVector moniGenerations;
	NumericVector moniPopSize;
	LogicalMatrix IDstops;

	void compute_freq_haplo()
	{
		if (dioecy) {
			freqHaplo["female"] = MATRIX_PRODUCT(freqGeno["female"], meiosisMat);
			freqHaplo["male"] = MATRIX_PRODUCT(freqGeno["male"], meiosisMat);
		}
		else
		{
			freqHaplo["female"] = MATRIX_PRODUCT(freqGeno["ind"], meiosisMat);
			freqHaplo["male"] = MATRIX_PRODUCT(freqGeno["ind"], meiosisMat);
		}
	}

	void compute_mean_fitness()
	{
		if (dioecy) {
			indMeanFitness = COMPUTE_MEAN_FITNESS_LIST(List::create(Named("female") = freqGeno["female"], Named("male") = freqGeno["male"]),
				List::create(indFit["female"], indFit["male"]));
			gamProdMeanFitness = COMPUTE_MEAN_FITNESS_LIST(List::create(Named("female") = freqGeno["female"], Named("male") = freqGeno["male"]), gamProdFit);
		}
		else {
			indMeanFitness = COMPUTE_MEAN_FITNESS_LIST(freqGeno, indFit);
			gamProdMeanFitness = COMPUTE_MEAN_FITNESS_LIST(List::create(Named("female") = freqGeno["ind"], Named("male") = freqGeno["ind"]), gamProdFit);
		}
		gamMeanFitness = COMPUTE_MEAN_FITNESS_LIST(freqHaplo, gamFit);
	}

	void outFunct()
	{
		Function f(nameOutFunct);
		List pop = List::create(
			Named("customOutput") = customOut,
			Named("gen") = gen,
			Named("freqGeno") = freqGeno,
			Named("freqHaplo") = freqHaplo,
			Named("freqAlleles") = freqAlleles
		);
		List res = f(pop);
		if (res[0]) { customOut.push_back(res[1]); }
	}


	void reset()
	{
		// Genotypic and allelic frequencies
		if (dioecy) {
			freqGeno = List::create(Named("ind") = initGenoFreq, Named("female") = initGenoFreq, Named("male") = initGenoFreq);
			NumericMatrix FA = MATRIX_PRODUCT(freqGeno["ind"], alleleFreqMat);
			freqAlleles = List::create(Named("ind") = FA, Named("female") = FA, Named("male") = FA);
			NumericMatrix FH = MATRIX_PRODUCT(freqGeno["ind"], meiosisMat);
			freqHaplo = List::create(Named("female") = FH, Named("male") = FH);
			moniFreqGeno = List::create(Named("ind") = freqGeno["ind"], Named("female") = freqGeno["female"], Named("male") = freqGeno["male"]);
			moniFreqAlleles = List::create(Named("ind") = freqAlleles["ind"], Named("female") = freqAlleles["female"], Named("male") = freqAlleles["male"]);
		}
		else {
			freqGeno = List::create(Named("ind") = initGenoFreq);
			freqAlleles = List::create(Named("ind") = MATRIX_PRODUCT(freqGeno["ind"], alleleFreqMat));
			NumericMatrix FH = MATRIX_PRODUCT(freqGeno["ind"], meiosisMat);
			freqHaplo = List::create(Named("female") = FH, Named("male") = FH);
			moniFreqGeno = List::create(Named("ind") = freqGeno["ind"]);
			moniFreqAlleles = List::create(Named("ind") = freqAlleles["ind"]);
		}

		// Mean fitnesses
		compute_mean_fitness();
		moniIndMeanFitness = indMeanFitness;
		moniGamProdMeanFitness = gamProdMeanFitness;
		moniGamMeanFitness = gamMeanFitness;

		// Population parameters
		gen = 0;
		moniGenerations = 0;
		moniPopSize = initPopSize;

		// Custom output
		customOut = List::create();
		outFunct();
	}

	Population(CharacterVector ID, bool RECORDING, int RECORD_GEN_GAP, bool DRIFT, int NB_HAPLO, int NB_GENO, CharacterVector ID_GENO, int NB_ALLELES, CharacterVector ID_ALLELES, NumericMatrix INIT_GENO_FREQ, NumericMatrix MEIOSIS_MAT,NumericMatrix GAMETOGENESIS_MAT, int POP_SIZE, int THRESHOLD, bool DIOECY,
		double SELF_RATE, List STOP_CONDITION, CharacterVector ID_STOP_CONDITION, NumericMatrix HAPLO_CROSS_MAT, NumericMatrix ALLELE_FREQ_MAT, List GAM_FIT, List IND_FIT, List GAM_PROD_FIT, bool DEMOGRAPHY, double GROWTH_RATE, int INIT_POP_SIZE, String NAME_OUT_FUNCT)
	{
		// Definition of attributes
		id = ID;
		nbHaplo = NB_HAPLO;
		nbGeno = NB_GENO;
		idGeno = ID_GENO;
		nbAlleles = NB_ALLELES;
		idAlleles = ID_ALLELES;
		popSize = POP_SIZE;
		popSize_current = INIT_POP_SIZE;
		dioecy = DIOECY;
		selfRate = SELF_RATE;
		meiosisMat = MEIOSIS_MAT;
		gametogenesisMat = GAMETOGENESIS_MAT;
		haploCrossMat = HAPLO_CROSS_MAT;
		alleleFreqMat = ALLELE_FREQ_MAT;
		gamFit = GAM_FIT;
		indFit = IND_FIT;
		gamProdFit = GAM_PROD_FIT;
		threshold = THRESHOLD;
		recording = RECORDING;
		recordGenGap = RECORD_GEN_GAP;
		drift = DRIFT;
		stopCondition = STOP_CONDITION;
		IDstopCondition = ID_STOP_CONDITION;
		IDstopCondition.push_back("unstopped");
		initGenoFreq = INIT_GENO_FREQ;
		demography = DEMOGRAPHY;
		growthRate = GROWTH_RATE;
		initPopSize = INIT_POP_SIZE;
		nameOutFunct = NAME_OUT_FUNCT;

		reset();
	}

	void demography_draw()
	{
		popSize_current = popSize_current + growthRate * popSize_current * (1 - popSize_current/popSize);
		popSize_current = rpois_simul(popSize_current);
	}

	void genetic_drift()
	{
		if (dioecy) {
			freqGeno["female"] = DRIFT(freqGeno["female"], (int)popSize_current / 2);
			freqGeno["male"] = DRIFT(freqGeno["male"], (int)popSize_current / 2);
			freqGeno["ind"] = MEAN_MATS(freqGeno["male"], freqGeno["female"]);
		}
		else {
			freqGeno["ind"] = DRIFT(freqGeno["ind"], (int)popSize_current);
		}
	}

	void reproduction()
	{
		if (dioecy) {
			NumericMatrix freqGenoOffspring(OUTCROSSING(nbHaplo, nbGeno, freqGeno["female"], freqGeno["male"], gametogenesisMat, haploCrossMat, gamProdFit["female"], gamProdFit["male"], gamFit["female"], gamFit["male"]));
			freqGeno["ind"] = freqGenoOffspring;
			freqGeno["female"] = freqGenoOffspring;
			freqGeno["male"] = freqGenoOffspring;
		}
		else {
			if (selfRate == 1) {
				NumericMatrix freqGenoOffspring(SELFING(nbHaplo, nbGeno, freqGeno["ind"], gametogenesisMat, haploCrossMat, gamProdFit["female"], gamProdFit["male"], gamFit["female"], gamFit["male"]));
				freqGeno["ind"] = freqGenoOffspring;
			}
			else if (selfRate == 0) {
				NumericMatrix freqGenoOffspring(OUTCROSSING(nbHaplo, nbGeno, freqGeno["ind"], freqGeno["ind"], gametogenesisMat, haploCrossMat, gamProdFit["female"], gamProdFit["male"], gamFit["female"], gamFit["male"]));
				freqGeno["ind"] = freqGenoOffspring;
			}
			else {
				NumericMatrix freqGenoSelfOffspring(SELFING(nbHaplo, nbGeno, freqGeno["ind"], gametogenesisMat, haploCrossMat, gamProdFit["female"], gamProdFit["male"], gamFit["female"], gamFit["male"]));
				NumericMatrix freqGenoOutOffspring(OUTCROSSING(nbHaplo, nbGeno, freqGeno["ind"], freqGeno["ind"], gametogenesisMat, haploCrossMat, gamProdFit["female"], gamProdFit["male"], gamFit["female"], gamFit["male"]));
				NumericMatrix freqGenoOffspring(1, nbGeno);
				for (int i = 0; i < nbGeno; i++) { freqGenoOffspring(0, i) = selfRate * freqGenoSelfOffspring(0, i) + (1 - selfRate) * freqGenoOutOffspring(0, i); }
				freqGeno["ind"] = freqGenoOffspring;
				freqGeno["female"] = freqGenoOffspring;
				freqGeno["male"] = freqGenoOffspring;
			}
		}
	}

	void individual_selection()
	{
		if (dioecy) {
			freqGeno["female"] = INDIVIDUAL_SELECTION(nbGeno, freqGeno["female"], indFit["female"]);
			freqGeno["male"] = INDIVIDUAL_SELECTION(nbGeno, freqGeno["male"], indFit["male"]);
			freqGeno["ind"] = MEAN_MATS(freqGeno["male"], freqGeno["female"]);
		}
		else {
			freqGeno["ind"] = INDIVIDUAL_SELECTION(nbGeno, freqGeno["ind"], indFit["ind"]);
		}
	}

	void compute_allele_freq()
	{
		if (dioecy) {
			freqAlleles["ind"] = MATRIX_PRODUCT(STANDARDISATION(freqGeno["ind"]), alleleFreqMat);
			freqAlleles["female"] = MATRIX_PRODUCT(STANDARDISATION(freqGeno["female"]), alleleFreqMat);
			freqAlleles["male"] = MATRIX_PRODUCT(STANDARDISATION(freqGeno["male"]), alleleFreqMat);
		}
		else {
			freqAlleles["ind"] = MATRIX_PRODUCT(STANDARDISATION(freqGeno["ind"]), alleleFreqMat);
		}
	}

	void standardisation()
	{
		freqGeno["ind"] = STANDARDISATION(freqGeno["ind"]);
		if (dioecy) {
			freqGeno["female"] = STANDARDISATION(freqGeno["female"]);
			freqGeno["male"] = STANDARDISATION(freqGeno["male"]);
		}
	}

	void recordings()
	{
		standardisation();

		moniFreqGeno = APPEND_RECORD_MATRIX_LIST(moniFreqGeno, freqGeno);
		moniFreqAlleles = APPEND_RECORD_MATRIX_LIST(moniFreqAlleles, freqAlleles);
		moniGenerations.push_back(gen);
		moniPopSize.push_back(popSize_current);

		compute_mean_fitness();
		moniIndMeanFitness = APPEND_RECORD_VECTOR_LIST(moniIndMeanFitness, indMeanFitness);
		moniGamProdMeanFitness = APPEND_RECORD_VECTOR_LIST(moniGamProdMeanFitness, gamProdMeanFitness);
		moniGamMeanFitness = APPEND_RECORD_VECTOR_LIST(moniGamMeanFitness, gamMeanFitness);
	}

	void next_generation()
	{
		// Generation incrementing
		gen++;

		// Reproduction (selection on gamete production + gametogenesis + mutation + selection on gametes + syngamy)
		reproduction();

		// Selection on individuals
		individual_selection();

		// Are there any living individuals left?
		if (SUM_MAT(freqGeno["ind"]) == 0.0) { popSize_current = 0; }

		// Demography
		if (demography) { demography_draw(); }

		standardisation();

		// Drift
		if (drift) { genetic_drift(); }

		// Computing the allelic and haplotypic frequencies
		compute_allele_freq();
		compute_freq_haplo();

		// Recording of frequencies (allelic and genotypic)
		if (recording && gen > 0 && (gen % recordGenGap) == 0) { recordings(); }

		// Custom output
		outFunct();
	}

	void simulate()
	{
		while (gen < threshold && !HAVE_TO_STOP(freqAlleles["ind"], stopCondition))
		{
			next_generation();
		}
		IDstops = WHICH_STOP(freqAlleles["ind"], stopCondition);
		if (moniGenerations(moniGenerations.size() - 1) != gen) { recordings();	}
	}

	List output(bool generations = true, bool populationSize = true, bool genotypicFreq = true, bool allelicFreq = true, bool stop = true, bool meanFitness = true, 
		bool recordGenoFreq = true, bool recordAlleleFreq = true, bool recordGen = true, bool recordPopSize = true, bool recordMeanFitness = true)
	{
		// Results list
		List res;
		if (generations) { 
			NumericMatrix genMat(1, 1);
			genMat(0, 0) = gen;
			res.push_back(genMat, "gen");
		}
		if (populationSize) { 
			NumericMatrix popSizeMat(1, 1);
			popSizeMat(0, 0) = popSize_current;
			res.push_back(popSizeMat, "popSize");
		}
		if (genotypicFreq) {
			colnames(freqGeno["ind"]) = idGeno;
			res.push_back(freqGeno["ind"]);
			if (dioecy) {
				colnames(freqGeno["female"]) = idGeno;
				colnames(freqGeno["male"]) = idGeno;
				res.push_back(freqGeno["female"], "(fem)");
				res.push_back(freqGeno["male"], "(mal)");
			}
		}
		if (allelicFreq) {
			colnames(freqAlleles["ind"]) = idAlleles;
			res.push_back(freqAlleles["ind"]);
			if (dioecy) {
				colnames(freqAlleles["female"]) = idAlleles;
				colnames(freqAlleles["male"]) = idAlleles;
				res.push_back(freqAlleles["female"], "(fem)");
				res.push_back(freqAlleles["male"], "(mal)");
			}
		}
		if (stop) {
			colnames(IDstops) = IDstopCondition;
			res.push_back(IDstops); 
		}
		if (meanFitness) {
			CharacterVector cn;
			if (dioecy) {
				cn = { "(fem)indMeanFit", "(mal)indMeanFit" };
				NumericMatrix indMeanFitnessMat(1, 2);
				indMeanFitnessMat(0, 0) = indMeanFitness["female"];
				indMeanFitnessMat(0, 1) = indMeanFitness["male"];
				colnames(indMeanFitnessMat) = cn;
				res.push_back(indMeanFitnessMat);
			}
			else {
				NumericMatrix indMeanFitnessMat(1, 1);
				indMeanFitnessMat(0, 0) = indMeanFitness["ind"];
				colnames(indMeanFitnessMat) = (CharacterVector)"indMeanFit";
				res.push_back(indMeanFitnessMat);
			}

			cn = { "(fem)gamProdMeanFit", "(mal)gamProdMeanFit" };
			NumericMatrix gamProdMeanFitnessMat(1, 2);
			gamProdMeanFitnessMat(0, 0) = gamProdMeanFitness["female"];
			gamProdMeanFitnessMat(0, 1) = gamProdMeanFitness["male"];
			colnames(gamProdMeanFitnessMat) = cn;
			res.push_back(gamProdMeanFitnessMat);

			cn = { "(fem)gamMeanFit", "(mal)gamMeanFit" };
			NumericMatrix gamMeanFitnessMat(1, 2);
			gamMeanFitnessMat(0, 0) = gamMeanFitness["female"];
			gamMeanFitnessMat(0, 1) = gamMeanFitness["male"];
			colnames(gamMeanFitnessMat) = cn;
			res.push_back(gamMeanFitnessMat);
		}

		// Record list
		List rec;
		if (recording) {

			if (recordGenoFreq) {
				colnames(moniFreqGeno["ind"]) = idGeno;
				rec.push_back(moniFreqGeno["ind"]);
				if (dioecy) {
					colnames(moniFreqGeno["female"]) = idGeno;
					colnames(moniFreqGeno["male"]) = idGeno;
					rec.push_back(moniFreqGeno["female"], "(fem)");
					rec.push_back(moniFreqGeno["male"], "(mal)");
				}
			}
			if (recordAlleleFreq) {
				colnames(moniFreqAlleles["ind"]) = idAlleles;
				rec.push_back(moniFreqAlleles["ind"]);
				if (dioecy) {
					colnames(moniFreqAlleles["female"]) = idAlleles;
					colnames(moniFreqAlleles["male"]) = idAlleles;
					rec.push_back(moniFreqAlleles["female"], "(fem)");
					rec.push_back(moniFreqAlleles["male"], "(mal)");
				}
			}
			if (recordGen) { 
				moniGenerations.attr("dim") = Dimension(moniGenerations.length(), 1);
				NumericMatrix moniGenerationsMat = as<NumericMatrix>(moniGenerations);
				rec.push_back(moniGenerationsMat, "gen");
			}
			if (recordPopSize) {
				moniPopSize.attr("dim") = Dimension(moniPopSize.length(), 1);
				NumericMatrix moniPopSizeMat = as<NumericMatrix>(moniPopSize);
				rec.push_back(moniPopSizeMat, "popSize");
			}
			if (recordMeanFitness) {
				NumericVector outVect;
				NumericMatrix outMat;
				if (dioecy)
				{
					outVect = moniIndMeanFitness["female"];
					outVect.attr("dim") = Dimension(moniGenerations.length(), 1);
					outMat = as<NumericMatrix>(outVect);
					rec.push_back(outMat, "(fem)indMeanFit");

					outVect = moniIndMeanFitness["male"];
					outVect.attr("dim") = Dimension(moniGenerations.length(), 1);
					outMat = as<NumericMatrix>(outVect);
					rec.push_back(outMat, "(mal)indMeanFit");
				}
				else {
					outVect = moniIndMeanFitness["ind"];
					outVect.attr("dim") = Dimension(moniGenerations.length(), 1);
					outMat = as<NumericMatrix>(outVect);
					rec.push_back(outMat, "indMeanFit");
				}

				outVect = moniGamProdMeanFitness["male"];
				outVect.attr("dim") = Dimension(moniGenerations.length(), 1);
				outMat = as<NumericMatrix>(outVect);
				rec.push_back(outMat, "(mal)gamProdMeanFit");

				outVect = moniGamProdMeanFitness["male"];
				outVect.attr("dim") = Dimension(moniGenerations.length(), 1);
				outMat = as<NumericMatrix>(outVect);
				rec.push_back(outMat, "(mal)gamProdMeanFit");

				outVect = moniGamMeanFitness["male"];
				outVect.attr("dim") = Dimension(moniGenerations.length(), 1);
				outMat = as<NumericMatrix>(outVect);
				rec.push_back(outMat, "(mal)gamMeanFit");

				outVect = moniGamMeanFitness["male"];
				outVect.attr("dim") = Dimension(moniGenerations.length(), 1);
				outMat = as<NumericMatrix>(outVect);
				rec.push_back(outMat, "(mal)gamMeanFit");
			}
		}
			
		return List::create(Named("results") = res, Named("records") = rec, Named("custom") = customOut);
	}
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Metapopulation

void GENOTYPE_MIGRATION(list<Population>& populations, NumericMatrix migMat, int nbGeno, int nbPop, int indTypeIndex)
{
	NumericMatrix freqGenoMat(nbGeno, nbPop);
	int i = 0;
	for (Population& pop : populations) {
		if (pop.popSize_current > 0) {}
		NumericMatrix fG = pop.freqGeno[indTypeIndex];
		for (int j = 0; j < nbGeno; j++)
		{
			if (pop.popSize_current > 0)
			{
				freqGenoMat(j, i) = fG(0, j);
			}
			else
			{
				freqGenoMat(j, i) = 0;
			}
		}
		i++;
	}
	NumericMatrix newFreqGenoMat = MATRIX_PRODUCT(freqGenoMat, migMat);
	i = 0;
	for (Population& pop : populations) {
		NumericMatrix fG = pop.freqGeno[indTypeIndex];
		for (int j = 0; j < nbGeno; j++)
		{
			fG(0, j) = newFreqGenoMat(j, i);
		}
		pop.freqGeno["ind"] = STANDARDISATION(fG);
		i++;
	}
}


class Metapopulation
{
	// Access specifier
public:

	list<Population> populations;
	int nbPop;
	NumericMatrix migMat;
	List ids;
	int threshold;
	int nbAlleles;
	int nbLoci;
	int nbGeno;
	int gen;
	bool dioecy;
	bool recording;

	void reset()
	{
		for (Population& pop : populations) {
			pop.reset();
		}
	}

	Metapopulation(int NB_POP, NumericMatrix MIG_MAT, List IDs, bool RECORDING, int RECORD_GEN_GAP, bool DRIFT, int NB_HAPLO, int NB_GENO, CharacterVector ID_GENO, int NB_ALLELES, CharacterVector ID_ALLELES, int NB_LOCI,
		List INIT_GENO_FREQs, NumericMatrix MEIOSIS_MAT, NumericMatrix GAMETOGENESIS_MAT, List POP_SIZEs, int THRESHOLD, bool DIOECY, List DEMOGRAPHIES, List GROWTH_RATES,
		List SELF_RATEs, List STOP_CONDITION, CharacterVector ID_STOP_CONDITION, NumericMatrix HAPLO_CROSS_MAT, NumericMatrix ALLELE_FREQ_MAT, List GAM_FIT, List IND_FIT, List GAM_PROD_FIT,
		List INIT_POP_SIZEs, String NAME_OUT_FUNCT)
	{
		// Definition of attributes
		dioecy = DIOECY;
		gen = 0;
		nbPop = NB_POP;
		migMat = MIG_MAT;
		nbAlleles = NB_ALLELES;
		nbGeno = NB_GENO;
		nbLoci = NB_LOCI;
		threshold = THRESHOLD;
		ids = IDs;
		recording = RECORDING;

		list<Population> pops;
		for (int i = 0; i < NB_POP; i++)
		{
			Population pop(IDs(i), RECORDING, RECORD_GEN_GAP, DRIFT, NB_HAPLO, NB_GENO, ID_GENO, NB_ALLELES, ID_ALLELES, INIT_GENO_FREQs(i), MEIOSIS_MAT, GAMETOGENESIS_MAT,
				POP_SIZEs(i), THRESHOLD, DIOECY, SELF_RATEs(i), STOP_CONDITION, ID_STOP_CONDITION, HAPLO_CROSS_MAT, ALLELE_FREQ_MAT, GAM_FIT(i), IND_FIT(i),
				GAM_PROD_FIT(i), DEMOGRAPHIES(i), GROWTH_RATES(i), INIT_POP_SIZEs(i), NAME_OUT_FUNCT);
			pops.push_back(pop);
		}
		populations = pops;

		reset();
	}

	void migration()
	{
		if (nbPop > 1)
		{
			// Migration of individuals
			NumericMatrix popSizes(1, nbPop);
			int i = 0;
			for (Population& pop : populations) {
				popSizes(0, i) = pop.popSize_current;
				i++;
			}
			popSizes = MATRIX_PRODUCT(popSizes, migMat);
			i = 0;
			for (Population& pop : populations) {
				if (pop.demography) { pop.popSize_current = popSizes(0, i); }
				i++;
			}

			if (dioecy)
			{
				GENOTYPE_MIGRATION(populations, migMat, nbGeno, nbPop, 0);
				GENOTYPE_MIGRATION(populations, migMat, nbGeno, nbPop, 1);
				GENOTYPE_MIGRATION(populations, migMat, nbGeno, nbPop, 2);
			}
			else
			{
				GENOTYPE_MIGRATION(populations, migMat, nbGeno, nbPop, 0);
			}
		}
	}

	void simulate()
	{
		gen = 0;
		bool hts = false;
		while (gen < threshold && !hts)
		{
			migration();
			hts = true;
			for (Population& pop : populations) {
				pop.next_generation();
				hts = hts && HAVE_TO_STOP(pop.freqAlleles["ind"], pop.stopCondition);
			}
			gen++;
		}
		for (Population& pop : populations) {
			pop.IDstops = WHICH_STOP(pop.freqAlleles["ind"], pop.stopCondition);
			LogicalMatrix tba(1, 1);
			if (gen == threshold) {	tba(0, 0) = true; }
			else { tba(0, 0) = false; }
			pop.IDstops = BOOL_COL_BIND(pop.IDstops, tba);
		}
		list<Population>::iterator pop = populations.begin();
		NumericVector moniGenerations = pop->moniGenerations;
		if (moniGenerations(moniGenerations.size() - 1) != gen) {
			for (Population& pop : populations) {
				pop.recordings();
			}
		}
	}

	List output(bool generations = true, bool populationSize = true, bool genotypicFreq = true, bool allelicFreq = true, bool stop = true, bool meanFitness = true, 
		bool recordGenoFreq = true, bool recordAlleleFreq = true, bool recordGen = true, bool recordPopSize = true, bool recordMeanFitness = true)
	{
		List res;
		for (Population& pop : populations) {
			if (recording) { res.push_back(pop.output(generations, populationSize, genotypicFreq, allelicFreq, stop, meanFitness, recordGenoFreq, recordAlleleFreq, recordGen, recordPopSize, recordMeanFitness)); }
			else { res.push_back(pop.output(generations, populationSize, genotypicFreq, allelicFreq, stop, meanFitness, false, false, false, false, false)); }
		}
		res.names() = ids;
		return res;
	}


};


//' Simulation of a metapopulation
//' 
//' @param nbPop number of populations in the metapopulation
//' @param ids population IDs
//' @param migMat migration matrix
//' @param nsim number of simulations
//' @param verbose boolean determining if the progress of the simulations should be displayed or not (useful in case of many simulations)
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
//' @param idGeno genotypes ID
//' @param nbAlleles number of alleles for each loci
//' @param idAlleles alleles ID
//' @param nbLoci number of loci
//' @param initGenoFreq list of initial genotype frequencies in the populations
//' @param meiosisMat meiosis matrix
//' @param gametogenesisMat gametogenesis matrix
//' @param popSize list population sizes
//' @param threshold threshold for simulations
//' @param dioecy whether the population(s) is dioecious or not (hermaphrodism)
//' @param selfRate list of the selfing rate in populations (only for hermaphroditic population)
//' @param stopCondition list of stop conditions
//' @param IDstopCondition vector of stop condition ID
//' @param haploCrossMat haplotypes crossing matrix
//' @param alleleFreqMat matrix for calculating allelic frequencies
//' @param gamFit fitness of gametes
//' @param indFit fitness of individuals
//' @param gamProdFit fitness for gamete production
//' @param demography list of population demographies
//' @param growthRate list of population growth rates
//' @param initPopSize list of initial population
//' @param nameOutFunct name of the custom output function
//'
//' @author Ehouarn Le Faou
//'
// [[Rcpp::export]]
List METAPOP_SIMULATION(int nbPop, List ids, NumericMatrix migMat, int nsim, bool verbose, bool recording, int recordGenGap, bool drift, 
	int nbHaplo, int nbGeno, CharacterVector idGeno, int nbAlleles, CharacterVector idAlleles, int nbLoci, List initGenoFreq, NumericMatrix meiosisMat, NumericMatrix gametogenesisMat, List popSize, int threshold, bool dioecy,
	List selfRate, List stopCondition, CharacterVector IDstopCondition, NumericMatrix haploCrossMat, NumericMatrix alleleFreqMat, List gamFit, List indFit, List gamProdFit, List demography, List growthRate, List initPopSize,
	String nameOutFunct) {
	
	Metapopulation metapop(nbPop, migMat, ids, recording, recordGenGap, drift, nbHaplo, nbGeno, idGeno, nbAlleles, idAlleles, nbLoci,
		initGenoFreq, meiosisMat,gametogenesisMat, popSize, threshold, dioecy, demography, growthRate,
		selfRate, stopCondition, IDstopCondition, haploCrossMat, alleleFreqMat, gamFit, indFit, gamProdFit,
		initPopSize, nameOutFunct);

	List res;
	Progress prog(nsim, verbose);
	for (int i = 0; i < nsim; i++)
	{
		metapop.reset();
		metapop.simulate();
		res.push_back(metapop.output(), "s" + to_string(i));
		prog.increment();
	}

	return res;
}
