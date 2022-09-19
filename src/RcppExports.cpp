// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SIMULATION_MULTIPLE
List SIMULATION_MULTIPLE(int nsim, bool recording, int recordGenGap, bool drift, int nbHaplo, int nbGeno, int nbAlleles, NumericMatrix initGenoFreq, NumericMatrix gametogenesisMat, int N, int threshold, bool dioecy, const double selfRate, List stopCondition, NumericMatrix haploCrossMat, NumericMatrix alleleFreqMat, NumericVector femgamFit, NumericVector malegamFit, NumericVector femindFit, NumericVector maleindFit, NumericVector indFit, NumericVector femProdFit, NumericVector maleProdFit, bool verbose);
RcppExport SEXP _Ease_SIMULATION_MULTIPLE(SEXP nsimSEXP, SEXP recordingSEXP, SEXP recordGenGapSEXP, SEXP driftSEXP, SEXP nbHaploSEXP, SEXP nbGenoSEXP, SEXP nbAllelesSEXP, SEXP initGenoFreqSEXP, SEXP gametogenesisMatSEXP, SEXP NSEXP, SEXP thresholdSEXP, SEXP dioecySEXP, SEXP selfRateSEXP, SEXP stopConditionSEXP, SEXP haploCrossMatSEXP, SEXP alleleFreqMatSEXP, SEXP femgamFitSEXP, SEXP malegamFitSEXP, SEXP femindFitSEXP, SEXP maleindFitSEXP, SEXP indFitSEXP, SEXP femProdFitSEXP, SEXP maleProdFitSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< bool >::type recording(recordingSEXP);
    Rcpp::traits::input_parameter< int >::type recordGenGap(recordGenGapSEXP);
    Rcpp::traits::input_parameter< bool >::type drift(driftSEXP);
    Rcpp::traits::input_parameter< int >::type nbHaplo(nbHaploSEXP);
    Rcpp::traits::input_parameter< int >::type nbGeno(nbGenoSEXP);
    Rcpp::traits::input_parameter< int >::type nbAlleles(nbAllelesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initGenoFreq(initGenoFreqSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gametogenesisMat(gametogenesisMatSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type dioecy(dioecySEXP);
    Rcpp::traits::input_parameter< const double >::type selfRate(selfRateSEXP);
    Rcpp::traits::input_parameter< List >::type stopCondition(stopConditionSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type haploCrossMat(haploCrossMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alleleFreqMat(alleleFreqMatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type femgamFit(femgamFitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type malegamFit(malegamFitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type femindFit(femindFitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type maleindFit(maleindFitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type indFit(indFitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type femProdFit(femProdFitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type maleProdFit(maleProdFitSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(SIMULATION_MULTIPLE(nsim, recording, recordGenGap, drift, nbHaplo, nbGeno, nbAlleles, initGenoFreq, gametogenesisMat, N, threshold, dioecy, selfRate, stopCondition, haploCrossMat, alleleFreqMat, femgamFit, malegamFit, femindFit, maleindFit, indFit, femProdFit, maleProdFit, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Ease_SIMULATION_MULTIPLE", (DL_FUNC) &_Ease_SIMULATION_MULTIPLE, 24},
    {NULL, NULL, 0}
};

RcppExport void R_init_Ease(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}