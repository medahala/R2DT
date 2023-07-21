// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// UdetC
NumericVector UdetC(NumericVector ep, double er, double el, double eag, double eal, NumericVector tp, double tr, double tl, double tag, double tal, double k1, double k2);
RcppExport SEXP _doseDet_UdetC(SEXP epSEXP, SEXP erSEXP, SEXP elSEXP, SEXP eagSEXP, SEXP ealSEXP, SEXP tpSEXP, SEXP trSEXP, SEXP tlSEXP, SEXP tagSEXP, SEXP talSEXP, SEXP k1SEXP, SEXP k2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ep(epSEXP);
    Rcpp::traits::input_parameter< double >::type er(erSEXP);
    Rcpp::traits::input_parameter< double >::type el(elSEXP);
    Rcpp::traits::input_parameter< double >::type eag(eagSEXP);
    Rcpp::traits::input_parameter< double >::type eal(ealSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tp(tpSEXP);
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< double >::type tl(tlSEXP);
    Rcpp::traits::input_parameter< double >::type tag(tagSEXP);
    Rcpp::traits::input_parameter< double >::type tal(talSEXP);
    Rcpp::traits::input_parameter< double >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< double >::type k2(k2SEXP);
    rcpp_result_gen = Rcpp::wrap(UdetC(ep, er, el, eag, eal, tp, tr, tl, tag, tal, k1, k2));
    return rcpp_result_gen;
END_RCPP
}
// efftox_add
int efftox_add(NumericMatrix ep, NumericMatrix tp, double r, double pie, double pit, double ptt, double pet, double ptc, double pec, int h_treated);
RcppExport SEXP _doseDet_efftox_add(SEXP epSEXP, SEXP tpSEXP, SEXP rSEXP, SEXP pieSEXP, SEXP pitSEXP, SEXP pttSEXP, SEXP petSEXP, SEXP ptcSEXP, SEXP pecSEXP, SEXP h_treatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ep(epSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp(tpSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type pie(pieSEXP);
    Rcpp::traits::input_parameter< double >::type pit(pitSEXP);
    Rcpp::traits::input_parameter< double >::type ptt(pttSEXP);
    Rcpp::traits::input_parameter< double >::type pet(petSEXP);
    Rcpp::traits::input_parameter< double >::type ptc(ptcSEXP);
    Rcpp::traits::input_parameter< double >::type pec(pecSEXP);
    Rcpp::traits::input_parameter< int >::type h_treated(h_treatedSEXP);
    rcpp_result_gen = Rcpp::wrap(efftox_add(ep, tp, r, pie, pit, ptt, pet, ptc, pec, h_treated));
    return rcpp_result_gen;
END_RCPP
}
// prospect_powC
int prospect_powC(NumericMatrix ep, double er, double el, double eag, double eal, NumericMatrix tp, double tr, double tl, double tag, double tal, double k1, double k2, double ustop, int h_treated);
RcppExport SEXP _doseDet_prospect_powC(SEXP epSEXP, SEXP erSEXP, SEXP elSEXP, SEXP eagSEXP, SEXP ealSEXP, SEXP tpSEXP, SEXP trSEXP, SEXP tlSEXP, SEXP tagSEXP, SEXP talSEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP ustopSEXP, SEXP h_treatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ep(epSEXP);
    Rcpp::traits::input_parameter< double >::type er(erSEXP);
    Rcpp::traits::input_parameter< double >::type el(elSEXP);
    Rcpp::traits::input_parameter< double >::type eag(eagSEXP);
    Rcpp::traits::input_parameter< double >::type eal(ealSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp(tpSEXP);
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< double >::type tl(tlSEXP);
    Rcpp::traits::input_parameter< double >::type tag(tagSEXP);
    Rcpp::traits::input_parameter< double >::type tal(talSEXP);
    Rcpp::traits::input_parameter< double >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< double >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< double >::type ustop(ustopSEXP);
    Rcpp::traits::input_parameter< int >::type h_treated(h_treatedSEXP);
    rcpp_result_gen = Rcpp::wrap(prospect_powC(ep, er, el, eag, eal, tp, tr, tl, tag, tal, k1, k2, ustop, h_treated));
    return rcpp_result_gen;
END_RCPP
}
// prospect_pow_addC
int prospect_pow_addC(NumericMatrix ep, double er, double el, double eag, double eal, NumericMatrix tp, double tr, double tl, double tag, double tal, double k1, double k2, double ptt, double pet, double ptc, double pec, int h_treated);
RcppExport SEXP _doseDet_prospect_pow_addC(SEXP epSEXP, SEXP erSEXP, SEXP elSEXP, SEXP eagSEXP, SEXP ealSEXP, SEXP tpSEXP, SEXP trSEXP, SEXP tlSEXP, SEXP tagSEXP, SEXP talSEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP pttSEXP, SEXP petSEXP, SEXP ptcSEXP, SEXP pecSEXP, SEXP h_treatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ep(epSEXP);
    Rcpp::traits::input_parameter< double >::type er(erSEXP);
    Rcpp::traits::input_parameter< double >::type el(elSEXP);
    Rcpp::traits::input_parameter< double >::type eag(eagSEXP);
    Rcpp::traits::input_parameter< double >::type eal(ealSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp(tpSEXP);
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< double >::type tl(tlSEXP);
    Rcpp::traits::input_parameter< double >::type tag(tagSEXP);
    Rcpp::traits::input_parameter< double >::type tal(talSEXP);
    Rcpp::traits::input_parameter< double >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< double >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< double >::type ptt(pttSEXP);
    Rcpp::traits::input_parameter< double >::type pet(petSEXP);
    Rcpp::traits::input_parameter< double >::type ptc(ptcSEXP);
    Rcpp::traits::input_parameter< double >::type pec(pecSEXP);
    Rcpp::traits::input_parameter< int >::type h_treated(h_treatedSEXP);
    rcpp_result_gen = Rcpp::wrap(prospect_pow_addC(ep, er, el, eag, eal, tp, tr, tl, tag, tal, k1, k2, ptt, pet, ptc, pec, h_treated));
    return rcpp_result_gen;
END_RCPP
}
// prospect_pow_add_u_C
int prospect_pow_add_u_C(NumericMatrix ep, double er, double el, double eag, double eal, NumericMatrix tp, double tr, double tl, double tag, double tal, double k1, double k2, double ustop, double pstop, int h_treated);
RcppExport SEXP _doseDet_prospect_pow_add_u_C(SEXP epSEXP, SEXP erSEXP, SEXP elSEXP, SEXP eagSEXP, SEXP ealSEXP, SEXP tpSEXP, SEXP trSEXP, SEXP tlSEXP, SEXP tagSEXP, SEXP talSEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP ustopSEXP, SEXP pstopSEXP, SEXP h_treatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ep(epSEXP);
    Rcpp::traits::input_parameter< double >::type er(erSEXP);
    Rcpp::traits::input_parameter< double >::type el(elSEXP);
    Rcpp::traits::input_parameter< double >::type eag(eagSEXP);
    Rcpp::traits::input_parameter< double >::type eal(ealSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp(tpSEXP);
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< double >::type tl(tlSEXP);
    Rcpp::traits::input_parameter< double >::type tag(tagSEXP);
    Rcpp::traits::input_parameter< double >::type tal(talSEXP);
    Rcpp::traits::input_parameter< double >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< double >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< double >::type ustop(ustopSEXP);
    Rcpp::traits::input_parameter< double >::type pstop(pstopSEXP);
    Rcpp::traits::input_parameter< int >::type h_treated(h_treatedSEXP);
    rcpp_result_gen = Rcpp::wrap(prospect_pow_add_u_C(ep, er, el, eag, eal, tp, tr, tl, tag, tal, k1, k2, ustop, pstop, h_treated));
    return rcpp_result_gen;
END_RCPP
}
// prospect_pow_u_C
int prospect_pow_u_C(NumericMatrix ep, double er, double el, double eag, double eal, NumericMatrix tp, double tr, double tl, double tag, double tal, double k1, double k2, double ustop, double pstop, int h_treated);
RcppExport SEXP _doseDet_prospect_pow_u_C(SEXP epSEXP, SEXP erSEXP, SEXP elSEXP, SEXP eagSEXP, SEXP ealSEXP, SEXP tpSEXP, SEXP trSEXP, SEXP tlSEXP, SEXP tagSEXP, SEXP talSEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP ustopSEXP, SEXP pstopSEXP, SEXP h_treatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ep(epSEXP);
    Rcpp::traits::input_parameter< double >::type er(erSEXP);
    Rcpp::traits::input_parameter< double >::type el(elSEXP);
    Rcpp::traits::input_parameter< double >::type eag(eagSEXP);
    Rcpp::traits::input_parameter< double >::type eal(ealSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp(tpSEXP);
    Rcpp::traits::input_parameter< double >::type tr(trSEXP);
    Rcpp::traits::input_parameter< double >::type tl(tlSEXP);
    Rcpp::traits::input_parameter< double >::type tag(tagSEXP);
    Rcpp::traits::input_parameter< double >::type tal(talSEXP);
    Rcpp::traits::input_parameter< double >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< double >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< double >::type ustop(ustopSEXP);
    Rcpp::traits::input_parameter< double >::type pstop(pstopSEXP);
    Rcpp::traits::input_parameter< int >::type h_treated(h_treatedSEXP);
    rcpp_result_gen = Rcpp::wrap(prospect_pow_u_C(ep, er, el, eag, eal, tp, tr, tl, tag, tal, k1, k2, ustop, pstop, h_treated));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _doseDet_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_doseDet_UdetC", (DL_FUNC) &_doseDet_UdetC, 12},
    {"_doseDet_efftox_add", (DL_FUNC) &_doseDet_efftox_add, 10},
    {"_doseDet_prospect_powC", (DL_FUNC) &_doseDet_prospect_powC, 14},
    {"_doseDet_prospect_pow_addC", (DL_FUNC) &_doseDet_prospect_pow_addC, 17},
    {"_doseDet_prospect_pow_add_u_C", (DL_FUNC) &_doseDet_prospect_pow_add_u_C, 15},
    {"_doseDet_prospect_pow_u_C", (DL_FUNC) &_doseDet_prospect_pow_u_C, 15},
    {"_doseDet_rcpp_hello_world", (DL_FUNC) &_doseDet_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_doseDet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}