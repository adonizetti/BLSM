#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP _BLSM_alpha_up(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BLSM_dst(SEXP);
extern SEXP _BLSM_lpY(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BLSM_lpYNODE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BLSM_lpz_dist(SEXP);
extern SEXP _BLSM_lpz_distNODE(SEXP, SEXP, SEXP);
extern SEXP _BLSM_mlpY(SEXP, SEXP, SEXP);
extern SEXP _BLSM_Z_up(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_BLSM_alpha_up",     (DL_FUNC) &_BLSM_alpha_up,     7},
  {"_BLSM_dst",          (DL_FUNC) &_BLSM_dst,          1},
  {"_BLSM_lpY",          (DL_FUNC) &_BLSM_lpY,          4},
  {"_BLSM_lpYNODE",      (DL_FUNC) &_BLSM_lpYNODE,      6},
  {"_BLSM_lpz_dist",     (DL_FUNC) &_BLSM_lpz_dist,     1},
  {"_BLSM_lpz_distNODE", (DL_FUNC) &_BLSM_lpz_distNODE, 3},
  {"_BLSM_mlpY",         (DL_FUNC) &_BLSM_mlpY,         3},
  {"_BLSM_Z_up",         (DL_FUNC) &_BLSM_Z_up,         7},
  {NULL, NULL, 0}
};

void R_init_BLSM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}