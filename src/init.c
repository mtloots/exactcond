#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern void C_ec_rc_null(const int*, const int*, const int*, const int*, const int*, const int*,
                         const double*, const double*, const int*, double*);
extern void C_ec_rc_test(const int*, const int*, const int*, const int*, const int*,
                         const double*, const double*, const int*, double*);
extern void C_ec_rc_two(const int*, const int*, const int*, const int*, const int*, const int*,
                        const int*, const double*, const double*, const double*, const int*,
                        double*);
extern void C_ec_sphere_maxquad(const double*, const double*, const int*, const int*,
                                const int*, const int*, double*);
extern void C_ec_pvalue(const double*, const double*, const int*, double*);

static const R_CMethodDef CEntries[] = {
    {"C_ec_rc_null",        (DL_FUNC) &C_ec_rc_null,        10},
    {"C_ec_rc_test",        (DL_FUNC) &C_ec_rc_test,         9},
    {"C_ec_rc_two",         (DL_FUNC) &C_ec_rc_two,         12},
    {"C_ec_sphere_maxquad", (DL_FUNC) &C_ec_sphere_maxquad,  7},
    {"C_ec_pvalue",         (DL_FUNC) &C_ec_pvalue,          4},
    {NULL, NULL, 0}
};

void R_init_exactcond(DllInfo *dll){
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
