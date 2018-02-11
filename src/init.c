#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void R_gamlr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void R_gamlr_cleanup();

static const R_CMethodDef CEntries[] = {
    {"R_gamlr",         (DL_FUNC) &R_gamlr,         32},
    {"R_gamlr_cleanup", (DL_FUNC) &R_gamlr_cleanup,  0},
    {NULL, NULL, 0}
};

void R_init_gamlr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
