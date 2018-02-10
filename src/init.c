#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gamlr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gamlr_cleanup();

static const R_CMethodDef CEntries[] = {
    {"gamlr",         (DL_FUNC) &gamlr,         32},
    {"gamlr_cleanup", (DL_FUNC) &gamlr_cleanup,  0},
    {NULL, NULL, 0}
};

void R_init_gamlr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}