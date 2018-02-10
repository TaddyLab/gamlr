#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void _gamlr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void _gamlr_cleanup();

static const R_CMethodDef CEntries[] = {
    {"_gamlr",         (DL_FUNC) &_gamlr,         32},
    {"_gamlr_cleanup", (DL_FUNC) &_gamlr_cleanup,  0},
    {NULL, NULL, 0}
};

void R_init_gamlr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}