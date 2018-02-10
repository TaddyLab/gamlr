#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void gamlr_cleanup();
extern void gamlr_inc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"gamlr_cleanup", (DL_FUNC) &gamlr_cleanup,  0},
    {"gamlr_inc",     (DL_FUNC) &gamlr_inc,     32},
    {NULL, NULL, 0}
};

void R_init_gamlr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}