#ifndef __RHELP_H__
#define __RHELP_H__

#include <stdio.h>
#include <time.h>

/* this is now covered by -D RPRINT flags in Makevars */
/*#define RPRINT*/

extern FILE *mystdout, *mystderr;

#ifndef RPRINT
void warning(const char *str, ...);
void error(const char *str, ...);
#else
#include <R_ext/Utils.h>
#include <R.h>
#endif

void R_FlushConsole(void); /* R < 2.3 does not have this in R.h (in Rinterface.h) */
void myprintf(FILE *outfile, const char *str, ...);
void myflush(FILE *outfile);
time_t my_r_process_events(time_t itime);
void myassert(int cnd);

#endif
