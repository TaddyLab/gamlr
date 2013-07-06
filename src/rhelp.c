/********************************************************************************
 *
 * Adapted from R.B. Gramacy's:
 * Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
 * Copyright (C) 2005, University of California
 *
 ********************************************************************************/

#include "rhelp.h"
#ifdef RPRINT
#include <R_ext/Print.h>
#include <R.h>
FILE *mystdout = (FILE*) 0;
FILE *mystderr = (FILE*) 1;
#endif
#include <stdarg.h>
#include <time.h>
#include <assert.h>

/*
 * myprintf:
 *
 * a function many different types of printing-- in particular, using
 * the Rprintf if the code happens to be compiled with RPRINT,
 * othersie fprintf (takes the same arguments as fprintf)
 */

void myprintf(FILE *outfile, const char *str, ...)
{
  va_list argp;
  va_start(argp, str);

  #ifdef RPRINT
  if(outfile == mystdout) Rvprintf(str, argp);
  else if(outfile == mystderr) REvprintf(str, argp);
  else vfprintf(outfile, str, argp);
  #else
  vfprintf(outfile, str, argp);
  #endif

  va_end(argp);
}


#ifndef RPRINT
/*
 * error:
 *
 * printf style function that reports errors to stderr
 */

void error(const char *str, ...)
{
  va_list argp;
  va_start(argp, str);

  myprintf(stderr, "ERROR: ");
  vfprintf(stderr, str, argp);

  va_end(argp);
  myflush(stderr);

  /* add a final newline */
  myprintf(stderr, "\n");

  /* kill the code */
  assert(0);
}


/*
 * warning:
 *
 * printf style function that reports warnings to stderr
 */

void warning(const char *str, ...)
{
  va_list argp;
  va_start(argp, str);

  myprintf(stderr, "WARNING: ");
  vfprintf(stderr, str, argp);

  va_end(argp);
  myflush(stderr);

  /* add a final newline */
  myprintf(stderr, "\n");
}
#endif

/*
 * myflush:
 *
 * a function for many different types of flushing--  in particular,
 * using * the R_FlushConsole the code happens to be compiled with
 * RPRINT, otherwise fflush
 */

void myflush(FILE *outfile)
{
#ifdef RPRINT
  R_FlushConsole();
#else
  fflush(outfile);
#endif
}


/*
 * my_r_process_events:
 *
 * at least every 1 second(s) pass control back to
 * R so that it can check for interrupts and/or
 * process other R-gui events
 */

time_t my_r_process_events(time_t itime)
{
#ifdef RPRINT
  time_t ntime = time(NULL);

  if(ntime - itime > 1) {
    R_FlushConsole();
    R_CheckUserInterrupt();
#if  (defined(HAVE_AQUA) || defined(Win32) || defined(Win64))
    R_ProcessEvents();
#endif
    itime = ntime;
  }
#endif
  return itime;
}

/* assert that uses proper R errors under RPRINT */

void myassert(int cnd){
  if(!cnd)
	error("failed assertion."); }
