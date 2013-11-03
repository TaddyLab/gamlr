// IO tools

#include <stdarg.h>
#ifdef RPRINT
#include <R.h>
#endif
#include "gui.h"

// print and warning/error wrappers

void speak(const char *str, ...)
{
  va_list argp;
  va_start(argp, str);

#ifdef RPRINT
  Rvprintf(str, argp);
#else
  vprintf(str, argp);
#endif

  va_end(argp);
}

void shout(const char *str, ...)
{
  va_list argp;
  va_start(argp, str);

#ifdef RPRINT
  REvprintf(str, argp);
#else
  vfprintf(stderr, str, argp);
#endif

  va_end(argp);
}

// pass control up to check for interrupts, etc.

time_t interact(time_t itime)
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
