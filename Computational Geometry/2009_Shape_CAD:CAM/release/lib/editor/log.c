/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* log.c */

/* CloseLog()
   GetLogDebug()
   GetLogKey()
   OpenLog()
   PrintLog()
   SetLogDebug()
   SetLogKey()
*/

#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include "editor.h"

int cftime(char *, char *, time_t *);

static FILE *fp = NULL;
static short debug = 0, key = LOG_FILE_STDOUT;
static int tps;

void OpenLog(char *logFile, char *type)   /* open log file */
{
  fp = fopen(logFile, type);

  tps = sysconf(_SC_CLK_TCK);
}

void SetLogKey(short k)       /* set print key */
{
  key = k;
}

void PrintLog(char *msg)      /* print to log file and standard output */
{
  struct tms buffer;
  time_t clck;
  long hh, mm, ss;
  char cpu[80], line[256], local[10];

  if (debug) {
    time(&clck);
    cftime(local, "%T", &clck);
    if (debug == 1)
      sprintf(line, "%s %s", local, msg);
    else {
      times(&buffer);                     /* get clock ticks */
      ss = buffer.tms_utime/tps;         /* convert to total CPU seconds */
      hh = ss/3600;                       /* get CPU hours */
      ss -= (hh*3600);
      mm = ss/60;                         /* get remaining CPU minutes */
      ss -= (mm*60);                      /* get remaining CPU seconds */
      sprintf(cpu, "%02d:%02d:%02d", hh, mm, ss);
      sprintf(line, "%s %s %s", local, cpu, msg);
    }
  }
  else
    sprintf(line, "%s", msg);

  if (key == LOG_FILE || key == LOG_FILE_STDOUT)
    if (fp) {
      fprintf(fp,"%s", line);  /* write to log file */
      fflush(fp);
    }
  if (key == LOG_STDOUT || key == LOG_FILE_STDOUT) {
    printf("%s", line);        /* write to standard output */
    fflush(stdout);
  }
}

void CloseLog(void)           /* close log file */
{
  fclose(fp);
}

short GetLogKey(void)         /* return print key */
{
  return (key);
}

void SetLogDebug(int dbg)
{
  debug = dbg;
}

int GetLogDebug(void)
{
  return debug;
}
