/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* path.c */

/* GetFileName()
   GetPathName()
*/

#include <stdio.h>
#include <string.h>
 char *strdup(const char *);

char *GetFileName(char *file)
{
  short i;

  for (i = strlen(file)-1; i > -2; i--) {
    if (i < 0)
      break;
    if (file[i] == '/')
      break;
  }

  return (strdup(&file[i+1]));
}

char *GetPathName(char *file)
{
  char *path;
  short i;

  path = strdup(file);
  for (i = strlen(path)-1; i > -2; i--) {
    if (i < 0)
      break;
    if (path[i] == '/')
      break;
  }
  path[i] = '\0';

  return (path);
}
