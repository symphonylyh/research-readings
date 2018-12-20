/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* MakeFileName()
 */

#include <ctype.h>
#include <string.h>
#include "editor.h"

void MakeFileName(const char *pattern, const char *ext, char *name)
{
  int i, j;

  for (i=j=0; pattern[i]; i++) {
    if (pattern[i] == ' ')          /* replace blank with underscore */
      name[j++] = '_';
    else if (isalnum(pattern[i]))   /* only copy alphanumerics */
      name[j++] = pattern[i];
  }
  if (!j)
    strcpy(name, "file");           /* default file name */
  else
    name[j] = '\0';

  strcat(name, ".");                /* separate name from extension */
  strcat(name, ext);                /* add file extension */
}
