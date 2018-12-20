/* Copyright (C) 1997 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* check.c */

/* CheckCurvOpen()
*/

#include "editor.h"

int CheckCurvOpen(ParCurv *egeom)
{
  int i, open = 1;

  for (i=1; i<egeom->order && open; i++)          /* check for starting */
    if (egeom->knots[i] != egeom->knots[i-1])     /* knot multiplicity */
      open = 0;

  for (i=egeom->ncontpts+1; i<egeom->order+egeom->ncontpts && open; i++)
    if (egeom->knots[i] != egeom->knots[i-1])     /* check end multiplicity */
      open = 0;

  return open;
}

int CheckSurfOpen(ParSurf *fgeom, int *uOpen, int *vOpen)
{
  int i, open = 1;

  *uOpen = 1;
  for (i=1; i<fgeom->uorder && *uOpen; i++)       /* check for starting */
    if (fgeom->uknots[i] != fgeom->uknots[i-1]) { /* u knot multiplicity */
      open   = 0;
      *uOpen = 0;
    }

  for (i=fgeom->ucontpts+1; i<fgeom->uorder+fgeom->ucontpts && *uOpen; i++)
    if (fgeom->uknots[i] != fgeom->uknots[i-1]) { /* check end multiplicity */
      open   = 0;
      *uOpen = 0;
    }

  *vOpen = 1;
  for (i=1; i<fgeom->vorder && *vOpen; i++)       /* check for starting */
    if (fgeom->vknots[i] != fgeom->vknots[i-1]) { /* v knot multiplicity */
      open   = 0;
      *vOpen = 0;
    }

  for (i=fgeom->vcontpts+1; i<fgeom->vorder+fgeom->vcontpts && *vOpen; i++)
    if (fgeom->vknots[i] != fgeom->vknots[i-1]) { /* check end multiplicity */
      open   = 0;
      *vOpen = 0;
    }

  return open;
}
