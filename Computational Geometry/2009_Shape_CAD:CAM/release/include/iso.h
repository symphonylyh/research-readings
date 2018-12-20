/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* iso.h */

#ifndef ISO_H
#define ISO_H

#ifdef __cplusplus
extern "C" {
#endif

void	CalcIsophotes(ParSurf *, short, short, double ***, double *);
void	DrawIsophotes(ParSurf *, double *, short, short, short,
		double ***, FILE *);
double	fcn(ParSurf *, double *, double, double, double);
int	main(unsigned, char **);
double	SolvIsophotes(ParSurf *, double *, double, double, double,
		 double, short);

#ifdef __cplusplus
}
#endif

#endif
