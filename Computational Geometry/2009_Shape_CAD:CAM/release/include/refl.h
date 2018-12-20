/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* refl.h */

#ifndef REFL_H
#define REFL_H

#ifdef __cplusplus
extern "C" {
#endif

void	CalcReflectionLines(ParSurf *, short, short, double ***, double *);
void	DrawReflectionLines(ParSurf *, double *, short, short, short,
		double ***, FILE *);
double	fcn(ParSurf *, double *, double, double, double);
int	main(unsigned argc, char **);
double	slv(ParSurf *, double *, double, double, double, double, short);

#ifdef __cplusplus
}
#endif

#endif
