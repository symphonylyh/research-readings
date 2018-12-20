/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* geod.h */

#ifndef GEOD_H
#define GEOD_H

#ifdef __cplusplus
extern "C" {
#endif

void	DrawGeodesics(ParSurf *, short, short, double, FILE *);
void	find_init(ParSurf *, double *);
void	fn(ParSurf *, double, double *, double *);
int	main(unsigned, char **);
void	rk4(ParSurf *, double *, double *, short, double, double, double *,
		void (*)(ParSurf *, double, double *, double *));
void	trn(double *, double *);

#ifdef __cplusplus
extern "C" {
#endif

#endif
