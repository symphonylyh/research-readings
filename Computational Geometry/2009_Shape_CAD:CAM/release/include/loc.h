/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* loc.h */

#ifndef LOC_H
#define LOC_H

#ifdef __cplusplus
extern "C" {
#endif

void	clear_vector(vector *);
void	DrawLinesOfCurvature(ParSurf *, short, short, double, FILE *);
void	find_dir(ParSurf *, vector *, vector *, double *, double *);
void	fn(ParSurf *, double, double *, double *);
int	main(unsigned, char *[]);
void	rk4(ParSurf *, double *, double *, short, double, double,
		double *, void (*)(ParSurf *, double, double *, double *));
void	trn(ParSurf *, double *, double *);

#ifdef __cplusplus
}
#endif

#endif
