#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

char *tempnam(const char *, const char *);

static short nBot, nLft, nRht, nTop;

int FindBot(double **pts, int nPts, double ***edge, int goal)
{
  FILE *fp;
  double **temp;
  short i;
  char *tmpfil;

  tmpfil = tempnam("/usr/tmp", "cont");
  if (fp = fopen(tmpfil, "w")) {
    nBot = 0;
    SearchBot(fp, pts, nPts, 0.0, 1.0, 0.0, 1.0, goal);
    fclose(fp);

    (*edge) = dbl_array2(nBot, 3);
    if (fp = fopen(tmpfil, "r")) {
      for (i=0; i<nBot; i++)
	fscanf(fp, "%lf %lf %lf", &((*edge)[i][0]), &((*edge)[i][1]),
	       &((*edge)[i][2]));
      fclose(fp);
    }
    unlink(tmpfil);
    free(tmpfil);
  }

  return nBot;
}

int FindLft(double **pts, int nPts, double ***edge, int goal)
{
  FILE *fp;
  double **temp;
  short i;
  char *tmpfil;

  tmpfil = tempnam("/usr/tmp", "cont");
  if (fp = fopen(tmpfil, "w")) {
    nLft = 0;
    SearchLft(fp, pts, nPts, 0.0, 1.0, 0.0, 1.0, goal);
    fclose(fp);

    (*edge) = dbl_array2(nLft, 3);
    if (fp = fopen(tmpfil, "r")) {
      for (i=0; i<nLft; i++)
	fscanf(fp, "%lf %lf %lf", &((*edge)[i][0]), &((*edge)[i][1]),
	       &((*edge)[i][2]));
      fclose(fp);
    }
    unlink(tmpfil);
    free(tmpfil);
  }

  return nLft;
}

int FindRht(double **pts, int nPts, double ***edge, int goal)
{
  FILE *fp;
  double **temp;
  short i;
  char *tmpfil;

  tmpfil = tempnam("/usr/tmp", "cont");
  if (fp = fopen(tmpfil, "w")) {
    nRht = 0;
    SearchRht(fp, pts, nPts, 0.0, 1.0, 0.0, 1.0, goal);
    fclose(fp);

    (*edge) = dbl_array2(nRht, 3);
    if (fp = fopen(tmpfil, "r")) {
      for (i=0; i<nRht; i++)
	fscanf(fp, "%lf %lf %lf", &((*edge)[i][0]), &((*edge)[i][1]),
	       &((*edge)[i][2]));
      fclose(fp);
    }
    unlink(tmpfil);
    free(tmpfil);
  }

  return nRht;
}

int FindTop(double **pts, int nPts, double ***edge, int goal)
{
  FILE *fp;
  double **temp;
  short i;
  char *tmpfil;

  tmpfil = tempnam("/usr/tmp", "cont");
  if (fp = fopen(tmpfil, "w")) {
    nTop = 0;
    SearchTop(fp, pts, nPts, 0.0, 1.0, 0.0, 1.0, goal);
    fclose(fp);

    (*edge) = dbl_array2(nTop, 3);
    if (fp = fopen(tmpfil, "r")) {
      for (i=0; i<nTop; i++)
	fscanf(fp, "%lf %lf %lf", &((*edge)[i][0]), &((*edge)[i][1]),
	       &((*edge)[i][2]));
      fclose(fp);
    }
    unlink(tmpfil);
    free(tmpfil);
  }

  return nTop;
}

void SearchTop(FILE *fp, double **pts, int nPts, double umin, double umax,
	       double vmin, double vmax, int goal)
{
  double umid, vmid;
  short i, nin, nout;

  vmid = (vmin + vmax)/2.0;

  for (i=nin=nout=0; i<nPts; i++)
    if (umin <= pts[i][0] && pts[i][0] <= umax && vmin <= pts[i][1]) {
      if (vmid <= pts[i][1] && pts[i][1] <= vmax)
	nout++;
      else
	nin++;
    }

  if (nin || nout) {
    umid = (umin + umax)/2.0;
    if (nout) {
      SearchTop(fp, pts, nPts, umin, umid, vmid, vmax, goal);
      SearchTop(fp, pts, nPts, umid, umax, vmid, vmax, goal);
    }
    else if (nin > goal) {
      SearchTop(fp, pts, nPts, umin, umid, vmin, vmid, goal);
      SearchTop(fp, pts, nPts, umid, umax, vmin, vmid, goal);
    }
    else {
      for (i=0; i<nPts; i++)
	if (umin <= pts[i][0] && pts[i][0] <= umax &&
	    vmin <= pts[i][1] && pts[i][1] <= vmid) {
	  fprintf(fp, "%f %f %f\n", pts[i][0], pts[i][1], pts[i][2]);
	  nTop++;
	}
    }
  }
}

void SearchLft(FILE *fp, double **pts, int nPts, double umin, double umax,
	       double vmin, double vmax, int goal)
{
  double umid, vmid;
  short i, nin, nout;

  umid = (umin + umax)/2.0;

  for (i=nin=nout=0; i<nPts; i++)
    if (pts[i][0] <= umax && vmin <= pts[i][1] && pts[i][1] <= vmax) {
      if (umin <= pts[i][0] && pts[i][0] <= umid)
	nout++;
      else
	nin++;
    }

  if (nin || nout) {
    vmid = (vmin + vmax)/2.0;
    if (nout) {
      SearchLft(fp, pts, nPts, umin, umid, vmin, vmid, goal);
      SearchLft(fp, pts, nPts, umin, umid, vmid, vmax, goal);
    }
    else if (nin > goal) {
      SearchLft(fp, pts, nPts, umid, umax, vmin, vmid, goal);
      SearchLft(fp, pts, nPts, umid, umax, vmid, vmax, goal);
    }
    else {
      for (i=0; i<nPts; i++)
	if (umid <= pts[i][0] && pts[i][0] <= umax &&
	    vmin <= pts[i][1] && pts[i][1] <= vmax) {
	  fprintf(fp, "%f %f %f\n", pts[i][0], pts[i][1], pts[i][2]);
	  nLft++;
	}
    }
  }
}

void SearchBot(FILE *fp, double **pts, int nPts, double umin, double umax,
	       double vmin, double vmax, int goal)
{
  double umid, vmid;
  short i, nin, nout;

  vmid = (vmin + vmax)/2.0;

  for (i=nin=nout=0; i<nPts; i++)
    if (umin <= pts[i][0] && pts[i][0] <= umax && pts[i][1] <= vmax) {
      if (vmin <= pts[i][1] && pts[i][1] <= vmid)
	nout++;
      else
	nin++;
    }

  if (nin || nout) {
    umid = (umin + umax)/2.0;
    if (nout) {
      SearchBot(fp, pts, nPts, umin, umid, vmin, vmid, goal);
      SearchBot(fp, pts, nPts, umid, umax, vmin, vmid, goal);
    }
    else if (nin > goal) {
      SearchBot(fp, pts, nPts, umin, umid, vmid, vmax, goal);
      SearchBot(fp, pts, nPts, umid, umax, vmid, vmax, goal);
    }
    else {
      for (i=0; i<nPts; i++)
	if (umin <= pts[i][0] && pts[i][0] <= umax &&
	    vmid <= pts[i][1] && pts[i][1] <= vmax) {
	  fprintf(fp, "%f %f %f\n", pts[i][0], pts[i][1], pts[i][2]);
	  nBot++;
	}
    }
  }
}

void SearchRht(FILE *fp, double **pts, int nPts, double umin, double umax,
	       double vmin, double vmax, int goal)
{
  double umid, vmid;
  short i, nin, nout;

  umid = (umin + umax)/2.0;

  for (i=nin=nout=0; i<nPts; i++)
    if (umin <= pts[i][0] && vmin <= pts[i][1] && pts[i][1] <= vmax) {
      if (umid <= pts[i][0] && pts[i][0] <= umax)
	nout++;
      else
	nin++;
    }

  if (nin || nout) {
    vmid = (vmin + vmax)/2.0;
    if (nout) {
      SearchRht(fp, pts, nPts, umid, umax, vmin, vmid, goal);
      SearchRht(fp, pts, nPts, umid, umax, vmid, vmax, goal);
    }
    else if (nin > goal) {
      SearchRht(fp, pts, nPts, umin, umid, vmin, vmid, goal);
      SearchRht(fp, pts, nPts, umin, umid, vmid, vmax, goal);
    }
    else {
      for (i=0; i<nPts; i++)
	if (umin <= pts[i][0] && pts[i][0] <= umid &&
	    vmin <= pts[i][1] && pts[i][1] <= vmax) {
	  fprintf(fp, "%f %f %f\n", pts[i][0], pts[i][1], pts[i][2]);
	  nRht++;
	}
    }
  }
}
