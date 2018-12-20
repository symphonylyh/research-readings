/*
   PostScript library based on SGI GL library 
   S. L. Abrams 
   Copyright (C) 1992 Massachusetts Institute of Technology 
   All rights reserved 
*/

#ifndef PS_H
#define PS_H

#include <stdio.h>
/*
   Cho beg
   #include <gl/gl.h> 
*/

#include <GL/gl.h>
/* Cho end */

#include "gen.h"
#include "editor.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PS_XPMAX 0
#define PS_YPMAX 1

#define PS_LANDSCAPE 1
#define PS_PORTRAIT 0

#define PS_MVIEWING 0
#define PS_MPROJECTION 1

#define PS_DEG_TO_RAD 0.017453293

#define PS_BGNCLOSEDLINE 0
#define PS_BGNLINE 1
#define PS_BGNPOINT 2
#define PS_BGNPOLYGON 3
#define PS_CHARSTR 4
#define PS_CLOSEDLINE 5
#define PS_CMOV 6
#define PS_DRAW 7
#define PS_LINE 8
#define PS_MOVE 9
#define PS_NOOP 10
#define PS_POINT 11
#define PS_POLYGON 12

#define PS_INSIDE 0
#define PS_TO_LEFT 1
#define PS_TO_RIGHT 2
#define PS_BELOW 4
#define PS_ABOVE 8

#define PS_MAXV 50

#define PS_PTSIZE 10

#define PS_DASH 1
#define PS_SOLID 0

void ps_arc(Coord x, Coord y, Coord radius, Angle start, Angle end);
void ps_arci(Icoord x, Icoord y, Icoord radius, Angle start, Angle end);
void ps_arcs(Scoord x, Scoord y, Scoord radius, Angle start, Angle end);
void ps_bgnclosedline(void);
void ps_bgnline(void);
void ps_bgnpoint(void);
void ps_bgnpolygon(void);
void ps_bounding_box(int, int);
void ps_charstr(String);
void ps_charstr_center(String);
void ps_charstr_right(String);
void ps_circ(Coord, Coord, Coord);
void ps_circf(Coord, Coord, Coord);
void ps_circfi(Icoord, Icoord, Icoord);
void ps_circi(Icoord, Icoord, Icoord);
int ps_clip_line(int,  int, int, int, int *, int *, int *, int *);
void ps_close(void);
void ps_cmov(Coord, Coord, Coord);
void ps_cmov2(Coord, Coord);
void ps_cmov2i(Icoord, Icoord);
void ps_colorimage(int, int, unsigned long *);
void ps_comment(String);
void ps_cpack(unsigned long);
void ps_dec2hex(unsigned char);
void ps_draw(Coord, Coord, Coord);
void ps_draw2(Coord, Coord);
void ps_draw2i(Icoord, Icoord);
void ps_drawi(Icoord, Icoord, Icoord);
void ps_endclosedline(void);
void ps_endline(void);
void ps_endpoint(void);
void ps_endpolygon(void);
void ps_error(String, String);
void ps_gconfig(int);
int  ps_getcmode(void);
long ps_getgdesc(long);
void ps_getsize(long *, long *);
void ps_ginit(String, String);
void ps_GridSurf(GridSurf *);
void ps_header(String);
void ps_identity(void);
void ps_linewidth(short);
void ps_ListCurv(ListCurv *);
void ps_lookat(Coord, Coord, Coord, Coord, Coord, Coord, Angle);
void ps_move(Coord, Coord, Coord);
void ps_move2(Coord, Coord);
void ps_move2i(Icoord, Icoord);
void ps_movei(Icoord, Icoord, Icoord);
void ps_newpage(void);
void ps_newpath(void);
void ps_object_to_screen(Coord, Coord, Coord, int *, int *);
FILE *ps_open(String, String);
void ps_ortho(Coord, Coord, Coord, Coord, Coord, Coord);
void ps_ortho2(Coord, Coord, Coord, Coord);
int  ps_outcode(int, int);
void ps_ParCurv(ParCurv *, int);
void ps_ParSurf(ParSurf *, int, int);
void ps_ParSurf2(ParSurf *, int, int, int, int);
void ps_perspective(Angle, float, Coord, Coord);
void ps_pnt(Coord, Coord, Coord);
void ps_pnt2(Coord, Coord);
void ps_pnt2i(Icoord, Icoord);
void ps_pnti(Icoord, Icoord, Icoord);
void ps_polarview(Coord, Angle, Angle, Angle);
void ps_prologue(void);
void ps_praxiteles(void);
void ps_rotate(Angle, char);
void ps_rotate_lookat(double, double, char);
void ps_scale(float x, float, float);
void ps_setcmode(void);
void ps_setgmode(void);
void ps_setgray(float);
void ps_setlinestyle(short);
void ps_setrgbcolor(float, float, float);
void ps_showpage(void);
void ps_stroke(void);
void ps_trailer(void);
void ps_translate(Coord, Coord, Coord);
void ps_v2d(double *);
void ps_v2f(float *);
void ps_v2i(long *);
void ps_v2s(short *);
void ps_v3d(double *);
void ps_v3f(float *);
void ps_v3i(long *);
void ps_v3s(short *);
void ps_viewport(Screencoord, Screencoord, Screencoord, Screencoord);
void ps_window(Coord, Coord, Coord, Coord, Coord, Coord);

#ifdef __cplusplus
}
#endif

#endif
