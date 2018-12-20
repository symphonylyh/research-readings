/* Copyright (C) 1994 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* countdown.h
*/

#ifndef PRAX_COUNTDOWN_H
#define PRAX_COUNTDOWN_H

#include <Xm/Xm.h>

#ifdef __cplusplus
extern "C" {
#endif

void	Countdown(float);
void	CountdownExposeCB(Widget, XtPointer, XtPointer);
void	CountdownGinitCB(Widget, XtPointer, XtPointer);
void	Counter(short);
void	CounterExposeCB(Widget, XtPointer, XtPointer);
void	CounterGinitCB(Widget, XtPointer, XtPointer);
void	EndCountdown(void);
void	EndCounter(void);
void	StartCountdown(Widget, char *);
void	StartCounter(Widget, char *, short, short);
void	UpdateCountdown(float);
void	UpdateCounter(short);

#ifdef __cplusplus
}
#endif

#endif
