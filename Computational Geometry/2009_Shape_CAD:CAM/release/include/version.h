/* Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA */
/* All rights reserved */

/* version.h */

#ifndef VERSION_H
#define VERSION_H

#ifdef __cplusplus
extern "C" {
#endif

char	*GetAppDate(void);
char	*GetAppName(void);
char	*GetAppVers(void);
char	*GetCopyrightYear(void);
char    *GetLinkTimestamp(void);

#ifdef __cplusplus
}
#endif

#endif
