/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* prax_cursor.h */

#ifndef PRAX_CURSOR_H
#define PRAX_CURSOS_H

#ifdef __cplusplus
extern "C" {
#endif

#define PRAX_CURSOR_WIDTH 16
#define PRAX_CURSOR_HEIGHT 16
#define PRAX_CURSOR_X_HOT 7
#define PRAX_CURSOR_Y_HOT 7

/* RotX cursor */

static char rotx_bits[] = {
   0x80, 0x01, 0xc0, 0x03, 0xe0, 0x07, 0xf0, 0x0f, 0x80, 0x01, 0x80, 0x01,
   0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x80, 0x01, 0x80, 0x01,
   0xf0, 0x0f, 0xe0, 0x07, 0xc0, 0x03, 0x80, 0x01};
static char rotx_mask[] = {
   0xc0, 0x03, 0xe0, 0x07, 0xf0, 0x0f, 0xf8, 0x1f, 0xf8, 0x1f, 0xc0, 0x03,
   0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xc0, 0x03, 0xf8, 0x1f,
   0xf8, 0x1f, 0xf0, 0x0f, 0xe0, 0x07, 0xc0, 0x03};

/* RotY cursor */

static char roty_bits[] = {
   0x80, 0x01, 0x80, 0x01, 0x80, 0x01, 0x80, 0x01, 0x88, 0x11, 0x8c, 0x31,
   0x8e, 0x71, 0xbf, 0xfd, 0xbf, 0xfd, 0x8e, 0x71, 0x8c, 0x31, 0x88, 0x11,
   0x80, 0x01, 0x80, 0x01, 0x80, 0x01, 0x80, 0x01};
static char roty_mask[] = {
   0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xd8, 0x1b, 0xdc, 0x3b, 0xde, 0x7b,
   0xff, 0xfb, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xde, 0x7b, 0xdc, 0x3b,
   0xd8, 0x1b, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03};

/* RotZ cursor */

static char rotz_bits[] = {
   0x00, 0x00, 0x80, 0x07, 0x60, 0x18, 0x10, 0x20, 0x08, 0x40, 0x7f, 0x40,
   0x3e, 0x80, 0x3e, 0x80, 0x1c, 0x80, 0x1c, 0x80, 0x08, 0x40, 0x00, 0x40,
   0x00, 0x20, 0x40, 0x18, 0x80, 0x07, 0x00, 0x00};
static char rotz_mask[] = {
   0x80, 0x07, 0xe0, 0x1f, 0xf0, 0x3f, 0x3c, 0x78, 0xff, 0xe0, 0xff, 0xe0,
   0x7f, 0xc0, 0x7f, 0xc0, 0x3e, 0xc0, 0x3e, 0xc0, 0x1c, 0xe0, 0x08, 0xe0,
   0x40, 0x70, 0xe0, 0x3c, 0xc0, 0x1f, 0x80, 0x07};

/* TrXy cursor */

static char trxy_bits[] = {
   0x80, 0x01, 0xc0, 0x03, 0xe0, 0x07, 0x80, 0x01, 0x80, 0x01, 0x84, 0x21,
   0x86, 0x61, 0xff, 0xff, 0xff, 0xff, 0x86, 0x61, 0x84, 0x21, 0x80, 0x01,
   0x80, 0x01, 0xe0, 0x07, 0xc0, 0x03, 0x80, 0x01};
static char trxy_mask[] = {
   0xc0, 0x03, 0xe0, 0x07, 0xf0, 0x0f, 0xf0, 0x0f, 0xcc, 0x33, 0xce, 0x73,
   0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xf3, 0xce, 0x73, 0xcc, 0x33,
   0xf0, 0x0f, 0xf0, 0x0f, 0xe0, 0x07, 0xc0, 0x03};

/* TrZ cursor */

static char trz_bits[] = {
   0x80, 0x01, 0xc0, 0x03, 0xf0, 0x0f, 0x80, 0x01, 0x80, 0x01, 0xc0, 0x03,
   0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xe0, 0x07, 0xfe, 0x7f, 0xfc, 0x3f,
   0xf0, 0x0f, 0xe0, 0x07, 0xc0, 0x03, 0x80, 0x01};
static char trz_mask[] = {
   0xc0, 0x03, 0xe0, 0x0f, 0xf8, 0x1f, 0xfc, 0x3f, 0xc0, 0x03, 0xe0, 0x07,
   0xe0, 0x07, 0xe0, 0x07, 0xe0, 0x07, 0xff, 0xff, 0xff, 0xff, 0xfe, 0x7f,
   0xfc, 0x3f, 0xf0, 0x0f, 0xe0, 0x07, 0xc0, 0x03};

/* Dist cursor */

static char dist_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x01, 0x80, 0x01, 0x80, 0x11, 0x88, 0x19, 0x98,
   0x1d, 0xb8, 0xff, 0xff, 0x1d, 0xb8, 0x19, 0x98, 0x11, 0x88, 0x01, 0x80,
   0x01, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
static char dist_mask[] = {
   0x00, 0x00, 0x03, 0xc0, 0x03, 0xc0, 0x33, 0xcc, 0x3b, 0xdc, 0x3f, 0xfc,
   0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x3f, 0xfc, 0x3b, 0xdc, 0x33, 0xcc,
   0x03, 0xc0, 0x03, 0xc0, 0x00, 0x00, 0x00, 0x00};

/* Fov cursor */

static char fov_bits[] = {
   0x00, 0x30, 0x00, 0x0c, 0x00, 0x03, 0x80, 0x1c, 0x60, 0x0c, 0x18, 0x14,
   0x04, 0x20, 0x02, 0x20, 0x04, 0x20, 0x18, 0x14, 0x60, 0x0c, 0x80, 0x1c,
   0x00, 0x03, 0x00, 0x0c, 0x00, 0x30, 0x00, 0x00};
static char fov_mask[] = {
   0x00, 0x3c, 0x00, 0x3f, 0x80, 0x0f, 0xe0, 0x1f, 0xf8, 0x3f, 0x7c, 0x3e,
   0x1e, 0x7e, 0x07, 0x70, 0x1e, 0x70, 0x7c, 0x7e, 0xf8, 0x3e, 0xe0, 0x3f,
   0x80, 0x1f, 0x00, 0x0f, 0x00, 0x3c, 0x00, 0x30};

/* Near cursor */

static char near_bits[] = {
   0xf0, 0xff, 0x10, 0x80, 0x00, 0x80, 0xff, 0x9f, 0xff, 0x9f, 0x03, 0x98,
   0x13, 0x98, 0x13, 0x98, 0x13, 0x98, 0x13, 0x98, 0x13, 0x98, 0xf3, 0xdb,
   0x03, 0x18, 0x03, 0x18, 0xff, 0x1f, 0xff, 0x1f};
static char near_mask[] = {
   0xf8, 0xff, 0xf8, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
   0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0xff, 0xff, 0xff, 0xff,
   0xff, 0xff, 0xff, 0x1f, 0xff, 0x1f, 0xff, 0x1f};

/* Far Cursor */

static char far_bits[] = {
   0xf8, 0xff, 0xf8, 0xff, 0x18, 0xc0, 0x00, 0xc0, 0xff, 0xcf, 0x01, 0xc8,
   0x19, 0xc8, 0x19, 0xc8, 0x19, 0xc8, 0x19, 0xc8, 0x19, 0xc8, 0xf9, 0xeb,
   0xf9, 0xeb, 0x01, 0x08, 0x01, 0x08, 0xff, 0x0f};
static char far_mask[] = {
   0xfc, 0xff, 0xfc, 0xff, 0xfc, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
   0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0xff, 0xff, 0xff, 0xff,
   0xff, 0xff, 0xff, 0xff, 0xff, 0x1f, 0xff, 0x1f};

/* Reset bitmap */

static char reset_bits[] = {
   0x00, 0x00, 0x40, 0x40, 0x60, 0x60, 0x70, 0x70, 0x78, 0x78, 0x7c, 0x7c,
   0x7e, 0x7e, 0x7f, 0x7f, 0x7f, 0x7f, 0x7e, 0x7e, 0x7c, 0x7c, 0x78, 0x78,
   0x70, 0x70, 0x60, 0x60, 0x40, 0x40, 0x00, 0x00};

#ifdef __cplusplus
}
#endif

#endif