#ifndef CPGDOBLE_H
#define CPGDOBLE_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  void cpgline_d(int n, const double *xpts, const double *ypts);
  void cpggray_d(const double *a, int idim, int jdim, int i1, int i2, int j1, int j2, double fg, double bg, const double *tr);
  void cpgcons_d(const double *a, int idim, int jdim, int i1, int i2, int j1, int j2, const double *c, int nc, const double *tr);
  void cpghist_d(int n, const double *data, double datmin, double datmax, int nbin, int pgflag);
  void cpgpt_d(int n, const double *xpts, const double *ypts, int symbol);
  void cpgerry_d(int n, const double *x, const double *y1, const double *y2, double t);
  int cpgcurs_d(double *x, double *y, char *ch_scalar);
  int cpgband_d(int mode, int posn, double xref, double yref, double *x, double *y, char *ch_scalar);
  void cpgpoly_d(int n, const double *xpts, const double *ypts);
  void pgLimits_d(int n, double *x, float *wmin, float *wmax);


#ifdef __cplusplus
}
#endif


#endif
