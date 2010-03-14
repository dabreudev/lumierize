#ifndef STHISTO_H
#define STHISTO_H

#ifdef __cplusplus
extern "C" {
#endif

  int *StHisto2_d(int n, double *a, int nbin, double *amin, double *amax);
  int *StHisto_d(int n, double *a, int nbin, double *amin, double *amax);
  int **StHisto2D_d(int n, double *x, double *y, int nbinx, double *xmin, double *xmax, int nbiny, double *ymin, double *ymax);
  int **StHisto2DFF_d(int n, double *x, double *y, int nbinx, double *xmin, double *xmax, double yff);


#ifdef __cplusplus
}
#endif


#endif
