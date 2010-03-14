#ifndef STHISTO_H
#define STHISTO_H

#ifdef __cplusplus
extern "C" {
#endif

  int *StHisto2_d(int n, double *a, int nbin, double *amin, double *amax);
  int *StHisto_d(int n, double *a, int nbin, double *amin, double *amax);


#ifdef __cplusplus
}
#endif


#endif
