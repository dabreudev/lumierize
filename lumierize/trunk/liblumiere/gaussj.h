#ifndef GAUSSJ_H
#define GAUSSJ_H

#ifdef __cplusplus
extern "C" {
#endif

  int gaussj_d(double **a, int n, double **b, int m);
  int gaussj(float **a, int n, float **b, int m);

#ifdef __cplusplus
}
#endif


#endif
