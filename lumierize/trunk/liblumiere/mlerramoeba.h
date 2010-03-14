#ifndef MLERRAMOEBA_H
#define MLERRAMOEBA_H

#ifdef __cplusplus
extern "C" {
#endif
  
  int mlerr_amo(int npt, float *x, float *y, int npar, float *par, float *sigpar,
                float ftol, int niter,
                float (*amofunc_main)(int, float *, float *, float *),
                int nconfl, float **covarpar);

  int mlerr_amo_d(int npt, double *x, double *y, int npar, double *par,
                  double *sigpar, double ftol, int niter,
                  double (*amofunc_main)(int, double *, double *, double *),
                  int nconfl, double **covarpar);


#ifdef __cplusplus
}
#endif


#endif
