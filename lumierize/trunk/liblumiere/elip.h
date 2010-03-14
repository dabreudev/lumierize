#ifndef ELIP_H
#define ELIP_H

#ifdef __cplusplus
extern "C" {
#endif

  int  MCElipN_d(int npt, int ndim, double **x, double **C);
  int  MCElip_d(int n, double *x, double *y,
          double *a, double *b, double *c, double *d, double *f, double *e);
  int ElipPar_d(double a, double b, double c, double d, double f, double e,
         double *x0, double *y0, double *semax, double *semin,
         double *t);



#ifdef __cplusplus
}
#endif


#endif
