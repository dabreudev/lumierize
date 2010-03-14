#ifndef MRQMIN_H
#define MRQMIN_H

#ifdef __cplusplus
extern "C" {
#endif
  
  int Mrq(float  x[],float y[],float sig[],int ndata,float a[],int ia[],int ma, float **covar, float *chisq, void (*funcs)(float, float [], float *, float [],int));
  int Mrq_d(double  x[],double y[],double sig[],int ndata,double a[],int ia[],int ma, double **covar, double *chisq, void (*funcs)(double, double [], double *, double [],int));


#ifdef __cplusplus
}
#endif


#endif
