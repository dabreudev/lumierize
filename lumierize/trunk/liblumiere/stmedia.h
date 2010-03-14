#ifndef STMEDIA_H
#define STMEDIA_H

#ifdef __cplusplus
extern "C" {
#endif

  float StMedia(int n, float *x, float *sigma);
  float StWeightMedia(int n, float *x, float *w,float *sigma);

  double StMedia_d(int n, double *x, double *sigma);
  double StWeightMedia_d(int n, double *x, double *w,double *sigma);
	double StErrWeightMedia_d(int n, double *x, double *err,double *sigma);

	double StSuma1_d(int n, double *a, int i);

#ifdef __cplusplus
}
#endif


#endif
