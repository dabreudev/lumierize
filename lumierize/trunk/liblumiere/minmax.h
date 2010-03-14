#ifndef MINMAX_H
#define MINMAX_H

#ifdef __cplusplus
extern "C" {
#endif

	float minf(float x1,float x2);
	float maxf(float x1,float x2);
	int mini(int   x1,int   x2);
	int   maxi(int   x1,int   x2);
	void MinMax_d(int n, double *a, double *min, double *max);
	void MinMax(int n, float *a, float *min, float *max);

#ifdef __cplusplus
}
#endif


#endif
