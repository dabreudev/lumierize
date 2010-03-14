#ifndef AMOEBA_H
#define AMOEBA_H

#ifdef __cplusplus
extern "C" {
#endif

	/*  Amoeba */
	/*  ------ */
	int Amoe_NR(int ndata, float *xdata, float *ydata,
			float *p, float *y, int ndim, float ftol,
			int itmax,float (*amofunc)(int, float *, float *, float *));
	void Amoe_Ini(int ndata, float *xdata, float *ydata,
			int ndim, float *p0, float *sigp0, float *p, float *y,
			float (*amofunc)(int, float *, float *, float *));
	int Amoeba(int npt, float *xp, float *yp,int ndim, float *p0, float *sig0,
			float ftol, int itmax,
			float (*amofunc)(int, float *, float *, float *));
	int Amoe_NR_d(int ndata, double *xdata, double *ydata,
			double *p, double *y, int ndim, double ftol,
			int itmax,double (*amofunc)(int, double *, double *, double *));
	void Amoe_Ini_d(int ndata, double *xdata, double *ydata,
			int ndim, double *p0, double *sigp0, double *p, double *y,
			double (*amofunc)(int, double *, double *, double *));
	int Amoeba_d(int npt,double *xp,double *yp,int ndim, double *p0, double *sig0,
			double  ftol, int itmax,
			double (*amofunc)(int, double *, double *, double *));

#ifdef __cplusplus
}
#endif


#endif
