#ifndef MLHIST_H
#define MLHIST_H

#ifdef __cplusplus
extern "C" {
#endif


	int  ML_g_g_d(int n,double *x,double *errx,double *mean,double *sigma,double *errmean,double *errsigma,double *covarmeansigma);
	int  ML_g_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma);
	int  MLA_g_g_d(int n,double *x,double *errx,double *mean,double *sigma,double *errmean,double *errsigma,double *covarmeansigma);
	int  MLA_g_g_f_d(int n,double *x,double *errx,double xfermi, double Tfermi,double *mean,double *sigma,double *errmean,double *errsigma,double *covarmeansigma);
	int  MLA_g_n_f_d(int n,double *x,double xfermi, double Tfermi, double *mean,double *sigma,double *errmean,double *errsigma,double *covarmeansigma);

	int  MLA_g_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma);
	int  ML_g_g_corr(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma);
	int  ML_g_g_corr_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma);
	int  MLA_h_g_d(int n,double *x,double *errx, int k, double *xk, double *Pk, double **covPk);
	int  MLA_h_g_f_d(int n,double *x,double *errx,double xfermi, double Tfermi, int k, double *xk, double *Pk, double **covPk);
	int  MLA_hh_gg_d(int n,double *x,double *errx, double *y, double *erry, int kx, int ky, double *xk, double *yk, double *Pk, double **covPk);
	int  MLA_ff_gg_d(int n,double *x,double *errx, double *y, double *erry, int kx, double *xk, double yff, double yinter, double *Pk, double **covPk);


#ifdef __cplusplus
}
#endif


#endif
