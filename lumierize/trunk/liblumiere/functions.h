#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

	int Fact(int n);
	double gammln(double xx);
	double lndergamm(double x);
	double gaussian(double x,double xmean,double sigma);
	double lngaussian(double x,double xmean,double sigma);
	double poidist(double x, double mean);
	double intgaussian(double x1, double x2, double xmean,double sigma);
	double int2dgaussian(double x1, double x2, double xmean, double y1, double y2, double ymean, double sigma);
	double gammq(double a,double x);
	double gammp(double a,double x);
	double incom(double a,double x);
	double erfcc(double x);
	void gcf( double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	double  Fermi(double x,double mu,double T);

#ifdef __cplusplus
}
#endif


#endif
