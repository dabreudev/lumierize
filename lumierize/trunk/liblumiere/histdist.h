#ifndef HISTDIST_H
#define HISTDIST_H

#ifdef __cplusplus
extern "C" {
#endif


	/* Estructuras con definici<F3>n de data de distribuciones */

	struct Histdist {
		double *xk;       /* Dimensionado k+1 */
		double *errxk;    /* Dimensionado k+1 */
		double *Pk;       /* Dimensionado k*/
		double *errPk;    /* Dimensionado k*/
		double **covarPk; /* Dimensionado k*k */
		int k;
		/* Pk[i]  es el valor de la distribuci<F3>n en el intervalo
	       x[i] - x[i+1]. El array x est<E1> en orden ascendente */
		/* errPk[i]=sqrt(covarPk[i][i]) debe cumplirse */
	};

	double Histfunc(double x, struct Histdist hd);

#ifdef __cplusplus
}
#endif


#endif
