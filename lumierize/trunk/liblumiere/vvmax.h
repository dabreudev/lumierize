#ifndef VVMAX_H
#define VVMAX_H

#include "schechter.h"
#include "cosmology.h"

#ifdef __cplusplus
extern "C" {
#endif


	/* Calculo de la funcion de Lum mediante V/Vmax */
	int  VVmax_M(int n,double *mag_sel, double *mag_cal,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  VVmax_L(int n,double *flux_sel, double *flux_cal,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
#ifdef __cplusplus
}
#endif


#endif
