#ifndef MLSWML_H
#define MLSWML_H

#include "mlprocess.h"
#include "cosmology.h"
#include "schechter.h"
#include "histdist.h"
#include "fermisel.h"
#include "surveydb.h"
#include "poselfunc.h"

#ifdef __cplusplus
extern "C" {
#endif
  


	/* Funciones ML para LF */
	int  MLA_SWML_M    (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_L    (int n,double *flux,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
	int  MLA_SWML_PO   (int n,double *flux,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, double mlim, double ewlim, struct Histdist ewd, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
	/* Con normalizacion de Poisson */
	int  MLA_SWML_p_M  (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_p_L  (int n,double *flux,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
	int  MLA_SWML_p_PO (int n,double *flux,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, double mlim, double ewlim, struct Histdist ewd, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
	/* Con Funcion de seleccion en flujo/magnitud y normalizacion de Poisson */

	int  MLA_SWML_p_f_M  (int n,double *magn,double *z,struct fermifsel_M fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_p_f_L  (int n,double *flux,double *z,struct fermifsel_L fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
	int  MLA_SWML_p_f_PO (int n,double *magn,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, double strrad, struct cosmo_param cosmo,struct Schlf_L *lf);

	/* Con errores en flujos o magnitudes*/
	int  MLA_SWML_gm_M  (int n,double *magn,double *errmag ,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_gf_L  (int n,double *flux,double *errflux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
	int  MLA_SWML_gm_p_M(int n,double *magn,double *errmag ,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_gf_p_L(int n,double *flux,double *errflux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);

	/* Con errores en z */

	int  MLA_SWML_gz_M  (int n,double *magn,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_gz_L  (int n,double *flux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
	int  MLA_SWML_gz_p_M(int n,double *magn,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_gz_p_L(int n,double *flux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);

	/* Con errores en ambos */

	int  MLA_SWML_gmz_M  (int n,double *magn,double *errmag ,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_gfz_L  (int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
	int  MLA_SWML_gmz_p_M(int n,double *magn,double *errmag ,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
	int  MLA_SWML_gfz_p_L(int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);

#ifdef __cplusplus
}
#endif


#endif
