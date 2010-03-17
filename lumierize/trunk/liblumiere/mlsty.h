#ifndef MLSTY_H
#define MLSTY_H

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
	int  MLA_STY_M     (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo);
	int  MLA_STY_M2    (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf);
	int  MLA_STY_L     (int n,double *flux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf, struct MLProcessInfo *mlinfo);
	int  MLA_STY_PO    (int n,double *flux,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, double mlim, double ewlim, struct Histdist ewd, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
	/* Con normalizacion de Poisson */
	int  MLA_STY_p_M   (int n,double *magSeln,double *magDistn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo);
	int  MLA_STY_p_M_wC (int n,double *magSel, double *magDist, double color_mean, double color_stddev, double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo);
	int  MLA_STY_p_L   (int n,double *flux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf, struct MLProcessInfo *mlinfo);
	int  MLA_STY_p_PO  (int n,double *flux,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, double mlim, double ewlim, struct Histdist ewd, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
	/* Con Funcion de seleccion en flujo/magnitud y normalizacion de Poisson */

	int  MLA_STY_p_f_M   (int n,double *magn,double *z,struct fermifsel_M fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf);
	int  MLA_STY_p_f_L   (int n,double *flux,double *z,struct fermifsel_L fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
	int  MLA_STY_p_f_PO  (int n,double *magn,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, double strrad, struct cosmo_param cosmo,struct Schlf_L *lf);

	/* Con errores en color */
	int MLA_STY_gc_p_M_wC(int n,double *magSeln, double *magDistn, double color_mean, double color_stddev, double *errColorn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo);


	/* Con errores en flujos o magnitudes*/

	int  MLA_STY_gm_M   (int n,double *magn,double *errmag ,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf);
	int  MLA_STY_gf_L   (int n,double *flux,double *errflux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
	int  MLA_STY_gm_p_M (int n,double *magn,double *errmag ,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo);
	int  MLA_STY_gf_p_L (int n,double *flux,double *errflux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);

	/* Con errores en z */

	int  MLA_STY_gz_M   (int n,double *magn,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf);
	int  MLA_STY_gz_L   (int n,double *flux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
	int  MLA_STY_gz_p_M (int n,double *magn,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf);
	int  MLA_STY_gz_p_L (int n,double *flux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);

	/* Con errores en ambos */

	int  MLA_STY_gmz_M   (int n,double *magn,double *errmag ,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf);
	int  MLA_STY_gfz_L   (int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
	int  MLA_STY_gmz_p_M (int n,double *magn,double *errmag ,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf);
	int  MLA_STY_gfz_p_L (int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
  int  MLA_STY_g_p_L   (int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);

	/* Con funciones de selecci�n y varias im�genes */

	int  MLA_STY_s_p_f_PO  (int n,double *magn,double *ew,double *z,int *isurvey,char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyDB sdb, double ewlim, struct Histdist ewd, struct cosmo_param cosmo,struct Schlf_L *lf);

	/* Con funciones de selecci�n y varias im�genes */

	int  MLA_STY_e_s_p_f_PO  (int n ,double *magn ,double *ew ,double *z ,double *ebv , double *cocnii, int *isurvey,char photband[51], float gamma, float delta, float Kcoc, float Aext, struct poselfunc fsel, struct SurveyDB sdb, double ewlim, struct Histdist ewd, struct cosmo_param cosmo,struct Schlf_L *lf);

	/* Simplemente de densidad */

	int  MLA_rho(int n,double *ra, double *dec, int *ima, struct SurveyDB sdb, double *probdetec, double *dens, double *errdens);

	/* Con casi todo lo que puedas imaginar (dabreu) ;) */
	int  MLA_STY_gmz_p_f_M_wC(int n,double *magSeln, double *magDistn, double *errmagDistn, double color_mean, double color_stddev, double *z, double *errz, struct fermifsel_M fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo, struct Schlf_M *lf, struct MLProcessInfo *mlinfo);
#ifdef __cplusplus
}
#endif


#endif
