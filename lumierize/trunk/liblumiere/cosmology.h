#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#ifdef __cplusplus
extern "C" {
#endif

	/* cosmology */
	struct cosmo_param
	{
		double H0; /* Hubble constant */
		double OM; /* Omega matter */
		double OL; /* Omega lambda */
		double Ok; /* Omega k */
		double DH; /* Hubble distance: c/H0 */
		struct cached_log10dlum_dVdz cached;
	};

	void cosmo_init(struct cosmo_param *cosmo, double H0, double OM, double OL);
	void cosmo_free(struct cosmo_param *cosmo);
	double Efunction(double z, struct cosmo_param cosmo);
	double Vol(double z, struct cosmo_param cosmo);
	double dVdz(double z, struct cosmo_param cosmo);
	double Lum(double z, double flux, struct cosmo_param cosmo);
	double dLumdflux(double z, struct cosmo_param cosmo);
	double Flux(double z, double Lum, struct cosmo_param cosmo);
	double mag(double z, double M, struct cosmo_param cosmo);
	double Mag(double z, double m, struct cosmo_param cosmo) ;
	void Mag_err(double z, double errz, double m, double errm,struct cosmo_param cosmo, double *Mag, double *errMag);
	double Z_l(double flux, double L, struct cosmo_param cosmo);
	double Z_m(double m, double M, struct cosmo_param cosmo);
	/* only valid for Omega lambda = 0 */
	double Vol_OmegaLambda0(double z, struct cosmo_param cosmo);
	double dVdz_OmegaLambda0(double z, struct cosmo_param cosmo);
	double mag_OmegaLambda0(double z, double M, struct cosmo_param cosmo);
	double Mag_OmegaLambda0(double z, double m, struct cosmo_param cosmo) ;
	double Z_m_OmegaLambda0(double m, double M, struct cosmo_param cosmo);
	//  double D_co(double z, struct cosmo_param cosmo);
	//  double D_lum(double z, struct cosmo_param cosmo);
	//  double D_ang(double z, struct cosmo_param cosmo);
	double Flux_ew_mag(double ew, double mag, char photband[51], float gamma, float delta, float Kcoc);
	void   Flux_ew_mag_err(double ew, double errew, double mag, double errmag, char photband[51], float gamma, float delta, float Kcoc, double *fluxline, double *errfluxline);
	double mag_ew_flux(double ew, double flux, char photband[51], float gamma, float delta, float Kcoc);
#ifdef __cplusplus
}
#endif


#endif
