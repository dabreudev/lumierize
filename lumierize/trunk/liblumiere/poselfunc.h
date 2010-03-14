#ifndef POSELFUNC_H
#define POSELFUNC_H

#ifdef __cplusplus
extern "C" {
#endif


	struct poselfunc {   /* Para prisma objetivo */
		/* Galaxy parameters */
		double *ewbin;
		int nEW;
		double *zbin;
		int nz;
		double *magbin;
		int nmag;

		/* Observational parameters */
		double *seeing;
		int nseeing;
		double *sky;
		int nsky;
		double *transparency;
		int ntransparency;
		char **instsetup;
		int ninstsetup;
		double *******p;
		double *******errp;
		struct fermifsel_M  ******pfermi;
		struct fermifsel_M  ******errpfermi;
	};

	void writeposelfunc(struct poselfunc SF, char selfile[101]);
	void readposelfunc(struct poselfunc *SF, char selfile[101]);
	void projectposf(struct poselfunc SF, double **p, double **x, int *nbin, int iproj);
	int  projectposf_fixonerot(struct poselfunc SF, double **p, double **x, int *nbin, int iproj, int ifix);
	void projectposf_fixone(struct poselfunc SF, double **p, double **x, int *nbin, int iproj, int ifix, double valuefix);
	void projectposf_fixtwo(struct poselfunc SF, double **p, double **x, int *nbin, int iproj, int ifix1, double valuefix1, int ifix2, double valufix2);
	void project2posf(struct poselfunc SF, double ***p, double **x, double **y, int *nx, int *ny, int iprojx, int iprojy);
	double prob_poselfunc_scale(struct poselfunc fsel, struct SurveyItem si, double ew, double z, double magn);
	double errprob_poselfunc(struct poselfunc fsel, struct SurveyItem si, double ew, double z, double magn);
	int sf_selectbin(double value, double *bins, int nbin);
	int sf_selectbin_char(char value[], char **bins, int nbin);

	void writediselfunc(struct poselfunc SF, char selfile[101]);
	void readdiselfunc(struct poselfunc *SF, char selfile[101]);
#ifdef __cplusplus
}
#endif


#endif
