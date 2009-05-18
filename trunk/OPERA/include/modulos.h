#ifndef OPERA_H
#define OPERA_H

#ifdef __cplusplus
extern "C" {
#endif
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <glob.h>
  /* #include <readline/readline.h> */
  /* #include <readline/history.h> */
  
  
  
#include "cpgplot.h"
#include "fitsio.h"
#include "cbutton.h"
#include "wcs.h"
  /*#include "lwcs.h"*/
#include "fitsfile.h"
#include "fitshead.h"
#include "imio.h"
#include "wcslib.h"
#include "wcscat.h"
  
  
  /* Funciones de allocateo */
  float *******tensor7_f(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
  float ******tensor6_f(int n1, int n2, int n3, int n4, int n5, int n6);
  float *****tensor5_f(int n1, int n2, int n3, int n4, int n5);
  float ****tensor4_f(int n1, int n2, int n3, int n4);
  float ***tensor_f(int n1, int n2, int n3);
  float **matrix_f(int n1, int n2);
  float *vector_f(int n);
  float **vector_pf(int n);
  float ***vector_ppf(int n);
  float ****vector_pppf(int n);
  float *****vector_ppppf(int n);
  float ******vector_pppppf(int n);
  float *******vector_ppppppf(int n);
  double *******tensor7_d(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
  double ******tensor6_d(int n1, int n2, int n3, int n4, int n5, int n6);
  double *****tensor5_d(int n1, int n2, int n3, int n4, int n5);
  double ****tensor4_d(int n1, int n2, int n3, int n4);
  double ***tensor_d(int n1, int n2, int n3);
  double **matrix_d(int n1, int n2);
  double *vector_d(int n);
  double **vector_pd(int n);
  double ***vector_ppd(int n);
  double ****vector_pppd(int n);
  double *****vector_ppppd(int n);
  double ******vector_pppppd(int n);
  double *******vector_ppppppd(int n);
  int *******tensor7_i(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
  int ******tensor6_i(int n1, int n2, int n3, int n4, int n5, int n6);
  int *****tensor5_i(int n1, int n2, int n3, int n4, int n5);
  int ****tensor4_i(int n1, int n2, int n3, int n4);
  int ***tensor_i(int n1, int n2, int n3);
  int **matrix_i(int n1, int n2);
  int *vector_i(int n);
  int **vector_pi(int n);
  int ***vector_ppi(int n);
  int ****vector_pppi(int n);
  int *****vector_ppppi(int n);
  int ******vector_pppppi(int n);
  int *******vector_ppppppi(int n);
  char **vector_s(int n, int nchar);
  char **vector_pps(int n);
  char *alloc_s(int nchar);
  void free_tensor7_f(float  *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7);
  void free_tensor7_d(double *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7);
  void free_tensor7_i(int    *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7);
  void free_tensor6_f(float  ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6);
  void free_tensor6_d(double ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6);
  void free_tensor6_i(int    ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6);
  void free_tensor5_f(float  *****tensor5, int n1, int n2, int n3, int n4, int n5);
  void free_tensor5_d(double *****tensor5, int n1, int n2, int n3, int n4, int n5);
  void free_tensor5_i(int    *****tensor5, int n1, int n2, int n3, int n4, int n5);
  void free_tensor4_f(float  ****tensor4, int n1, int n2,int n3, int n4);
  void free_tensor4_d(double ****tensor4, int n1, int n2,int n3, int n4);
  void free_tensor4_i(int    ****tensor4, int n1, int n2,int n3, int n4);
  void free_tensor_f(float  ***tensor, int n1, int n2, int n3);
  void free_tensor_d(double ***tensor, int n1, int n2, int n3);
  void free_tensor_i(int    ***tensor, int n1, int n2, int n3);
  void free_matrix_f(float  **mat, int n1, int n2);
  void free_matrix_d(double **mat, int n1, int n2);
  void free_matrix_i(int    **mat, int n1, int n2);
  void free_vector_s(char **v, int n, int nchar);
  void DeleteRecord_f(float *vec, int n, int irecord);
  void DeleteRecord_i(int *vec, int n, int irecord);
  void DeleteRecord_d(double *vec, int n, int irecord);
  void DeleteRecord_s(char *vec, int n, int nrec, int irecord);
  

/* Funciones inline */

#define SWAP(a,b) {float swap_temp=(a);(a)=(b);(b)=swap_temp;}
#define SWAP_D(a,b) {double swap_temp=(a);(a)=(b);(b)=swap_temp;}
  float minf(float x1,float x2);
  float maxf(float x1,float x2);
  int   mini(int   x1,int   x2);
  int   maxi(int   x1,int   x2);
  
  
  
  /* Funciones ML */

  struct MLProcessInfo
  {
    unsigned int nIter;
    unsigned int ntries;
    double MLmax;
  };
  
  int mlerr_amo(int npt, float *x, float *y, int npar, float *par, float *sigpar, float ftol, int niter, float (*amofunc_main)(int, float *, float *, float *), int nconfl, float **covarpar);
  int mlerr_amo_d(int npt, double *x, double *y, int npar, double *par, double *sigpar, double ftol, int niter, double (*amofunc_main)(int, double *, double *, double *), int nconfl, double **covarpar);


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


  
  
  struct  plt_cte {
    
    float ap,bp,cp,dp,ep,fp,gp,hp,pp,qp;
    float ae,be,ce,de,ee,fe,ge,he,pe,qe;
    float amdx1,amdx2,amdx3,amdx4,amdx5,amdx6,amdx7,amdx8,amdx9;
    float amdx10,amdx11,amdx12,amdx13,amdx14,amdx15,amdx16,amdx17,amdx18,amdx19,amdx20;
    float amdy1,amdy2,amdy3,amdy4,amdy5,amdy6,amdy7,amdy8,amdy9;
    float amdy10,amdy11,amdy12,amdy13,amdy14,amdy15,amdy16,amdy17,amdy18,amdy19,amdy20;
    float xpixelsz,ypixelsz;
  };
  
  float  circ_aper(float *ima, int nx, int ny, float x,float y,float r,int *npix);
  float  circ_aper_err(float *ima, float *errima, int nx, int ny, float x,float y,float r,int *npix,float *errflux);
  float  circ_aper_bak(float *ima, int nx, int ny, float x,float y,float r,int *npix);
  float  circ_aper_err_bak(float *ima, float *errima, int nx, int ny, float x,float y,float r,int *npix,float *errflux);
  float  elip_aper(float *ima, int nx, int ny, float x,float y,float a, float b, float pa,int *npix);
  
  
  float *SplNodX(int ndat, float *x, float *y, int nodos, float *xnodo, float *ynodo, float ftol);
  int Tridiag(float *a,float *b,float *c,float *r,float *u,int n);

  void LeeWord(char a[],int nw,char word[]);
  int ReadNumcol(char file[],int col,float *vector,int *lvec,int *nlin);
  int ReadDoublecol(char file[],int col,double *vector,int *lvec,int *nlin);
  int ReadHourcol(char file[],int col,float *vector,int *lvec,int *nlin);
  int ReadCharcol(char file[],int col,char *vector,int *lvec,int charlen,int *nlin);
  int ReadWCScol(char file[],int col,double *vector,int *lvec,int *nlin);
  int Disper( float cosx, float cosy, float cocin, float alfa, float *cosxs, float *cosys );
  
  int getline(char s[],int lim);
  int fgetline(FILE *fp,char s[],int lim);
  void kbdpause(void);
  int readi(int n);
  float readf(float n);
  char readc(char c);
  double readd(double n);
  void reads(const char* def,char *outstr);
  
  float   Lagr2(float *x, float *y, int n, float x0);
  float   Lagr4(float *x, float *y, int n, float x0);
  double   Lagr2_d(double *x, double *y, int n, double x0);
  double   Lagr4_d(double *x, double *y, int n, double x0);
  
  /* Integracion de funciones */
  void gauher_d(double x[], double w[], int n);
  double gaussinther_d(double (*funk)(double), double offset, double scale, int n);
  void gauher(float  x[], float  w[], int n);
  float  gaussinther(float  (*funk)(float ), float offset, float  scale, int n);
  void gauleg_d(double x1, double x2, double  x[], double  w[], int n);
  double gaussintleg_d(double (*funk)(double), double x1, double x2, int n);
  void gauleg(float  x1, float  x2, float   x[], float   w[], int n);
  float  gaussintleg(float  (*funk)(float ), float  x1, float  x2, int n);
  void gaulag_d(double x[], double w[], int n, double alf);
  double gaussintlag_d(double (*funk)(double), double scale, double alfa, int n);
  void gaulag(float  x[], float  w[], int n, float alf);
  float  gaussintlag(float  (*funk)(float ), float  scale, float alfa, int n);
  float intimapix_lin(float *ima, int nx, int ny, float xp, float yp, float nullval);
  float intspecpix_lin(float *spec, int nx, float xp, float nullval);

  float sumspecpix_lin(float *spec,int nx,float pix1,float pix2);
  float sumspecpix_vista(float *spec,int nx,float pix1,float pix2);
  

  /* Funciones de numeros aleatorios */

  int    Pdev(double prob);
  double Poidev(double xm);
  double Expdev(void);
  double Gasdev(void);
  double Bnldev(double pp,int n);
  double Histdev(int k,double *Pk,double *xk);
  double Powdev(double xmin, double xmax,double alfa);
  double Constdev(double xmin, double xmax);
  float ran1(long *idum);
  float ran2(long *idum);

  /* Estructuras con definición de data de distribuciones */

  struct Histdist {
    double *xk;       /* Dimensionado k+1 */
    double *errxk;    /* Dimensionado k+1 */
    double *Pk;       /* Dimensionado k*/
    double *errPk;    /* Dimensionado k*/
    double **covarPk; /* Dimensionado k*k */
    int k;
    /* Pk[i]  es el valor de la distribución en el intervalo
       x[i] - x[i+1]. El array x está en orden ascendente */
    /* errPk[i]=sqrt(covarPk[i][i]) debe cumplirse */
  };
  
  
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
  double Histfunc(double x, struct Histdist hd);

  
  
  
  void Gausfilter2D(float *ima_input, float *ima_output,int nx, int ny, float sigma);
  /* void Sort(int n,float *x,float *y); */
  /* void sort(int n,float ra[]); */
  /* int  sort2(int *n,float ra[]); */
  void indexx(unsigned int n, double *arr, unsigned int *indx);
  void Quartil(int n,float *x,float *first, float *median,float *third);
  void Quartil_d(int n,double *x,double *first, double *median,double *third);
  float poly(float x,int g, float *c);
  double poly_d(double x,int g, double *c);
  void covcor(float *a,int k,int l,float *med,float *b,float *c);
  void eigen(float *m,int l,float *eigenval,float *eigenvec);
  
  
  /* Funciones para trabajar con imagenes FITS */
  /* ******************************************* */
  
  struct image {
    char file[200];
    long naxes[2];
    long npixels;
    fitsfile *filefits;
    float datamin, datamax;
    float *array;
    int nx,ny;
    int aloc_flag;
  };
  struct wcsimage {
    struct image image;
    struct WorldCoor *wcs;
    int wcsset;
  };
  void ReadImage(struct image *im); 
  void ReadWCSImage(struct wcsimage *im); 
  void SaveImage(struct image im);
  void SaveImage_16(struct image im);
  void PlotImage(struct image *im);
  void PlotImage_sec(struct image *im, int x1, int x2, int y1, int y2);
  void PlotImage_sec_cuts(struct image *im, int x1, int x2, int y1, int y2, float bg, float fg);

  
  /* Funciones para trabajar con espectros FITS */
  /* ******************************************* */


  struct spectrum {
    char file[300];
    long naxes[2];
    long npixels;
    fitsfile *filefits;
    float datamin, datamax;
    float *spec;
    float *ldo;
    int nx;
    float ldomin;
    float deltaldo;
    int aloc_flag;
    int alocldo_flag;
  };
  
  void ReadSpec(struct spectrum *sp);
  void CloseSpec(struct spectrum *sp);
  void ReadSpec_buf(struct spectrum *sp);
  void SaveSpec(struct spectrum sp);
  void FillSpec(struct spectrum *sp, int npix, float *spec, float ldomin, float ldomax);
  void CopySpec(struct spectrum *spdest, struct spectrum sporig);
  void PlotSpec(struct spectrum sp);
  void PlotSpec_zoom(struct spectrum sp, float x1, float x2, float y1, float y2);
  void PlotSpec_ov(struct spectrum sp);
  void PlotSpec_err(struct spectrum sp, struct spectrum errsp);
  void PlotSpec_err_zoom(struct spectrum sp, struct spectrum errsp, float x1, float x2, float y1, float y2);
  void PlotSpec_pix(struct spectrum sp);
  void PlotSpec_pix_ov(struct spectrum sp);
  void PlotSpec_pix_err(struct spectrum sp, struct spectrum errsp);

  /* Funciones de FFT  y correlacion cruzada */
  /*******************************************/

  void  fftcorrel(int n, float *data1, float *data2, float *xcorr, float *fcorr);
  float offsetpix(int n, float *spec1, float *spec2);
  void  cfft(int n, float *xr, float *xi, int imode);
  int   fft2power(int n);
  float cosbell(int i, int n, float fl);

  /* Funciones para transformar espectros */
  /***************************************/

  void  shiftspec(float xshift, int n, float *specorig, float *newspec);

  /* Funciones para trabajar con funciones de dispersion de prismas */
  struct disper_prism {
    float A,B,C;
    float tampix;
  };
  
  float pr_ldo2pix(float ldo,struct disper_prism DP);
  float pr_pix2ldo(float x,struct disper_prism DP);
  float pr_dpdl(float ldo,struct disper_prism DP);
  float pr_dldp(float x,struct disper_prism DP);

  void SP_ldo2pix(struct spectrum *SP_pix, struct spectrum SP_lam, struct disper_prism DP);
  void SP_pix2ldo(struct spectrum *SP_lam, struct spectrum SP_pix, struct disper_prism DP);
    
  /* Funciones y estrucutras para trabajar con bases de datos de Surveys */
  /* ******************************************************************* */

  struct SurveyItem {
    
    /*Astrometric information */
    double alfac, deltac;
    float xdim, ydim;
    float rot;
    float epoch;
    float stdx, stdy;
    
    /*Flux calibration variables. mag = a + b log10(flux)*/
    float a,b;
    float erra,errb,covab;
    float rms;
    
    /*Limiting magnitude*/
    float mmaximum;
    float mmode;
    float mfermicut;
    float deltafermi;
    float gammapowerlaw;
    float errmfermicut;
    float errdeltafermi;
    float errgammapowerlaw;
    
    /*Seeing variables */
    float seeing;
    float stdseeing;
    
    /*Image variables */
    char image[101];

    /*File cross variables */
    char calfile[101];
    int nobj;

    /* Sky britghness variables */
    float sky;
    float sigsky;

    /* Spectra variables */
    char specfile[101];
    char respfile[101];

   /* Transparency */
   float transparency;   /* con respecto a la teorica de 1, por ejemplo, o por lo menos algo proporcional a la transparencia  */

   /* Instrumental setup */
   char instsetup[51];
    
  };
  
  
  struct SurveyDB {
    
    struct SurveyItem *si;
    int nitems;
    
  };
  
  void ReadSurDB(char *dbfile, struct SurveyDB  *sdb);
  void SaveSurDB(char *dbfile, struct SurveyDB sdb);
  int  whithinimage(double ra, double dec, struct SurveyItem sitem);
  float Surveyrad(struct SurveyDB sdb);
  void PlotSurDB(struct SurveyDB sdb, int labflag);
  void PlotSurDB_zoom(struct SurveyDB sdb, float ramin, float ramax, float decmin, float decmax, int labflag);


  
  /* Para trabajar con funciones de seleccion */
  struct fermifsel_M  {
    double  magcut;   
    double  deltamag; 
  };
  struct fermifsel_L  {
    double  fluxcut;   
    double  deltaflux; 
  };
  struct diselfunc {   /* Para imagen directa */
    float *sky;
    int nsky;
    float *transparency;
    int ntransparency;
    char *instsetup;
    int ninstsetup;
    struct fermifsel  *p;
    struct fermifsel  *errp;
  };

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
  double prob_poselfunc_scale(struct poselfunc fsel, struct SurveyItem si, double ew, double z, double magn);
  double errprob_poselfunc(struct poselfunc fsel, struct SurveyItem si, double ew, double z, double magn);
  double errprob_poselfunc(struct poselfunc fsel, struct SurveyItem si, double ew, double z, double magn);
  int sf_selectbin(double value, double *bins, int nbin);
  int sf_selectbin_char(char value[], char **bins, int nbin);

  




  void writediselfunc(struct poselfunc SF, char selfile[101]);
  void readdiselfunc(struct poselfunc *SF, char selfile[101]);


  
  /* Funciones para trabajar con funciones de luminosidad y cosmologia */
  /* #################################################### */

  struct cached_log10dlum_dVdz /* log10dlum and dVdz */
  {
    int N; /* Number of points */
    double *log10dlum; /* log10 dlum values */
    double *dVdz; /* dVdz values */
    double *z;  /* z for the values */
    double *invlog10dlum; /* dlum values sampled */
    double *invz;    /* z values for inversedlum */
    double maxlog10dlum; /* max of inversedlum */
    double minlog10dlum; /* min of inversedlum */
    double *DM; /* comoving distance (transverse) */
  };

  struct cosmo_param
  {
    double H0; /* Hubble constant */
    double OM; /* Omega matter */
    double OL; /* Omega lambda */
    double Ok; /* Omega k */
    double DH; /* Hubble distance: c/H0 */
    struct cached_log10dlum_dVdz cached;
  };

  struct Schlf_M {
    double alfa;
    double erralfa;
    double phistar;
    double errphistar;
    double Mstar;
    double errMstar;

    double covaralfaphistar;
    double covaralfaMstar;
    double covarphistarMstar;
  };

  struct Schlf_L {
    double alfa;
    double erralfa;
    double phistar;
    double errphistar;
    double Lstar;  /*Esta debe ir SIN logaritmos */
    double errLstar;

    double covaralfaphistar;
    double covaralfaLstar;
    double covarphistarLstar;
  };

  struct Steplf_M 
  {
    int nbin;
    double *magni;
    double *errmagni;
    double *lnlf;     /* Estan en log(Phi_i) neperiano */
                      /* Usamos el criterio de que cuando la FL vale 0 lnlf 
                         vale menos infinito (-1/0.) */
    double *errlnlf;  /* Es el error de lf (expresado como antes) */
    double **covarlnlf;
    double *lf;
    double *errlf;
    /* Cuidado!! El bin j va desde magni[j] a magni[j+1], 
       luego magni esta dimensionado
       de nbin+1 y lf de nbin. 
       De esta manera, el intervalo es magni[j+1]-magni[j] */
    /* magni[0] contiene el valor mas negativo, 
       es decir, la magnitud mas brillante.
       Esta por tanto en orden ascendente */
  };

  struct Steplf_L 
  {
    int nbin; 
    double *lumi;    /* Estas luminosidades estan en log(Lum_i) neperiano */
    double *errlumi; /* Es el error en luminosidades (expresado en log(Lum)) */
    double *lnlf;      /* Estan en log(Phi_i) neperiano */
                       /* Usamos el criterio de que cuando la FL vale 0 lnlf 
                         vale menos infinito (-1/0.) */
    double *errlnlf;   /* Es el error de lf (expresado como antes) */

    double **covarlnlf;
    /* Esta funcion representa la funcion de lum en unidades de 
       (#/Mpc3/W) para el intervalo de luminosidades
       exp(lumi[i]) - exp(lumi[i+1]), con lo que si pintas la funcion
       de luminosidad en luminosidades te sale constante en todo ese 
       intervalo. Otra cosa es que luego  tu definas los intervalos 
       para que sean logaritmicos y luego pintes el eje X de Lumi
       en luminosidades */
  };

  double Schechter_M(double M, struct Schlf_M lf);
  double Schechter_L(double L, struct Schlf_L lf);
  double Schechterdev_M(struct Schlf_M lf, double Mlow,double Mup);
  double zSchdev_M(struct Schlf_M lf, double zlow,double zup,double mlow,double mup, struct cosmo_param cosmo);
  double Schechterdev_L(struct Schlf_L lf, double Llow,double Lup);
  double zSchdev_L(struct Schlf_L lf, double zlow,double zup,double ffaint,double fbright, struct cosmo_param cosmo);
  double Int_sch_M(struct Schlf_M lf, double zlow,double zup,double mlim,struct cosmo_param cosmo);
  double Int_sch_L(struct Schlf_L lf, double zlow,double zup,double fluxlim,struct cosmo_param cosmo);
  double Int_sch_M_wC(struct Schlf_M lf, double zlow,double zup, double color_mean, double color_stddev,double mlim,struct cosmo_param cosmo);
  double Int_sch_L_wC(struct Schlf_L lf, double zlow,double zup, double color_mean, double color_stddev,double fluxlim,struct cosmo_param cosmo);
  double Int_sch_f_M(struct Schlf_M lf, double zlow,double zup,struct fermifsel_M fsel, struct cosmo_param cosmo);
  double Int_sch_f_M_wC(struct Schlf_M lf, double zlow,double zup, double color_mean, double color_stddev, struct fermifsel_M fsel, struct cosmo_param cosmo);
  double Int_sch_f_L(struct Schlf_L lf, double zlow,double zup,struct fermifsel_L fsel,struct cosmo_param cosmo);
  double Int_sch_PO(struct Schlf_L lf, char photband[51], double gamma, double delta, double Kcoc, double zlow,double zup,double mlim, double ewlim, struct Histdist ewd, struct cosmo_param cosmo);
  double Int_sch_f_PO(struct Schlf_L lf, char photband[51], double gamma, double delta, double Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, struct cosmo_param cosmo);
  double Int_sch_s_f_PO(struct Schlf_L lf, char photband[51], double gamma, double delta, double Kcoc, struct poselfunc fsel, struct SurveyDB sdb , double ewlim, struct Histdist ewd, struct cosmo_param cosmo);
  double Int_sch_e_s_f_PO(struct Schlf_L lf, char photband[51], double gamma, double delta, double Kcoc, double A_ext, double meanebv, double meancocnii, struct poselfunc fsel, struct SurveyDB sdb , double ewlim, struct Histdist ewd, struct cosmo_param cosmo);


  double Sch_rhoz_L(double z, struct Schlf_L lf, double zlow,double zup,double fluxlim,struct cosmo_param cosmo);
  double Sch_rhoz_M(double z, struct Schlf_L lf, double zlow,double zup,double mlim,struct cosmo_param cosmo);

  double Step_M(double M, struct Steplf_M lf);
  double Step_L(double L, struct Steplf_L lf);
  double Steplfdev_M(struct Steplf_M lf, double Mlow,double Mup);
  double zSteplfdev_M(struct Steplf_M lf, double zlow,double zup,double mlow,double mup, struct cosmo_param cosmo);
  double Steplfdev_L(struct Steplf_L lf, double Llow,double Lup);
  double zSteplfdev_L(struct Steplf_L lf, double zlow,double zup,double ffaint,double fbright, struct cosmo_param cosmo);
  double Int_step_M(struct Steplf_M lf, double zlow,double zup,double mlim,struct cosmo_param cosmo);
  double Int_step_L(struct Steplf_L lf, double zlow,double zup,double fluxlim, struct cosmo_param cosmo);
  double Int_step_f_M(struct Steplf_M lf, double zlow,double zup,struct fermifsel_M fsel,struct cosmo_param cosmo);
  double Int_step_f_L(struct Steplf_L lf, double zlow,double zup,struct fermifsel_L fsel, struct cosmo_param cosmo);

  int  FitSch2StepLF_L(struct Steplf_L lfstep, struct Schlf_L *lfsch, double *chisq);
  int  FitSch2StepLF_M(struct Steplf_M lfstep, struct Schlf_M *lfsch, double *chisq);

  void PlotStepLF_L(struct Steplf_L lf);
  void PlotStepLF_M(struct Steplf_M lf);
  void PlotStepLF_L_ov(struct Steplf_L lf);
  void PlotStepLF_M_ov(struct Steplf_M lf);
  void PlotSchLF_L(struct Schlf_L lf);
  void PlotSchLF_M(struct Schlf_M lf);
  void PlotSchLF_L_ov(struct Schlf_L lf);
  void PlotSchLF_M_ov(struct Schlf_M lf);
  void PlotStepSchLF_L(struct Steplf_L lfstep, struct Schlf_L lfsch);
  void PlotStepSchLF_M(struct Steplf_M lfstep, struct Schlf_M lfsch);

  void PrintStepLF_M(struct Steplf_M lf);
  void PrintStepLF_L(struct Steplf_L lf);

  /* cosmology */
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



  /* Calculo de la funcion de Lum mediante V/Vmax */
  int  VVmax_M(int n,double *mag_sel, double *mag_cal,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  VVmax_L(int n,double *flux_sel, double *flux_cal,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);


  /* Funciones ML para LF */
  int  MLA_STY_M     (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo);
  int  MLA_STY_M2    (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf); 
  int  MLA_STY_L     (int n,double *flux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf, struct MLProcessInfo *mlinfo); 
  int  MLA_STY_PO    (int n,double *flux,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, double mlim, double ewlim, struct Histdist ewd, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
  int  MLA_SWML_M    (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  MLA_SWML_L    (int n,double *flux,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
  int  MLA_SWML_PO   (int n,double *flux,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, double mlim, double ewlim, struct Histdist ewd, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
  /* Con normalizacion de Poisson */
  int  MLA_STY_p_M   (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo); 
  int  MLA_STY_p_M_wC (int n,double *magSel, double *magDist, double color_mean, double color_stddev, double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo);
  int  MLA_STY_p_L   (int n,double *flux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf, struct MLProcessInfo *mlinfo);
  int  MLA_STY_p_PO  (int n,double *flux,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, double mlim, double ewlim, struct Histdist ewd, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
  int  MLA_SWML_p_M  (int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  MLA_SWML_p_L  (int n,double *flux,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
  int  MLA_SWML_p_PO (int n,double *flux,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, double mlim, double ewlim, struct Histdist ewd, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf);
  /* Con Funcion de seleccion en flujo/magnitud y normalizacion de Poisson */

  int  MLA_STY_p_f_M   (int n,double *magn,double *z,struct fermifsel_M fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf); 
  int  MLA_STY_p_f_L   (int n,double *flux,double *z,struct fermifsel_L fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf); 
  int  MLA_STY_p_f_PO  (int n,double *magn,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, double strrad, struct cosmo_param cosmo,struct Schlf_L *lf);
  int  MLA_SWML_p_f_M  (int n,double *magn,double *z,struct fermifsel_M fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf); 
  int  MLA_SWML_p_f_L  (int n,double *flux,double *z,struct fermifsel_L fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf); 
  int  MLA_SWML_p_f_PO (int n,double *magn,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, double strrad, struct cosmo_param cosmo,struct Schlf_L *lf);

  /* Con errores en color */

  int MLA_STY_gc_p_M_wC(int n,double *magSeln, double *magDistn, double color_mean, double color_stddev, double *errColorn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo);
 
 
  /* Con errores en flujos o magnitudes*/

  int  MLA_STY_gm_M   (int n,double *magn,double *errmag ,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf); 
  int  MLA_STY_gf_L   (int n,double *flux,double *errflux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf); 
  int  MLA_STY_gm_p_M (int n,double *magn,double *errmag ,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo); 
  int  MLA_STY_gf_p_L (int n,double *flux,double *errflux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf); 
  int  MLA_SWML_gm_M  (int n,double *magn,double *errmag ,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  MLA_SWML_gf_L  (int n,double *flux,double *errflux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
  int  MLA_SWML_gm_p_M(int n,double *magn,double *errmag ,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  MLA_SWML_gf_p_L(int n,double *flux,double *errflux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);

  /* Con errores en z */

  int  MLA_STY_gz_M   (int n,double *magn,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf); 
  int  MLA_STY_gz_L   (int n,double *flux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf); 
  int  MLA_STY_gz_p_M (int n,double *magn,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf); 
  int  MLA_STY_gz_p_L (int n,double *flux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf); 
  int  MLA_SWML_gz_M  (int n,double *magn,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  MLA_SWML_gz_L  (int n,double *flux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
  int  MLA_SWML_gz_p_M(int n,double *magn,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  MLA_SWML_gz_p_L(int n,double *flux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
  
  /* Con errores en ambos */

  int  MLA_STY_gmz_M   (int n,double *magn,double *errmag ,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf); 
  int  MLA_STY_gfz_L   (int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf); 
  int  MLA_STY_gmz_p_M (int n,double *magn,double *errmag ,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf); 
  int  MLA_STY_gfz_p_L (int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf); 
  int  MLA_SWML_gmz_M  (int n,double *magn,double *errmag ,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  MLA_SWML_gfz_L  (int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);
  int  MLA_SWML_gmz_p_M(int n,double *magn,double *errmag ,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf);
  int  MLA_SWML_gfz_p_L(int n,double *flux,double *errflux,double *z, double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf);

  /* Con funciones de selección y varias imágenes */

  int  MLA_STY_s_p_f_PO  (int n,double *magn,double *ew,double *z,int *isurvey,char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyDB sdb, double ewlim, struct Histdist ewd, struct cosmo_param cosmo,struct Schlf_L *lf);
  
  /* Con funciones de selección y varias imágenes */

  int  MLA_STY_e_s_p_f_PO  (int n ,double *magn ,double *ew ,double *z ,double *ebv , double *cocnii, int *isurvey,char photband[51], float gamma, float delta, float Kcoc, float Aext, struct poselfunc fsel, struct SurveyDB sdb, double ewlim, struct Histdist ewd, struct cosmo_param cosmo,struct Schlf_L *lf);

  /* Simplemente de densidad */

  int  MLA_rho(int n,double *ra, double *dec, int *ima, struct SurveyDB sdb, double *probdetec, double *dens, double *errdens);

  /* Con casi todo lo que puedas imaginar (dabreu) ;) */
  int  MLA_STY_gmz_p_f_M_wC(int n,double *magSeln, double *magDistn, double *errmagDistn, double color_mean, double color_stddev, double *z, double *errz, struct fermifsel_M fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo, struct Schlf_M *lf, struct MLProcessInfo *mlinfo);

  
  int FileNLin(char *file);


  void Plac2Ecu(double x, double y, double ac, double dc, double *a, double *d);
  void Ecu2Plac(double a, double d, double ac, double dc, double *x, double *y);

  void PlateCtes_3(int ndat, float *xdat, float *ydat, float *alfa,float *delta,float alfac,float deltac,struct plt_cte *constantes);
  int apmread(
	      double  cra,            /* Search center J2000 right ascension in degrees */
	      double  cdec,           /* Search center J2000 declination in degrees */
	      double  dra,            /* Search half width in right ascension in degrees */
	      double  ddec,           /* Search half-width in declination in degrees */
	      double  drad,           /* Limiting separation in degrees (ignore if 0) */
	      int     sysout,         /* Search coordinate system */
	      double  eqout,          /* Search coordinate equinox */
	      double  epout,          /* Proper motion epoch (0.0 for no proper motion) */
	      double  mag1,
	      double  mag2,      /* Limiting magnitudes (none if equal) */
	      int     classd,         /* Desired object class (-1=all, 0=stars, 3=nonstars) */
	      int     nstarmax,       /* Maximum number of stars to be returned */
	      double  *gnum,          /* Array of APM numbers (returned) */
	      double  *gra,           /* Array of right ascensions (returned) */
	      double  *gdec,          /* Array of declinations (returned) */
	      double  *gmag,          /* Array of magnitudes (returned) */
	      int     *gtype,         /* Array of object types (returned) */
	      int     nlog);           /* 1 for diagnostics */
  
  int apmread_file(
		   char    file[51],
		   double  cra,       /* Search center J2000 right ascension in degrees */
		   double  cdec,      /* Search center J2000 declination in degrees */
		   double  dra,       /* Search half width in right ascension in degrees */
		   double  ddec,      /* Search half-width in declination in degrees */
		   double  drad,      /* Limiting separation in degrees (ignore if 0) */
		   int     sysout,    /* Search coordinate system */
		   double  eqout,     /* Search coordinate equinox */
		   double  epout,     /* Proper motion epoch (0.0 for no proper motion) */
		   double  mag1,
		   double  mag2,      /* Limiting magnitudes (none if equal) */
		   int     classd,    /* Desired object class (-1=all, 0=stars, 3=nonstars) */
		   int     nstarmax,  /* Maximum number of stars to be returned */
		   double  *gnum,     /* Array of APM numbers (returned) */
		   double  *gra,      /* Array of right ascensions (returned) */
		   double  *gdec,     /* Array of declinations (returned) */
		   double  *gmag,     /* Array of magnitudes (returned) */
		   int     *gtype,    /* Array of object types (returned) */
		   int     nlog);     /* 1 for diagnostics */
  
  
  int openparfile(fitsfile **fptr,  /* O - FITS file pointer                 */
	      const char *filename, /* I - name of file to open              */
	      int readwrite,        /* I - 0 = open readonly; 1 = read/write */
	      int *status);         /* IO - error status                     */
  
  /* Para leer del DSS */
 int retrieveDSS(struct wcsimage *dssimage, double ra, double dec, float xsize, float ysize); 

  
  int f_kfls(FILE *fp,char *pclave, char *var);
  
  /*  */
  /* 	util.h */
  
  /* 	Version mayo	1996		Oscar Alonso */
  /* 					Depto. Astrofisica */
  /* 					Facultad CC Fisicas */
  
  /* --------------------------------------------------------------------------- */
  
  /*   Utilidades variadas */
  
  /* --------------------------------------------------------------------------- */
  /* */ 
  
  
  struct tabla {
    int   ndat;
    float *x;
    float *y;
  };
  
  
  /*  Utilidades para lectura por teclado */
  /*  ----------------------------------- */
  int LeeInt(char *texto, int def);
  float LeeFloat(char *texto, float def);
  
  
  /*  Utilidades de numeros aleatorios */
  /*  -------------------------------- */
  void RNDInit(void);
  float RNDGauss(float s);
  int RNDPoiss(int l);
  
  /*  Utilidades para lectura de ficheros */
  /*  ----------------------------------- */
  int FileSize(char *file);
  int FileLLin(FILE *fp, int ncol, int nin, int *colin, float *colval);
  int FileSkL(FILE *fp, int nlin);
  FILE *FileLog(char *logfile, char *none);
  
  /*  Utilidades para trabajo la estructura tabla */
  /*  ------------------------------------------- */
  int  TabDimen(int n, struct tabla *t);
  int  TabRFile(char *file, int skip, struct tabla *t);
  void TabFree(struct tabla *t);
  
  /*  Interpolacion e integracion de tablas */
  /*  ------------------------------------- */
  int Locate(int n, float *x, float xx);
  float Interp1(int n, float *x, float *y, float xx);
  float Simpson(int n, float *x, float *y, float xmin, float xmax);
  float *Spline(int n, float *x, float *y, float yp1, float ypn);
  float SplInt(int n, float *xa, float *ya, float *y2a, float x);
  
  /*  Procesado de imagenes */
  /*  --------------------- */
  int   I_Intp1p(int nx, int ny, float *ima, float pixx, float pixy, 
		 float *value);
  int I_Shift(float *ima, int nx, int ny, float shiftx, float shifty);
  float *I_Binn(float *ima,int nx, int ny, int bin, int *nxbin, int *nybin);
  int I_FlipV(float *ima, int nx, int ny);
  int I_FlipH(float *ima, int nx, int ny);
  
  
  /*  Utilidades de PGPLOT */
  /*  -------------------- */
  void cpgelip(float xc,float yc,float a,float b,float t);
  void cpgqvpXY(float lxmin,  float lxmax,  float lymin,  float lymax,
		float xmin,   float xmax,   float ymin,   float ymax,
		float *vpxmin,float *vpxmax,float *vpymin,float *vpymax);
  void cpgsvpXY(float lxmin,  float lxmax,  float lymin,  float lymax,
		float xmin,   float xmax,   float ymin,   float ymax);
  void pgLimits_d(int n, double *x, float *wmin, float *wmax);
  void pgLimits(int n, float *x, float *wmin, float *wmax);
  int pgLimitY(int n, float *x, float *y, float xmin, float xmax,
	       float *wymin, float *wymax);
  void cpghist3d(float *ima, int nx, int ny,  int i1, int i2, int j1, int j2);
  int cpgdraw3d(float *arrayen, int kx, int ny, float size, float angle);
  int cpgdodraw3d(double *arrayen, int kx, int ny, float size, float angle);
  void cpgline_d(int n, const double *xpts, const double *ypts);
  void cpggray_d(const double *a, int idim, int jdim, int i1, int i2, int j1, int j2, double fg, double bg, const double *tr);
  void cpgcons_d(const double *a, int idim, int jdim, int i1, int i2, int j1, int j2, const double *c, int nc, const double *tr);
  void cpghist_d(int n, const double *data, double datmin, double datmax, int nbin, int pgflag);
  void cpgpt_d(int n, const double *xpts, const double *ypts, int symbol);
  void cpgerry_d(int n, const double *x, const double *y1, const double *y2, double t);
  int cpgcurs_d(double *x, double *y, char *ch_scalar); 
  int cpgband_d(int mode, int posn, double xref, double yref, double *x, double *y, char *ch_scalar);
  void cpgpoly_d(int n, const double *xpts, const double *ypts);

  
  /*  Utilidades de Estadistica */
  /*  ------------------------- */
  float StSuma1(int n, float *a, int i);
  float StSuma2(int n, float *a, int i, float *b, int j);
  float StSuma3(int n, float *a, int i, float *b, int j,float *c, int k);
  float StMedia(int n, float *x, float *sigma);
  float StWeightMedia(int n, float *x, float *w, float *sigma);
  float StErrWeightMedia(int n, float *x, float *err, float *sigma);
  float StModa(int n, float *a, int nbin, float *amin, float *amax);
  int  *StHisto(int n, float *a, int nbin, float *amin, float *amax);
  int  *StHisto2(int n, float *a, int nbin, float *amin, float *amax);
  int  *StCumHisto(int n, float *a, int nbin, float *amin, float *amax);
  int  *StCumHisto_d(int n, double *a, int nbin, double *amin, double *amax);
  int  **StHisto2D_d(int n, double *x, double *y, int nbinx, double *xmin, double *xmax, int nbiny, double *ymin, double *ymax);
  float StCovar(int n, float *a, float *b, float *sigmaa, float *sigmab);
  double StCovar_d(int n, double *a, double *b, double *sigmaa, double *sigmab);
  
  /* Para numeres double */
  /*---------------------*/
  double StSuma1_d(int n, double *a, int i);
  double StSuma2_d(int n, double *a, int i, double *b, int j);
  double StSuma3_d(int n, double *a, int i, double *b, int j,double *c, int k);
  double StMedia_d(int n, double *x, double *sigma);
  double StWeightMedia_d(int n, double *x, double *w, double *sigma);
  double StErrWeightMedia_d(int n, double *x, double *err, double *sigma);
  double StModa_d(int n, double *a, int nbin, double *amin, double *amax);
  int  *StHisto_d(int n, double *a, int nbin, double *amin, double *amax);
  int  *StHisto2_d(int n, double *a, int nbin, double *amin, double *amax);
  
  /*  Utilidades de ajustes por Minimos Cuadrados */
  /*  ------------------------------------------- */
  float MCP1(int n, float *x, float *y, float *a, float *b);
  double MCP1_d(int n, double *x, double *y, double *a, double *b);
  void  MCLine(int n, float *x, float *y, float *a, float *erra);
  void  MCP1Weight_err(int n, float *x, float *y, float *sig,float *a, float *b, float *erra, float *errb,float *covab, float *chi2, float *q);
  void  MCP1Weight_err_d(int n, double *x, double *y, double *sig,double *a, double *b, double *erra, double *errb,double *covab, double *chi2, float *q);
  void  MCP1_err(int n, float *x, float *y, float *a, float *b, float *erra, float *errb,float *covab, float *chi2, float *q);
  void  MCP1_err_d(int n, double *x, double *y,double *a, double *b, double *erra, double *errb,double *covab, double *chi2, float *q);
  float MCPN(int n, float *x, float *y, int g, float *c);
  float MCPm1x(int n, float *x, float *y, float *p, float *q, float *x0);
  float MCPm1(int n, float *x, float *y, float *p, float *q);
  float MCPm2(int n, float *x, float *y, float *p, float *q);
  float MCPm4(int n, float *x, float *y, float *p, float *q, float *s);
  int   MCElip(int n, float *x, float *y,
	       float *a, float *b, float *c, float *d, float *f, float *e);
  int   MCElip_d(int n, double *x, double *y,
	      double *a, double *b, double *c, double *d, double *f, double *e);
  int  MCElipN(int npt, int ndim, float **x, float **C);
  int  MCElipN_d(int npt, int ndim, double **x, double  **C);
  struct Pol_deg2 {
     float a;
     float b;
     float c;
  };
  void FitPol_deg2(float x1, float y1,float x2,float y2,float x3,float y3,struct Pol_deg2 *pol);
  int ResolvEcu_deg2(struct Pol_deg2 pol, float *x1, float *x2);

  
  void svdfit( float *x, float *y, float *sig, int n, float *a, int ma, float **cvm, float *chisq,void (*funcs)(float,float [],int));
  void svdvar(float **v,int ma, float w[],float **cvm);
  void svdcmp(float **a,int m,int n,float w[],float **v); 
  float pythag(float a,float b); 
  void svbksb(float **u,float w[],float **v,int m,int n,float b[],float x[]);
  void svdfit_d( double *x, double *y, double *sig, int n, double *a, int ma, double **cvm, double *chisq,void (*funcs)(double,double [],int));
  void svdvar_d(double **v,int ma, double w[],double **cvm);
  void svdcmp_d(double **a,int m,int n,double w[],double **v); 
  double pythag_d(double a,double b); 
  void svbksb_d(double **u,double w[],double **v,int m,int n,double b[],double x[]);
  
  int Mrq(float  x[],float y[],float sig[],int ndata,float a[],int ia[],int ma, float **covar, float *chisq, void (*funcs)(float, float [], float *, float [],int));
  void mrqmin(float  x[],float y[],float sig[],int ndata,float a[],int ia[],int ma, float **covar,float **alpha, float *chisq, void (*funcs)(float, float [], float *, float [],int), float *alamda);
  void mrqcof(float x[],float y[],float sig[],int ndata,float a[],int ia[],int ma,float **alpha,float beta[],float *chisq,void (*funcs)(float ,float [],float *,float [],int));
  void covsrt(float **covar,int ma,int ia[],int mfit);
  
  
  int Mrq_d(double  x[],double y[],double sig[],int ndata,double a[],int ia[],int ma, double **covar, double *chisq, void (*funcs)(double, double [], double *, double [],int));
  void mrqmin_d(double  x[],double y[],double sig[],int ndata,double a[],int ia[],int ma, double **covar,double **alpha, double *chisq, void (*funcs)(double, double [], double *, double [],int), double *alamda);
  void mrqcof_d(double x[],double y[],double sig[],int ndata,double a[],int ia[],int ma,double **alpha,double beta[],double *chisq,void (*funcs)(double ,double [],double *,double [],int));
  void covsrt_d(double **covar,int ma,int ia[],int mfit);
  
   
  
  /*  Utilidades de Vectores */
  /*  ---------------------- */
  float VDist(int nc, float *a, float *b);
  void  VSuma(int nc, float *a, float *b, float *c);
  float VSumaCmp(int nc, float *a);
  float VModulo(int nc, float *a);
  void  VNorma1(int nc, float *a);
  void  VNormaIC(int nc, float *a);
  void  VPrintf(FILE *stream, int nc, float *a, char *ident);
  float VProdEsc(int nc, float *a, float *b);
  
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
  
  /* Genetic Algorithms */
  /* ------------------ */

  struct ga_set {
    unsigned int npop;
    unsigned int ngeneration;
    unsigned int ncod;
    double pcross;
    int imut;
    double pmut;
    double pmutmin;
    double pmutmax;
    double fdif;
    int irep;
    int ielite;
    int verbose;
  };

  void ga_setdefault(struct ga_set *gs);
  void ga_rnkpop(unsigned int n, double *arrin, unsigned int  *indx, unsigned int *rank);
  void ga_select(unsigned int np, unsigned int *jfit, double fdif, int *idad);
  void ga_adjmut(unsigned int np, double *fitns, unsigned int *ifit, double pmutmn, double pmutmx, double *pmut);
  void ga_encode(unsigned int ndim, unsigned int ncod, double *ph, int *gen);
  void ga_decode(unsigned int ndim, unsigned int ncod, double *ph, int *gen);
  void ga_mutate(unsigned int ndim, unsigned int ncod, double pmut, int *gen);
  void ga_cross(unsigned int n, unsigned int nd, double pcross, int *gen1, int *gen2);
  void ga_genrep(unsigned int n , unsigned int ip, double **ph, double **newph);
  void ga_newpop(double (*gafunc)(int, double *), double *pscale, double *pzero, int ielite, unsigned int ndim, unsigned int np, double **oldph, double **newph, unsigned int *ifit, unsigned int *jfit, double *fitns, unsigned int *nnew);

  
  /*  Varios */
  /*  ------ */
  void SlpTime(FILE *stream);
  int Baraja(int n, int *a);
  void MinMax(int n, float *a, float *min, float *max);
  void MinMax_d(int n, double *a, double *min, double *max);
  void MinMax_l(int n, long  *a, long  *min, long  *max);

  int   SELGauss(int n, float *a, float *c, float *x);
  int   SELGauss_d(int n, double *a, double *c, double *x);
  int gaussj(float **a, int n, float **b, int m);
  int gaussj_d(double **a, int n, double **b, int m);
  int   TrSCoo(int n, float *a, float *b, float *x, float *y, float *coeff,int g,
	       float *sigx, float *sigy);
  int ElipPar(float a, float b, float c, float d, float f, float e,
	      float *x0, float *y0, float *semax, float *semin,
	      float *t);
  int ElipPar_d(double a, double b, double c, double d, double f, double e,
	      double *x0, double *y0, double *semax, double *semin,
	      double *t);
  int TestDiv0(float a, float b, float *c, float tol);			   
  int TestDiv0_d(double a, double b, double *c, double tol);			   
  
  /* 	lectorfits.h          */
  
  /* 	Version febrero 1996			Dpto. Astrofisica (UCM) */
  
  /* 						Juan Carlos Vallejo */
  /* 						Oscar Alonso */
  
  /* ---------------------------------------------------------------------------/ */
  
  /*   Ayuda y definiciones de las funciones para la lectura de imagnes FITS */
  
  /* --------------------------------------------------------------------------- */
  
  
#define MAX_KEYLAB  10
  
  struct headfitsmosaic {         /* KEYWORDS FITS MOSAICOS                    */
    int    nmosx;           /* No. de objetos del mosaico en X           */
    int    nmosy;           /* No. de objetos del mosaico en Y           */
    int    cajax;           /* Tamano de la caja X de los objetos (pix)  */
    int    cajay;           /* Tamano de la caja Y de los objetos (pix)  */
  };
  
  struct headfits {                /* KEYWORDS STANDARD FITS                   */
    char   simple_c;           /* Formato fits ?                           */
    int    simple;           /* Formato fits ?                           */
    int    bitpix;           /* No. de bits por pixel                    */
    int    naxis;            /* No. de ejes de la imagen                 */
    int    naxis1;           /* No. de pixels del eje X                  */
    int    naxis2;           /* No. de pixels del eje Y                  */
    float  crval1;           /* Coordenada WORLD X del primer pixel      */
    float  crval2;           /* Coordenada WORLD Y del primer pixel      */
    float  cdelt1;           /* Tamano fisico del pixel X (incremento)   */
    float  cdelt2;           /* Tamano fisico del pixel Y (incremento)   */
    char   ctype1[51];       /* Tipo de coordenadas X                    */
    char   ctype2[51];       /* Tipo de coordenadas Y                    */
    float  bscale;           /* Z = Zfits * bscale + bzero               */
    float  bzero;            /*                                          */
    float  datamax;          /* Valor maximo de la imagen                */
    float  datamin;          /* Valor minimo de la imagen                */
    char   bunit[51];        /* Unidades de los pixels                   */
    char   object[51];       /* Objeto                                   */
    char   filename[51];     /* Origen de la imagen                      */
    char   origin[51];       /* Origen de la imagen                      */
    char   date[51];         /* Fecha                                    */
    char   instrume[51];     /* Instrumento                              */
    /* Lineas de Historia y de Comentatrio      */
    char   history[MAX_KEYLAB][71];
    char   comment[MAX_KEYLAB][71];
    struct headfitsmosaic    /* Estructura con datos si la imagen es     */
    mosaic;       /* un mosaico (INSTRUME = MOSAICO)          */
    
  };
  
  
  /* --------------------------------------------------------------------------- */
  int     fitsbhn   (char *ima);
  int     fitsrh    (char *ima,struct headfits *h);
  int     fitsinfo  (char *ima,FILE *stream);
  float   *fitsrf   (char *ima,struct headfits *h);
  float   *fitsrfb  (char *ima,struct headfits *h,
		     int *x1, int *x2, int *y1, int *y2);
  float   *fitsrfm  (char *name,struct headfits *h, int nmos);
  float   *fitsrmm  (char *comodin, struct headfits *h, int nobj);
  int     fitswf    (char *name,struct headfits *h,float *ima,char autoescala);
  /* ------------------------------------------------------------------------- */
  int     f_bhn     (FILE *fp);
  int     f_fill    (FILE *fp);
  int     f_rh      (FILE *fp,struct headfits *h);
  void    f_rf      (FILE *fp,struct headfits *h,float *ima);
  void    f_rfb     (FILE *fp,struct headfits *h, float *ima,
		     int x1, int x2, int y1, int y2);
  int     f_testh   (struct headfits *h);
  void    f_wflin   (FILE *fp,struct headfits *h,float *linea);
  /* ------------------------------------------------------------------------- */
  int     fitskfi  (char *ima,char *pclave, int    *var);
  int     fitskff  (char *ima,char *pclave, float  *var);
  int     fitskfd  (char *ima,char *pclave, double *var);
  int     fitskfs  (char *ima,char *pclave, char   *var);
  int     fitskfc  (char *ima,char *pclave, char   *var);
  int     fitskfl  (char *ima,char *pclave, char label[][71]);
  
  int     f_kfi    (FILE *fp,char *pclave, int    *var);
  int     f_kff    (FILE *fp,char *pclave, float  *var);
  int     f_kfd    (FILE *fp,char *pclave, double *var);
  int     f_kfs    (FILE *fp,char *pclave, char   *var);
  int     f_kfc    (FILE *fp,char *pclave, char   *var);
  int     f_kfl    (FILE *fp,char *pclave, char label[][71]);
  int     f_actlab (char label[][71], char *var);
  void    f_clrh   (struct headfits *h);
  void    f_clrl   (char label[][71]);
  
  
  int     fitskmi  (char *ima,char *pclave,int   var);
  int     fitskmf  (char *ima,char *pclave,float var);
  int     fitskmc  (char *ima,char *pclave,char  var);
  int     fitskms  (char *ima,char *pclave,char *var);
  
  int     f_kmi    (FILE *fp,char *pclave,int   var);
  int     f_kmf    (FILE *fp,char *pclave,float var);
  int     f_kmc    (FILE *fp,char *pclave,char  var);
  int     f_kms    (FILE *fp,char *pclave,char *var);
  
  int     f_kwi    (FILE *fp,char *pclave,int   var, char *comment);
  int     f_kwf    (FILE *fp,char *pclave,float var, char *comment);
  int     f_kwc    (FILE *fp,char *pclave,char  var, char *comment);
  int     f_kws    (FILE *fp,char *pclave,char *var, char *comment);
  
  int     f_wh     (FILE *fp, struct headfits *h);
  /* ------------------------------------------------------------------------- */
  int     fits_pv  (float *ima, struct headfits *h, float x, float y, float *pv);
  
  
  /************************************************************************/
  /*									*/
  /*               	SUBRUTINAS DE ASTROMETRIA			*/
  /*		   							*/
  /*    	       		    ­­­­­ EN C !!!!!				*/
  /*									*/
  /*									*/
  /*   Oscar Alonso                                         Junio 1994	*/
  /*                              Madrid					*/
  /*									*/
  /************************************************************************/
  
  
  
  float hms2r(int ah, int am, float as);
  float gms2r(char dsig, int dg, int dm, float ds);
  void r2hms(double r, int *ah, int *am, float *as);
  void r2gms(double r, char *dsig, int *dg, int *dm, float *ds);
  char *r2hmsch(float r);
  char *r2gmsch(float r);
  
  void Precesa(double *p, float ai, float di, float *af, float *df);
  double *MPrec(float datei, float datef);
  
  float Readhms(void);
  float Readgms(void);
  
  void Ecu2Cart(float ar, float dec, float *v);
  
  void stdastr(char *ar,char *dec,char *stdformat);
  void strtoast(char *c, char *ident, float *ar, float *dec);
  
#ifdef __cplusplus
}
#endif


#endif
