#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <cpgplot.h>
#include "alloc.h"
#include "amoeba.h"
#include "cosmology.h"
#include "schechter.h"
#include "functions.h"
#include "mlprocess.h"
#include "mlsty.h"
#include "mlswml.h"
#include "stmedia.h"
#include "minmax.h"
#include "readkbd.h"
#include "sthisto.h"
#include "lagr.h"
#include "step.h"
#include "random.h"
#include "vvmax.h"

#define ZMIN 0.00001
#define NSTEP_LF 200
#define NSTEP_Z  500
#define NSTEP_MAG 500
/* #define NSTEP_Z  200 */
/* #define NSTEP_MAG 200 */
#define FTOL 1.e-10
#define DEBUG 0
#define DEBUG2 0
#define USE_PGPLOT 1

/*  Las formulas de cosmologia estan cogidas del articulo de Hogg */
/*    astroph/9905116  */

struct sample_data {
  int ngalax;
  double *z;
  double *lum;
  double *mag;
  double *ew;
  char   **image;
};

/* struct to contain selection and distribution magnitudes */
struct sample_data_sel_dist 
{
  int ngalax;
  double *z;
  double *lumSel;
  double *lumDist;
  double *magSel;
  double *magDist;
  char   **image;
};

/* dabreu */
struct sample_data_mag_err 
{
  int ngalax;
  double *z;
  double *mag;
  double *mag_err;
};

struct sample_data_wC_errColor
{
  int ngalax;
  double *z;
  double *magSel;
  double *magDist;
  double *errColor;
};

struct sample_data_gmz_wC
{
  int ngalax;
  double *z;
  double *errz;
  double *magSel;
  double *magDist;
  double *errMagDist;
};

struct lf_param
{
  double mag_min;
  double mag_max;
  double lum_min;
  double lum_max;
  double zlow;
  double zup;
  double area;
  int islum;
  int npoint;
};

struct lum_func_ceg {
  struct poselfunc fsel;
  struct SurveyDB sdb;
  struct Schlf_L lf;
  struct Histdist ewd;
  double ewlim;
  char photband[51];
  double gamma,delta,Kcoc;

};


void get_sample(struct sample_data *sample);
void get_sample_sel_dist(struct sample_data_sel_dist *sample);
void get_sample_mag_err(struct sample_data_mag_err *sample);
void get_sample_wC_errColor(struct sample_data_wC_errColor *sample);
void get_sample_gmz_wC(struct sample_data_gmz_wC *sample);
void get_sample_ceg(struct sample_data *sample);
void set_lf_parameters(struct lf_param *lf);
void set_lf_ceg(struct lum_func_ceg *lf);
void set_cosmology(void);
void VVmax(void);
void STY(void);
void STY_MAG_ERR(void);
void STY_wC(void);
void STY_wC_errColor(void);
void STY_gmz_f_M_wC(void);
void SWML(void);
void CEG(void);
void Calc_Num(void);
void Calc_Num_wC(void);
void Generate_Cat_M(void);
void Generate_Cat_M_wC(void);
void Generate_Cat_M_C(void);
void Generate_Cat_L(void);
void set_Schechter_M(void);
void set_Schechter_L(void);



double c=299792.46; /* // En km/s */
/* //double alfa,phistar,Mstar; */
/* dabreu */
//const double kk=144543977.0745928; /* pow(10.,-0.4*-20.4) */
struct Schlf_M schlf_M={-1.3,0.,0.0033,0.,-20.4,0.,0.,0.,0.};
struct cosmo_param cosmo;
/* utilizamos la struct de modulos.h */
struct Schlf_L schlf_L={-1.3,0.0,0.0033,0.0,144543977.0745928,0.0,0.0,0.0,0.0};

double fglobal;

int main(void)
{
  char opt='E';

  do{
    printf("\n C Compute number of galaxies with a given luminosity function \n"); 
    printf(" D Compute number of galaxies with a given luminosity function using color distribution.\n"); 
    printf(" V Calculate luminosity function by V/Vmax method\n"); 
    
    printf(" L Calculate luminosity function by CEG method\n"); 
    printf(" M Calculate luminosity function by STY maximum likelihood method\n"); 
    /* dabreu */
    printf(" N Calculate luminosity function by STY method with gaussian errors in magnitude\n");
    printf(" O Calculate luminosity function by STY_wC method\n");
    printf(" P Calculate luminosity function by STY_wC with gaussian errors in color\n");
    printf(" R Calculate luminosity function by STY_wC with gaussian errors in mag and z\n");
    printf(" W Calculate luminosity function by SWML maximum likelihood method\n"); 
/*    printf(" O Calculate luminosity function by C- method\n");  */
    printf(" G Generate a random catalogue for a given LF by Montecarlo simulations in magnitudes\n");
    printf(" H Generate a random catalogue for a given LF by Montecarlo simulations in luminosities\n");
    printf(" I Generate a random catalogue for a given LF by Montecarlo simulations in magnitudes using color distribution.\n");
    printf(" E Exit\n");
    opt=readc(opt);
    switch (opt) { 
    case 'C':
    case 'c':
      set_cosmology();
      Calc_Num();
      cosmo_free(&cosmo);
      break;
    case 'D':
    case 'd':
      set_cosmology();
      Calc_Num_wC();
      cosmo_free(&cosmo);
      break;
    case 'V':
    case 'v':
      set_cosmology();
      VVmax();
      cosmo_free(&cosmo);
      break;
    case 'M':
    case 'm':
      set_cosmology();
      STY();
      cosmo_free(&cosmo);
      break;
    case 'N':
    case 'n':
      set_cosmology();
      STY_MAG_ERR();
      cosmo_free(&cosmo);
      break;
    case 'O':
    case 'o':
      set_cosmology();
      STY_wC();
      cosmo_free(&cosmo);
      break;
    case 'P':
    case 'p':
      set_cosmology();
      STY_wC_errColor();
      cosmo_free(&cosmo);
      break;
    case 'R':
    case 'r':
      set_cosmology();
      STY_gmz_f_M_wC();
      cosmo_free(&cosmo);
      break;
    case 'W':
    case 'w':
      set_cosmology();
      SWML();
      cosmo_free(&cosmo);
      break;
    case 'L':
    case 'l':
      set_cosmology();
      CEG();
      cosmo_free(&cosmo);
      break;
    case 'G':
    case 'g':
      Generate_Cat_M();
      break;
    case 'H':
    case 'h':
      Generate_Cat_L();
      break;
    case 'I':
    case 'i':
      Generate_Cat_M_wC();
      break;
    }

  }while(opt!='E' && opt!='e');

  return(0);
}

void STY(void)
{
  static struct lf_param sty;
  struct sample_data_sel_dist sample;
  static double mlim=0;
  static double flim=0;
  double zlow;

  int status;
  
  struct Schlf_M lfsch_M;
  struct Schlf_L lfsch_L;

  static int poissonflag=0;

  /* dabreu */
  static int plots=0; /* para no hacer las gr�ficas */
  FILE *fout;
  static char resultFileName[300]="";
  
  /* Information about the ML process */
  struct MLProcessInfo mlprocess;
  
  get_sample_sel_dist(&sample);  
  set_lf_parameters(&sty);
  zlow=sty.zlow;
  zlow = (zlow < ZMIN ? ZMIN : zlow);
  
  printf(" The following question allows to use an improvement of STY that allows to \ncompute");
  printf(" normalization using Maximum Likelihood. It assumes that number of \ndetected sources");
  printf(" follows a Poisson event. It takes into account Poisson \nerrors both in Mstar, alfa");
  printf(" and Phi_star. These errors are used to be \nlarger than for traditional method, wich");
  printf(" does not take into account the \ninfluence of Poisson detections in Mstar, alfa and Phi_star\n");
  printf(" Use modified method to compute normalization at ML time (0=no/1=yes)? \n");
  poissonflag=readi(poissonflag);
 
  if(sty.islum) 
  {
    printf(" Input the limiting flux: ");
    flim=readd(flim);
    if(poissonflag) status=MLA_STY_p_L(sample.ngalax,sample.lumDist,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    else            status=  MLA_STY_L(sample.ngalax,sample.lumDist,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    printf("\n Solutioner exited with error status %d\n",status);
    printf("\n   Solution found at iteration %d\n",mlprocess.nIter);
    printf("\n   Likelihood function: %g\n",mlprocess.MLmax);
    printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lfsch_L.Lstar),lfsch_L.alfa,lfsch_L.phistar,log10(lfsch_L.phistar));
    printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_L.errLstar/lfsch_L.Lstar/log(10.),lfsch_L.erralfa,lfsch_L.errphistar,lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g Covar(log(Lstar),alpha): %g\n",lfsch_L.covaralfaLstar,lfsch_L.covaralfaLstar/lfsch_L.Lstar/log(10.));
    printf(" Covar(Lstar,Phistar): %g Covar(log(Lstar),log(Phistar)): %g\n",lfsch_L.covarphistarLstar,lfsch_L.covarphistarLstar/lfsch_L.Lstar/log(10.)/lfsch_L.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_L.covaralfaphistar,lfsch_L.covaralfaphistar/lfsch_L.phistar/log(10.));
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD logLstar alpha Phistar log\n");
    printf("#FL_HEAD E_logLstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",log10(lfsch_L.Lstar),lfsch_L.alfa,lfsch_L.phistar,log10(lfsch_L.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_L.errLstar/lfsch_L.Lstar/log(10.),lfsch_L.erralfa,lfsch_L.errphistar,lfsch_L.errphistar/lfsch_L.phistar/log(10.));

    /* dabreu */
    /* fichero para los resultados */
    printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 LOG_L_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_LOG_L_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "# 9 ML_MAX\n");
    fprintf(fout, "# 10 N_ITER\n");
    fprintf(fout, "%g\t", log10(lfsch_L.Lstar));
    fprintf(fout, "%g\t", lfsch_L.alfa);
    fprintf(fout, "%g\t", lfsch_L.phistar);
    fprintf(fout, "%g\t", log10(lfsch_L.phistar));
    fprintf(fout, "%g\t", lfsch_L.errLstar/lfsch_L.Lstar/log(10.));
    fprintf(fout, "%g\t", lfsch_L.erralfa);
    fprintf(fout, "%g\t", lfsch_L.errphistar);
    fprintf(fout, "%g\t", lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    fprintf(fout, "%g\t", mlprocess.MLmax);
    fprintf(fout, "%i\n", mlprocess.nIter);
    fclose(fout);
  }
  else 
  {
    printf(" Input the limiting magnitude: ");
    mlim=readd(mlim);
    if(poissonflag) status=MLA_STY_p_M(sample.ngalax,sample.magSel,sample.magDist,sample.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess);
    else            status=  MLA_STY_M(sample.ngalax,sample.magDist,sample.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess);
    printf("\n Solutioner exited with error status %d\n",status);
    printf("\n   Solution found at iteration %d\n",mlprocess.nIter);
    printf("\n   Likelihood function: %g\n",mlprocess.MLmax);
    printf(" Mstar:  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf(" E_Mstar:    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Lstar,Phistar): %g Covar(Lstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD Mstar alpha Phistar log\n");
    printf("#FL_HEAD E_Mstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
/*    printf(" Covar(Lstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Lstar,Phistar): %g Covar(Lstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    */

    /* dabreu */
    /* fichero para los resultados */
    printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 M_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_M_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "# 9 COVAR_ALPHA_M_STAR\n");
    fprintf(fout, "# 10 COVAR_ALPHA_PHISTAR\n");
    fprintf(fout, "# 11 COVAR_ALPHA_LOGPHISTAR\n");
    fprintf(fout, "# 12 COVAR_PHISTAR_M_STAR\n");
    fprintf(fout, "# 13 COVAR_LOGPHISTAR_M_STAR\n");
    fprintf(fout, "# 14 ML_MAX\n");
    fprintf(fout, "# 15 N_ITER\n");
    fprintf(fout, "%g\t", lfsch_M.Mstar);
    fprintf(fout, "%g\t", lfsch_M.alfa);
    fprintf(fout, "%g\t", lfsch_M.phistar);
    fprintf(fout, "%g\t", log10(lfsch_M.phistar));
    fprintf(fout, "%g\t", lfsch_M.errMstar);
    fprintf(fout, "%g\t", lfsch_M.erralfa);
    fprintf(fout, "%g\t", lfsch_M.errphistar);
    fprintf(fout, "%g\t", lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", lfsch_M.covaralfaMstar);
    fprintf(fout, "%g\t", lfsch_M.covaralfaphistar);
    fprintf(fout, "%g\t", lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", lfsch_M.covarphistarMstar);
    fprintf(fout, "%g\t", lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", mlprocess.MLmax);
    fprintf(fout, "%i\n", mlprocess.nIter);
    fclose(fout);
  }

  /* dabreu */
  printf(" Do you want plots of the LF (1=yes, 0=no)?\n");
  plots=readf(plots);
  if(plots) {
     cpgopen("?");
     cpgask(0);
     cpgclos();
  }
}

void STY_wC(void)
{
  static struct lf_param sty;
  struct sample_data_sel_dist sample;
  static double mlim=0;
  double zlow;

  double color_mean=0.;
  double color_stddev=0.;

  int status=0;
  
  struct Schlf_M lfsch_M;

  static int poissonflag=0;

  static int plots=0; /* para no hacer las gr�ficas */
  FILE *fout;
  static char resultFileName[200]="";
  
  /* Information about the ML process */
  struct MLProcessInfo mlprocess;
  
  get_sample_sel_dist(&sample);
  printf("Input the color mean:\n");
  color_mean=readd(color_mean);
  printf("Input the color stddev:\n");
  color_stddev=readd(color_stddev);
  set_lf_parameters(&sty);
  zlow=sty.zlow;
  zlow = (zlow < ZMIN ? ZMIN : zlow);
  
  printf(" The following question allows to use an improvement of STY that allows to \ncompute");
  printf(" normalization using Maximum Likelihood. It assumes that number of \ndetected sources");
  printf(" follows a Poisson event. It takes into account Poisson \nerrors both in Mstar, alfa");
  printf(" and Phi_star. These errors are used to be \nlarger than for traditional method, wich");
  printf(" does not take into account the \ninfluence of Poisson detections in Mstar, alfa and Phi_star\n");
  printf(" Use modified method to compute normalization at ML time (0=no/1=yes)? \n");
  poissonflag=readi(poissonflag);
 
  if(sty.islum) 
  {
    printf(" This method is not (yet) available for luminosities (only for magnitudes).\n");
    /* printf(" Input the limiting flux: ");
    flim=readd(flim);
    if(poissonflag) status=MLA_STY_p_L(sample.ngalax,sample.lum,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    else            status=  MLA_STY_L(sample.ngalax,sample.lum,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    printf("\n Solutioner exited with error status %d\n",status);
    printf("\n   Solution found at iteration %d\n",mlprocess.nIter);
    printf("\n   Likelihood function: %g\n",mlprocess.MLmax);
    printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lfsch_L.Lstar),lfsch_L.alfa,lfsch_L.phistar,log10(lfsch_L.phistar));
    printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_L.errLstar/lfsch_L.Lstar/log(10.),lfsch_L.erralfa,lfsch_L.errphistar,lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g Covar(log(Lstar),alpha): %g\n",lfsch_L.covaralfaLstar,lfsch_L.covaralfaLstar/lfsch_L.Lstar/log(10.));
    printf(" Covar(Lstar,Phistar): %g Covar(log(Lstar),log(Phistar)): %g\n",lfsch_L.covarphistarLstar,lfsch_L.covarphistarLstar/lfsch_L.Lstar/log(10.)/lfsch_L.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_L.covaralfaphistar,lfsch_L.covaralfaphistar/lfsch_L.phistar/log(10.));
    printf("##################################################\n");
    printf("#FL_HEAD logLstar alpha Phistar log\n");
    printf("#FL_HEAD E_logLstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",log10(lfsch_L.Lstar),lfsch_L.alfa,lfsch_L.phistar,log10(lfsch_L.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_L.errLstar/lfsch_L.Lstar/log(10.),lfsch_L.erralfa,lfsch_L.errphistar,lfsch_L.errphistar/lfsch_L.phistar/log(10.)); */

    /* dabreu */
    /* fichero para los resultados */
    /* printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 LOG_L_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_LOG_L_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "%g\t", log10(lfsch_L.Lstar));
    fprintf(fout, "%g\t", lfsch_L.alfa);
    fprintf(fout, "%g\t", lfsch_L.phistar);
    fprintf(fout, "%g\t", log10(lfsch_L.phistar));
    fprintf(fout, "%g\t", lfsch_L.errLstar/lfsch_L.Lstar/log(10.));
    fprintf(fout, "%g\t", lfsch_L.erralfa);
    fprintf(fout, "%g\t", lfsch_L.errphistar);
    fprintf(fout, "%g\n", lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    fclose(fout); */
  }
  else 
  {
    printf(" Input the limiting magnitude: ");
    mlim=readd(mlim);
    if(poissonflag) status=MLA_STY_p_M_wC(sample.ngalax,sample.magSel,sample.magDist,color_mean,color_stddev,sample.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess);
    /* else            status=  MLA_STY_M(sample.ngalax,sample.mag,sample.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess); */
    else printf("This method is not (yet) available (only poisson)\n");
    printf("\n Solutioner exited with error status %d\n",status);
    printf("\n   Solution found at iteration %d\n",mlprocess.nIter);
    printf("\n   Likelihood function: %g\n",mlprocess.MLmax);
    printf(" Mstar:  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf(" E_Mstar:    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    printf(" Covar(Mstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Mstar,Phistar): %g Covar(Mstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD Mstar alpha Phistar log\n");
    printf("#FL_HEAD E_Mstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
/*    printf(" Covar(Lstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Lstar,Phistar): %g Covar(Lstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    */

    /* fichero para los resultados */
    printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 M_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_M_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "# 9 COVAR_ALPHA_M_STAR\n");
    fprintf(fout, "# 10 COVAR_ALPHA_PHISTAR\n");
    fprintf(fout, "# 11 COVAR_ALPHA_LOGPHISTAR\n");
    fprintf(fout, "# 12 COVAR_PHISTAR_M_STAR\n");
    fprintf(fout, "# 13 COVAR_LOGPHISTAR_M_STAR\n");
    fprintf(fout, "# 14 ML_MAX\n");
    fprintf(fout, "# 15 N_ITER\n");
    fprintf(fout, "%g\t", lfsch_M.Mstar);
    fprintf(fout, "%g\t", lfsch_M.alfa);
    fprintf(fout, "%g\t", lfsch_M.phistar);
    fprintf(fout, "%g\t", log10(lfsch_M.phistar));
    fprintf(fout, "%g\t", lfsch_M.errMstar);
    fprintf(fout, "%g\t", lfsch_M.erralfa);
    fprintf(fout, "%g\t", lfsch_M.errphistar);
    fprintf(fout, "%g\t", lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", lfsch_M.covaralfaMstar);
    fprintf(fout, "%g\t", lfsch_M.covaralfaphistar);
    fprintf(fout, "%g\t", lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", lfsch_M.covarphistarMstar);
    fprintf(fout, "%g\t", lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", mlprocess.MLmax);
    fprintf(fout, "%i\n", mlprocess.nIter);
    fclose(fout);
  }

  printf(" Do you want plots of the LF (1=yes, 0=no)?\n");
  plots=readf(plots);
  if(plots) {
     cpgopen("?");
     cpgask(0);
     cpgclos();
  }
}

/* Estamos trabajando en ellou */
void STY_wC_errColor(void)
{
  static struct lf_param sty;
  struct sample_data_wC_errColor sample;
  static double mlim=0;
  double zlow;

  double color_mean=0.;
  double color_stddev=0.;

  int status=0;
  
  struct Schlf_M lfsch_M;

  static int poissonflag=0;

  static int plots=0; /* para no hacer las gr�ficas */
  FILE *fout;
  static char resultFileName[200]="";
  
  /* Information about the ML process */
  struct MLProcessInfo mlprocess;
  
  get_sample_wC_errColor(&sample);
  printf("Input the color mean:\n");
  color_mean=readd(color_mean);
  printf("Input the color stddev:\n");
  color_stddev=readd(color_stddev);
  set_lf_parameters(&sty);
  zlow=sty.zlow;
  zlow = (zlow < ZMIN ? ZMIN : zlow);
  
  printf(" The following question allows to use an improvement of STY that allows to \ncompute");
  printf(" normalization using Maximum Likelihood. It assumes that number of \ndetected sources");
  printf(" follows a Poisson event. It takes into account Poisson \nerrors both in Mstar, alfa");
  printf(" and Phi_star. These errors are used to be \nlarger than for traditional method, wich");
  printf(" does not take into account the \ninfluence of Poisson detections in Mstar, alfa and Phi_star\n");
  printf(" Use modified method to compute normalization at ML time (0=no/1=yes)? \n");
  poissonflag=readi(poissonflag);
 
  if(sty.islum) 
  {
    printf(" This method is not (yet) available for luminosities (only for magnitudes).\n");
    /* printf(" Input the limiting flux: ");
    flim=readd(flim);
    if(poissonflag) status=MLA_STY_p_L(sample.ngalax,sample.lum,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    else            status=  MLA_STY_L(sample.ngalax,sample.lum,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    printf("\n Solutioner exited with error status %d\n",status);
    printf("\n   Solution found at iteration %d\n",mlprocess.nIter);
    printf("\n   Likelihood function: %g\n",mlprocess.MLmax);
    printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lfsch_L.Lstar),lfsch_L.alfa,lfsch_L.phistar,log10(lfsch_L.phistar));
    printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_L.errLstar/lfsch_L.Lstar/log(10.),lfsch_L.erralfa,lfsch_L.errphistar,lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g Covar(log(Lstar),alpha): %g\n",lfsch_L.covaralfaLstar,lfsch_L.covaralfaLstar/lfsch_L.Lstar/log(10.));
    printf(" Covar(Lstar,Phistar): %g Covar(log(Lstar),log(Phistar)): %g\n",lfsch_L.covarphistarLstar,lfsch_L.covarphistarLstar/lfsch_L.Lstar/log(10.)/lfsch_L.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_L.covaralfaphistar,lfsch_L.covaralfaphistar/lfsch_L.phistar/log(10.));
    printf("##################################################\n");
    printf("#FL_HEAD logLstar alpha Phistar log\n");
    printf("#FL_HEAD E_logLstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",log10(lfsch_L.Lstar),lfsch_L.alfa,lfsch_L.phistar,log10(lfsch_L.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_L.errLstar/lfsch_L.Lstar/log(10.),lfsch_L.erralfa,lfsch_L.errphistar,lfsch_L.errphistar/lfsch_L.phistar/log(10.)); */

    /* dabreu */
    /* fichero para los resultados */
    /* printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 LOG_L_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_LOG_L_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "%g\t", log10(lfsch_L.Lstar));
    fprintf(fout, "%g\t", lfsch_L.alfa);
    fprintf(fout, "%g\t", lfsch_L.phistar);
    fprintf(fout, "%g\t", log10(lfsch_L.phistar));
    fprintf(fout, "%g\t", lfsch_L.errLstar/lfsch_L.Lstar/log(10.));
    fprintf(fout, "%g\t", lfsch_L.erralfa);
    fprintf(fout, "%g\t", lfsch_L.errphistar);
    fprintf(fout, "%g\n", lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    fclose(fout); */
  }
  else 
  {
    printf(" Input the limiting magnitude: ");
    mlim=readd(mlim);
    /* printf("ngalax %i\n", sample.ngalax); */
    if(poissonflag) status=MLA_STY_gc_p_M_wC(sample.ngalax,sample.magSel,sample.magDist,color_mean,color_stddev,sample.errColor,sample.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess);
    /* else            status=  MLA_STY_M(sample.ngalax,sample.mag,sample.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess); */
    else printf("This method is not (yet) available (only poisson)\n");
    printf("\n Solutioner exited with error status %d\n",status);
    printf("\n   Solution found at iteration %d\n",mlprocess.nIter);
    printf("\n   Likelihood function: %g\n",mlprocess.MLmax);
    printf(" Mstar:  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf(" E_Mstar:    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Lstar,Phistar): %g Covar(Lstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD Mstar alpha Phistar log\n");
    printf("#FL_HEAD E_Mstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
/*    printf(" Covar(Lstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Lstar,Phistar): %g Covar(Lstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    */

    /* fichero para los resultados */
    printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 M_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_M_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "# 9 ML_MAX\n");
    fprintf(fout, "# 10 N_ITER\n");
    fprintf(fout, "%g\t", lfsch_M.Mstar);
    fprintf(fout, "%g\t", lfsch_M.alfa);
    fprintf(fout, "%g\t", lfsch_M.phistar);
    fprintf(fout, "%g\t", log10(lfsch_M.phistar));
    fprintf(fout, "%g\t", lfsch_M.errMstar);
    fprintf(fout, "%g\t", lfsch_M.erralfa);
    fprintf(fout, "%g\t", lfsch_M.errphistar);
    fprintf(fout, "%g\t", lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", mlprocess.MLmax);
    fprintf(fout, "%i\n", mlprocess.nIter);
    fclose(fout);
  }

  printf(" Do you want plots of the LF (1=yes, 0=no)?\n");
  plots=readf(plots);
  if(plots) {
     cpgopen("?");
     cpgask(0);
     cpgclos();
  }
}

/* Ya tenemos un nuevo miembro en la casa */
void STY_gmz_f_M_wC(void)
{
  static struct lf_param sty;
  struct sample_data_gmz_wC sample;
  struct fermifsel_M fsel;
  double zlow;

  double color_mean=0.;
  double color_stddev=0.;

  int status=0;
  
  struct Schlf_M lfsch_M;

  static int plots=0; /* para no hacer las gr�ficas */
  FILE *fout;
  static char resultFileName[200]="";
  
  /* Information about the ML process */
  struct MLProcessInfo mlprocess;
  
  get_sample_gmz_wC(&sample);
  printf("Input the color mean:\n");
  color_mean=readd(color_mean);
  printf("Input the color stddev:\n");
  color_stddev=readd(color_stddev);
  set_lf_parameters(&sty);
  zlow=sty.zlow;
  zlow = (zlow < ZMIN ? ZMIN : zlow);
  
  if(sty.islum) 
  {
    printf(" This method is not (yet) available for luminosities (only for magnitudes).\n");
    /* printf(" Input the limiting flux: ");
    flim=readd(flim);
    if(poissonflag) status=MLA_STY_p_L(sample.ngalax,sample.lum,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    else            status=  MLA_STY_L(sample.ngalax,sample.lum,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    printf("\n Solutioner exited with error status %d\n",status);
    printf("\n   Solution found at iteration %d\n",mlprocess.nIter);
    printf("\n   Likelihood function: %g\n",mlprocess.MLmax);
    printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lfsch_L.Lstar),lfsch_L.alfa,lfsch_L.phistar,log10(lfsch_L.phistar));
    printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_L.errLstar/lfsch_L.Lstar/log(10.),lfsch_L.erralfa,lfsch_L.errphistar,lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g Covar(log(Lstar),alpha): %g\n",lfsch_L.covaralfaLstar,lfsch_L.covaralfaLstar/lfsch_L.Lstar/log(10.));
    printf(" Covar(Lstar,Phistar): %g Covar(log(Lstar),log(Phistar)): %g\n",lfsch_L.covarphistarLstar,lfsch_L.covarphistarLstar/lfsch_L.Lstar/log(10.)/lfsch_L.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_L.covaralfaphistar,lfsch_L.covaralfaphistar/lfsch_L.phistar/log(10.));
    printf("##################################################\n");
    printf("#FL_HEAD logLstar alpha Phistar log\n");
    printf("#FL_HEAD E_logLstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",log10(lfsch_L.Lstar),lfsch_L.alfa,lfsch_L.phistar,log10(lfsch_L.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_L.errLstar/lfsch_L.Lstar/log(10.),lfsch_L.erralfa,lfsch_L.errphistar,lfsch_L.errphistar/lfsch_L.phistar/log(10.)); */

    /* dabreu */
    /* fichero para los resultados */
    /* printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 LOG_L_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_LOG_L_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "%g\t", log10(lfsch_L.Lstar));
    fprintf(fout, "%g\t", lfsch_L.alfa);
    fprintf(fout, "%g\t", lfsch_L.phistar);
    fprintf(fout, "%g\t", log10(lfsch_L.phistar));
    fprintf(fout, "%g\t", lfsch_L.errLstar/lfsch_L.Lstar/log(10.));
    fprintf(fout, "%g\t", lfsch_L.erralfa);
    fprintf(fout, "%g\t", lfsch_L.errphistar);
    fprintf(fout, "%g\n", lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    fclose(fout); */
  }
  else 
  {
    printf(" Input limiting magnitude (.5 probability of detection): ");
    fsel.magcut=readd(fsel.magcut);
    printf(" Input selection function sharpness (0 for a step function): ");
    fsel.deltamag=readd(fsel.deltamag); 
    /* Actually calling the method */
    status=MLA_STY_gmz_p_f_M_wC(sample.ngalax,sample.magSel,sample.magDist, sample.errMagDist,color_mean,color_stddev,sample.z,sample.errz,fsel,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess);
    printf("\n Solutioner exited with error status %d\n",status);
    printf("\n   Solution found at iteration %d\n",mlprocess.nIter);
    printf("\n   Likelihood function: %g\n",mlprocess.MLmax);
    printf(" Mstar:  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf(" E_Mstar:    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Lstar,Phistar): %g Covar(Lstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD Mstar alpha Phistar log\n");
    printf("#FL_HEAD E_Mstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
/*    printf(" Covar(Lstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Lstar,Phistar): %g Covar(Lstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    */

    /* fichero para los resultados */
    printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 M_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_M_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "# 9 COVAR_ALPHA_M_STAR\n");
    fprintf(fout, "# 10 COVAR_ALPHA_PHISTAR\n");
    fprintf(fout, "# 11 COVAR_ALPHA_LOGPHISTAR\n");
    fprintf(fout, "# 12 COVAR_PHISTAR_M_STAR\n");
    fprintf(fout, "# 13 COVAR_LOGPHISTAR_M_STAR\n");
    fprintf(fout, "# 14 ML_MAX\n");
    fprintf(fout, "# 15 N_ITER\n");
    fprintf(fout, "%g\t", lfsch_M.Mstar);
    fprintf(fout, "%g\t", lfsch_M.alfa);
    fprintf(fout, "%g\t", lfsch_M.phistar);
    fprintf(fout, "%g\t", log10(lfsch_M.phistar));
    fprintf(fout, "%g\t", lfsch_M.errMstar);
    fprintf(fout, "%g\t", lfsch_M.erralfa);
    fprintf(fout, "%g\t", lfsch_M.errphistar);
    fprintf(fout, "%g\t", lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", lfsch_M.covaralfaMstar);
    fprintf(fout, "%g\t", lfsch_M.covaralfaphistar);
    fprintf(fout, "%g\t", lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", lfsch_M.covarphistarMstar);
    fprintf(fout, "%g\t", lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", mlprocess.MLmax);
    fprintf(fout, "%i\n", mlprocess.nIter);
    fclose(fout);
  }

  printf(" Do you want plots of the LF (1=yes, 0=no)?\n");
  plots=readf(plots);
  if(plots) {
     cpgopen("?");
     cpgask(0);
     cpgclos();
  }
}

void CEG(void)
{
/*   struct lum_func sty; */

  int iter;
  int i,j;
  int *isurvey;
  
  static struct lum_func_ceg ceg;
  struct sample_data sample;
  double *pew,*xew;
  int nbinew;
  double ewminh,ewmaxh;

  int *ewhisto;
  int nbinhist=12;

  
  get_sample_ceg(&sample);
  set_lf_ceg(&ceg); 
  isurvey=vector_i(sample.ngalax);
  

  for(i=0;i<ceg.sdb.nitems;i++) {
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q1C1")) 
      ceg.sdb.si[i].transparency=0.188*pow(10.,4.6)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q1C2"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.8)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q1C3"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.65)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q1C4"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.65)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q2C1"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.7)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q2C2"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.65)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q2C3"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.7)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q2C4"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.65)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q3C2"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.75)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q3C3"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.75)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q3C4"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.7)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q4C1"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.85)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q4C2"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.85)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q4C3"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.85)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
    if(!strcmp(ceg.sdb.si[i].instsetup,"Q4C4"))
       ceg.sdb.si[i].transparency=0.188*pow(10.,4.40)/pow(10.,(16-ceg.sdb.si[i].a)/ceg.sdb.si[i].b);
       printf(" Image %d %s has transparency %f sky %f seeing %f with setup <<%s>>\n",i,ceg.sdb.si[i].image,ceg.sdb.si[i].transparency,ceg.sdb.si[i].sky,ceg.sdb.si[i].seeing,ceg.sdb.si[i].instsetup);
  }
 
  for(j=0;j<sample.ngalax;j++) {
    printf(" This object is %s\n",sample.image[j]);
    isurvey[j]=-1;
    for(i=0;i<ceg.sdb.nitems;i++) {
      if(!strcmp(sample.image[j],ceg.sdb.si[i].image)) {
	isurvey[j]=i;
	printf(" Object %d has image %d\n",j,i);
	break;
      }
    }
    if(isurvey[j]==-1) {
      printf(" ERROR: Image %s not found in database. Exiting\n",sample.image[j]);
      exit(1); 
    }

    printf(" obj %d mag %f ew %f z %f\n",j,sample.mag[j],sample.ew[j],sample.z[j]);  
  }

  projectposf(ceg.fsel,&pew,&xew,&nbinew,5);
    
  ewminh=ceg.ewlim;
  ewmaxh=ceg.fsel.ewbin[ceg.fsel.nEW-1];
  ewhisto=StHisto2_d(sample.ngalax,sample.ew,nbinhist,&ewminh,&ewmaxh);
  printf(" ewmin %f ewmax %f\n",ewminh,ewmaxh);
  for(i=0;i<nbinhist;i++) {
    printf(" bin %d his %d\n",i,ewhisto[i]);
  }  
  ceg.ewd.k=nbinhist;
  ceg.ewd.xk=vector_d(nbinhist+1);
  ceg.ewd.errxk=vector_d(nbinhist+1);
  ceg.ewd.Pk=vector_d(nbinhist);
  ceg.ewd.errPk=vector_d(nbinhist);
  ceg.ewd.covarPk=matrix_d(nbinhist,nbinhist);

  for(i=0;i<nbinhist;i++) {
    ceg.ewd.xk[i]=ewminh+i*(ewmaxh-ewminh)/(nbinhist+0.);
    ceg.ewd.Pk[i]=(float)ewhisto[i]/Lagr2_d(xew,pew,nbinew,(ceg.ewd.xk[i]+ceg.ewd.xk[i+1])/2.);
    printf(" A partir de %f  vale %f (fd %f)\n",ceg.ewd.xk[i],ceg.ewd.Pk[i],Lagr2_d(xew,pew,nbinew,(ceg.ewd.xk[i]+ceg.ewd.xk[i+1])/2.));
  } 
  ceg.ewd.xk[nbinhist]=ewmaxh;

  iter=MLA_STY_s_p_f_PO(sample.ngalax,sample.mag,sample.ew,sample.z,isurvey,ceg.photband, ceg.gamma,ceg.delta,ceg.Kcoc,ceg.fsel,ceg.sdb,ceg.ewlim,ceg.ewd,cosmo,&ceg.lf);

  printf("\n Solution found at iteration %d\n",iter);
  printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(ceg.lf.Lstar),ceg.lf.alfa,ceg.lf.phistar,log10(ceg.lf.phistar));
  printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",ceg.lf.errLstar/ceg.lf.Lstar/log(10.),ceg.lf.erralfa,ceg.lf.errphistar,ceg.lf.errphistar/ceg.lf.phistar/log(10.));
  printf(" Covar(Lstar,alpha): %g Covar(log(Lstar),alpha): %g\n",ceg.lf.covaralfaLstar,ceg.lf.covaralfaLstar/ceg.lf.Lstar/log(10.));
  printf(" Covar(Lstar,Phistar): %g Covar(log(Lstar),log(Phistar)): %g\n",ceg.lf.covarphistarLstar,ceg.lf.covarphistarLstar/ceg.lf.Lstar/log(10.)/ceg.lf.phistar/log(10.));
  printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",ceg.lf.covaralfaphistar,ceg.lf.covaralfaphistar/ceg.lf.phistar/log(10.));

  cpgopen("?");
  cpgask(0);

  cpgclos();
   
}

void STY_MAG_ERR(void)
{
/* dabreu */
/* STY utilizando los errores en la magnitud */
  struct sample_data_mag_err sample_mag_err;
  static double mlim=0;
/*  static double flim=0; */
  double zlow;

  int iter;
  
  static struct lf_param sty;
  struct Schlf_M lfsch_M;
/*  struct Schlf_L lfsch_L; */

  static int poissonflag=0;
  static int plots=0;
  FILE *fout;
  static char resultFileName[300]="";

  /* Information about ML process */
  struct MLProcessInfo mlprocess;
  
  get_sample_mag_err(&sample_mag_err);  
  set_lf_parameters(&sty); 
  zlow=sty.zlow;
  zlow = (zlow < ZMIN ? ZMIN : zlow);
  
  printf(" The following question allows to use an improvement of STY that allows to \ncompute");
  printf(" normalization using Maximum Likelihood. It assumes that number of \ndetected sources");
  printf(" follows a Poisson event. It takes into account Poisson \nerrors both in Mstar, alfa");
  printf(" and Phi_star. These errors are used to be \nlarger than for traditional method, wich");
  printf(" does not take into account the \ninfluence of Poisson detections in Mstar, alfa and Phi_star\n");
  printf(" Use modified method to compute normalization at ML time (0=no/1=yes)? \n");
  poissonflag=readi(poissonflag);
 
  if(sty.islum) {
    printf(" This method is not available: ");
    return;
  }
  else
  { 
    printf(" Input the limiting magnitude: ");
    mlim=readd(mlim);
    if(poissonflag) iter=MLA_STY_gm_p_M(sample_mag_err.ngalax,sample_mag_err.mag,sample_mag_err.mag_err,sample_mag_err.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess);
    else    
    {
      printf(" This method is not available: ");
      return;
    }
    printf("\n Solution found at iteration %d\n",iter);
    printf(" Mstar:  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf(" E_Mstar:    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g\n",lfsch_M.covaralfaMstar);
    printf(" Covar(Lstar,Phistar): %g Covar(Lstar,log(Phistar)): %g\n",lfsch_M.covarphistarMstar,lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lfsch_M.covaralfaphistar,lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD Mstar alpha Phistar log\n");
    printf("#FL_HEAD E_Mstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",lfsch_M.Mstar,lfsch_M.alfa,lfsch_M.phistar,log10(lfsch_M.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfsch_M.errMstar,lfsch_M.erralfa,lfsch_M.errphistar,lfsch_M.errphistar/lfsch_M.phistar/log(10.));

    /* fichero para los resultados */
    printf(" Input file name to write results: ");
    reads(resultFileName,resultFileName);
    if((fout=fopen(resultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",resultFileName);
      return;
    }
    fprintf(fout, "# 1 M_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_M_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "# 9 COVAR_ALPHA_M_STAR\n");
    fprintf(fout, "# 10 COVAR_ALPHA_PHISTAR\n");
    fprintf(fout, "# 11 COVAR_ALPHA_LOGPHISTAR\n");
    fprintf(fout, "# 12 COVAR_PHISTAR_M_STAR\n");
    fprintf(fout, "# 13 COVAR_LOGPHISTAR_M_STAR\n");
    fprintf(fout, "# 14 ML_MAX\n");
    fprintf(fout, "# 15 N_ITER\n");
    fprintf(fout, "%g\t", lfsch_M.Mstar);
    fprintf(fout, "%g\t", lfsch_M.alfa);
    fprintf(fout, "%g\t", lfsch_M.phistar);
    fprintf(fout, "%g\t", log10(lfsch_M.phistar));
    fprintf(fout, "%g\t", lfsch_M.errMstar);
    fprintf(fout, "%g\t", lfsch_M.erralfa);
    fprintf(fout, "%g\t", lfsch_M.errphistar);
    fprintf(fout, "%g\t", lfsch_M.errphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", lfsch_M.covaralfaMstar);
    fprintf(fout, "%g\t", lfsch_M.covaralfaphistar);
    fprintf(fout, "%g\t", lfsch_M.covaralfaphistar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", lfsch_M.covarphistarMstar);
    fprintf(fout, "%g\t", lfsch_M.covarphistarMstar/lfsch_M.phistar/log(10.));
    fprintf(fout, "%g\t", mlprocess.MLmax);
    fprintf(fout, "%i\n", mlprocess.nIter);
    fclose(fout);
  }

  /* dabreu */
  printf(" Do you want plots of the LF (1=yes, 0=no)?\n");
  plots=readf(plots);
  if(plots) {
     cpgopen("?");
     cpgask(0);
     cpgclos();
  } 
}


void SWML(void)
{
/*   struct lum_func sty; */
  struct sample_data sample;
  static double mlim=0;
  static double flim=0;
  double zlow;

  int iter;
  
  static struct lf_param swml;
  struct Steplf_M lfstep_M;
  struct Schlf_M  lfschfit_M;

  struct Steplf_L lfstep_L;
  struct Schlf_L  lfschfit_L;
  double chisq;

  static int poissonflag=0;


  int j;
  
  get_sample(&sample);  
  set_lf_parameters(&swml); 
  zlow=swml.zlow;
  zlow = (zlow < ZMIN ? ZMIN : zlow);
  lfstep_M.covarlnlf=matrix_d(swml.npoint,swml.npoint);
  lfstep_L.covarlnlf=matrix_d(swml.npoint,swml.npoint);
  
  printf(" The following question allows to use an improvement of STY that allows to \ncompute");
  printf(" normalization using Maximum Likelihood. It assumes that number of \ndetected sources");
  printf(" follows a Poisson event. It takes into account Poisson \nerrors both in Mstar, alfa");
  printf(" and Phi_star. These errors are used to be \nlarger than for traditional method, wich");
  printf(" does not take into account the \ninfluence of Poisson detections in Mstar, alfa and Phi_star\n");
  printf(" Use modified method to compute normalization at ML time (0=no/1=yes)? \n");
  poissonflag=readi(poissonflag);


  if(swml.islum) {
    printf(" Input the limiting flux: ");
    flim=readf(flim);
    cpgopen("?");
    if(poissonflag)    cpgswin((float)swml.lum_min,(float)swml.lum_max,-45.,-30.);
    else               cpgswin((float)swml.lum_min,(float)swml.lum_max,-5.,1.);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
    lfstep_L.nbin    =swml.npoint;
    lfstep_L.lumi    =vector_d(swml.npoint+1);
    lfstep_L.lnlf    =vector_d(swml.npoint);
    lfstep_L.errlumi =vector_d(swml.npoint+1);
    lfstep_L.errlnlf =vector_d(swml.npoint);

    for(j=0;j<swml.npoint;j++) {
      lfstep_L.lumi[j]=(swml.lum_min+(swml.lum_max-swml.lum_min)*j/swml.npoint)*log(10);
    }
    lfstep_L.lumi[swml.npoint]=swml.lum_max*log(10);


    if(poissonflag) iter=MLA_SWML_p_L(sample.ngalax,sample.lum,sample.z,flim,swml.area,zlow,swml.zup,cosmo,&lfstep_L);  
    else            iter=MLA_SWML_L(sample.ngalax,sample.lum,sample.z,flim,swml.area,zlow,swml.zup,cosmo,&lfstep_L);  
    printf("\n Solution found at iteration %d\n",iter);
    for(j=0;j<swml.npoint;j++) {
      printf(" Lum %11g - %11g    LF %11g (log=%7g) Err_LF %11g (log=%7g)\n",lfstep_L.lumi[j]/log(10),lfstep_L.lumi[j+1]/log(10),exp(lfstep_L.lnlf[j]),log10(exp(lfstep_L.lnlf[j])),exp(lfstep_L.lnlf[j])*lfstep_L.errlnlf[j],lfstep_L.errlnlf[j]/log(10.));
    }
    PlotStepLF_L(lfstep_L);

  }
  else {
    printf(" Input the limiting magnitude: ");
    mlim=readf(mlim);
    cpgopen("?");
    cpgswin((float)swml.mag_min,(float)swml.mag_max,-5.,1.);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
 
    lfstep_M.nbin     =swml.npoint;
    lfstep_M.magni    =vector_d(swml.npoint+1);
    lfstep_M.lnlf     =vector_d(swml.npoint);
    lfstep_M.errmagni =vector_d(swml.npoint+1);
    lfstep_M.errlnlf  =vector_d(swml.npoint);
    for(j=0;j<swml.npoint;j++) {
      lfstep_M.magni[j]=swml.mag_min+(swml.mag_max-swml.mag_min)*j/swml.npoint;
    }
    lfstep_M.magni[swml.npoint]=swml.mag_max;
    if(poissonflag)  iter=MLA_SWML_p_M(sample.ngalax,sample.mag,sample.z,mlim,swml.area,zlow,swml.zup,cosmo,&lfstep_M);
    else             iter=MLA_SWML_M(sample.ngalax,sample.mag,sample.z,mlim,swml.area,zlow,swml.zup,cosmo,&lfstep_M);
    printf("\n Solution found at iteration %d\n",iter);
    for(j=0;j<swml.npoint;j++) printf(" Mag %11g - %11g    LF %11g (log=%7g) Err_LF %11g (log=%7g)\n",lfstep_M.magni[j],lfstep_M.magni[j+1],exp(lfstep_M.lnlf[j]),log10(exp(lfstep_M.lnlf[j])),exp(lfstep_M.lnlf[j])*lfstep_M.errlnlf[j],lfstep_M.errlnlf[j]/log(10.));
    printf(" Plotting function\n");
    PlotStepLF_M(lfstep_M);
    printf(" After plotting\n");   
   
  }

  cpgopen("?");  
  cpgask(0);


  if(swml.islum) 
  {
    FitSch2StepLF_L(lfstep_L, &lfschfit_L, &chisq); 
    PlotStepSchLF_L(lfstep_L,lfschfit_L);   
    free(lfstep_L.lumi);
    free(lfstep_L.lnlf);
    free(lfstep_L.errlumi);
    free(lfstep_L.errlnlf);
  }
  else           {
    FitSch2StepLF_M(lfstep_M, &lfschfit_M, &chisq); 
    PlotStepSchLF_M(lfstep_M,lfschfit_M);
    free(lfstep_M.magni);
    free(lfstep_M.lnlf);
    free(lfstep_M.errmagni); 
    free(lfstep_M.errlnlf);
  }
  cpgclos();

  free_matrix_d(lfstep_M.covarlnlf,swml.npoint,swml.npoint);
  free_matrix_d(lfstep_L.covarlnlf,swml.npoint,swml.npoint);
}


void VVmax(void)
{
  int i,j;
  
  /* Variables de la funcion de lum */
  struct Steplf_M lf_vvmax_M;
  struct Steplf_L lf_vvmax_L;
  static struct lf_param vvmax_param={-40,40,0,0,0,10,41252,1,10};
  struct sample_data_sel_dist sample;
  struct Schlf_M  lfschfit_M;
  struct Schlf_L  lfschfit_L;

  static double mlim=1.;
  static double llim=1.;
  double M=0;
  double L=0;
/*   int vol; */
  double zmax;
  //double *x,*y,*sigy,*yfit;
  double dm=0;
  double dl=0;

  /* Variables para el test V/Vmax */
  double *vvmaxtest;
  int   *ngaltest;           /* n�mero de galaxias en cada bin del test */
  double mmin_t=0,mmax_t=0;       /* magnitudes m�nimas y m�ximas desde las cuales se hara el test V/Vmax. */
  double lmin_t=0,lmax_t=0;       /* luminosidades m�nimas y m�ximas desde las cuales se hara el test V/Vmax. */
  int   n_t;                 /* n�mero de puntos para el test */
  double *mhist=NULL;              /* array para las magnitudes del histograma */
  double *loglhist=NULL;           /* array para las luminosidades del histograma */
  double *mlim_t=NULL;
  double *llim_t=NULL;
  double *logllim_t=NULL;
  char  cnul;
  double fnul;
  double zlow;
  double zup;

  /* Variables para el ajuste a Schechter */
  double chisq;
  //double dum1,dum2,dum3;
  //int nfit;

  /* Variables para los plots */
  int pg1=0,pg2;
  //double ymin,ymax;
  /* dabreu */
  static int plots=0; /* para no hacer gr�ficas */

  /* Variables para escribir los resultados en un fichero */
  FILE *fout = NULL;
  static char schfitResultFileName[300]="";
  static char vvmaxResultFileName[300]="";

  /* long results file (for test) */
  static int allVVmax=0; /* do not want this file */
  static char allVVmaxResultsFileName[300]="";
  
  /* Obtengo los datos y los parametros de la LF */
  get_sample_sel_dist(&sample);
  set_lf_parameters(&vvmax_param);

  /* Allocateo la LF */
  lf_vvmax_M.nbin      = vvmax_param.npoint;
  lf_vvmax_M.magni     =vector_d(lf_vvmax_M.nbin+1);
  lf_vvmax_M.errmagni  =vector_d(lf_vvmax_M.nbin+1);
  lf_vvmax_M.lnlf      =vector_d(lf_vvmax_M.nbin);
  lf_vvmax_M.errlnlf   =vector_d(lf_vvmax_M.nbin);
  lf_vvmax_M.lf        =vector_d(lf_vvmax_M.nbin);
  lf_vvmax_M.errlf     =vector_d(lf_vvmax_M.nbin);
  lf_vvmax_M.covarlnlf =matrix_d(lf_vvmax_M.nbin,lf_vvmax_M.nbin);
  lf_vvmax_M.ngalbin   =vector_i(lf_vvmax_M.nbin);
  lf_vvmax_L.nbin      =vvmax_param.npoint;
  lf_vvmax_L.lumi      =vector_d(lf_vvmax_L.nbin+1);
  lf_vvmax_L.errlumi   =vector_d(lf_vvmax_L.nbin+1);
  lf_vvmax_L.lnlf      =vector_d(lf_vvmax_L.nbin);
  lf_vvmax_L.errlnlf   =vector_d(lf_vvmax_L.nbin);
  lf_vvmax_L.covarlnlf =matrix_d(lf_vvmax_L.nbin,lf_vvmax_L.nbin);

  printf(" Performing the V/Vmax test\n");
  if(vvmax_param.islum) 
  {
    MinMax_d(sample.ngalax,sample.lumSel,&lmin_t,&lmax_t);
    printf(" Fluxes range: %g - %g\n",lmin_t,lmax_t);
    lmin_t*=0.3;lmax_t*=1.2;
  }
  else 
  {
    MinMax_d(sample.ngalax,sample.magSel,&mmin_t,&mmax_t);
    printf(" Apparent magnitude range: %g - %g\n",mmin_t,mmax_t);
    mmin_t-=1.;mmax_t+=4.;
  }
  n_t=100;
  vvmaxtest   =malloc(n_t*sizeof(double));
  ngaltest    =malloc(n_t*sizeof(int));
  if(vvmax_param.islum) {
    loglhist      =malloc(sample.ngalax*sizeof(double));
    logllim_t     =malloc(n_t*sizeof(double));
    llim_t        =malloc(n_t*sizeof(double));
  }
  else {
    mhist       =malloc(sample.ngalax*sizeof(double));
    mlim_t      =malloc(n_t*sizeof(double));
  }
/*   printf(" Input lower redshift: "); */
  zlow=vvmax_param.zlow;
  zup=vvmax_param.zup;
  zlow = (zlow < ZMIN ? ZMIN : zlow);
  if(DEBUG) printf(" zlow %f zup %f\n",zlow,zup);

/*   printf(" Input limiting redshift (0=none, magnitude limited sample): "); */
/*   zup=readf(0.); */
  /* Aqui comienza el Test VVmax */
  for(i=0;i<n_t;i++)
  {
    if(vvmax_param.islum)
    {
      llim_t[i]=pow(10.,i*(log10(lmax_t)-log10(lmin_t))/(n_t-1.)+log10(lmin_t));
      logllim_t[i]=i*(log10(lmax_t)-log10(lmin_t))/(n_t-1.)+log10(lmin_t);
      if(DEBUG2) printf(" llim %g\n",llim_t[i]);
    }
    else
    {
      mlim_t[i]=i*(mmax_t-mmin_t)/(n_t-1.)+mmin_t;
      if(DEBUG) printf(" mlim %f\n",mlim_t[i]);
    }
    vvmaxtest[i]=0;
    ngaltest[i]=0;
    for(j=0;j<sample.ngalax;j++)
    {
      if(vvmax_param.islum)
      {
	loglhist[j]=log10(sample.lumSel[j]);
	L=Lum(sample.z[j],sample.lumSel[j],cosmo);
	if(DEBUG2)        printf(" Galax %d z %f flux %g Lum %g mag %f Mag %f\n",j,sample.z[j],sample.lumSel[j],L,sample.magSel[j],Mag(sample.z[j],sample.magSel[j],cosmo));
	if(L>1e36) 
 	printf(" Galax %d z %f flux %g Lum %g mag %f Mag %f\n",j,sample.z[j],sample.lumSel[j],L,sample.magSel[j],Mag(sample.z[j],sample.magSel[j],cosmo)); 

      }
      else
      {
	mhist[j]=sample.magSel[j];
	M=Mag(sample.z[j],sample.magSel[j],cosmo);
	if(DEBUG2)        printf(" Galax %d z %f mag %f Mag %f\n",j,sample.z[j],sample.magSel[j],M);
      }
      /*       if(sample.mag[j]<mlim_t[i]) {   Esta es la antigua, cuando no podia hacer en luminosidades */
      /* Cuidado, sample.lum[j] (o mag[j]) puede no estar alocateado, 
	 por eso pongo antes vvmax.islum, que comprueba antes y si no lo cumple 
	 no pasa a la siguiente */
/*       if(( !vvmax.islum && sample.mag[j]<mlim_t[i])  || (vvmax.islum && sample.lum[j]<llim_t[i]) ) { */
      if(vvmax_param.islum)
      {
	if(sample.lumSel[j]>llim_t[i])
        {
	  if(DEBUG2) printf(" Entro la flux=%f con llim=%f\n",sample.lumSel[j],llim_t[i]);
	  zmax=Z_l(llim_t[i],L,cosmo);
	  if(zmax==0) zmax=ZMIN;
	  if(zmax>zup && zup!=0) zmax=zup;
	  if(zmax<zlow && zlow!=0) {printf(" zlow %f zmax %f\n ERROR!!\n",zlow,zmax);exit(1);}
	  if(DEBUG2) printf(" zmax %f %g %g %g\n",zmax, L, llim,sample.z[j]); 
	  /* 		vvmaxtest[i]+=Vol_LF(sample.z[j])/Vol_LF(zmax);  VVmax tradicional */
	  vvmaxtest[i]+=(Vol(sample.z[j],cosmo)-Vol(zlow,cosmo))/(Vol(zmax,cosmo)-Vol(zlow,cosmo));  /* Cambiado para un zlow!=0 */
	  ngaltest[i]++;
	}
      }
      else
      {
	if(sample.magSel[j]<mlim_t[i])
        {
	  if(DEBUG2) printf(" Entro la m=%f con mlim=%f\n",sample.magSel[j],mlim_t[i]);
	  zmax=Z_m(mlim_t[i],M,cosmo);
          if(DEBUG2) printf("Z_m: m %g M %g\n",mlim_t[i],M);
	  if(zmax==0) zmax=ZMIN;
	  if(zmax>zup && zup!=0) zmax=zup;
	  if(zmax<zlow && zlow!=0) {printf(" zlow %f zmax %f\n ERROR!!\n",zlow,zmax);exit(1);}
	  if(DEBUG2) printf(" zmax %f %f %f %f\n",zmax, M, mlim,sample.z[j]); 
	  /* 		vvmaxtest[i]+=Vol_LF(sample.z[j])/Vol_LF(zmax);  VVmax tradicional */
	  vvmaxtest[i]+=(Vol(sample.z[j],cosmo)-Vol(zlow,cosmo))/(Vol(zmax,cosmo)-Vol(zlow, cosmo));  /* Cambiado para un zlow!=0 */
	  ngaltest[i]++;
	}
	if(DEBUG2) printf(" No entro la m=%f con mlim=%f\n",sample.magSel[j],mlim_t[i]);
      }
    }
    vvmaxtest[i]/=ngaltest[i];
    if(vvmax_param.islum)    printf(" Flux %g V/Vmax %f Ngal %d\n",llim_t[i],vvmaxtest[i],ngaltest[i]);
    else               printf(" Mag %g V/Vmax %f Ngal %d\n",mlim_t[i],vvmaxtest[i],ngaltest[i]);
  }
  /* dabreu */
  if(DEBUG2) printf("Despues de 'Mag g V/Vmax f Ngal'\n");
  printf(" Do you want plots(1=yes, 0=no)?\n");
  plots=readf(plots);
  if(plots) 
  {
    pg1=cpgopen("/xserve");
    if(vvmax_param.islum)
    {
        cpgswin(log10(lmin_t),log10(lmax_t),0.,sample.ngalax/2.);
        cpghist_d(sample.ngalax,loglhist,log10(lmin_t),log10(lmax_t),15,1);
        cpglab("Apparent fluxes (W/m2)","<V/Vmax>","V/Vmax test");
        cpgbox("BCTNSL",0,0,"BTNS",0,0);
    }
    else
    {
        cpgswin(mmin_t,mmax_t,0.,sample.ngalax/2.);
      /*     cpghist_d(sample.ngalax,mhist,mmin_t,mmax_t,15,1); */
        cpglab("Magnitude","<V/Vmax>","V/Vmax test");
        cpgbox("BCTNS",0,0,"BTNS",0,0);
    }
    if(vvmax_param.islum)
    {
      cpgswin(log10(lmin_t),log10(lmax_t),0.,1.);
      cpgbox("BCTNSL",0,0,"CTMS",0,0);
      cpgpt_d(n_t,logllim_t,vvmaxtest,1); 
      printf(" ANTES abd ad\n");
      printf(" The V/Vmax method requires to know the limiting flux.\n Click with the left button to select the limiting luminosity with the mouse or right button to select it by keyboard\n");
      cpgband_d(6,1,0.,0.,&llim,&fnul,&cnul); 
      if(DEBUG) printf(" Pasae band\n");
      llim=pow(10.,llim);
      if(DEBUG) printf("KK Limiting flux: %g\n",llim);
    }
    else
    {
        cpgswin(mmin_t,mmax_t,0.,1.);
        cpgbox("BCTNS",0,0,"CTMS",0,0);
        cpgpt_d(n_t,mlim_t,vvmaxtest,1); 
      printf(" The V/Vmax method requires to know the limiting magnitude.\n Click with the left button to select the limiting magnitude with the mouse or right button to select it by keyboard\n");
        cpgband_d(6,1,0.,0.,&mlim,&fnul,&cnul); 
    }
    /*   //printf(" OPION %c\n",cnul); */
    if(cnul=='X')
    {
      if(vvmax_param.islum)
      {
        printf(" Input the limiting flux: ");
        llim=readd(llim);
      }
      else
      {
        printf(" Input the limiting magnitude: ");
        mlim=readd(mlim);
      }
    }
  /* dabreu */
  }   /* final del if(plots) */ 
  else 
  {
    if(vvmax_param.islum) 
    {
      printf(" Input the limiting flux: ");
      llim=readd(llim);
    }
    else 
    {
      printf(" Input the limiting magnitude: ");
      mlim=readd(mlim);
    }
  }


  /* Calculo los intervalos de la LF (en luminosidades o mangitudes)   */
  if(vvmax_param.islum)
  {
    dl=(vvmax_param.lum_max-vvmax_param.lum_min)/lf_vvmax_M.nbin;
    for(i=0;i<=lf_vvmax_L.nbin;i++) 
      lf_vvmax_L.lumi[i]=vvmax_param.lum_min+i*dl;
    printf(" Limiting luminosity: %g. Lum bin %g From %g to %g\n",llim,dl,vvmax_param.lum_max,vvmax_param.lum_min);
  }
  else 
  {
    dm=(vvmax_param.mag_max-vvmax_param.mag_min)/lf_vvmax_M.nbin;
    for(i=0;i<=lf_vvmax_M.nbin;i++) 
      lf_vvmax_M.magni[i]=vvmax_param.mag_min+i*dm;
    printf(" Limiting magnitude: %f\n",mlim);
  }

  if(DEBUG) printf(" Intervalo %g\n",dm);


  /* Calculo la FL */
  if(vvmax_param.islum)
  {
    VVmax_L(sample.ngalax, sample.lumSel, sample.lumDist, sample.z, llim, vvmax_param.area, vvmax_param.zlow, vvmax_param.zup, cosmo, &lf_vvmax_L);
    PrintStepLF_L(lf_vvmax_L);
  }
  else
  {
    VVmax_M(sample.ngalax, sample.magSel, sample.magDist, sample.z, mlim, vvmax_param.area, vvmax_param.zlow, vvmax_param.zup, cosmo, &lf_vvmax_M);
    PrintStepLF_M(lf_vvmax_M);
  }

  /* C�lculo del ajuste a Schechter */
  if(vvmax_param.islum)
  {
    FitSch2StepLF_L(lf_vvmax_L, &lfschfit_L, &chisq);
    printf("##################################################\n");
    printf("#FL_HEAD logLstar alpha Phistar log\n");
    printf("#FL_HEAD E_logLstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",log10(lfschfit_L.Lstar),lfschfit_L.alfa,lfschfit_L.phistar,log10(lfschfit_L.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfschfit_L.errLstar/lfschfit_L.Lstar/log(10),lfschfit_L.erralfa,lfschfit_L.errphistar,lfschfit_L.errphistar/lfschfit_L.phistar/log(10));

    /* dabreu */
    /* fichero para los resultados */
    printf(" Input file name to write schfit results: ");
    reads(schfitResultFileName,schfitResultFileName);
    if((fout=fopen(schfitResultFileName,"w")) ==NULL) 
    {
      printf(" Couldn't open %s for writing\n",schfitResultFileName);
      return;
    }
    fprintf(fout, "# 1 LOG_L_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_LOG_L_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "%g\t", log10(lfschfit_L.Lstar));
    fprintf(fout, "%g\t", lfschfit_L.alfa);
    fprintf(fout, "%g\t", lfschfit_L.phistar);
    fprintf(fout, "%g\t", log10(lfschfit_L.phistar));
    fprintf(fout, "%g\t", lfschfit_L.errLstar/lfschfit_L.Lstar/log(10));
    fprintf(fout, "%g\t", lfschfit_L.erralfa);
    fprintf(fout, "%g\t", lfschfit_L.errphistar);
    fprintf(fout, "%g\n", lfschfit_L.errphistar/lfschfit_L.phistar/log(10));
    fclose(fout); 
  }
  else
  {
    FitSch2StepLF_M(lf_vvmax_M, &lfschfit_M, &chisq);
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD Mstar alpha Phistar log\n");
    printf("#FL_HEAD E_Mstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",lfschfit_M.Mstar,lfschfit_M.alfa,lfschfit_M.phistar,log10(lfschfit_M.phistar));
    printf("#FL_ERR %g %g %g %g\n",lfschfit_M.errMstar,lfschfit_M.erralfa,lfschfit_M.errphistar,lfschfit_M.errphistar/lfschfit_L.phistar/log(10));

    /* dabreu */
    /* fichero para los resultados */
    printf(" Input file name to write schfit results: ");
    reads(schfitResultFileName,schfitResultFileName);
    if((fout=fopen(schfitResultFileName,"w")) ==NULL) {
      printf(" Couldn't open %s for writing\n",schfitResultFileName);
      return;
    }
    fprintf(fout, "# 1 M_STAR\n");
    fprintf(fout, "# 2 ALPHA\n");
    fprintf(fout, "# 3 PHISTAR\n");
    fprintf(fout, "# 4 LOG_PHISTAR\n");
    fprintf(fout, "# 5 ERR_M_STAR\n");
    fprintf(fout, "# 6 ERR_ALPHA\n");
    fprintf(fout, "# 7 ERR_PHISTAR\n");
    fprintf(fout, "# 8 ERR_LOG_PHISTAR\n");
    fprintf(fout, "# 9 COVAR_ALPHA_M_STAR\n");
    fprintf(fout, "# 10 COVAR_ALPHA_PHISTAR\n");
    fprintf(fout, "# 11 COVAR_ALPHA_LOGPHISTAR\n");
    fprintf(fout, "# 12 COVAR_PHISTAR_M_STAR\n");
    fprintf(fout, "# 13 COVAR_LOGPHISTAR_M_STAR\n");
    fprintf(fout, "# 14 ML_MAX\n");
    fprintf(fout, "# 15 N_ITER\n");
    fprintf(fout, "%g\t", lfschfit_M.Mstar);
    fprintf(fout, "%g\t", lfschfit_M.alfa);
    fprintf(fout, "%g\t", lfschfit_M.phistar);
    fprintf(fout, "%g\t", log10(lfschfit_M.phistar));
    fprintf(fout, "%g\t", lfschfit_M.errMstar);
    fprintf(fout, "%g\t", lfschfit_M.erralfa);
    fprintf(fout, "%g\t", lfschfit_M.errphistar);
    fprintf(fout, "%g\t", lfschfit_M.errphistar/lfschfit_M.phistar/log(10));
    fprintf(fout, "INDEF\t");
    fprintf(fout, "INDEF\t");
    fprintf(fout, "INDEF\t");
    fprintf(fout, "INDEF\t");
    fprintf(fout, "INDEF\t");
    fprintf(fout, "INDEF\t");
    fprintf(fout, "INDEF\n");
    fclose(fout);

    /* to store values of LF in each bin */
    printf("Do you want VVmax LF values for each bin? (1=yes, 0=no):\n");
    allVVmax=readi(allVVmax);
    if(allVVmax)
    {
      printf("Give a file name:\n");
      reads(allVVmaxResultsFileName, allVVmaxResultsFileName);
      if ((fout=fopen(allVVmaxResultsFileName, "w")) == NULL)
      {
         printf("Couldn't open %s for writing\n",allVVmaxResultsFileName);
         return;
      }
      fprintf(fout, "# 1 MAG\n");
      fprintf(fout, "# 2 LN_VVMAX\n");
      fprintf(fout, "# 3 VVMAX\n");
      fprintf(fout, "# 4 ERR_LN_VVMAX\n");
      fprintf(fout, "# 5 NGALBIN\n");
      for(i=0;i<lf_vvmax_M.nbin;i++)
      {
         fprintf(fout, "%g\t%g\t", lf_vvmax_M.magni[i], lf_vvmax_M.lnlf[i]);
         fprintf(fout, "%g\t%g\t",exp(lf_vvmax_M.lnlf[i]), lf_vvmax_M.errlnlf[i]);
         fprintf(fout, "%d\n",lf_vvmax_M.ngalbin[i]);
      }
      fclose(fout);
    }

  }

  /* dabreu */
  if(plots)
  {
    pg2=cpgopen("?");
    
    if(vvmax_param.islum)
    {
      PlotStepSchLF_L(lf_vvmax_L, lfschfit_L);
    }
    else
    {
      PlotStepSchLF_M(lf_vvmax_M, lfschfit_M);
    }

    /* dabreu */
    /* fichero para los resultados */
    printf(" Input file name to write vvmax results: ");
    reads(vvmaxResultFileName,vvmaxResultFileName);
    if((fout=fopen(vvmaxResultFileName,"w")) ==NULL) 
    {
      printf(" Couldn't open %s for writing\n",vvmaxResultFileName);
      return;
    }
    
    if(vvmax_param.islum)
    {
      fprintf(fout, "# 1 FLUX\n");
      fprintf(fout, "# 2 LN_VVMAX\n");
    
      for(i=0;i<lf_vvmax_L.nbin;i++)
      {
        fprintf(fout, "%g\t", lf_vvmax_L.lumi[i]);
        fprintf(fout, "%g\t", lf_vvmax_L.lnlf[i]);
        fprintf(fout, "\n");
      }
      fclose(fout);    
    }
    
    else
    {
      fprintf(fout, "# 1 MAG\n");
      fprintf(fout, "# 2 LN_VVMAX\n");
      fprintf(fout, "# 3 VVMAX\n");
      fprintf(fout, "# 4 ERR_LN_VVMAX\n");
    
      for(i=0;i<lf_vvmax_M.nbin;i++)
      {
        fprintf(fout, "%g\t", lf_vvmax_M.magni[i]);
        fprintf(fout, "%g\t", lf_vvmax_M.lnlf[i]);
        fprintf(fout, "%g\t%g",exp(lf_vvmax_M.lnlf[i]), lf_vvmax_M.errlnlf[i]);
        fprintf(fout, "\n");
      }
      fclose(fout);
    }

/*    ymin=1e38;
    ymax=-1e38;
    for(i=0;i<vvmax_param.npoint;i++) 
    {
      if(lf_vvmax_M.lnlf[i]!=-1/0.) 
      {
        if(log10(exp(lf_vvmax_M.lnlf[i]))>ymax) ymax=(double)log10(exp(lf_vvmax_M.lnlf[i]));
        if(log10(exp(lf_vvmax_M.lnlf[i]))<ymin) ymin=(double)log10(exp(lf_vvmax_M.lnlf[i]));
      }
    }
    ymin=ymin-2.5;
    ymax=ymax+2.5;

    printf(" Cuts %f %f \n",ymin,ymax);

    if(vvmax_param.islum)
    {
      cpgswin(vvmax_param.lum_max+0.1,vvmax_param.lum_min-0.1,ymin,ymax);
      cpgbox("BCTNSL",0,0,"BCTNS",0,0);
    }
    else
    {
        cpgswin(vvmax_param.mag_max+1,vvmax_param.mag_min-1,ymin,ymax);
        cpgbox("BCTNS",0,0,"BCTNS",0,0);
    }
    cpgsci(1);
    cpgsch(1.2);

    printf(" Antes del for\n");

    for(i=0;i<nfit;i++) 
    {
      printf("Valor de i= %i \n",i);
      printf("x %g \n",x[i]);
      printf("y %g \n",y[i]);
      cpgpt1(x[i],y[i],17);
      dum1=x[i];
      dum2=y[i]-sigy[i];
      dum3=y[i]+sigy[i];
      printf(" dum2 %f dum3 %f \n",dum2,dum3);
      cpgerry_d(1,&dum1,&dum2,&dum3,0.35); 
    }
    cpgsci(4);
    cpgline_d(nfit,x,yfit);
    cpgsci(1);
    cpgsch(1.);
  
    if(vvmax_param.islum) cpglab("Luminosity (W)","log(Phi\\dlum\\u) (#/Mpc3/W) ","Luminosity function by the V/Vmax method computed with fluxes");
    else cpglab("Mag","log(Phi\\dmag\\u) (#/Mpc3/Mag) ","Luminosity function by the V/Vmax method computed with magnitudes");
    
*/  
    cpgclos();
    cpgslct(pg1); 
    cpgclos();
  }  /* final del if del 2� plot */


  free(lf_vvmax_M.magni);
  free(lf_vvmax_M.errmagni);
  free(lf_vvmax_M.lnlf);
  free(lf_vvmax_M.errlnlf);
  free(lf_vvmax_M.lf);
  free(lf_vvmax_M.errlf);
  free(lf_vvmax_M.ngalbin);
  free_matrix_d(lf_vvmax_M.covarlnlf, lf_vvmax_M.nbin,lf_vvmax_M.nbin);
  free(lf_vvmax_L.lumi);
  free(lf_vvmax_L.errlumi);
  free(lf_vvmax_L.lnlf);
  free(lf_vvmax_L.errlnlf);
  free_matrix_d(lf_vvmax_L.covarlnlf, lf_vvmax_L.nbin,lf_vvmax_L.nbin);
}

void Calc_Num(void)
{ 


  static double zup=0.5;
  static double zlow=0;
  double ztmp=0;
  static double mlimite=20;
  double mlim;
  double Mabs;
  
  double Ntot;
  double Nbin;
  static char opt='M';

  double Mabsmin,Mabsmax;
  double mlimmin,mlimmax;

  /* dabreu */
  static int plots=0; /* para no hacer las gr�ficas */
  
  printf(" Do you want plots of the LF (1=yes, 0=no)?\n");
  plots=readf(plots);
  
  if(plots) {
     cpgbeg(0,"?",2,1);
  }

  printf("\n L Give LF in luminosity\n"); 
  printf(" M Give LF in magnitudes\n");
  opt=readc(opt); 
  
  switch (opt) { 
  case 'L':
  case 'l':
    printf(" Introduce apparent flux limit (Watt/m2 , 1 erg/s/cm2 = 10^(-3) Watt/m2): ");
    mlimite=-2.5*log10(readf(pow(10.,-0.4*mlimite)));
    mlim=mlimite;
    set_Schechter_L();
    break;
  default:
    opt='M';
    printf(" Introduce apparent magnitude limit: ");
    mlimite=readf(mlimite);
    mlim=mlimite;
    set_Schechter_M();
    break;
  }

  printf(" Introduce lower redshift limit: ");
  if(zlow<=ZMIN) zlow=0;
  zlow=readf(zlow);
/*   //scanf("%f",&zlow); */
  if(zlow==0) zlow=ZMIN;
  printf(" Introduce upper redshift limit (0=none, magnitude limited survey: ");
  zup=readf(zup);
  /*   //scanf("%f",&zup); */
  
  if(zup==0) zup=10.;

/*   //  printf(" Phi star %f\n",phistar); */

   printf(" Surveyed volume %e Mpc3\n",(Vol(zup,cosmo)-Vol(zlow,cosmo))/1.e18);

  /* Ahora hago un plot de las galaxias que se van a detectar en funcion 
     de la magnmitud absoluta */
  if(plots){
     cpgpanl(2,1);
     cpgscf(2);
     cpgswin(-10,-24.,0.,1000000./41252.);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpgsch(1.5);
     cpglab("Absolute magnitude","Number/degree","Number of objects/degree per 0.5 bin in magnitude");
     cpgsch(3.);
  }

  if(opt=='M') {
    printf("#Objects per absolute magnitud interval\n");
    printf("#Mmin   Mmax    Zmin  Zmax     N_sky     N/deg\n");
    printf("############################################\n");
    Mabsmin=-7;
    Mabsmax=-24;
  }
  else {
    printf("#Objects per luminosity interval logLstar %f\n",-0.4*(schlf_L.Lstar-.25));
    printf("#log(Lmin)   log10(Lmax)   Zmin  Zmax     N_sky     N/deg\n");
    printf("############################################\n");
    Mabsmin=-80;
    Mabsmax=-120;
  }
  for(Mabs=Mabsmin;Mabs>Mabsmax;Mabs-=.5) {
    if(opt=='M')ztmp=Z_m(mlim,Mabs,cosmo);
    else        ztmp=Z_l(pow(10.,-0.4*mlim),pow(10.,-0.4*Mabs),cosmo);
    if(ztmp<zlow) ztmp=zlow;
    if(ztmp>zup) ztmp=zup;
/*     //printf(" Absolute magnitude: %f Upper redshift limit: %f  %f\n",Mabs,ztmp,Mag_LF(ztmp,mlim)); */
   Nbin=0.5*Schechter_M(Mabs,schlf_M)*(Vol(ztmp,cosmo)-Vol(zlow,cosmo))/1.e18;
    if(opt=='M')     printf(" %6.2f %6.2f   %5.3f %5.3f  %10g  %10g/deg\n",Mabs-.25,Mabs+.25,zlow,ztmp,Nbin,Nbin/41252.);
    else             printf(" %6.2f %6.2f   %5.3f %5.3f  %10g  %10g/deg\n",-0.4*(Mabs-.25),-0.4*(Mabs+.25),zlow,ztmp,Nbin,Nbin/41252.);
      if(plots) cpgpt1(Mabs,Nbin/41252.,3);
  }

  if(opt=='M') {
    printf("#Objects per absolute magnitud interval\n");
    printf("#Mmin   Mmax    Zmin  Zmax     N_sky     N/deg\n");
    printf("############################################\n");
    Mabsmin=-9;
    Mabsmax=-24;
  }
  else {
    printf("#Objects per luminosity interval\n");
    printf("#log(Lmin)   log10(Lmax)   Zmin  Zmax     N_sky     N/deg\n");
    printf("############################################\n");
    Mabsmin=-80;
    Mabsmax=-120;
  }
  for(Mabs=Mabsmin;Mabs>Mabsmax;Mabs-=1.) {
    if(opt=='M')ztmp=Z_m(mlim,Mabs,cosmo);
    else        ztmp=Z_l(pow(10.,-0.4*mlim),pow(10.,-0.4*Mabs),cosmo);
    if(ztmp>zup) ztmp=zup;
    if(ztmp<zlow) ztmp=zlow;
    /*     //printf(" Absolute magnitude: %f Upper redshift limit: %f  %f\n",Mabs,ztmp,Mag_LF(ztmp,mlim)); */
       Nbin=Schechter_M(Mabs,schlf_M)*(Vol(ztmp,cosmo)-Vol(zlow,cosmo))/1.e18;
    if(opt=='M')     printf(" %6.2f %6.2f   %5.3f %5.3f  %10g  %10g/deg\n",Mabs-.25,Mabs+.25,zlow,ztmp,Nbin,Nbin/41252.);
    else             printf(" %6.2g %6.2g   %5.3f %5.3f  %10g  %10g/deg\n",-0.4*(Mabs-.25),-0.4*(Mabs+.25),zlow,ztmp,Nbin,Nbin/41252.);
      if(plots) cpgpt1(Mabs,Nbin/41252.,3);
  }
  Nbin=0;

  for(Mabs=-18.75;Mabs>-24.25;Mabs-=.5) {
    if(opt=='M')ztmp=Z_m(mlim,Mabs,cosmo);
    else        ztmp=Z_l(pow(10.,-0.4*mlim),pow(10.,-0.4*Mabs),cosmo);
    if(ztmp>zup) ztmp=zup;
    if(ztmp<zlow) ztmp=zlow;
    /*     printf(" Absolute magnitude: %f Upper redshift limit: %f  %f\n",Mabs,ztmp,Mag_LF(ztmp,mlim)); */
    Nbin+=0.5*Schechter_M(Mabs,schlf_M)*(Vol(ztmp,cosmo)-Vol(zlow,cosmo))/1.e18;
      if(plots) cpgpt1(Mabs,Nbin/41252.,3);
  }
/*   printf(" -24.25 -18.75   %5.3f %5.3f  %10g  %10g/deg\n",zlow,ztmp,Nbin,Nbin/41252.); */


/*   //  Lo mismo pero en magnitud aparente  */
  if(opt=='M') {
    printf("#Objects per apparent magnitud interval\n");
    printf("#mmin   mmax       N_sky     N/deg\n");
    printf("############################################\n");
    mlimmin=15;
    mlimmax=28;
  }
  else {
    printf("#Objects per flux interval\n");
    printf("#log(Fmin)   log10(Fmax)   N_sky     N/deg\n");
    printf("############################################\n");
    mlimmin=35;
    mlimmax=50;
  }

  if(plots) cpgpanl(1,1);
  
  if(opt=='M')   Ntot=Int_sch_M(schlf_M,zlow,zup,mlim,cosmo);
  else           Ntot=Int_sch_L(schlf_L,zlow,zup,pow(10.,-0.4*mlim),cosmo);
      
  if(plots){
    cpgscf(2);
    cpgswin(13.,28.,0.,40000000./41252.);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
    cpgsch(1.5);
    cpglab("Limiting magnitude","Number/degree","Number of objects/degree per 0.5 bin in magnitude");
    cpgsch(3.);
  }
  for(mlim=mlimmin;mlim<mlimmax;mlim+=.25) {
    if(opt=='M')  Nbin=-Int_sch_M(schlf_M,zlow,zup,mlim-.125,cosmo)+Int_sch_M(schlf_M,zlow,zup,mlim+.125,cosmo); 
    else          Nbin=-Int_sch_L(schlf_L,zlow,zup,pow(10.,-0.4*(mlim-.125)),cosmo)+Int_sch_L(schlf_L,zlow,zup,pow(10.,-0.4*(mlim+.125)),cosmo); 
    if(opt=='M')    printf(" %8.4f %8.4f  %10g  %10g/deg\n",mlim-.125,mlim+.125,Nbin,Nbin/41252.);
    else            printf(" %8.4f %8.4f  %10g  %10g/deg\n",-0.4*(mlim-.125),-0.4*(mlim+.125),Nbin,Nbin/41252.);
      if(plots) cpgpt1(mlim,Nbin/41252.,3);
  }
  
  printf(" Total number of objects: %e\n\n",Ntot);
  printf("-----------------------------\n");
  printf(" Objects/degree: %e\n",Ntot/41252.);
  printf("-----------------------------\n");
}


void Calc_Num_wC(void)
{ 
  /* like Calc_Num but using Int_sch_M_wC instead of Int_sch_M */
  static double zup=0.5;
  static double zlow=0;
  double ztmp=0;
  static double mlimite=20;
  double mlim;
  double color_mean=0., color_stddev=0.;
  double Mabs;
  
  double Ntot;
  double Nbin;
  static char opt='M';

  double Mabsmin,Mabsmax;
  double mlimmin,mlimmax;

  static int plots=0; /* para no hacer las gr�ficas */
  
  printf(" Do you want plots of the LF (1=yes, 0=no)?\n");
  plots=readf(plots);
  
  if(plots) {
     cpgbeg(0,"?",2,1);
  }

  printf("\n L Give LF in luminosity\n"); 
  printf(" M Give LF in magnitudes\n");
  opt=readc(opt); 
  
  switch (opt) { 
  case 'L':
  case 'l':
    printf(" Introduce apparent flux limit (Watt/m2 , 1 erg/s/cm2 = 10^(-3) Watt/m2): ");
    mlimite=-2.5*log10(readf(pow(10.,-0.4*mlimite)));
    mlim=mlimite;
    set_Schechter_L();
    break;
  default:
    opt='M';
    printf(" Introduce apparent magnitude limit: ");
    mlimite=readf(mlimite);
    mlim=mlimite;
    set_Schechter_M();
    printf("Introduce mean color:");
    color_mean=readf(color_mean);
    printf("Introduce stddev for the color:");
    color_stddev=readf(color_stddev);
    break;
  }

  printf(" Introduce lower redshift limit: ");
  if(zlow<=ZMIN) zlow=0;
  zlow=readf(zlow);
/*   //scanf("%f",&zlow); */
  if(zlow==0) zlow=ZMIN;
  printf(" Introduce upper redshift limit (0=none, magnitude limited survey: ");
  zup=readf(zup);
  /*   //scanf("%f",&zup); */
  
  if(zup==0) zup=10.;

/*   //  printf(" Phi star %f\n",phistar); */

   printf(" Surveyed volume %e Mpc3\n",(Vol(zup,cosmo)-Vol(zlow,cosmo))/1.e18);

  /* Ahora hago un plot de las galaxias que se van a detectar en funcion 
     de la magnmitud absoluta */
  if(plots){
     cpgpanl(2,1);
     cpgscf(2);
     cpgswin(-10,-24.,0.,1000000./41252.);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpgsch(1.5);
     cpglab("Absolute magnitude","Number/degree","Number of objects/degree per 0.5 bin in magnitude");
     cpgsch(3.);
  }

  if(opt=='M') {
    printf("#Objects per absolute magnitud interval\n");
    printf("#Mmin   Mmax    Zmin  Zmax     N_sky     N/deg\n");
    printf("############################################\n");
    Mabsmin=-7;
    Mabsmax=-24;
  }
  else {
    printf("#Objects per luminosity interval logLstar %f\n",-0.4*(schlf_L.Lstar-.25));
    printf("#log(Lmin)   log10(Lmax)   Zmin  Zmax     N_sky     N/deg\n");
    printf("############################################\n");
    Mabsmin=-80;
    Mabsmax=-120;
  }
  for(Mabs=Mabsmin;Mabs>Mabsmax;Mabs-=.5) {
    if(opt=='M')ztmp=Z_m(mlim,Mabs,cosmo);
    else        ztmp=Z_l(pow(10.,-0.4*mlim),pow(10.,-0.4*Mabs),cosmo);
    if(ztmp<zlow) ztmp=zlow;
    if(ztmp>zup) ztmp=zup;
/*     //printf(" Absolute magnitude: %f Upper redshift limit: %f  %f\n",Mabs,ztmp,Mag_LF(ztmp,mlim)); */
   Nbin=0.5*Schechter_M(Mabs,schlf_M)*(Vol(ztmp,cosmo)-Vol(zlow,cosmo))/1.e18;
    if(opt=='M')     printf(" %6.2f %6.2f   %5.3f %5.3f  %10g  %10g/deg\n",Mabs-.25,Mabs+.25,zlow,ztmp,Nbin,Nbin/41252.);
    else             printf(" %6.2f %6.2f   %5.3f %5.3f  %10g  %10g/deg\n",-0.4*(Mabs-.25),-0.4*(Mabs+.25),zlow,ztmp,Nbin,Nbin/41252.);
      if(plots) cpgpt1(Mabs,Nbin/41252.,3);
  }

  if(opt=='M') {
    printf("#Objects per absolute magnitud interval\n");
    printf("#Mmin   Mmax    Zmin  Zmax     N_sky     N/deg\n");
    printf("############################################\n");
    Mabsmin=-9;
    Mabsmax=-24;
  }
  else {
    printf("#Objects per luminosity interval\n");
    printf("#log(Lmin)   log10(Lmax)   Zmin  Zmax     N_sky     N/deg\n");
    printf("############################################\n");
    Mabsmin=-80;
    Mabsmax=-120;
  }
  for(Mabs=Mabsmin;Mabs>Mabsmax;Mabs-=1.) {
    if(opt=='M')ztmp=Z_m(mlim,Mabs,cosmo);
    else        ztmp=Z_l(pow(10.,-0.4*mlim),pow(10.,-0.4*Mabs),cosmo);
    if(ztmp>zup) ztmp=zup;
    if(ztmp<zlow) ztmp=zlow;
    /*     //printf(" Absolute magnitude: %f Upper redshift limit: %f  %f\n",Mabs,ztmp,Mag_LF(ztmp,mlim)); */
       Nbin=Schechter_M(Mabs,schlf_M)*(Vol(ztmp,cosmo)-Vol(zlow,cosmo))/1.e18;
    if(opt=='M')     printf(" %6.2f %6.2f   %5.3f %5.3f  %10g  %10g/deg\n",Mabs-.25,Mabs+.25,zlow,ztmp,Nbin,Nbin/41252.);
    else             printf(" %6.2g %6.2g   %5.3f %5.3f  %10g  %10g/deg\n",-0.4*(Mabs-.25),-0.4*(Mabs+.25),zlow,ztmp,Nbin,Nbin/41252.);
      if(plots) cpgpt1(Mabs,Nbin/41252.,3);
  }
  Nbin=0;

  for(Mabs=-18.75;Mabs>-24.25;Mabs-=.5) {
    if(opt=='M')ztmp=Z_m(mlim,Mabs,cosmo);
    else        ztmp=Z_l(pow(10.,-0.4*mlim),pow(10.,-0.4*Mabs),cosmo);
    if(ztmp>zup) ztmp=zup;
    if(ztmp<zlow) ztmp=zlow;
    /*     printf(" Absolute magnitude: %f Upper redshift limit: %f  %f\n",Mabs,ztmp,Mag_LF(ztmp,mlim)); */
    Nbin+=0.5*Schechter_M(Mabs,schlf_M)*(Vol(ztmp,cosmo)-Vol(zlow,cosmo))/1.e18;
      if(plots) cpgpt1(Mabs,Nbin/41252.,3);
  }
/*   printf(" -24.25 -18.75   %5.3f %5.3f  %10g  %10g/deg\n",zlow,ztmp,Nbin,Nbin/41252.); */


/*   //  Lo mismo pero en magnitud aparente  */
  if(opt=='M') {
    printf("#Objects per apparent magnitud interval\n");
    printf("#mmin   mmax       N_sky     N/deg\n");
    printf("############################################\n");
    mlimmin=15;
    mlimmax=28;
  }
  else {
    printf("#Objects per flux interval\n");
    printf("#log(Fmin)   log10(Fmax)   N_sky     N/deg\n");
    printf("############################################\n");
    mlimmin=35;
    mlimmax=50;
  }

  if(plots) cpgpanl(1,1);
  
  if(opt=='M')   Ntot=Int_sch_M_wC(schlf_M,zlow,zup,color_mean,color_stddev,mlim,cosmo);
  else           Ntot=Int_sch_L(schlf_L,zlow,zup,pow(10.,-0.4*mlim),cosmo); /* FIXME: A�adir L_wC */
      
  if(plots){
    cpgscf(2);
    cpgswin(13.,28.,0.,40000000./41252.);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
    cpgsch(1.5);
    cpglab("Limiting magnitude","Number/degree","Number of objects/degree per 0.5 bin in magnitude");
    cpgsch(3.);
  }
  for(mlim=mlimmin;mlim<mlimmax;mlim+=.25) {
    if(opt=='M')  Nbin=-Int_sch_M_wC(schlf_M,zlow,zup,color_mean,color_stddev,mlim-.125,cosmo)+Int_sch_M_wC(schlf_M,zlow,zup,color_mean,color_stddev,mlim+.125,cosmo); 
    else          Nbin=-Int_sch_L(schlf_L,zlow,zup,pow(10.,-0.4*(mlim-.125)),cosmo)+Int_sch_L(schlf_L,zlow,zup,pow(10.,-0.4*(mlim+.125)),cosmo); 
    if(opt=='M')    printf(" %8.4f %8.4f  %10g  %10g/deg\n",mlim-.125,mlim+.125,Nbin,Nbin/41252.);
    else            printf(" %8.4f %8.4f  %10g  %10g/deg\n",-0.4*(mlim-.125),-0.4*(mlim+.125),Nbin,Nbin/41252.);
      if(plots) cpgpt1(mlim,Nbin/41252.,3);
  }
  
  printf(" Total number of objects: %e\n\n",Ntot);
  printf("-----------------------------\n");
  printf(" Objects/degree: %e\n",Ntot/41252.);
  printf("-----------------------------\n");
}


void set_lf_parameters(struct lf_param *lf) 
{
  
  static char opt='M';
  double pi=3.1415926535897932384;
  
  printf("\n L Compute LF in luminosity\n"); 
  printf(" M Compute LF in magnitudes\n");
  opt=readc(opt); 
  switch (opt) { 
  case 'L':
  case 'l':
    lf->islum=1;
    printf(" Input minimum log10(luminosity): ");
    lf->lum_min=(double)readf(lf->lum_min);
    printf(" Input maximum log10(luminosity): ");
    lf->lum_max=(double)readf(lf->lum_max);
    
    break;
  default:
    printf(" Input minimum magnitude : ");
    lf->mag_min=(double)readf(lf->mag_min);
    printf(" Input maximum magnitude : ");
    lf->mag_max=(double)readf(lf->mag_max);
    lf->islum=0;
    break;
  }
  
  printf(" Input number of bins    : ");
  lf->npoint=readi(lf->npoint);
  printf(" Input lower redshift: ");
  lf->zlow=(double)readf(lf->zlow); 
  printf(" Input limiting redshift (0=none, magnitude limited sample): ");
  lf->zup=(double)readf(lf->zup); 


  printf(" Input area covered by the survey (square degrees): ");
  lf->area=readf(lf->area*41252/4/pi);

  lf->area=lf->area/41252*4*pi;
}

void set_lf_ceg(struct lum_func_ceg *lf)
{
  
  static char fsfile[101]="";
  static char surveyfile[101]="";

  printf(" Name of file with selection function: ");
  reads(fsfile,fsfile);
  printf(" Name of survey database: ");
  reads(surveyfile,surveyfile);

  readposelfunc(&(lf->fsel),fsfile);
  ReadSurDB(surveyfile,&(lf->sdb));

  printf(" Limiting EW: ");
  lf->ewlim=readd(lf->ewlim);

  printf(" Photometric band: ");
  reads(lf->photband,lf->photband);
  printf(" QE of broad passband at emission line wavelenght: ");
  lf->gamma=readd(lf->gamma);
  printf(" Delta of passband for broad band: ");
  lf->delta=readd(lf->delta);
  printf(" Ratio of continuum ar emission line and broad band: ");
  lf->Kcoc=readd(lf->Kcoc);

  
}



void get_sample(struct sample_data *sample)
{

  char opt;
  static char datafile[300]="";
  static int colz=2,collum=3,colmag=4;
  double *z,*magni,*lum;
  int *logic1,*logic2,*logic3;
  int ndat;
  int i,j=0;

  printf("\n F Get sample data by file\n"); 
  printf(" K Get sample data by keyboard\n");
  opt=readc('F');  
  switch (opt) { 
  case 'F':
  case 'f':
    printf(" Input file name with data: ");
    reads(datafile,datafile);
    printf(" Input column with z data: ");
    colz=readi(colz);
    /* if(sample.islum) { */
    printf(" Input column with apparent luminosity (fluxes in W/m2)  data: ");
    collum=readi(collum);
    /*     }else{ */
    printf(" Input column with apparent magnitude data: ");
    colmag=readi(colmag);
    /*     } */
    ndat=FileNLin(datafile);
    z  =malloc(ndat*sizeof(double));
    lum=malloc(ndat*sizeof(double));
    magni=malloc(ndat*sizeof(double));
    logic1=malloc(ndat*sizeof(int));
    logic2=malloc(ndat*sizeof(int));
    logic3=malloc(ndat*sizeof(int));
    sample->mag=malloc(ndat*sizeof(double));
    sample->lum=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    ReadDoublecol(datafile,colz  ,z  ,logic1,&ndat);
    ReadDoublecol(datafile,colmag,magni,logic2,&ndat);
    ReadDoublecol(datafile,collum,lum,logic3,&ndat);
    for(i=0;i<ndat;i++) {
/*       //printf(" z %f log %d\n",z[i],log[i]); */
      if(logic1[i] && logic2[i] && logic3[i]) {
/* 	//printf(" YESSS\n"); */
	sample->mag[j]=(double)magni[i];
	sample->lum[j]=(double)lum[i];
	sample->z[j]  =(double)z[i];
/* 	//printf(" Por aqui %d %f %f \n",j,sample->mag[j],sample->z[j]); */
	j++;
      }
    }
    sample->ngalax=j;
    break;
  case 'K':
  case 'k':
    j=0;
    sample->z  =malloc(sizeof(double));
    sample->mag=malloc(sizeof(double));
    sample->lum=malloc(sizeof(double));
    opt='y';
    while(opt=='y' ){
      printf(" Input redshift  : ");
      sample->z[j]=(double)readf(0.);
      printf(" Input magnitude : ");
      sample->mag[j]=(double)readf(0.);
      printf(" Input luminosity (W): ");
      sample->lum[j]=(double)readf(0.);
      j++;
      sample->z  =realloc(sample->z  ,(j+1)*sizeof(double));
      sample->mag=realloc(sample->mag,(j+1)*sizeof(double));
      sample->lum=realloc(sample->lum,(j+1)*sizeof(double));
      
      printf(" Other?: ");
      opt=readc('y');
	
    }
    sample->ngalax=j+1;
    break;
  }

}

void get_sample_sel_dist(struct sample_data_sel_dist *sample) 
{
  char opt;
  static char datafile[300]="";
  static int colz=3,collum_sel=1,collum_dist=2,colmag_sel=1,colmag_dist=2;
  double *z,*magSel,*magDist,*lumSel,*lumDist;
  int *logic1,*logic2,*logic3,*logic4,*logic5;
  int ndat;
  int i,j=0;

  printf("\n F Get sample data by file\n"); 
  printf(" K Get sample data by keyboard\n");
  opt=readc('F');  
  switch (opt) { 
  case 'F':
  case 'f':
    printf(" Input file name with data: ");
    reads(datafile,datafile);
    printf(" Input column with z data: ");
    colz=readi(colz);
    /* if(sample.islum) { */
    printf(" Input column with selection apparent luminosity (fluxes in W/m2) data: ");
    collum_sel=readi(collum_sel);
    /*     }else{ */
    printf(" Input column with selection apparent magnitude data: ");
    colmag_sel=readi(colmag_sel);
    printf(" Input column with apparent distribution luminosity (fluxes in W/m2)  data: ");
    collum_dist=readi(collum_dist);
    /*     }else{ */
    printf(" Input column with apparent distribution magnitude data: ");
    colmag_dist=readi(colmag_dist);
    /*     } */
    ndat=FileNLin(datafile);
    z  =malloc(ndat*sizeof(double));
    lumSel=malloc(ndat*sizeof(double));
    magSel=malloc(ndat*sizeof(double));
    lumDist=malloc(ndat*sizeof(double));
    magDist=malloc(ndat*sizeof(double));
    logic1=malloc(ndat*sizeof(int));
    logic2=malloc(ndat*sizeof(int));
    logic3=malloc(ndat*sizeof(int));
    logic4=malloc(ndat*sizeof(int));
    logic5=malloc(ndat*sizeof(int));
    sample->magSel=malloc(ndat*sizeof(double));
    sample->lumSel=malloc(ndat*sizeof(double));
    sample->magDist=malloc(ndat*sizeof(double));
    sample->lumDist=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    ReadDoublecol(datafile,colz  ,z  ,logic1,&ndat);
    ReadDoublecol(datafile,colmag_sel,magSel,logic2,&ndat);
    ReadDoublecol(datafile,collum_sel,lumSel,logic3,&ndat);
    ReadDoublecol(datafile,colmag_dist,magDist,logic4,&ndat);
    ReadDoublecol(datafile,collum_dist,lumDist,logic5,&ndat);
    for(i=0;i<ndat;i++) { 
      /*       //printf(" z %f log %d\n",z[i],log[i]); */
      if(logic1[i] && logic2[i] && logic3[i] && logic4[i] && logic5[i]) {
        /* 	//printf(" YESSS\n"); */
        sample->magSel[j]=(double)magSel[i];
        sample->lumSel[j]=(double)lumSel[i];
        sample->magDist[j]=(double)magDist[i];
        sample->lumDist[j]=(double)lumDist[i];
        sample->z[j]  =(double)z[i];
        /* 	//printf(" Por aqui %d %f %f \n",j,sample->mag[j],sample->z[j]); */
        j++;
      }
    }
    sample->ngalax=j;
    break;
  case 'K':
  case 'k':
    j=0;
    sample->z  =malloc(sizeof(double));
    sample->magSel=malloc(sizeof(double));
    sample->lumSel=malloc(sizeof(double));
    sample->magDist=malloc(sizeof(double));
    sample->lumDist=malloc(sizeof(double));
    opt='y';
    while(opt=='y' ){
      printf(" Input redshift  : ");
      sample->z[j]=(double)readf(0.);
      printf(" Input selection magnitude : ");
      sample->magSel[j]=(double)readf(0.);
      printf(" Input selection luminosity (W): ");
      sample->lumSel[j]=(double)readf(0.);
      printf(" Input calculation magnitude : ");
      sample->magDist[j]=(double)readf(0.);
      printf(" Input calculation luminosity (W): ");
      sample->lumDist[j]=(double)readf(0.);
      j++;
      sample->z  =realloc(sample->z  ,(j+1)*sizeof(double));
      sample->magSel=realloc(sample->magSel,(j+1)*sizeof(double));
      sample->lumSel=realloc(sample->lumSel,(j+1)*sizeof(double));
      sample->magDist=realloc(sample->magDist,(j+1)*sizeof(double));
      sample->lumDist=realloc(sample->lumDist,(j+1)*sizeof(double));
      
      printf(" Other?: ");
      opt=readc('y');
	
    }
    sample->ngalax=j+1;
    break;
  }

}

void get_sample_wC_errColor(struct sample_data_wC_errColor *sample) 
{
  char opt;
  static char datafile[300]="";
  static int colz=2,colmag_sel=1,colmag_dist=2,colerrColor=3;
  double *z,*magSel,*magDist,*errColor;
  int *logic1,*logic2,*logic3,*logic4;

  int ndat;
  int i,j=0;

  printf("\n F Get sample data by file\n"); 
  printf(" K Get sample data by keyboard\n");
  opt=readc('F');  
  switch (opt) { 
  case 'F':
  case 'f':
    printf(" Input file name with data: ");
    reads(datafile,datafile);
    printf(" Input column with z data: ");
    colz=readi(colz);
    printf(" Input column with selection apparent magnitude data: ");
    colmag_sel=readi(colmag_sel);
    printf(" Input column with apparent calculation magnitude data: ");
    colmag_dist=readi(colmag_dist);
    printf(" Input column with error in color data: ");
    colerrColor=readi(colerrColor);
    ndat=FileNLin(datafile);
    z  =malloc(ndat*sizeof(double));
    magSel=malloc(ndat*sizeof(double));
    magDist=malloc(ndat*sizeof(double));
    errColor=malloc(ndat*sizeof(double));
    logic1=malloc(ndat*sizeof(int));
    logic2=malloc(ndat*sizeof(int));
    logic3=malloc(ndat*sizeof(int));
    logic4=malloc(ndat*sizeof(int));
    sample->magSel  =malloc(ndat*sizeof(double));
    sample->magDist =malloc(ndat*sizeof(double));
    sample->z       =malloc(ndat*sizeof(double));
    sample->errColor=malloc(ndat*sizeof(double));
    ReadDoublecol(datafile,colz  ,z  ,logic1,&ndat);
    ReadDoublecol(datafile,colmag_sel,magSel,logic2,&ndat);
    ReadDoublecol(datafile,colmag_dist,magDist,logic3,&ndat);
    ReadDoublecol(datafile,colerrColor,errColor,logic4,&ndat);
    for(i=0;i<ndat;i++) { 
/*       //printf(" z %f log %d\n",z[i],log[i]); */
      if(logic1[i] && logic2[i] && logic3[i] && logic4[i]) {
/* 	//printf(" YESSS\n"); */
	sample->magSel[j]  =(double)magSel[i];
	sample->magDist[j] =(double)magDist[i];
	sample->z[j]       =(double)z[i];
        sample->errColor[j]=(double)errColor[i];
/* 	//printf(" Por aqui %d %f %f \n",j,sample->mag[j],sample->z[j]); */
	j++;
      }
    }
    /* printf("ngalax dentro de get_sample: %i",j); */
    sample->ngalax=j;
    break;
  case 'K':
  case 'k':
    printf(" Not available.\n");
    break;
  }

}

void get_sample_mag_err(struct sample_data_mag_err *sample)
{
/* dabreu */
/* para leer cat�logo con magnitud y error en la magnitud */
  char opt;
  static char datafile[300]="";
  static int colz=2,colmag=3,colmagerr=4;
  double *z,*magni,*magerr;
  int *logic1,*logic2,*logic3;
  int ndat;
  int i,j=0;

  printf("\n F Get sample data by file\n"); 
  printf(" K Get sample data by keyboard\n");
  opt=readc('F');  
  switch (opt) { 
  case 'F':
  case 'f':
    printf(" Input file name with data: ");
    reads(datafile,datafile);
    printf(" Input column with z data: ");
    colz=readi(colz);
    printf(" Input column with apparent magnitude data: ");
    colmag=readi(colmag);
    printf(" Input column with magnitude error data: ");
    colmagerr=readi(colmagerr);
    
    ndat=FileNLin(datafile);
    z  =malloc(ndat*sizeof(double));
    magni=malloc(ndat*sizeof(double));
    magerr=malloc(ndat*sizeof(double));
    logic1=malloc(ndat*sizeof(int));
    logic2=malloc(ndat*sizeof(int));
    logic3=malloc(ndat*sizeof(int));
    sample->mag=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    sample->mag_err=malloc(ndat*sizeof(double));
    ReadDoublecol(datafile,colz  ,z  ,logic1,&ndat);
    ReadDoublecol(datafile,colmag,magni,logic2,&ndat);
    ReadDoublecol(datafile,colmagerr,magerr,logic3,&ndat);
    for(i=0;i<ndat;i++) {
      if(logic1[i] && logic2[i] && logic3[i]) {
	sample->mag[j]=(double)magni[i];
	sample->z[j]  =(double)z[i];
	sample->mag_err[j]  =(double)magerr[i];
	j++;
      }
    }
    sample->ngalax=j;
    break;
  case 'K':
  case 'k':
    j=0;
    sample->z  =malloc(sizeof(double));
    sample->mag=malloc(sizeof(double));
    sample->mag_err=malloc(sizeof(double));
    opt='y';
    while(opt=='y' ){
      printf(" Input redshift  : ");
      sample->z[j]=(double)readf(0.);
      printf(" Input magnitude : ");
      sample->mag[j]=(double)readf(0.);
      printf(" Input magnitude error : ");
      sample->mag_err[j]=(double)readf(0.);
      j++;
      sample->z  =realloc(sample->z  ,(j+1)*sizeof(double));
      sample->mag=realloc(sample->mag,(j+1)*sizeof(double));
      sample->mag_err=realloc(sample->mag_err,(j+1)*sizeof(double));
      
      printf(" Other?: ");
      opt=readc('y');
	
    }
    sample->ngalax=j+1;
    break;
  }

}

void get_sample_gmz_wC (struct sample_data_gmz_wC *sample)
{
  /* to read catalog with magsel, magdist, errmagdist, z, errz */
  char opt;
  static char datafile[300]="";
  static int colz=2, colerrz=3,colmagSel=4,colmagDist=5,colerrMagDist=6;
  double *z,*errz,*magSel,*magDist,*errMagDist;
  int *logic1,*logic2,*logic3,*logic4,*logic5;

  int ndat;
  int i,j=0;

  printf("\n F Get sample data by file\n"); 
  printf(" K Get sample data by keyboard\n");
  opt=readc('F');  
  switch (opt)
  { 
  case 'F':
  case 'f':
    printf(" Input file name with data: ");
    reads(datafile,datafile);
    printf(" Input column with z data: ");
    colz=readi(colz);
    printf(" Input column with errz data: ");
    colerrz=readi(colerrz);
    printf(" Input column with selection magnitude data: ");
    colmagSel=readi(colmagSel);
    printf(" Input column with distribution magnitude data: ");
    colmagDist=readi(colmagDist);
    printf(" Input column with distribution magnitude error data: ");
    colerrMagDist=readi(colerrMagDist);
    
    ndat=FileNLin(datafile);
    z          =malloc(ndat*sizeof(double));
    errz       =malloc(ndat*sizeof(double));
    magSel     =malloc(ndat*sizeof(double));
    magDist    =malloc(ndat*sizeof(double));
    errMagDist =malloc(ndat*sizeof(double));
    logic1=malloc(ndat*sizeof(int));
    logic2=malloc(ndat*sizeof(int));
    logic3=malloc(ndat*sizeof(int));
    logic4=malloc(ndat*sizeof(int));
    logic5=malloc(ndat*sizeof(int));
    sample->magSel=malloc(ndat*sizeof(double));
    sample->magDist=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    sample->errz  =malloc(ndat*sizeof(double));
    sample->errMagDist=malloc(ndat*sizeof(double));
    ReadDoublecol(datafile,colz  ,z  ,logic1,&ndat);
    ReadDoublecol(datafile,colerrz,errz,logic2,&ndat);
    ReadDoublecol(datafile,colmagSel,magSel,logic3,&ndat);
    ReadDoublecol(datafile,colmagDist,magDist,logic4,&ndat);
    ReadDoublecol(datafile,colerrMagDist,errMagDist,logic5,&ndat);
    for(i=0;i<ndat;i++)
    {
      if(logic1[i] && logic2[i] && logic3[i] && logic4[i] && logic5[i])
      {
        sample->z[j]=(double)z[i];
        sample->errz[j]=(double)errz[i];
	sample->magSel[j]=(double)magSel[i];
	sample->magDist[j]  =(double)magDist[i];
	sample->errMagDist[j]  =(double)errMagDist[i];
	j++;
      }
    }
    sample->ngalax=j;
    break;
  case 'K':
  case 'k':
    //j=0;
    printf("You are really crazy...\n");
    printf("...try with other method with less numbers.\n");
    /*sample->z  =malloc(sizeof(double));
    sample->mag=malloc(sizeof(double));
    sample->mag_err=malloc(sizeof(double));
    opt='y';
    while(opt=='y' ){
      printf(" Input redshift  : ");
      sample->z[j]=(double)readf(0.);
      printf(" Input magnitude : ");
      sample->mag[j]=(double)readf(0.);
      printf(" Input magnitude error : ");
      sample->mag_err[j]=(double)readf(0.);
      j++;
      sample->z  =realloc(sample->z  ,(j+1)*sizeof(double));
      sample->mag=realloc(sample->mag,(j+1)*sizeof(double));
      sample->mag_err=realloc(sample->mag_err,(j+1)*sizeof(double));
      
      printf(" Other?: ");
      opt=readc('y');
	
    }
    sample->ngalax=j+1; */
    break;
  }
}

void get_sample_ceg(struct sample_data *sample)
{

  char opt;
  static char datafile[100]="";
  static int colz=2,collum=3,colmag=4,colew=5,colimage=6;
  double *z,*magni,*lum,*ew;
  char *image;
  int *logic1,*logic2,*logic3,*logic4,*logic5;
  int ndat;
  int i,j=0;

  printf("\n F Get sample data by file\n");
  printf(" K Get sample data by keyboard\n");
  opt=readc('F');  
  switch (opt) { 
  case 'F':
  case 'f':
    printf(" Input file name with data: ");
    reads(datafile,datafile);
    printf(" Input column with z data: ");
    colz=readi(colz);
    printf(" Input column with apparent magnitude data: ");
    colmag=readi(colmag);
    printf(" Input column with equivalent width: ");
    colew=readi(colew);
    printf(" Input column with image: ");
    colimage=readi(colimage);
    /*     } */
    ndat=FileNLin(datafile);
    z    =malloc(ndat*sizeof(double));
    lum  =malloc(ndat*sizeof(double));
    magni  =malloc(ndat*sizeof(double));
    ew   =malloc(ndat*sizeof(double));
    image=malloc(ndat*51*sizeof(char));
    logic1=malloc(ndat*sizeof(int));
    logic2=malloc(ndat*sizeof(int));
    logic3=malloc(ndat*sizeof(int));
    logic4=malloc(ndat*sizeof(int));
    logic5=malloc(ndat*sizeof(int));
    sample->mag=malloc(ndat*sizeof(double));
    sample->lum=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    sample->ew =malloc(ndat*sizeof(double));
    sample->image=vector_s(ndat,51);
    ReadDoublecol(datafile,colz  ,z  ,logic1,&ndat);
    ReadDoublecol(datafile,colmag,mag,logic2,&ndat);
    ReadDoublecol(datafile,collum,lum,logic3,&ndat);
    ReadDoublecol(datafile,colew ,ew ,logic4,&ndat);
    ReadCharcol  (datafile,colimage ,image ,logic5,51,&ndat);
    for(i=0;i<ndat;i++) {
      /*       //printf(" z %f log %d\n",z[i],log[i]); */
      if(logic1[i] && logic2[i] && logic3[i] && logic4[i]) {
        /* 	//printf(" YESSS\n"); */
        sample->mag[j]=(double)magni[i];
        sample->lum[j]=(double)lum[i];
        sample->z[j]  =(double)z[i];
        sample->ew[j]  =(double)ew[i];
        strcpy(sample->image[j],image+i*51);
        /* 	//printf(" Por aqui %d %f %f \n",j,sample->mag[j],sample->z[j]); */
        j++;
      }
    }
    sample->ngalax=j;
    break;
  case 'K':
  case 'k':
    j=0;
    sample->z  =malloc(sizeof(double));
    sample->mag=malloc(sizeof(double));
    sample->lum=malloc(sizeof(double));
    sample->ew =malloc(sizeof(double));
    opt='y';
    while(opt=='y' ){
      printf(" Input redshift  : ");
      sample->z[j]=(double)readf(0.);
      printf(" Input magnitude : ");
      sample->mag[j]=(double)readf(0.);
      printf(" Input equivalente width: ");
      sample->ew[j]=(double)readf(0.);
      j++;
      sample->z  =realloc(sample->z  ,(j+1)*sizeof(double));
      sample->mag=realloc(sample->mag,(j+1)*sizeof(double));
      sample->lum=realloc(sample->lum,(j+1)*sizeof(double));
      sample->ew =realloc(sample->ew ,(j+1)*sizeof(double));
      
      printf(" Other?: ");
      opt=readc('y');
	
    }
    sample->ngalax=j+1;
    break;
  }

}


void set_cosmology(void)
{
  static double H0=70;
  static double OM=0.3;
  static double OL=0.7;
  printf(" Input H0: ");
  H0=readf(H0);
  printf(" Input Omega Matter: ");
  OM=readf(OM);
  printf(" Input Omega Lambda: ");
  OL=readf(OL);
  cosmo_init(&cosmo, H0, OM, OL);
}

void set_Schechter_M(void)
{
  printf(" Schechter parameters function: \n");
  printf(" M star: ");
  schlf_M.Mstar=readf(schlf_M.Mstar);
  printf(" Phi star: ");
  schlf_M.phistar=readf(schlf_M.phistar);
  printf(" alfa: ");
  schlf_M.alfa=readf(schlf_M.alfa);

}

void set_Schechter_L(void)
{

  printf(" Schechter parameters function: \n");
  printf(" L star (Watts): ");
  /* dabreu
  schlf.Mstar=-2.5*log10(readd(pow(10.,-0.4*schlf.Mstar))); */
  schlf_L.Lstar=readd(schlf_L.Lstar);
  printf(" Phi star: ");
  schlf_L.phistar=readf(schlf_L.phistar);
  printf(" alfa: ");
  schlf_L.alfa=readf(schlf_L.alfa);
}

void Generate_Cat_M(void)
{
  static double mlow=0,mup=0,merror_mean=0,merror_stddev=0;
  double Mlow,Mup;
  static double zlow=0,zup=0.5;
  static double area=41252.;
  /* dabreu */
  double zerror_mean=0, zerror_stddev=0;
  double *Msample,*zsample,*msample,*merror,*mobserved;
  double *zerror,*zobserved;
  static int plots=0; /* para no hacer las gr�ficas */

  unsigned int i;
/*   double xx[1000],yy[1000]; */
  double nobjAllSky, nobjMean;
  unsigned int nobj;
/*   char cnul; */
/*   double kkk=1.1; */
/*   double fnul; */
  char snul[1000];
/*   double *zgrid; */
/*   double *probgrid; */
  FILE *fout;
  static char filename[100]="";
  double mh1,mh2,Mh1,Mh2;

  /* Setup random number generator */
  gsl_rng* rng;
  struct timeval tv;  /* algo falta porque me da un error con el timeval */

  rng = gsl_rng_alloc(gsl_rng_default);
  gettimeofday(&tv, NULL);
  gsl_rng_set(rng,tv.tv_usec);

  /* Setup cosmology and LF parameters */
  set_cosmology();
  set_Schechter_M();
  printf(" Input lower limit redshift: ");
  zlow=readf(zlow);
  zlow = (zlow < ZMIN ? ZMIN : zlow);
/*   //  printf(" zlow %f\n",zlow); */
  printf(" Input upper limit redshift (0=none, magnitude limiting survey): ");
  zup= readf(zup);
  if(zup==0) zup=10.;

  /* dabreu */
  printf(" Input the mean error for the redshift: ");
  zerror_mean = readf(zerror_mean);
  printf(" Input the mean deviation for the error: ");
  zerror_stddev = readf(zerror_stddev);

  printf(" Input the limiting magnitude: ");
  mlow=readf(mlow);
  printf(" Input the maximum magnitude: ");
  mup =readf(mup);
  /* dabreu */
  printf(" Input the mean error for the magnitud: ");
  merror_mean = readf(merror_mean);
  printf(" Input the mean deviation for the error: ");
  merror_stddev = readf(merror_stddev);
    
  printf(" Input area covered by the survey (square degrees): ");
  area=readf(area);

  /* dabreu */
  printf(" Do you want plots of the catalog (1=yes, 0=no)?\n");
  plots=readf(plots);

  /*   printf(" aaaaaaaaaLIM mag %f %f\n",mlow,mup); */
  nobjAllSky=Int_sch_M(schlf_M, zlow, zup, mlow, cosmo); 
  printf("nobjAllSky %g\n", nobjAllSky);
  nobjMean=(nobjAllSky/41252.*area);
  nobj=gsl_ran_poisson(rng, nobjMean);
  
  zsample=malloc(nobj*sizeof(double));
  Msample=malloc(nobj*sizeof(double));
  msample=malloc(nobj*sizeof(double));
  /* dabreu */
  zerror=malloc(nobj*sizeof(double));
  zobserved=malloc(nobj*sizeof(double));
  merror=malloc(nobj*sizeof(double));
  mobserved=malloc(nobj*sizeof(double));

  printf(" Number of galaxies generated: %d\n",nobj);

/*   Teordist(schlf.Mstar,schlf.alfa,zlow,zup,mlow,mup) */

  printf("\n 000000000 / %9d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",nobj);
  for(i=0;i<nobj;i++) {
    printf("%9d\b\b\b\b\b\b\b\b\b",i);
/*     //printf(" %d\n",i); */
    zsample[i]= zSchdev_M(schlf_M, zlow, zup, mlow, mup, cosmo);
    Mlow=Mag(zsample[i],mlow,cosmo);
    Mup=Mag(zsample[i],mup,cosmo);
    Msample[i]= Schechterdev_M(schlf_M,Mlow,Mup);
    msample[i]= mag(zsample[i],Msample[i],cosmo);
    /* dabreu */
    while((merror[i] = merror_mean + Gasdev() * merror_stddev) < 0);
    mobserved[i] = msample[i] + Gasdev()*merror[i];
    while((zerror[i] = zerror_mean + Gasdev() * zerror_stddev) < 0);
    while((zobserved[i] = zsample[i] + Gasdev()*zerror[i]) < 0);
  }
  MinMax_d(nobj,msample,&mh1,&mh2);
  MinMax_d(nobj,Msample,&Mh1,&Mh2);
  
  if(plots) {
     printf(" Now, some nice plots about the catalogue generated\n");
     cpgbeg(0,"?",2,2);
     /*   cpgpanl(1,1); */
     cpgswin(zlow*.8,zup*1.2,0.,nobj/5);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,zsample,zlow*.8,zup*1.2,20,1); 
     cpglab("Redshift","Number of galaxies"," Redshift distribution");
     /*   cpgpage(); */
     cpgpanl(2,1);
     cpgswin(mh1-1.,mh2+1.,0.,nobj/2);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,msample,mh1-1.,mh2+1.,20,1); 
     cpglab("Apparent magnitude","Number of galaxies"," Apparente magnitude distribution");
     /*     cpgpage(); */
     cpgpanl(1,2);
     cpgswin(Mh1-1.,Mh2+1.,0.,nobj/2);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,Msample,Mh1-1.,Mh2+1.,20,1); 
     cpglab("Absolute magnitude","Number of galaxies"," Absolute magnitude distribution");
  }
     
  printf(" Input file name to write simulated catalog: ");
  reads(filename,filename);
  if((fout=fopen(filename,"w")) ==NULL) {
    printf(" Couldn't open %s for writing\n",filename);
    return;
  }
  fprintf(fout,"# Simulated catalog generated by Lumfunc with %d objects\n",nobj);
  fprintf(fout,"# with redshift limits: %f-%f. Apparent magnitude limits: %f-%f\n",zlow,zup,mlow,mup);
  fprintf(fout,"# Schechter function: Mstar %f Phistar %f Alfa %f\n",schlf_M.Mstar,schlf_M.phistar,schlf_M.alfa);
  fprintf(fout,"# Cosmology H0= %f  OM= %f  OL = %f\n",cosmo.H0,cosmo.OM,cosmo.OL);
  fprintf(fout,"# Square degrees covered = %f  (%f fraction of the sky)\n",area,area/41252.);
  /*fprintf(fout,"# Redshift   m_app   Mabs     merror   mobserved\n"); dabreu */
  fprintf(fout,"#\n"
          "# 1 Z (redshift without observational errors)\n"
	  "# 2 Z_OBS (observed redshift taking into account observational errors\n"
          "# 3 Z_ERR (error in redshift)\n"
          "# 4 M_APP (apparent magnitude without observational errors)\n"
          "# 5 M_ABS (absolute magnitude)\n"
	  "# 6 M_ERR (error in apparent magnitude)\n"
          "# 7 M_OBS (observed magnitudes taking into account observational errors\n"); 
  for(i=0;i<nobj;i++) {
   /* fprintf(fout," %8.6f  %8.3f %8.3f",zsample[i],msample[i],Msample[i]); */
    /* dabreu */
    fprintf(fout," %8.6f  %8.6f %8.6f",zsample[i],zobserved[i],zerror[i]);
    fprintf(fout," %8.6f %8.6f",msample[i],Msample[i]);
    fprintf(fout," %8.6f  %8.6f\n",merror[i],mobserved[i]);
  }
  fclose(fout);

/*   for(i=0;i<1000;i++) { */
/*     xx[i]=mlow+(mup-mlow)*i/999; */
/*     yy[i]=Schechter_M_LF(xx[i])*nobj/nbin*abs(mlow-mup)*kkk/schlf.phistar*fglobal; */
/*     //printf(" xx %f yy %g\n",xx[i],yy[i]); */
/*   } */
/*   for(i=0;i<nobj;i++) { */
/*     printf(" %d MAg %f\n",i,Msample[i]); */
/*   } */
  sprintf(snul,"Random galaxy distribution following a Schechter distribution with \\ga= %5.2f",schlf_M.alfa);
/*   cpgpage(); */
/*   cpgswin(mup,mlow+3,0.,nobj/nbin*abs(mup-mlow)/5.); */
/*   cpglab("Absolute Magnitude","Number of galaxies",snul); */
/*  cpgask(0); */ 
/*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
/* cpgline_d(1000,xx,yy); */ 
/*   cpghist_d(nobj,Msample,mup,mlow,nbin,1); */
  cpgend(); 
  cosmo_free(&cosmo);

}


/* es una copia de Generate_Cat_M ? */
void Generate_Cat_M_C(void)
{
  static double mlow=0,mup=0,merror_mean=0,merror_stddev=0;
  double Mlow,Mup;
  static double zlow=0,zup=0.5;
  static double area=41252.;
  /* dabreu */
  double zerror_mean=0, zerror_stddev=0;
  double *Msample,*zsample,*msample,*merror,*mobserved;
  double *zerror,*zobserved;
  static int plots=0; /* para no hacer las gr�ficas */

  unsigned int i;
/*   double xx[1000],yy[1000]; */
  double nobjAllSky, nobjMean;
  unsigned int nobj;
/*   char cnul; */
/*   double kkk=1.1; */
/*   double fnul; */
  char snul[1000];
/*   double *zgrid; */
/*   double *probgrid; */
  FILE *fout;
  static char filename[100]="";
  double mh1,mh2,Mh1,Mh2;

  /* Setup random number generator */
  gsl_rng* rng;
  struct timeval tv;  /* algo falta porque me da un error con el timeval */

  rng = gsl_rng_alloc(gsl_rng_default);
  gettimeofday(&tv, NULL);
  gsl_rng_set(rng,tv.tv_usec);
  
  set_cosmology();
  set_Schechter_M();
  printf(" Input lower limit redshift: ");
  zlow=readf(zlow);
  zlow = (zlow < ZMIN ? ZMIN : zlow);
/*   //  printf(" zlow %f\n",zlow); */
  printf(" Input upper limit redshift (0=none, magnitude limiting survey): ");
  zup= readf(zup);
  if(zup==0) zup=10.;

  /* dabreu */
  printf(" Input the mean error for the redshift: ");
  zerror_mean = readf(zerror_mean);
  printf(" Input the mean deviation for the error: ");
  zerror_stddev = readf(zerror_stddev);

  printf(" Input the limiting magnitude: ");
  mlow=readf(mlow);
  printf(" Input the maximum magnitude: ");
  mup =readf(mup);
  /* dabreu */
  printf(" Input the mean error for the magnitud: ");
  merror_mean = readf(merror_mean);
  printf(" Input the mean deviation for the error: ");
  merror_stddev = readf(merror_stddev);
    
  printf(" Input area covered by the survey (square degrees): ");
  area=readf(area);

  /* dabreu */
  printf(" Do you want plots of the catalog (1=yes, 0=no)?\n");
  plots=readf(plots);

/*  nobj=(int)Int_sch_M(schlf_M, zlow, zup, mlow, cosmo); 
  nobj=(int)((double)nobj/41252.*area);
  nobj=Poidev(nobj); */
  /* dabreu */
  nobjAllSky=Int_sch_M(schlf_M, zlow, zup, mlow, cosmo); 
  printf("nobjAllSky %g\n", nobjAllSky);
  nobjMean=(nobjAllSky/41252.*area);
  nobj=gsl_ran_poisson(rng, nobjMean);

  zsample=malloc(nobj*sizeof(double));
  Msample=malloc(nobj*sizeof(double));
  msample=malloc(nobj*sizeof(double));
  /* dabreu */
  zerror=malloc(nobj*sizeof(double));
  zobserved=malloc(nobj*sizeof(double));
  merror=malloc(nobj*sizeof(double));
  mobserved=malloc(nobj*sizeof(double));

  printf(" Number of galaxies generated: %d\n",nobj);

/*   Teordist(schlf.Mstar,schlf.alfa,zlow,zup,mlow,mup) */

  printf("\n 000000000 / %9d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",nobj);
  for(i=0;i<nobj;i++) {
    printf("%9d\b\b\b\b\b\b\b\b\b",i);
/*     //printf(" %d\n",i); */
    zsample[i]= zSchdev_M(schlf_M, zlow, zup, mlow, mup, cosmo);
    Mlow=Mag(zsample[i],mlow,cosmo);
    Mup=Mag(zsample[i],mup,cosmo);
    Msample[i]= Schechterdev_M(schlf_M,Mlow,Mup);
    msample[i]= mag(zsample[i],Msample[i],cosmo);
    /* dabreu */
    while((merror[i] = merror_mean + Gasdev() * merror_stddev) < 0);
    mobserved[i] = msample[i] + Gasdev()*merror[i];
    while((zerror[i] = zerror_mean + Gasdev() * zerror_stddev) < 0);
    while((zobserved[i] = zsample[i] + Gasdev()*zerror[i]) < 0);
  }
  MinMax_d(nobj,msample,&mh1,&mh2);
  MinMax_d(nobj,Msample,&Mh1,&Mh2);
  
  if(plots) {
     printf(" Now, some nice plots about the catalogue generated\n");
     cpgbeg(0,"?",2,2);
     /*   cpgpanl(1,1); */
     cpgswin(zlow*.8,zup*1.2,0.,nobj/5);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,zsample,zlow*.8,zup*1.2,20,1); 
     cpglab("Redshift","Number of galaxies"," Redshift distribution");
     /*   cpgpage(); */
     cpgpanl(2,1);
     cpgswin(mh1-1.,mh2+1.,0.,nobj/2);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,msample,mh1-1.,mh2+1.,20,1); 
     cpglab("Apparent magnitude","Number of galaxies"," Apparente magnitude distribution");
     /*     cpgpage(); */
     cpgpanl(1,2);
     cpgswin(Mh1-1.,Mh2+1.,0.,nobj/2);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,Msample,Mh1-1.,Mh2+1.,20,1); 
     cpglab("Absolute magnitude","Number of galaxies"," Absolute magnitude distribution");
  }
     
  printf(" Input file name to write simulated catalog: ");
  reads(filename,filename);
  if((fout=fopen(filename,"w")) ==NULL) {
    printf(" Couldn't open %s for writing\n",filename);
    return;
  }
  fprintf(fout,"# Simulated catalog generated by Lumfunc with %d objects\n",nobj);
  fprintf(fout,"# with redshift limits: %f-%f. Apparent magnitude limits: %f-%f\n",zlow,zup,mlow,mup);
  fprintf(fout,"# Schechter function: Mstar %f Phistar %f Alfa %f\n",schlf_M.Mstar,schlf_M.phistar,schlf_M.alfa);
  fprintf(fout,"# Cosmology H0= %f  OM= %f  OL = %f\n",cosmo.H0,cosmo.OM,cosmo.OL);
  fprintf(fout,"# Square degrees covered = %f  (%f fraction of the sky)\n",area,area/41252.);
  /*fprintf(fout,"# Redshift   m_app   Mabs     merror   mobserved\n"); dabreu */
  fprintf(fout,"#\n"
          "# 1 Z (redshift without observational errors)\n"
          "# 2 Z_ERR (error in redshift)\n"
	  "# 3 Z_OBS (observed redshift taking into account observational errors\n"
          "# 4 M_APP (apparent magnitude without observational errors)\n"
          "# 5 M_ABS (absolute magnitude)\n"
	  "# 6 M_ERR (error in apparent magnitude)\n"
          "# 7 M_OBS (observed magnitudes taking into account observational errors\n"); 
  for(i=0;i<nobj;i++) {
   /* fprintf(fout," %8.6f  %8.3f %8.3f",zsample[i],msample[i],Msample[i]); */
    /* dabreu */
    fprintf(fout," %8.6f  %8.6f %8.6f",zsample[i],zerror[i],zobserved[i]);
    fprintf(fout," %8.3f %8.3f",msample[i],Msample[i]);
    fprintf(fout," %8.6f  %8.6f\n",merror[i],mobserved[i]);
  }
  fclose(fout);

/*   for(i=0;i<1000;i++) { */
/*     xx[i]=mlow+(mup-mlow)*i/999; */
/*     yy[i]=Schechter_M_LF(xx[i])*nobj/nbin*abs(mlow-mup)*kkk/schlf.phistar*fglobal; */
/*     //printf(" xx %f yy %g\n",xx[i],yy[i]); */
/*   } */
/*   for(i=0;i<nobj;i++) { */
/*     printf(" %d MAg %f\n",i,Msample[i]); */
/*   } */
  sprintf(snul,"Random galaxy distribution following a Schechter distribution with \\ga= %5.2f",schlf_M.alfa);
/*   cpgpage(); */
/*   cpgswin(mup,mlow+3,0.,nobj/nbin*abs(mup-mlow)/5.); */
/*   cpglab("Absolute Magnitude","Number of galaxies",snul); */
/*  cpgask(0); */ 
/*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
/* cpgline_d(1000,xx,yy); */ 
/*   cpghist_d(nobj,Msample,mup,mlow,nbin,1); */
  cpgend(); 
  cosmo_free(&cosmo);
} 

/* dabreu */
/* toda la funci�n Generate_Cat_L() es una adaptaci�n de Generate_Cat() */

void Generate_Cat_L(void)
{
  static double fluxlow=0.,fluxup=0.,fluxerr_mean=0.,fluxerr_stddev=0.;
  double Flow,Fup;
  static double zlow=0,zup=0.5;
  static double area=41252.;
  double *Lsample,*zsample,*fluxsample,*fluxerror,*fluxobserved;
  static int plots=0; /* para no hacer las gr�ficas */

  unsigned int i;
  double nobjAllSky, nobjMean;
  unsigned int nobj;
  char snul[1000];
  FILE *fout;
  static char filename[100]="";
  double fh1,fh2,Fh1,Fh2;

  /* Setup random number generator */
  gsl_rng* rng;
  struct timeval tv;  /* algo falta porque me da un error con el timeval */

  rng = gsl_rng_alloc(gsl_rng_default);
  gettimeofday(&tv, NULL);
  gsl_rng_set(rng,tv.tv_usec);

  /* Setup cosmology and LF parameters */
  set_cosmology();
  set_Schechter_L();
  printf(" Input lower limit redshift: ");
  zlow=readf(zlow);
  zlow = (zlow < ZMIN ? ZMIN : zlow);
/*   //  printf(" zlow %f\n",zlow); */
  printf(" Input upper limit redshift (0=none, magnitude limiting survey): ");
  zup= readf(zup);
  if(zup==0) zup=10.;

  printf(" Input the limiting flux: ");
  fluxlow=readd(fluxlow);
  printf(" Input the maximum flux: ");
  fluxup =readd(fluxup);

  /* no funcionan los errores, creo que por error en lectura o algo as� */
  printf(" Input the mean error for the flux: ");
  fluxerr_mean = readd(fluxerr_mean);
  printf(" Input the mean deviation for the error: ");
  fluxerr_stddev = readd(fluxerr_stddev);
    
  printf(" Input area covered by the survey (square degrees): ");
  area=readf(area);

  printf(" Do you want plots of the catalog (1=yes, 0=no)?\n");
  plots=readf(plots);
  
  nobjAllSky=(int)Int_sch_L(schlf_L, zlow, zup, fluxlow, cosmo);
  printf("nobjAllSky %g\n", nobjAllSky);
  nobjMean=(nobjAllSky/41252.*area);
  nobj=gsl_ran_poisson(rng, nobjMean);

  zsample=malloc(nobj*sizeof(double));
  Lsample=malloc(nobj*sizeof(double));
  fluxsample=malloc(nobj*sizeof(double));

  fluxerror=malloc(nobj*sizeof(double));
  fluxobserved=malloc(nobj*sizeof(double));

  printf(" Number of galaxies generated: %d\n",nobj);

  printf("\n 000000000 / %9d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",nobj);
  for(i=0;i<nobj;i++) {
    printf("%9d\b\b\b\b\b\b\b\b\b",i);
    /* utilizamos zSchdev_L de modulos/schechterdev.c */
    zsample[i]= zSchdev_L(schlf_L,zlow,zup,fluxlow,fluxup,cosmo);
    Flow=Lum(zsample[i],fluxlow,cosmo);
    Fup=Lum(zsample[i],fluxup,cosmo);
    Lsample[i]= Schechterdev_L(schlf_L,Flow,Fup);
    fluxsample[i]=Flux(zsample[i],Lsample[i],cosmo);

/*    printf("fluxerr_mean = %g\n",fluxerr_mean); */

    while((fluxerror[i] = fluxerr_mean + Gasdev() * fluxerr_stddev) < 0);
    while((fluxobserved[i] = fluxsample[i] + Gasdev()*fluxerror[i]) < 0);
  }
  MinMax_d(nobj,fluxsample,&fh1,&fh2);
  MinMax_d(nobj,Lsample,&Fh1,&Fh2);
  
  if(plots) {
     printf(" Now, some nice plots about the catalogue generated\n");
     cpgbeg(0,"?",2,2);
     cpgswin(zlow*.8,zup*1.2,0.,nobj/5);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,zsample,zlow*.8,zup*1.2,20,1); 
     cpglab("Redshift","Number of galaxies"," Redshift distribution");
     cpgpanl(2,1);
     cpgswin(fh1-1.,fh2+1.,0.,nobj/2);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,fluxsample,fh1-1.,fh2+1.,20,1);
     cpglab("Apparent magnitude","Number of galaxies"," Apparente magnitude distribution");
     cpgpanl(1,2);
     cpgswin(Fh1-1.,Fh2+1.,0.,nobj/2);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,Lsample,Fh1-1.,Fh2+1.,20,1); 
     cpglab("Absolute magnitude","Number of galaxies"," Absolute magnitude distribution");
  }
     
  printf(" Input file name to write simulated catalog: ");
  reads(filename,filename);
  if((fout=fopen(filename,"w")) ==NULL) {
    printf(" Couldn't open %s for writing\n",filename);
    return;
  }
  fprintf(fout,"# Simulated catalog generated by Lumfunc with %d objects\n",nobj);
  fprintf(fout,"# with redshift limits: %f-%f. Apparent flux limits:%f-%f\n",zlow,zup,fluxlow,fluxup);
  fprintf(fout,"# Schechter function: Lstar %f Phistar %f Alfa %f\n",schlf_L.Lstar,schlf_L.phistar,schlf_L.alfa);
  fprintf(fout,"# Cosmology H0= %f OM= %f OL = %f\n",cosmo.H0,cosmo.OM,cosmo.OL);
  fprintf(fout,"# Square degrees covered = %f  (%f fraction of the sky)\n",area,area/41252.);
  /*fprintf(fout,"# Redshift   m_app   Mabs     merror   mobserved\n"); dabreu */
  fprintf(fout,"#\n# 1 Z  (redshift)\n# 2 F_APP (apparent flux without observational errors)\n# "
          "3 L_ABS (absolute luminosity)\n# 4 F_ERR (error in apparent flux)\n# "
          "5 F_OBS (observed flux taking into account observational errors\n"); 
  for(i=0;i<nobj;i++) {
    fprintf(fout," %8.6f  %10.5g %10.5g",zsample[i],fluxsample[i],Lsample[i]);
    fprintf(fout," %10.5g  %10.5g\n",fluxerror[i],fluxobserved[i]);
  }
  fclose(fout);

  sprintf(snul,"Random galaxy distribution following a Schechter distribution with \\ga= %5.2f",schlf_L.alfa);
  cpgend(); 
  cosmo_free(&cosmo);
}

/* Generate_Cat_M_wC -> generate catalogs in two bands for selection using Color distribution */

void Generate_Cat_M_wC(void)
{
  static double mSelUp=0;
  static double mSelError_mean=0, mSelError_stddev=0;
  static double mDistLow=0, mDistUp=0, MDistLow=0, MDistUp=0;
  static double color_mean=0, color_stddev=0;
  static double colorError_mean=0, colorError_stddev=0;
  //double Mlow,Mup;
  static double zlow=0,zup=0.5;
  static double area=41252.;
  static struct fermifsel_M fsel={21.,0.};
  double zerror_mean=0, zerror_stddev=0;
  double *Msample,*zsample,*msample,*merror,*mobserved;
  double *zerror,*zobserved;
  double *colorsample;
  double *MDist, *mDist, *mDistObserved, *mSelError, *mSelObserved, *mSel;
  double *colorError, *colorObserved;
  static int plots=0; /* not to do graphs */

  /* color = mDist - mSel */

  unsigned int iobj = 0;
  double nobjAllSky, nobjMean;
  unsigned int nobj;
  char snul[1000];
  FILE *fout;
  static char filename[100]="";
  double mh1,mh2,Mh1,Mh2;

  /* Setup random number generator */
  gsl_rng* rng;
  struct timeval tv;  /* algo falta porque me da un error con el timeval */

  rng = gsl_rng_alloc(gsl_rng_default);
  gettimeofday(&tv, NULL);
  gsl_rng_set(rng,tv.tv_usec);

  /* Setup cosmology and LF parameters */
  set_cosmology();
  set_Schechter_M();
  printf(" Input lower limit redshift: ");
  zlow=readf(zlow);
  zlow = (zlow < ZMIN ? ZMIN : zlow);
/*   //  printf(" zlow %f\n",zlow); */
  printf(" Input upper limit redshift (0=none, magnitude limited survey): ");
  zup= readf(zup);
  if(zup==0) zup=10.;

  printf(" Input the mean error for the redshift: ");
  zerror_mean = readf(zerror_mean);
  printf(" Input the mean deviation for the error: ");
  zerror_stddev = readf(zerror_stddev);

  printf(" Input the limiting magnitude \n"
         "(.5 of detect probability in the magnitude of selection): ");
  fsel.magcut=readf(fsel.magcut);
  printf(" Input selection function sharpness (0 for a step function): ");
  fsel.deltamag=readf(fsel.deltamag);
  printf(" Input the maximum magnitude (in the magnitude of selection): ");
  mSelUp =readf(mSelUp);
  printf(" Input the mean error for the selection magnitude: ");
  mSelError_mean = readf(mSelError_mean);
  printf(" Input the mean deviation for the error of selection magnitude: ");
  mSelError_stddev = readf(mSelError_stddev);
  printf(" Input the mean color (distributed magnitude - selection magnitude): ");
  color_mean = readf(color_mean);
  printf(" Input the deviation for the color distribution: ");
  color_stddev = readf(color_stddev);
  printf(" Input the mean color error: ");
  colorError_mean = readf(colorError_mean);
  printf(" Input the dispersion in errors of the color: ");
  colorError_stddev = readf(colorError_stddev);
    
  printf(" Input area covered by the survey (square degrees): ");
  area=readf(area);

  printf(" Do you want plots of the catalog (1=yes, 0=no)?\n");
  plots=readf(plots);

  /*   printf(" aaaaaaaaaLIM mag %f %f\n",mlow,mup); */
  if (color_stddev==0.) /* to avoid problems in Int_sch_M_wC */
  {
    printf("color_stddev=0 -> calling Int_sch_M instead of Int_sch_M_wC\n");
    nobjAllSky=Int_sch_f_M(schlf_M, zlow, zup, fsel, cosmo);
  }
  else
  {
    if(fsel.deltamag == 0)
      nobjAllSky=Int_sch_M_wC(schlf_M, zlow, zup,
                              color_mean, color_stddev, fsel.magcut,
                              cosmo);
    else
      nobjAllSky=Int_sch_f_M_wC(schlf_M, zlow, zup,
                              color_mean, color_stddev, fsel, cosmo);
  }
  printf("Total number of objects in the sky %g\n", nobjAllSky);
  nobjMean=(nobjAllSky/41252.*area);
  printf("Number of objects per degree %.10g\n", nobjAllSky/41252.);
  printf("Mean number of objects in the survey %g\n", nobjMean);
  if(DEBUG2) printf("nobjMean %g\n",nobjMean);
  if(DEBUG2) printf("Before gsl_ran_poisson\n");
  nobj=gsl_ran_poisson(rng, nobjMean);
  if(DEBUG2) printf("After gsl_ran_poisson\n");
  
  zsample=malloc(nobj*sizeof(double));
  Msample=malloc(nobj*sizeof(double));
  msample=malloc(nobj*sizeof(double));
  zerror=malloc(nobj*sizeof(double));
  zobserved=malloc(nobj*sizeof(double));
  merror=malloc(nobj*sizeof(double));
  mobserved=malloc(nobj*sizeof(double));
  colorsample=malloc(nobj*sizeof(double));
  MDist=malloc(nobj*sizeof(double));
  mDist=malloc(nobj*sizeof(double));
  mDistObserved=malloc(nobj*sizeof(double));
  mSelError=malloc(nobj*sizeof(double));
  mSelObserved=malloc(nobj*sizeof(double));
  mSel=malloc(nobj*sizeof(double));
  colorError=malloc(nobj*sizeof(double));
  colorObserved=malloc(nobj*sizeof(double));

  printf(" Number of galaxies generated (after poisson dist): %d\n",nobj);

  printf("\n 000000000 / %9d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",nobj);
  while(iobj<nobj)
  {
    double probdetec;
    printf("%9d\b\b\b\b\b\b\b\b\b",iobj);
/*     //printf(" %d\n",i); */
    colorsample[iobj] = color_mean + Gasdev() * color_stddev;
    mDistLow = fsel.magcut + colorsample[iobj] + 10 * fsel.deltamag;
    mDistUp = mSelUp + colorsample[iobj];
    zsample[iobj]= zSchdev_M(schlf_M, zlow, zup, mDistLow, mDistUp, cosmo);
    MDistLow=Mag(zsample[iobj],mDistLow,cosmo);
    MDistUp =Mag(zsample[iobj],mDistUp,cosmo);
    MDist[iobj]= Schechterdev_M(schlf_M,MDistLow,MDistUp);
    mDist[iobj]= mag(zsample[iobj],MDist[iobj],cosmo);
    mSel[iobj] = mDist[iobj] - colorsample[iobj]; /* color = dist - sel */
    probdetec = Fermi(mSel[iobj], fsel.magcut, fsel.deltamag); 
    if(gsl_ran_binomial(rng, probdetec, 1) == 0)
    {
       /* Object is not detected. The binomial dist gives the number
          of sucesses in n trials (here trials is 1) */
       continue;
    }
    while((mSelError[iobj] = mSelError_mean + Gasdev() * mSelError_stddev) < 0);
    mSelObserved[iobj] = mSel[iobj] + Gasdev()*mSelError[iobj];
    while((colorError[iobj] = colorError_mean + Gasdev() * colorError_stddev) < 0);
    colorObserved[iobj] = colorsample[iobj] + Gasdev()*colorError[iobj];
    while((zerror[iobj] = zerror_mean + Gasdev() * zerror_stddev) < 0);
    while((zobserved[iobj] = zsample[iobj] + Gasdev()*zerror[iobj]) < 0);
    mDistObserved[iobj] = mSelObserved[iobj] + colorObserved[iobj];
    if(DEBUG2) printf("z = %g\n", zsample[iobj]);
    ++iobj;
  }
  MinMax_d(nobj,msample,&mh1,&mh2);
  MinMax_d(nobj,Msample,&Mh1,&Mh2);
  
  if(plots) {
     printf(" Now, some nice plots about the catalogue generated\n");
     cpgbeg(0,"?",2,2);
     /*   cpgpanl(1,1); */
     cpgswin(zlow*.8,zup*1.2,0.,nobj/5);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,zsample,zlow*.8,zup*1.2,20,1); 
     cpglab("Redshift","Number of galaxies"," Redshift distribution");
     /*   cpgpage(); */
     cpgpanl(2,1);
     cpgswin(mh1-1.,mh2+1.,0.,nobj/2);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,msample,mh1-1.,mh2+1.,20,1); 
     cpglab("Apparent magnitude","Number of galaxies"," Apparent magnitude distribution");
     /*     cpgpage(); */
     cpgpanl(1,2);
     cpgswin(Mh1-1.,Mh2+1.,0.,nobj/2);
     cpgbox("BCTNS",0,0,"BCTNS",0,0);
     cpghist_d(nobj,Msample,Mh1-1.,Mh2+1.,20,1); 
     cpglab("Absolute magnitude","Number of galaxies"," Absolute magnitude distribution");
  }
     
  printf(" Input file name to write simulated catalog: ");
  reads(filename,filename);
  if((fout=fopen(filename,"w")) ==NULL) {
    printf(" Couldn't open %s for writing\n",filename);
    return;
  }
  fprintf(fout,"# Simulated catalog generated by Lumfunc with %d objects\n",nobj);
  fprintf(fout,"# Mean number of objects for this LF: %f (%f in all sky)\n",nobjMean, nobjAllSky);
  fprintf(fout,"# Redshift limits: %f-%f.\n",zlow, zup);
  fprintf(fout,"# Fermi selection function: selection magcut %f deltamag %f\n",fsel.magcut,fsel.deltamag);
  fprintf(fout,"# Brightest apparent selection magnitude %f\n",mSelUp);
  fprintf(fout,"# Schechter function: Mstar %f Phistar %f Alfa %f\n",schlf_M.Mstar,schlf_M.phistar,schlf_M.alfa);
  fprintf(fout,"# Mean color (dist - sel): %f, Stddev color: %f\n",color_mean, color_stddev);
  fprintf(fout,"# Cosmology H0= %f OM = %f OL = %f\n",cosmo.H0,cosmo.OM,cosmo.OL);
  fprintf(fout,"# Square degrees covered = %f  (%f fraction of the sky)\n",area,area/41252.);
  /*fprintf(fout,"# Redshift   m_app   Mabs     merror   mobserved\n"); dabreu */
  fprintf(fout,"#\n"
          "# 1 Z (redshift without observational errors)\n"
	  "# 2 Z_OBS (observed redshift taking into account observational errors)\n"
          "# 3 Z_ERR (error in redshift)\n"
          "# 4 M_DIST_APP (apparent distribution magnitude)\n"
          "# 5 M_DIST_ABS (absolute distribution magnitude)\n"
          "# 6 M_SEL_APP  (apparent selection magnitud)\n"
          "# 7 M_SEL_APP_OBS (observed selection magnitudes taking into account observational errors)\n"
	  "# 8 M_SEL_ERR (observational error in apparent selection magnitude)\n"
          "# 9 COLOR (color for the object = mDist - mSel)\n"
          "# 10 COLOR_OBS (observed color for the object taking into account obs. errors)\n"
          "# 11 COLOR_ERR (error in the observed color for the object = mDist - mSel\n"
          "# 12 M_DIST_APP_OBS (M_SEL_APP_OBS + COLOR_OBS)\n"
          "# 13 M_DIST_APP_ERR (M_SEL_APP_ERR + COLOR_ERR\n");
          
  for(iobj=0;iobj<nobj;iobj++) {
   /* fprintf(fout," %8.6f  %8.3f %8.3f",zsample[i],msample[i],Msample[i]); */
    fprintf(fout," %8.6f  %8.6f %8.6f",zsample[iobj],zobserved[iobj],zerror[iobj]);
    fprintf(fout," %8.6f %8.6f",mDist[iobj],MDist[iobj]);
    fprintf(fout," %8.6f %8.6f",mSel[iobj],mSelObserved[iobj]);
    fprintf(fout," %8.6f",mSelError[iobj]);
    fprintf(fout," %8.6f", colorsample[iobj]);
    fprintf(fout," %8.6f", colorObserved[iobj]);
    fprintf(fout," %8.6f", colorError[iobj]);
    fprintf(fout," %8.6f", mDistObserved[iobj]);
    fprintf(fout," %8.6f\n", mSelError[iobj]+colorError[iobj]);
  }
  fclose(fout);

/*   for(i=0;i<1000;i++) { */
/*     xx[i]=mlow+(mup-mlow)*i/999; */
/*     yy[i]=Schechter_M_LF(xx[i])*nobj/nbin*abs(mlow-mup)*kkk/schlf.phistar*fglobal; */
/*     //printf(" xx %f yy %g\n",xx[i],yy[i]); */
/*   } */
/*   for(i=0;i<nobj;i++) { */
/*     printf(" %d MAg %f\n",i,Msample[i]); */
/*   } */
  sprintf(snul,"Random galaxy distribution following a Schechter distribution with \\ga= %5.2f",schlf_M.alfa);
/*   cpgpage(); */
/*   cpgswin(mup,mlow+3,0.,nobj/nbin*abs(mup-mlow)/5.); */
/*   cpglab("Absolute Magnitude","Number of galaxies",snul); */
/*  cpgask(0); */ 
/*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
/* cpgline_d(1000,xx,yy); */ 
/*   cpghist_d(nobj,Msample,mup,mlow,nbin,1); */
  cpgend(); 
  cosmo_free(&cosmo);
}

