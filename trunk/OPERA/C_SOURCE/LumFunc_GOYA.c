#include "modulos.h"  
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

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

/* creamos una estrucura para a�adir una m de seleccion y una m de c�lculo
dabreu */
struct sample_data_sel_cal 
{
  int ngalax;
  double *z;
  double *lum_sel;
  double *lum_cal;
  double *mag_sel;
  double *mag_cal;
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
void get_sample_sel_cal(struct sample_data_sel_cal *sample);
void get_sample_mag_err(struct sample_data_mag_err *sample);
void get_sample_ceg(struct sample_data *sample);
void set_lf_parameters(struct lf_param *lf);
void set_lf_ceg(struct lum_func_ceg *lf);
void set_cosmology();
void VVmax();
void STY();
void STY_MAG_ERR();
void SWML();
void CEG();
void Calc_Num();
void Calc_Num_wC();
void Generate_Cat_M();
void Generate_Cat_M_wC();
void Generate_Cat_L();
void set_Schechter_M();
void set_Schechter_L();



double c=299792.46; /* // En km/s */
/* //double alfa,phistar,Mstar; */
/* dabreu */
const double kk=144543977.0745928; /* pow(10.,-0.4*-20.4) */
struct Schlf_M schlf_M={-1.3,0.,0.0033,0.,-20.4,0.,0.,0.,0.};
struct cosmo_param cosmo={75,0.5};
/* utilizamos la struct de modulos.h */
struct Schlf_L schlf_L={-1.3,0.0,0.0033,0.0,144543977.0745928,0.0,0.0,0.0,0.0};

double fglobal;

int main()
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
      break;
    case 'D':
    case 'd':
      set_cosmology();
      Calc_Num_wC();
      break;
    case 'V':
    case 'v':
      set_cosmology();
      VVmax();
      break;
    case 'M':
    case 'm':
      set_cosmology();
      STY();
      break;
    case 'N':
    case 'n':
       set_cosmology();
       STY_MAG_ERR();
       break;
    case 'W':
    case 'w':
      set_cosmology();
      SWML();
      break;
    case 'L':
    case 'l':
      set_cosmology();
      CEG();
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

void STY() 
{
  static struct lf_param sty;
  struct sample_data sample;
  static double mlim=0;
  static double flim=0;
  double zlow;

  int iter;
  
  struct Schlf_M lfsch_M;
  struct Schlf_L lfsch_L;

  static int poissonflag=0;

  /* dabreu */
  static int plots=0; /* para no hacer las gr�ficas */
  FILE *fout;
  static char resultFileName[200]="";
  
  /* Information about the ML process */
  struct MLProcessInfo mlprocess;
  
  get_sample(&sample);  
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
    if(poissonflag) iter=MLA_STY_p_L(sample.ngalax,sample.lum,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    else            iter=  MLA_STY_L(sample.ngalax,sample.lum,sample.z,flim,sty.area,zlow,sty.zup,cosmo,&lfsch_L,&mlprocess);
    printf("\n Solution found at iteration %d\n",iter);
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
    fprintf(fout, "%g\t", log10(lfsch_L.Lstar));
    fprintf(fout, "%g\t", lfsch_L.alfa);
    fprintf(fout, "%g\t", lfsch_L.phistar);
    fprintf(fout, "%g\t", log10(lfsch_L.phistar));
    fprintf(fout, "%g\t", lfsch_L.errLstar/lfsch_L.Lstar/log(10.));
    fprintf(fout, "%g\t", lfsch_L.erralfa);
    fprintf(fout, "%g\t", lfsch_L.errphistar);
    fprintf(fout, "%g\n", lfsch_L.errphistar/lfsch_L.phistar/log(10.));
    fclose(fout);
  }
  else 
  {
    printf(" Input the limiting magnitude: ");
    mlim=readd(mlim);
    if(poissonflag) iter=MLA_STY_p_M(sample.ngalax,sample.mag,sample.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess);
    else            iter=  MLA_STY_M(sample.ngalax,sample.mag,sample.z,mlim,sty.area,zlow,sty.zup,cosmo,&lfsch_M,&mlprocess);
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
    fprintf(fout, "%g\t", lfsch_M.Mstar);
    fprintf(fout, "%g\t", lfsch_M.alfa);
    fprintf(fout, "%g\t", lfsch_M.phistar);
    fprintf(fout, "%g\t", log10(lfsch_M.phistar));
    fprintf(fout, "%g\t", lfsch_M.errMstar);
    fprintf(fout, "%g\t", lfsch_M.erralfa);
    fprintf(fout, "%g\t", lfsch_M.errphistar);
    fprintf(fout, "%g\n", lfsch_M.errphistar/lfsch_M.phistar/log(10.));
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

void CEG() 
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

void STY_MAG_ERR() 
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
  else { 
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


void SWML() 
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


void VVmax()
{
  int i,j;
  
  /* Variables de la funcion de lum */
  struct Steplf_M lf_vvmax_M;
  struct Steplf_L lf_vvmax_L;
  static struct lf_param vvmax_param={-40,40,0,0,0,10,41252,1,10};
  struct sample_data_sel_cal sample;
  struct Schlf_M  lfschfit_M;
  struct Schlf_L  lfschfit_L;

  static double mlim=1.;
  static double llim=1.;
  double M=0;
  double L=0;
/*   int vol; */
  double zmax;
  double *x,*y,*sigy,*yfit;
  double dm=0;
  double dl=0;

  /* Variables para el test V/Vmax */
  double *vvmaxtest;
  int   *ngaltest;           /* numero de galaxias en cada bin del test */
  double mmin_t=0,mmax_t=0;       /* magnitudes minimas y maximas desde las cuales se hara el test V/Vmax. */
  double lmin_t=0,lmax_t=0;       /* luminosidades minimas y maximas desde las cuales se hara el test V/Vmax. */
  int   n_t;                 /* numero de puntos para el test */
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
  double dum1,dum2,dum3;
  int nfit;

  /* Variables para los plots */
  int pg1,pg2;
  double ymin,ymax;
  /* dabreu */
  static int plots=0; /* para no hacer gr�ficas */

  /* Variables para escribir los resultados en un fichero */
  FILE *fout = NULL;
  static char schfitResultFileName[200]="";
  static char vvmaxResultFileName[200]="";
  
  /* Obtengo los datos y los parametros de la LF */
  get_sample_sel_cal(&sample);
  set_lf_parameters(&vvmax_param);

  /* Allocateo la LF */
  lf_vvmax_M.nbin      = vvmax_param.npoint;
  lf_vvmax_M.magni     =vector_d(lf_vvmax_M.nbin+1);
  lf_vvmax_M.errmagni  =vector_d(lf_vvmax_M.nbin+1);
  lf_vvmax_M.lnlf      =vector_d(lf_vvmax_M.nbin);
  lf_vvmax_M.errlnlf   =vector_d(lf_vvmax_M.nbin);
  lf_vvmax_M.covarlnlf =matrix_d(lf_vvmax_M.nbin,lf_vvmax_M.nbin);
  lf_vvmax_L.nbin      =vvmax_param.npoint;
  lf_vvmax_L.lumi      =vector_d(lf_vvmax_L.nbin+1);
  lf_vvmax_L.errlumi   =vector_d(lf_vvmax_L.nbin+1);
  lf_vvmax_L.lnlf      =vector_d(lf_vvmax_L.nbin);
  lf_vvmax_L.errlnlf   =vector_d(lf_vvmax_L.nbin);
  lf_vvmax_L.covarlnlf =matrix_d(lf_vvmax_L.nbin,lf_vvmax_L.nbin);

  printf(" Performing the V/Vmax test\n");
  if(vvmax_param.islum) 
  {
    MinMax_d(sample.ngalax,sample.lum_sel,&lmin_t,&lmax_t);
    printf(" Fluxes range: %g - %g\n",lmin_t,lmax_t);
    lmin_t*=0.3;lmax_t*=1.2;
  }
  else 
  {
    MinMax_d(sample.ngalax,sample.mag_sel,&mmin_t,&mmax_t);
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
	loglhist[j]=log10(sample.lum_sel[j]);
	L=Lum(sample.z[j],sample.lum_sel[j],cosmo);
	if(DEBUG2)        printf(" Galax %d z %f flux %g Lum %g mag %f Mag %f\n",j,sample.z[j],sample.lum_sel[j],L,sample.mag_sel[j],Mag(sample.z[j],sample.mag_sel[j],cosmo));
	if(L>1e36) 
 	printf(" Galax %d z %f flux %g Lum %g mag %f Mag %f\n",j,sample.z[j],sample.lum_sel[j],L,sample.mag_sel[j],Mag(sample.z[j],sample.mag_sel[j],cosmo)); 

      }
      else
      {
	mhist[j]=sample.mag_sel[j];
	M=Mag(sample.z[j],sample.mag_sel[j],cosmo);
	if(DEBUG2)        printf(" Galax %d z %f mag %f Mag %f\n",j,sample.z[j],sample.mag_sel[j],M);
      }
      /*       if(sample.mag[j]<mlim_t[i]) {   Esta es la antigua, cuando no podia hacer en luminosidades */
      /* Cuidado, sample.lum[j] (o mag[j]) puede no estar alocateado, 
	 por eso pongo antes vvmax.islum, que comprueba antes y si no lo cumple 
	 no pasa a la siguiente */
/*       if(( !vvmax.islum && sample.mag[j]<mlim_t[i])  || (vvmax.islum && sample.lum[j]<llim_t[i]) ) { */
      if(vvmax_param.islum)
      {
	if(sample.lum_sel[j]>llim_t[i])
        {
	  if(DEBUG2) printf(" Entro la flux=%f con llim=%f\n",sample.lum_sel[j],llim_t[i]);
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
	if(sample.mag_sel[j]<mlim_t[i])
        {
	  if(DEBUG2) printf(" Entro la m=%f con mlim=%f\n",sample.mag_sel[j],mlim_t[i]);
	  zmax=Z_m(mlim_t[i],M,cosmo);
	  if(zmax==0) zmax=ZMIN;
	  if(zmax>zup && zup!=0) zmax=zup;
	  if(zmax<zlow && zlow!=0) {printf(" zlow %f zmax %f\n ERROR!!\n",zlow,zmax);exit(1);}
	  if(DEBUG2) printf(" zmax %f %f %f %f\n",zmax, M, mlim,sample.z[j]); 
	  /* 		vvmaxtest[i]+=Vol_LF(sample.z[j])/Vol_LF(zmax);  VVmax tradicional */
	  vvmaxtest[i]+=(Vol(sample.z[j],cosmo)-Vol(zlow,cosmo))/(Vol(zmax,cosmo)-Vol(zlow, cosmo));  /* Cambiado para un zlow!=0 */
	  ngaltest[i]++;
	}
	if(DEBUG2) printf(" No entro la m=%f con mlim=%f\n",sample.mag_sel[j],mlim_t[i]);
      }
    }
    vvmaxtest[i]/=ngaltest[i];
    if(vvmax_param.islum)    printf(" Flux %g V/Vmax %f Ngal %d\n",llim_t[i],vvmaxtest[i],ngaltest[i]);
    else               printf(" Mag %g V/Vmax %f Ngal %d\n",mlim_t[i],vvmaxtest[i],ngaltest[i]);
  }
  /* dabreu */
  printf("Despues de 'Mag g V/Vmax f Ngal'\n");
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
      if(DEBUG) printf("KKK Limiting flux: %g\n",llim);
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
    VVmax_L(sample.ngalax, sample.lum_sel, sample.lum_cal, sample.z, llim, vvmax_param.area, vvmax_param.zlow, vvmax_param.zup, cosmo, &lf_vvmax_L);
    PrintStepLF_L(lf_vvmax_L);
  }
  else
  {
    VVmax_M(sample.ngalax, sample.mag_sel, sample.mag_cal, sample.z, mlim, vvmax_param.area, vvmax_param.zlow, vvmax_param.zup, cosmo, &lf_vvmax_M);
    PrintStepLF_M(lf_vvmax_M);
  }

  /* Calculo del ajuste a Schechter */
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
    fprintf(fout, "%g\t", lfschfit_M.Mstar);
    fprintf(fout, "%g\t", lfschfit_M.alfa);
    fprintf(fout, "%g\t", lfschfit_M.phistar);
    fprintf(fout, "%g\t", log10(lfschfit_M.phistar));
    fprintf(fout, "%g\t", lfschfit_M.errMstar);
    fprintf(fout, "%g\t", lfschfit_M.erralfa);
    fprintf(fout, "%g\t", lfschfit_M.errphistar);
    fprintf(fout, "%g\n", lfschfit_M.errphistar/lfschfit_M.phistar/log(10));
    fclose(fout);
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
    
      for(i=0;i<lf_vvmax_M.nbin;i++)
      {
        fprintf(fout, "%g\t", lf_vvmax_M.magni[i]);
        fprintf(fout, "%g\t", lf_vvmax_M.lnlf[i]);
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
  free_matrix_d(lf_vvmax_M.covarlnlf, lf_vvmax_M.nbin,lf_vvmax_M.nbin);
  free(lf_vvmax_L.lumi);
  free(lf_vvmax_L.errlumi);
  free(lf_vvmax_L.lnlf);
  free(lf_vvmax_L.errlnlf);
  free_matrix_d(lf_vvmax_L.covarlnlf, lf_vvmax_L.nbin,lf_vvmax_L.nbin);
}

void Calc_Num()
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


void Calc_Num_wC()
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
  static char datafile[200]="";
  static int colz=2,collum=3,colmag=4;
  double *z,*mag,*lum;
  int *log1,*log2,*log3;
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
    mag=malloc(ndat*sizeof(double));
    log1=malloc(ndat*sizeof(int));
    log2=malloc(ndat*sizeof(int));
    log3=malloc(ndat*sizeof(int));
    sample->mag=malloc(ndat*sizeof(double));
    sample->lum=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    ReadDoublecol(datafile,colz  ,z  ,log1,&ndat);
    ReadDoublecol(datafile,colmag,mag,log2,&ndat);
    ReadDoublecol(datafile,collum,lum,log3,&ndat);
    for(i=0;i<ndat;i++) {
/*       //printf(" z %f log %d\n",z[i],log[i]); */
      if(log1[i] && log2[i] && log3[i]) {
/* 	//printf(" YESSS\n"); */
	sample->mag[j]=(double)mag[i];
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

void get_sample_sel_cal(struct sample_data_sel_cal *sample) 
{
  char opt;
  static char datafile[200]="";
  static int colz=3,collum_sel=1,collum_cal=2,colmag_sel=1,colmag_cal=2;
  double *z,*mag_sel,*mag_cal,*lum_sel,*lum_cal;
  int *log1,*log2,*log3,*log4,*log5;
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
    printf(" Input column with apparent calculation luminosity (fluxes in W/m2)  data: ");
    collum_cal=readi(collum_cal);
    /*     }else{ */
    printf(" Input column with apparent calculation magnitude data: ");
    colmag_cal=readi(colmag_cal);
    /*     } */
    ndat=FileNLin(datafile);
    z  =malloc(ndat*sizeof(double));
    lum_sel=malloc(ndat*sizeof(double));
    mag_sel=malloc(ndat*sizeof(double));
    lum_cal=malloc(ndat*sizeof(double));
    mag_cal=malloc(ndat*sizeof(double));
    log1=malloc(ndat*sizeof(int));
    log2=malloc(ndat*sizeof(int));
    log3=malloc(ndat*sizeof(int));
    log4=malloc(ndat*sizeof(int));
    log5=malloc(ndat*sizeof(int));
    sample->mag_sel=malloc(ndat*sizeof(double));
    sample->lum_sel=malloc(ndat*sizeof(double));
    sample->mag_cal=malloc(ndat*sizeof(double));
    sample->lum_cal=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    ReadDoublecol(datafile,colz  ,z  ,log1,&ndat);
    ReadDoublecol(datafile,colmag_sel,mag_sel,log2,&ndat);
    ReadDoublecol(datafile,collum_sel,lum_sel,log3,&ndat);
    ReadDoublecol(datafile,colmag_cal,mag_cal,log4,&ndat);
    ReadDoublecol(datafile,collum_cal,lum_cal,log5,&ndat);
    for(i=0;i<ndat;i++) { 
/*       //printf(" z %f log %d\n",z[i],log[i]); */
      if(log1[i] && log2[i] && log3[i] && log4[i] && log5[i]) {
/* 	//printf(" YESSS\n"); */
	sample->mag_sel[j]=(double)mag_sel[i];
	sample->lum_sel[j]=(double)lum_sel[i];
	sample->mag_cal[j]=(double)mag_cal[i];
	sample->lum_cal[j]=(double)lum_cal[i];
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
    sample->mag_sel=malloc(sizeof(double));
    sample->lum_sel=malloc(sizeof(double));
    sample->mag_cal=malloc(sizeof(double));
    sample->lum_cal=malloc(sizeof(double));
    opt='y';
    while(opt=='y' ){
      printf(" Input redshift  : ");
      sample->z[j]=(double)readf(0.);
      printf(" Input selection magnitude : ");
      sample->mag_sel[j]=(double)readf(0.);
      printf(" Input selection luminosity (W): ");
      sample->lum_sel[j]=(double)readf(0.);
      printf(" Input calculation magnitude : ");
      sample->mag_cal[j]=(double)readf(0.);
      printf(" Input calculation luminosity (W): ");
      sample->lum_cal[j]=(double)readf(0.);
      j++;
      sample->z  =realloc(sample->z  ,(j+1)*sizeof(double));
      sample->mag_sel=realloc(sample->mag_sel,(j+1)*sizeof(double));
      sample->lum_sel=realloc(sample->lum_sel,(j+1)*sizeof(double));
      sample->mag_cal=realloc(sample->mag_cal,(j+1)*sizeof(double));
      sample->lum_cal=realloc(sample->lum_cal,(j+1)*sizeof(double));
      
      printf(" Other?: ");
      opt=readc('y');
	
    }
    sample->ngalax=j+1;
    break;
  }

}

void get_sample_mag_err(struct sample_data_mag_err *sample)
{
/* dabreu */
/* para leer cat�logo con magnitud y error en la magnitud */
  char opt;
  static char datafile[100]="";
  static int colz=2,colmag=3,colmagerr=4;
  double *z,*mag,*magerr;
  int *log1,*log2,*log3;
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
    colmag=readi(colmagerr);
    
    ndat=FileNLin(datafile);
    z  =malloc(ndat*sizeof(double));
    mag=malloc(ndat*sizeof(double));
    magerr=malloc(ndat*sizeof(double));
    log1=malloc(ndat*sizeof(int));
    log2=malloc(ndat*sizeof(int));
    log3=malloc(ndat*sizeof(int));
    sample->mag=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    sample->mag_err=malloc(ndat*sizeof(double));
    ReadDoublecol(datafile,colz  ,z  ,log1,&ndat);
    ReadDoublecol(datafile,colmag,mag,log2,&ndat);
    ReadDoublecol(datafile,colmagerr,magerr,log3,&ndat);
    for(i=0;i<ndat;i++) {
      if(log1[i] && log2[i] && log3[i]) {
	sample->mag[j]=(double)mag[i];
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

void get_sample_ceg(struct sample_data *sample)
{

  char opt;
  static char datafile[100]="";
  static int colz=2,collum=3,colmag=4,colew=5,colimage=6;
  double *z,*mag,*lum,*ew;
  char *image;
  int *log1,*log2,*log3,*log4,*log5;
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
    mag  =malloc(ndat*sizeof(double));
    ew   =malloc(ndat*sizeof(double));
    image=malloc(ndat*51*sizeof(char));
    log1=malloc(ndat*sizeof(int));
    log2=malloc(ndat*sizeof(int));
    log3=malloc(ndat*sizeof(int));
    log4=malloc(ndat*sizeof(int));
    log5=malloc(ndat*sizeof(int));
    sample->mag=malloc(ndat*sizeof(double));
    sample->lum=malloc(ndat*sizeof(double));
    sample->z  =malloc(ndat*sizeof(double));
    sample->ew =malloc(ndat*sizeof(double));
    sample->image=vector_s(ndat,51);
    ReadDoublecol(datafile,colz  ,z  ,log1,&ndat);
    ReadDoublecol(datafile,colmag,mag,log2,&ndat);
    ReadDoublecol(datafile,collum,lum,log3,&ndat);
    ReadDoublecol(datafile,colew ,ew ,log4,&ndat);
    ReadCharcol  (datafile,colimage ,image ,log5,51,&ndat);
    for(i=0;i<ndat;i++) {
/*       //printf(" z %f log %d\n",z[i],log[i]); */
      if(log1[i] && log2[i] && log3[i] && log4[i]) {
/* 	//printf(" YESSS\n"); */
	sample->mag[j]=(double)mag[i]; 
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


void set_cosmology()
{
  printf(" Input H0: ");
  cosmo.H0=readf(cosmo.H0);
  printf(" Input q0: ");
  cosmo.q0=readf(cosmo.q0);
  
}

void set_Schechter_M() 
{
  printf(" Schechter parameters function: \n");
  printf(" M star: ");
  schlf_M.Mstar=readf(schlf_M.Mstar);
  printf(" Phi star: ");
  schlf_M.phistar=readf(schlf_M.phistar);
  printf(" alfa: ");
  schlf_M.alfa=readf(schlf_M.alfa);

}

void set_Schechter_L() 
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

void Generate_Cat_M()
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

  int i;
/*   double xx[1000],yy[1000]; */
  int nobj;
/*   char cnul; */
/*   double kkk=1.1; */
/*   double fnul; */
  char snul[1000];
/*   double *zgrid; */
/*   double *probgrid; */
  FILE *fout;
  static char filename[100]="";
  double mh1,mh2,Mh1,Mh2;

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
  nobj=1000;
  nobj=(int)Int_sch_M(schlf_M, zlow, zup, mlow, cosmo); 
  nobj=(int)((double)nobj/41252.*area);
  nobj=Poidev(nobj);
  
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
    zobserved[i] = zsample[i] + Gasdev()*zerror[i];
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
  fprintf(fout,"# Cosmology H0= %f  q0= %f  Lambda = 0\n",cosmo.H0,cosmo.q0);
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

}


/* es una copia de Generate_Cat_M ? */
void Generate_Cat_M_C()
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

  int i;
/*   double xx[1000],yy[1000]; */
  double nobj_mean;
  int nobj;
/*   char cnul; */
/*   double kkk=1.1; */
/*   double fnul; */
  char snul[1000];
/*   double *zgrid; */
/*   double *probgrid; */
  FILE *fout;
  static char filename[100]="";
  double mh1,mh2,Mh1,Mh2;

  /* dabreu */
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
  nobj_mean=Int_sch_M(schlf_M, zlow, zup, mlow, cosmo); 
  nobj_mean=nobj_mean/41252.*area;
  nobj=(int)gsl_ran_poisson(rng, nobj_mean);

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
    zobserved[i] = zsample[i] + Gasdev()*zerror[i];
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
  fprintf(fout,"# Cosmology H0= %f  q0= %f  Lambda = 0\n",cosmo.H0,cosmo.q0);
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

} 

/* dabreu */
/* toda la funci�n Generate_Cat_L() es una adaptaci�n de Generate_Cat() */

void Generate_Cat_L()
{
  static double fluxlow=0.,fluxup=0.,fluxerr_mean=0.,fluxerr_stddev=0.;
  double Flow,Fup;
  static double zlow=0,zup=0.5;
  static double area=41252.;
  double *Lsample,*zsample,*fsample,*ferror,*fobserved; 
  static int plots=0; /* para no hacer las gr�ficas */

  int i;
  int nobj;
  char snul[1000];
  FILE *fout;
  static char filename[100]="";
  double fh1,fh2,Fh1,Fh2;

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
  
  nobj=1000;

  nobj=(int)Int_sch_L(schlf_L, zlow, zup, fluxlow, cosmo);
  nobj=(int)((double)nobj/41252.*area);
  nobj=Poidev(nobj);

  zsample=malloc(nobj*sizeof(double));
  Lsample=malloc(nobj*sizeof(double));
  fsample=malloc(nobj*sizeof(double));

  ferror=malloc(nobj*sizeof(double));
  fobserved=malloc(nobj*sizeof(double));

  printf(" Number of galaxies generated: %d\n",nobj);

  printf("\n 000000000 / %9d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",nobj);
  for(i=0;i<nobj;i++) {
    printf("%9d\b\b\b\b\b\b\b\b\b",i);
    /* utilizamos zSchdev_L de modulos/schechterdev.c */
    zsample[i]= zSchdev_L(schlf_L,zlow,zup,fluxlow,fluxup,cosmo);
    Flow=Lum(zsample[i],fluxlow,cosmo);
    Fup=Lum(zsample[i],fluxup,cosmo);
    Lsample[i]= Schechterdev_L(schlf_L,Flow,Fup);
    fsample[i]=Flux(zsample[i],Lsample[i],cosmo);

/*    printf("fluxerr_mean = %g\n",fluxerr_mean); */

    while((ferror[i] = fluxerr_mean + Gasdev() * fluxerr_stddev) < 0);
    while((fobserved[i] = fsample[i] + Gasdev()*ferror[i]) < 0);
  }
  MinMax_d(nobj,fsample,&fh1,&fh2);
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
     cpghist_d(nobj,fsample,fh1-1.,fh2+1.,20,1); 
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
  fprintf(fout,"# Cosmology H0= %f  q0= %f  Lambda = 0\n",cosmo.H0,cosmo.q0);
  fprintf(fout,"# Square degrees covered = %f  (%f fraction of the sky)\n",area,area/41252.);
  /*fprintf(fout,"# Redshift   m_app   Mabs     merror   mobserved\n"); dabreu */
  fprintf(fout,"#\n# 1 Z  (redshift)\n# 2 F_APP (apparent flux without observational errors)\n# "
          "3 L_ABS (absolute luminosity)\n# 4 F_ERR (error in apparent flux)\n# "
          "5 F_OBS (observed flux taking into account observational errors\n"); 
  for(i=0;i<nobj;i++) {
    fprintf(fout," %8.6f  %10.5g %10.5g",zsample[i],fsample[i],Lsample[i]);
    fprintf(fout," %10.5g  %10.5g\n",ferror[i],fobserved[i]);
  }
  fclose(fout);

  sprintf(snul,"Random galaxy distribution following a Schechter distribution with \\ga= %5.2f",schlf_L.alfa);
  cpgend(); 
}

/* Generate_Cat_M_wC -> generate catalogs in two bands for selection using Color distribution */

void Generate_Cat_M_wC()
{
  static double mSelLow=0, mSelUp=0, mDistError_mean=0, mDistError_stddev=0;
  static double mDistLow=0, mDistUp=0, MDistLow=0, MDistUp=0;
  static double color_mean=0, color_stddev=0;
  double Mlow,Mup;
  static double zlow=0,zup=0.5;
  static double area=41252.;
  double zerror_mean=0, zerror_stddev=0;
  double *Msample,*zsample,*msample,*merror,*mobserved;
  double *zerror,*zobserved;
  double *colorsample;
  double *MDist, *mDist, *mDistError, *mDistObserved, *mSel;
  static int plots=0; /* not to do graphs */

  /* color = mDist - mSel */

  int i;
  int nobj;
  char snul[1000];
  FILE *fout;
  static char filename[100]="";
  double mh1,mh2,Mh1,Mh2;

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

  printf(" Input the limiting magnitude (in the magnitude of selection): ");
  mSelLow=readf(mSelLow);
  printf(" Input the maximum magnitude (in the magnitude of selection): ");
  mSelUp =readf(mSelUp);
  printf(" Input the mean error for the originally distributed magnitude: ");
  mDistError_mean = readf(mDistError_mean);
  printf(" Input the mean deviation for the error of magnitude: ");
  mDistError_stddev = readf(mDistError_stddev);
  printf(" Input the mean color (distributed magnitude - selection magnitude): ");
  color_mean = readf(color_mean);
  printf(" Input the mean deviation for the error of the color: ");
  color_stddev = readf(color_stddev);
    
  printf(" Input area covered by the survey (square degrees): ");
  area=readf(area);

  printf(" Do you want plots of the catalog (1=yes, 0=no)?\n");
  plots=readf(plots);

  /*   printf(" aaaaaaaaaLIM mag %f %f\n",mlow,mup); */
  nobj=1000;
  nobj=(int)Int_sch_M_wC(schlf_M, zlow, zup, color_mean, color_stddev, mSelLow, cosmo); 
  nobj=(int)((double)nobj/41252.*area);
  nobj=Poidev(nobj);
  
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
  mDistError=malloc(nobj*sizeof(double));
  mDistObserved=malloc(nobj*sizeof(double));
  mSel=malloc(nobj*sizeof(double));

  printf(" Number of galaxies generated: %d\n",nobj);

  /* Setting the approx limiting magnitude in the originally distributed mag */
  mDistLow = mSelLow + color_mean + 5 * color_stddev;
  mDistUp  = mSelUp  + color_mean - 5 * color_stddev;
  
  printf("\n 000000000 / %9d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",nobj);
  for(i=0;i<nobj;i++) {
    printf("%9d\b\b\b\b\b\b\b\b\b",i);
/*     //printf(" %d\n",i); */
    zsample[i]= zSchdev_M(schlf_M, zlow, zup, mDistLow, mDistUp, cosmo);  /* cambiar por zSchdev_M_wC */
    MDistLow=Mag(zsample[i],mDistLow,cosmo);
    MDistUp =Mag(zsample[i],mDistUp,cosmo);
    MDist[i]= Schechterdev_M(schlf_M,MDistLow,MDistUp);
    mDist[i]= mag(zsample[i],MDist[i],cosmo);
    while((mDistError[i] = mDistError_mean + Gasdev() * mDistError_stddev) < 0);
    mDistObserved[i] = mDist[i] + Gasdev()*mDistError[i];
    while((zerror[i] = zerror_mean + Gasdev() * zerror_stddev) < 0);
    zobserved[i] = zsample[i] + Gasdev()*zerror[i];
    colorsample[i] = color_mean + Gasdev() * color_stddev;
    mSel[i] = mDist[i] - colorsample[i]; /* color = dist - sel */
    if(mSel[i] > mSelLow) /* This galaxy does not pass the selection function */
      i--; //FIXME; Esto hay que ponerlo mejor, as� no mola.
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
  fprintf(fout,"# with redshift limits: %f-%f. Apparent mag limits: %f-%f\n",zlow,zup,mDistLow,mDistUp);
  fprintf(fout,"# Schechter function: Mstar %f Phistar %f Alfa %f\n",schlf_M.Mstar,schlf_M.phistar,schlf_M.alfa);
  fprintf(fout,"# Cosmology H0= %f  q0= %f  Lambda = 0\n",cosmo.H0,cosmo.q0);
  fprintf(fout,"# Square degrees covered = %f  (%f fraction of the sky)\n",area,area/41252.);
  /*fprintf(fout,"# Redshift   m_app   Mabs     merror   mobserved\n"); dabreu */
  fprintf(fout,"#\n"
          "# 1 Z (redshift without observational errors)\n"
          "# 2 Z_ERR (error in redshift)\n"
	  "# 3 Z_OBS (observed redshift taking into account observational errors)\n"
          "# 4 M_APP (apparent magnitude without observational errors)\n"
          "# 5 M_ABS (absolute magnitude)\n"
	  "# 6 M_ERR (error in apparent magnitude)\n"
          "# 7 M_OBS (observed magnitudes taking into account observational errors)\n"
          "# 8 M_SEL (selection magnitud)\n"
          "# 9 COLOR (color for the object = mDist - mSel)\n");
  for(i=0;i<nobj;i++) {
   /* fprintf(fout," %8.6f  %8.3f %8.3f",zsample[i],msample[i],Msample[i]); */
    fprintf(fout," %8.6f  %8.6f %8.6f",zsample[i],zerror[i],zobserved[i]);
    fprintf(fout," %8.3f %8.3f",mDist[i],MDist[i]);
    fprintf(fout," %8.6f  %8.6f",mDistError[i],mDistObserved[i]);
    fprintf(fout," %8.6f",mSel[i]);
    fprintf(fout," %8.6f\n", colorsample[i]);
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
}

