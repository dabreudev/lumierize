#include "modulos.h"  
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

/* creamos una estrucura para añadir una m de seleccion y una m de cálculo
dabreu */
struct sample_data_sel_cal {
  int ngalax;
  double *z;
  double *lum_sel;
  double *lum_cal;
  double *mag_sel;
  double *mag_cal;
  double *ew;
  char   **image;
};

/* dabreu */
struct sample_data_mag_err {
  int ngalax;
  double *z;
  double *mag;
  double *mag_err;
};

struct lum_func {
  double mag_min;
  double mag_max;
  double lum_min;
  double lum_max;
  double zlow;
  double zup;
  double area;
  int islum;
  int npoint;
  double *loglf;
  double *lf; 
  double *errlf;
  double *errloglf;
  double *cumlf;
  int    *ngalbin;
  double *mag;
  double *lum;
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
/* dabreu */
void get_sample_sel_cal(struct sample_data_sel_cal *sample);
void get_sample_mag_err(struct sample_data_mag_err *sample);
void get_sample_ceg(struct sample_data *sample);
void set_lf_parameters(struct lum_func *lf);
void set_lf_ceg(struct lum_func_ceg *lf);
void set_cosmology();
void VVmax();
void STY();
/* dabreu */
void STY_MAG_ERR();
void SWML();
void CEG();
void Calc_Num();
void Generate_Cat_M();
/* dabreu */
void Generate_Cat_L();
void set_Schechter_M();
void set_Schechter_L();
void FitSchechter_L(struct lum_func lf, struct Schlf_L *schlf, double *chisq);
void FitSchechter_M(struct lum_func lf, struct Schlf_M *schlf, double *chisq);
void FitSchechter_Step_L(struct Steplf_L lf, struct Schlf_L *schlf, double *chisq);
void FitSchechter_Step_M(struct Steplf_M lf, struct Schlf_M *schlf, double *chisq);

double Schechter_M_LF(double Mag,double Mstar,double alfa,double phistar);
double Schechter_L_LF(double L,double Lstar,double alfa,double phistar);
double Amoe_Funk_d(int n, double *x, double *y, double *p);
double amofunc_schechterfitmag_d(int n, double *x, double *y, double *p);
void  mrqfunc_schechterfitmag_d(double x,double *p,double *y,double *dyda,int n);
double amofunc_schechterfitlum_d(int n, double *x, double *y, double *p);
void  mrqfunc_schechterfitlum_d(double x,double *p,double *y,double *dyda,int n);


int amoeba_use=0;
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
/*   double vol1[1000],vol2[1000],zz[1000]; */
/*    int i=0;  */
/*    int j;  */
/*    double zz,zzz; */
   
/*    double mm; */
/*    cosmo.H0=50;cosmo.q0=0.5; */
/*    for(i=0;i<100;i++) { */
/*      zz=i/10.; */
/*      mm=Mag_LF(zz,15.); */
/*      zzz=Z(mm-15.); */
/*      printf(" LOS DOS zz %f zzz %f\n\n",zz,zzz); */
/*    } */


  /*ESto que viene aqui lo he hecho para comprobar que las funciones ..ter_M y ..ter_L estan bien */
/*   double m,l; */
/*   double fact; */
/*   double Mstar,alfa,phistar,Lstar; */
/*   int i; */
/*   int dm=1.; */
/*   fact=1e26; */
/*   alfa=-0.9; */
/*   phistar=0.0033; */
/*   Mstar=-20; */
/*   Lstar=fact*pow(10.,-0.4*Mstar); */
/*   for(i=0;i<(int)(10/dm);i++) { */
/*     m=-12.-i*dm; */
/*     l=fact*pow(10.,-0.4*m); */
/*     printf(" m      %g    l   %g  dl %g\n",m,l,l*log(10.)*0.4*0.2); */
/*     printf(" S_M    %g    S_L %g\n",Schechter_M(m,Mstar,alfa,phistar),Schechter_L_LF(l,Lstar,alfa,phistar)); */
/*     printf(" S_M_DM %g S_L_DL %g\n",dm*Schechter_M_LF(m,Mstar,alfa,phistar),l*log(10.)*0.4*dm*Schechter_L_LF(l,Lstar,alfa,phistar)); */
/*   } */

/*   exit(1); */

  do{
    printf("\n C Compute number of galaxies with a given luminosity function \n"); 
    printf(" V Calculate luminosity function by V/Vmax method\n"); 
    
    printf(" L Calculate luminosity function by CEG method\n"); 
    printf(" M Calculate luminosity function by STY maximum likelihood method\n"); 
    /* dabreu */
    printf(" N Calculate luminosity function by STY method with gaussian errors in magnitude\n");
    printf(" W Calculate luminosity function by SWML maximum likelihood method\n"); 
/*    printf(" O Calculate luminosity function by C- method\n");  */
    printf(" G Generate a random catalogue for a given LF by Montecarlo simulations in magnitudes\n");
    printf(" H Generate a random catalogue for a given LF by Montecarlo simulations in luminosities\n");
    printf(" E Exit\n");
    opt=readc(opt);  
    switch (opt) { 
    case 'C':
    case 'c':
      set_cosmology();
      Calc_Num();
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
      /* dabreu */
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
    }



  }while(opt!='E' && opt!='e');

  return(0);
}

void STY() 
{
  static struct lum_func sty;
  struct sample_data sample;
  static double mlim=0;
  static double flim=0;
  double zlow;

  int iter;
  
  struct Schlf_M lfsch_M;
  struct Schlf_L lfsch_L;

  static int poissonflag=0;

  /* dabreu */
  static int plots=0; /* para no hacer las gráficas */
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
  
  static struct lum_func sty;
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
  
  static struct lum_func swml;
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
  lfstep_M.covarlf=matrix_d(swml.npoint,swml.npoint);
  lfstep_L.covarlf=matrix_d(swml.npoint,swml.npoint);
  
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
    lfstep_L.nbin=swml.npoint;
    lfstep_L.lumi    =vector_d(swml.npoint+1);
    lfstep_L.lf      =vector_d(swml.npoint);
    lfstep_L.errlumi =vector_d(swml.npoint+1);
    lfstep_L.errlf   =vector_d(swml.npoint);

    for(j=0;j<swml.npoint;j++) {
      lfstep_L.lumi[j]=(swml.lum_min+(swml.lum_max-swml.lum_min)*j/swml.npoint)*log(10);
    }
    lfstep_L.lumi[swml.npoint]=swml.lum_max*log(10);


    if(poissonflag) iter=MLA_SWML_p_L(sample.ngalax,sample.lum,sample.z,flim,swml.area,zlow,swml.zup,cosmo,&lfstep_L);  
    else            iter=MLA_SWML_L(sample.ngalax,sample.lum,sample.z,flim,swml.area,zlow,swml.zup,cosmo,&lfstep_L);  
    printf("\n Solution found at iteration %d\n",iter);
    for(j=0;j<swml.npoint;j++) {
      printf(" Lum %11g - %11g    LF %11g (log=%7g) Err_LF %11g (log=%7g)\n",lfstep_L.lumi[j]/log(10),lfstep_L.lumi[j+1]/log(10),exp(lfstep_L.lf[j]),log10(exp(lfstep_L.lf[j])),exp(lfstep_L.lf[j])*lfstep_L.errlf[j],lfstep_L.errlf[j]/log(10.));
    }
    PlotStepLF_L(lfstep_L);

  }
  else {
    printf(" Input the limiting magnitude: ");
    mlim=readf(mlim);
    cpgopen("?");
    cpgswin((float)swml.mag_min,(float)swml.mag_max,-5.,1.);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
 
    lfstep_M.nbin=swml.npoint;
    lfstep_M.magni   =vector_d(swml.npoint+1);
    lfstep_M.lf      =vector_d(swml.npoint);
    lfstep_M.errmagni=vector_d(swml.npoint+1);
    lfstep_M.errlf   =vector_d(swml.npoint);
    for(j=0;j<swml.npoint;j++) {
      lfstep_M.magni[j]=swml.mag_min+(swml.mag_max-swml.mag_min)*j/swml.npoint;
    }
    lfstep_M.magni[swml.npoint]=swml.mag_max;
    if(poissonflag)  iter=MLA_SWML_p_M(sample.ngalax,sample.mag,sample.z,mlim,swml.area,zlow,swml.zup,cosmo,&lfstep_M);
    else             iter=MLA_SWML_M(sample.ngalax,sample.mag,sample.z,mlim,swml.area,zlow,swml.zup,cosmo,&lfstep_M);
    printf("\n Solution found at iteration %d\n",iter);
    for(j=0;j<swml.npoint;j++) printf(" Mag %11g - %11g    LF %11g (log=%7g) Err_LF %11g (log=%7g)\n",lfstep_M.magni[j],lfstep_M.magni[j+1],exp(lfstep_M.lf[j]),log10(exp(lfstep_M.lf[j])),exp(lfstep_M.lf[j])*lfstep_M.errlf[j],lfstep_M.errlf[j]/log(10.));
    printf(" Plotting function\n");
    PlotStepLF_M(lfstep_M);
    printf(" After plotting\n");   
   
  }

  cpgopen("?");  
  cpgask(0);


  if(swml.islum) {
    FitSchechter_Step_L(lfstep_L, &lfschfit_L, &chisq); 
    PlotStepSchLF_L(lfstep_L,lfschfit_L);   
    free(lfstep_L.lumi);
    free(lfstep_L.lf);
    free(lfstep_L.errlumi);
    free(lfstep_L.errlf);
  }
  else           {
    FitSchechter_Step_M(lfstep_M, &lfschfit_M, &chisq); 
    PlotStepSchLF_M(lfstep_M,lfschfit_M);
    free(lfstep_M.magni);
    free(lfstep_M.lf);
    free(lfstep_M.errmagni); 
    free(lfstep_M.errlf);
  }
  cpgclos();

  free_matrix_d(lfstep_M.covarlf,swml.npoint,swml.npoint);
  free_matrix_d(lfstep_L.covarlf,swml.npoint,swml.npoint);
}


void VVmax()
{
  int i,j;
  static struct lum_func vvmax={-40,40,0,0,0,10,41252,1,10};
/*   struct lum_func vvmax_dif; */
  struct sample_data_sel_cal sample;
/*   struct Schlf_M  lfschfit_M; */
/*   struct Schlf_L  lfschfit_L; */

  static double mlim=1.;
  static double llim=1.;
  double M=0;
  double L=0;
/*   int vol; */
  double zmax;
  double *x,*y,*sigy,*yfit;
  double dm=0;
  double dl=0;
/*   //Variables para el test V/Vmax */
  double *vvmaxtest;
  int   *ngaltest;           /* //numero de galaxias en cada bin del test */
  double mmin_t=0,mmax_t=0;       /* //magnitudes minimas y maximas desde las cuales se hara el test V/Vmax. */
  double lmin_t=0,lmax_t=0;       /* //luminosidades minimas y maximas desde las cuales se hara el test V/Vmax. */
  int   n_t;                 /* //numero de puntos para el test */
  double *mhist=NULL;              /* //array para las magnitudes del histograma */
  double *loglhist=NULL;           /* //array para las luminosidades del histograma */
  double *mlim_t=NULL;
  double *llim_t=NULL;
  double *logllim_t=NULL;
  char  cnul;
  double fnul;
  double zlow;
  double zup;
/*   int iter; */
  double par[3],sigpar[3];
  double **covarpar;
  int   ipar[3];
  double chisq;
  double dum1,dum2,dum3;
  int nfit;

  int pg1,pg2;
  int iter;
  double pi=3.1415926535897932384;

  double ymin,ymax;

  /* dabreu */
  static int plots=0; /* para no hacer gráficas */
  /* para escribir los resultados en un fichero */
  FILE *fout = NULL;
  static char resultFileName[200]="";

  covarpar=matrix_d(3,3);


  ipar[0]=1;ipar[1]=1;ipar[2]=1;

  get_sample_sel_cal(&sample);
  set_lf_parameters(&vvmax);
/*   //set_cosmology(); */
/*   //printf(" */
  for(i=0;i<sample.ngalax;i++) {
/*     //printf("i %d  z %f mag %f\n",i,sample.z[i],sample.mag[i]); */
  }
  
/*   printf(" Input number of bins in magnitude: "); */
/*   vvmax.npoint=(int)readf(10.); */
  
  vvmax.lf      =malloc(vvmax.npoint*sizeof(double));
  vvmax.loglf   =malloc(vvmax.npoint*sizeof(double));
  vvmax.errlf   =malloc(vvmax.npoint*sizeof(double));
  vvmax.errloglf=malloc(vvmax.npoint*sizeof(double));
  vvmax.cumlf   =malloc(vvmax.npoint*sizeof(double));
  if(vvmax.islum) vvmax.lum    =malloc(vvmax.npoint*sizeof(double));
  else            vvmax.mag    =malloc(vvmax.npoint*sizeof(double));
  vvmax.ngalbin=malloc(vvmax.npoint*sizeof(int));
  printf(" Performing the V/Vmax test\n");
  if(vvmax.islum) {
    MinMax_d(sample.ngalax,sample.lum,&lmin_t,&lmax_t);
    printf(" Fluxes range: %g - %g\n",lmin_t,lmax_t);
    lmin_t*=0.3;lmax_t*=1.2;
  }
  else {
    MinMax_d(sample.ngalax,sample.mag,&mmin_t,&mmax_t);
    printf(" Apparent magnitude range: %g - %g\n",mmin_t,mmax_t);
    mmin_t-=1.;mmax_t+=4.;
  }
  n_t=100;
  vvmaxtest   =malloc(n_t*sizeof(double));
  ngaltest    =malloc(n_t*sizeof(int));
  if(vvmax.islum) {
    loglhist      =malloc(sample.ngalax*sizeof(double));
    logllim_t     =malloc(n_t*sizeof(double));
    llim_t        =malloc(n_t*sizeof(double));
  }
  else {
    mhist       =malloc(sample.ngalax*sizeof(double));
    mlim_t      =malloc(n_t*sizeof(double));
  }
/*   printf(" Input lower redshift: "); */
  zlow=vvmax.zlow;
  zup=vvmax.zup;
  zlow = (zlow < ZMIN ? ZMIN : zlow);
  if(DEBUG) printf(" zlow %f zup %f\n",zlow,zup);

/*   printf(" Input limiting redshift (0=none, magnitude limited sample): "); */
/*   zup=readf(0.); */
  for(i=0;i<n_t;i++) {
    if(vvmax.islum) {
      llim_t[i]=pow(10.,i*(log10(lmax_t)-log10(lmin_t))/(n_t-1.)+log10(lmin_t));
      logllim_t[i]=i*(log10(lmax_t)-log10(lmin_t))/(n_t-1.)+log10(lmin_t);
      if(DEBUG2) printf(" llim %g\n",llim_t[i]);
    }
    else           {
      mlim_t[i]=i*(mmax_t-mmin_t)/(n_t-1.)+mmin_t;
      if(DEBUG) printf(" mlim %f\n",mlim_t[i]);
    }
    vvmaxtest[i]=0;
    ngaltest[i]=0;
    for(j=0;j<sample.ngalax;j++) {
      if(vvmax.islum)  {
	loglhist[j]=log10(sample.lum[j]);
	L=Lum(sample.z[j],sample.lum[j],cosmo);
	if(DEBUG2)        printf(" Galax %d z %f flux %g Lum %g mag %f Mag %f\n",j,sample.z[j],sample.lum[j],L,sample.mag[j],Mag(sample.z[j],sample.mag[j],cosmo));
	if(L>1e36) 
 	printf(" Galax %d z %f flux %g Lum %g mag %f Mag %f\n",j,sample.z[j],sample.lum[j],L,sample.mag[j],Mag(sample.z[j],sample.mag[j],cosmo)); 

      }
      else {
	mhist[j]=sample.mag[j];
	M=Mag(sample.z[j],sample.mag[j],cosmo);
	if(DEBUG2)        printf(" Galax %d z %f mag %f Mag %f\n",j,sample.z[j],sample.mag[j],M);
      }
      /*       if(sample.mag[j]<mlim_t[i]) {   Esta es la antigua, cuando no podia hacer en luminosidades */
      /* Cuidado, sample.lum[j] (o mag[j]) puede no estar alocateado, 
	 por eso pongo antes vvmax.islum, que comprueba antes y si no lo cumple 
	 no pasa a la siguiente */
/*       if(( !vvmax.islum && sample.mag[j]<mlim_t[i])  || (vvmax.islum && sample.lum[j]<llim_t[i]) ) { */
      if(vvmax.islum) {
	if(sample.lum[j]>llim_t[i]) {
	  if(DEBUG2) printf(" Entro la flux=%f con llim=%f\n",sample.lum[j],llim_t[i]);
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
      else {
	if(sample.mag[j]<mlim_t[i]) {
	  if(DEBUG2) printf(" Entro la m=%f con mlim=%f\n",sample.mag[j],mlim_t[i]);
	  zmax=Z_m(mlim_t[i],M,cosmo);
	  if(zmax==0) zmax=ZMIN;
	  if(zmax>zup && zup!=0) zmax=zup;
	  if(zmax<zlow && zlow!=0) {printf(" zlow %f zmax %f\n ERROR!!\n",zlow,zmax);exit(1);}
	  if(DEBUG2) printf(" zmax %f %f %f %f\n",zmax, M, mlim,sample.z[j]); 
	  /* 		vvmaxtest[i]+=Vol_LF(sample.z[j])/Vol_LF(zmax);  VVmax tradicional */
	  vvmaxtest[i]+=(Vol(sample.z[j],cosmo)-Vol(zlow,cosmo))/(Vol(zmax,cosmo)-Vol(zlow, cosmo));  /* Cambiado para un zlow!=0 */
	  ngaltest[i]++;
	}
	if(DEBUG2) printf(" No entro la m=%f con mlim=%f\n",sample.mag[j],mlim_t[i]);
      }
    }
    vvmaxtest[i]/=ngaltest[i];
    if(vvmax.islum)    printf(" Flux %g V/Vmax %f Ngal %d\n",llim_t[i],vvmaxtest[i],ngaltest[i]);
    else               printf(" Mag %g V/Vmax %f Ngal %d\n",mlim_t[i],vvmaxtest[i],ngaltest[i]);
  }
  /* dabreu */
  printf("Despues de 'Mag g V/Vmax f Ngal'\n");
  printf(" Do you want plots(1=yes, 0=no)?\n");
  plots=readf(plots);
  if(plots) 
  {
    pg1=cpgopen("/xserve");
    if(vvmax.islum) {
        cpgswin(log10(lmin_t),log10(lmax_t),0.,sample.ngalax/2.);
        cpghist_d(sample.ngalax,loglhist,log10(lmin_t),log10(lmax_t),15,1);
        cpglab("Apparent fluxes (W/m2)","<V/Vmax>","V/Vmax test");
        cpgbox("BCTNSL",0,0,"BTNS",0,0);
    }
    else {
        cpgswin(mmin_t,mmax_t,0.,sample.ngalax/2.);
      /*     cpghist_d(sample.ngalax,mhist,mmin_t,mmax_t,15,1); */
        cpglab("Magnitude","<V/Vmax>","V/Vmax test");
        cpgbox("BCTNS",0,0,"BTNS",0,0);
    }
    if(vvmax.islum) {
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
    else {
        cpgswin(mmin_t,mmax_t,0.,1.);
        cpgbox("BCTNS",0,0,"CTMS",0,0);
        cpgpt_d(n_t,mlim_t,vvmaxtest,1); 
      printf(" The V/Vmax method requires to know the limiting magnitude.\n Click with the left button to select the limiting magnitude with the mouse or right button to select it by keyboard\n");
        cpgband_d(6,1,0.,0.,&mlim,&fnul,&cnul); 
    }
    /*   //printf(" OPION %c\n",cnul); */
    if(cnul=='X') {
      if(vvmax.islum) {
        printf(" Input the limiting flux: ");
        llim=readd(llim);
      }
      else {
        printf(" Input the limiting magnitude: ");
        mlim=readd(mlim);
        }
    }
  /* dabreu */
  }   /* final del if(plots) */ 
  else 
  {
      if(vvmax.islum) {
        printf(" Input the limiting flux: ");
        llim=readd(llim);
      }
      else {
        printf(" Input the limiting magnitude: ");
        mlim=readd(mlim);
        }
  }
  
  if(vvmax.islum){
    dl=(vvmax.lum_max-vvmax.lum_min)/(vvmax.npoint-1.);
    printf(" Limiting luminosity: %g. Lum bin %g From %g to %g\n",llim,dl,vvmax.lum_max,vvmax.lum_min);
  }
  else {
    printf(" Limiting magnitude: %f\n",mlim);
    dm=(vvmax.mag_max-vvmax.mag_min)/(vvmax.npoint-1.);
  }

  if(DEBUG) printf(" Intervalo %g\n",dm);

/*   exit(1);  */
  for(i=0;i<vvmax.npoint;i++) 
  {
    /*     //Bucle principal con los bins en magnitudes */
    if(vvmax.islum)    vvmax.lum[i]  =pow(10.,vvmax.lum_min+i*dl);
    else               vvmax.mag[i]  =vvmax.mag_min+i*dm;
    /*     //printf(" NI ESTO \n"); */
    vvmax.lf[i]=0.;
    vvmax.ngalbin[i]=0;
/*     //printf(" AQUI\n"); */
    for(j=0;j<sample.ngalax;j++) 
    {
      /*       //Aqui lo mejor seria tener ordenado en magnitudes la muestra */
      /*       //Pero como, no, pues lo hago con un for y un if */

      /*       if(( !vvmax.islum && fabs(M-vvmax.mag[i])<dm/2.) || (vvmax.islum && fabs(M-vvmax.lum[i])<dl/2.)) { */
      if(vvmax.islum) { 
	L=Lum(sample.z[j],sample.lum[j],cosmo);
	if(DEBUG) printf(" L %g flux %g z %g L1 %g L2 %g llim %g\n",log10(L),sample.lum[j],sample.z[j],log10(vvmax.lum[i])-dl/2.,log10(vvmax.lum[i])+dl/2.,llim);
	if(fabs(log10(L)-log10(vvmax.lum[i]))<dl/2.) { 
	  zmax=Z_l(llim,L,cosmo);      /* Esto no es asi, tengo que comprobarlo */
	  if(zmax>zup && zup!=0) zmax=zup;
	  if(zmax>zlow) {
	    vvmax.lf[i]+=1./(Vol(zmax,cosmo)-Vol(zlow,cosmo));
	    vvmax.ngalbin[i]++;
	  }
	  if(DEBUG) printf(" La he metido LF %g 1/Vol  %g zmax %f zlow %f \n",vvmax.lf[i],1./(Vol(zmax,cosmo)-Vol(zlow,cosmo)),zmax,zlow);
	}
      }
      else {
	M=Mag(sample.z[j],sample.mag[j],cosmo);
	if(fabs(M-vvmax.mag[i])<dm/2.) {
	  zmax=Z_m(mlim,M,cosmo);       /* Esto no es asi, tengo que comprobarlo */
	  if(zmax>zup && zup!=0) zmax=zup;
	  if(zmax>zlow) {
	    vvmax.lf[i]+=1./(Vol(zmax,cosmo)-Vol(zlow,cosmo));
	    vvmax.ngalbin[i]++;
	  }
	}
      }
    }
    if(vvmax.islum) {
      if(DEBUG) printf(" lf %g fac %g dl %f\n",vvmax.lf[i],1/dl/vvmax.lum[i]/log(10),dl);
      vvmax.lf[i]=vvmax.lf[i]/vvmax.area*4*pi/dl/vvmax.lum[i]/log(10)*1.e18; /* El 1.e18 es para pasar de parsecs3 a megaparsecs3. El dl/vvmax.lum[i] es diferencial de L ( ya que dl es en realidad d(logL)) */
      vvmax.errlf[i]=vvmax.lf[i]/sqrt(vvmax.ngalbin[i]); /* Esto es como suponer que la funcion de lum es proporcional a una variable poissonian que es el numero de de galaxias en ese bin*/
    }
    else {
      vvmax.lf[i]=vvmax.lf[i]/vvmax.area*4*pi/dm*1.e18; /* El 1.e18 es para pasar de parsecs3 a megaparsecs3 */
      vvmax.errlf[i]=vvmax.lf[i]/sqrt(vvmax.ngalbin[i]); /* Esto es como suponer que la funcion de lum es proporcional a una variable poissonian que es el numero de de galaxias en ese bin*/
    }
    if(vvmax.lf[i]!=0) {
      vvmax.loglf[i]=log10(vvmax.lf[i]); 
      vvmax.errloglf[i]=vvmax.errlf[i]/vvmax.lf[i]/log(10.); 
    }
    else {
      vvmax.loglf[i]=0.; 
      vvmax.errloglf[i]=0.; 
    }

    if(DEBUG) printf(" Voy a escribir\n");
    /*     El dm es para dividir por el intervalo en magnitudes. */
    /*     if(i>0) vvmax.lf[i]=(vvmax.cumlf[i]-vvmax.cumlf[i-1])/(vvmax.mag[i]-vvmax.mag[i-1]); */
    /*     else vvmax.lf[i]=0; */
    /*     if(vvmax.islum) printf(" Lum %15g Func Lum %15g Diff LF %15g Err_LF %15f Ngal %d\n",vvmax.lum[i],vvmax.cumlf[i],vvmax.loglf[i],vvmax.errloglf[i],vvmax.ngalbin[i]); */
/*     if(vvmax.islum) printf(" Lum %15g lum-dl/2 %15g (%f) lum+dl/2 %15g (%f) LF %15f Err_LF %15f Ngal %d\n",vvmax.lum[i],pow(10.,log10(vvmax.lum[i])-dl/2.),log10(vvmax.lum[i])-dl/2.,pow(10.,log10(vvmax.lum[i])+dl/2.),log10(vvmax.lum[i])+dl/2.,vvmax.loglf[i],vvmax.errloglf[i],vvmax.ngalbin[i]); */
/*     else            printf(" Mag %15f mag-dl/2 %15g      mag+dl/2 %15g      LF %15g Err_LF %15f Ngal %d\n",vvmax.mag[i],vvmax.mag[i]-dm/2.,vvmax.mag[i]+dm/2.,vvmax.loglf[i],vvmax.errloglf[i],vvmax.ngalbin[i]); */
    if(vvmax.islum) printf(" Lum %11g (log=%7g)  LF %11g (log=%9g) Err_LF %11g (log=%9g) LF(log L)=%11g (log=%7g) Err_LF(logL)=%11g (log=%7g)  Ngal %d\n",vvmax.lum[i],log10(vvmax.lum[i]),vvmax.lf[i],vvmax.loglf[i],vvmax.errlf[i],vvmax.errloglf[i],vvmax.lf[i]*vvmax.lum[i],log10(vvmax.lf[i]*vvmax.lum[i]),vvmax.errlf[i]*vvmax.lum[i],vvmax.errloglf[i],vvmax.ngalbin[i]);
    else            printf(" Mag %11g    LF %11g (log=%7g) Err_LF %11g (log=%7g) Ngal %d\n",vvmax.mag[i],vvmax.lf[i],vvmax.loglf[i],vvmax.errlf[i],vvmax.errloglf[i],vvmax.ngalbin[i]);
    /*     else            printf(" Mag %15f Func Lum %15g Diff LF %15g Err_LF %15f Ngal %d\n",vvmax.mag[i],vvmax.cumlf[i],vvmax.loglf[i],vvmax.errloglf[i],vvmax.ngalbin[i]); */
  }

  
  x   =malloc(vvmax.npoint*sizeof(double));
  y   =malloc(vvmax.npoint*sizeof(double));
  sigy=malloc(vvmax.npoint*sizeof(double));
  yfit=malloc(vvmax.npoint*sizeof(double));
  par[2]=0.;
  nfit=0;
  for(i=0;i<vvmax.npoint;i++) {
    if(vvmax.lf[i]!=0) {
      if(vvmax.islum)     x[nfit]=log10((double)vvmax.lum[i]);
      else                x[nfit]=(double)vvmax.mag[i];
      y[nfit]   =(double)vvmax.loglf[i];
      sigy[nfit]=(double)vvmax.errloglf[i];
      par[2]+=vvmax.lf[i];
      nfit++;
    }
  }
  if(DEBUG) printf(" Pasa\n");
  /*   MinMax(vvmax.npoint,y,&fnul,par+2); */
  if(vvmax.islum)  par[2]=(vvmax.lum[vvmax.npoint]-vvmax.lum[0])*par[2]/vvmax.npoint;
  else             par[2]=(vvmax.mag[vvmax.npoint]-vvmax.mag[0])*par[2]/vvmax.npoint;
  sigpar[2]=par[2]/10.;
  par[0]=-20.0  ;sigpar[0]=.1;   /* Mstar */ 
  par[1]=-0.6    ;sigpar[1]=0.05;   /* alfa */ 
  par[2]=0.0033;sigpar[2]=par[2]/10;   /* Phistar */ 
  if(vvmax.islum) {
    par[0]=36;
    sigpar[0]=.5;
  }
  else {
    par[0]=-23.;
    sigpar[0]=.5;   /* Mstar */
  }
  if(DEBUG) printf(" A este \n");
  if(vvmax.islum) {
    if(DEBUG) printf(" Llamo amoeba\n");
    iter=Amoeba_d(nfit,x,y,3,par,sigpar,FTOL,10000,amofunc_schechterfitlum_d); 
/*     printf("Amoeba Schechter_M parameters fit: \n");  */
/*     printf(" log(Lstar (W)):  %g   alpha:  %g  Phistar:  %g\n",par[0],par[1],par[2]);  */
    iter=Mrq_d(x,y,sigy,nfit,par,ipar,3,covarpar,&chisq,mrqfunc_schechterfitlum_d);
    printf("Mrq Fit Schechter_L parameters fit: \n");
    printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",par[0],par[1],par[2],log10(par[2]));
    printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[2][2])/par[2]/log(10)); 
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD logLstar alpha Phistar log\n");
    printf("#FL_HEAD E_logLstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",par[0],par[1],par[2],log10(par[2]));
    printf("#FL_ERR %g %g %g %g\n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[2][2])/par[2]/log(10));

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
    fprintf(fout, "%g\t", par[0]);
    fprintf(fout, "%g\t", par[1]);
    fprintf(fout, "%g\t", par[2]);
    fprintf(fout, "%g\t", log10(par[2]));
    fprintf(fout, "%g\t", sqrt(covarpar[0][0]));
    fprintf(fout, "%g\t", sqrt(covarpar[1][1]));
    fprintf(fout, "%g\t", sqrt(covarpar[2][2]));
    fprintf(fout, "%g\n", sqrt(covarpar[2][2])/par[2]/log(10));
    fclose(fout);  
  }
  
  else   {
/*  iter=Amoeba_d(nfit,x,y,3,par,sigpar,FTOL,10000,amofunc_schechterfitmag_d);*/
/*     printf("Amoeba Schechter_M parameters fit: \n"); */
/*     printf(" Mstar:  %g   alpha:  %g  Phistar:  %g\n",par[0],par[1],par[2]); */
/*     for(i=0;i<nfit;i++) { */
/*       yfit[i]=log10(Schechter_M_LF(x[i],par[0],par[1],par[2])); */
/*     } */
/*     cpgsci(3); */
/*     cpgline_d(nfit,x,yfit); */

    /* Solo hago el ajuste con Mrq, sin usar Amoeba */
    
    iter=Mrq_d(x,y,sigy,nfit,par,ipar,3,covarpar,&chisq,mrqfunc_schechterfitmag_d);
    printf(" Mrq Fit Schechter_M parameters fit:\n");
    printf(" Mstar:    %g   alpha:    %g  Phistar  :  %g (log=%g) \n",par[0],par[1],par[2],log10(par[2]));
    printf(" E_Mstar:  %g   E_alpha:  %g  E_Phistar:  %g (log=%g)  \n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[2][2])/par[2]/log(10));
    /* dabreu */
    printf("##################################################\n");
    printf("#FL_HEAD Mstar alpha Phistar log\n");
    printf("#FL_HEAD E_Mstar E_alpha E_Phistar log\n");
    printf("#FL_DATA %g %g %g %g\n",par[0],par[1],par[2],log10(par[2]));
    printf("#FL_ERR %g %g %g %g\n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[2][2])/par[2]/log(10));

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
    fprintf(fout, "%g\t", par[0]);
    fprintf(fout, "%g\t", par[1]);
    fprintf(fout, "%g\t", par[2]);
    fprintf(fout, "%g\t", log10(par[2]));
    fprintf(fout, "%g\t", sqrt(covarpar[0][0]));
    fprintf(fout, "%g\t", sqrt(covarpar[1][1]));
    fprintf(fout, "%g\t", sqrt(covarpar[2][2]));
    fprintf(fout, "%g\n", sqrt(covarpar[2][2])/par[2]/log(10));
    fclose(fout);
  }
  
/*   par[0]=-33.5  ;sigpar[0]=.1;   //Mstar */
/*   par[1]=-.9    ;sigpar[1]=0.01;   //alfa */
/*   par[2]=0.004;sigpar[2]=par[2]/10;   //Phistar */

  for(i=0;i<nfit;i++) {
    if(vvmax.islum)    yfit[i]=Schechter_L_LF(pow(10.,x[i]),pow(10.,par[0]),par[1],par[2]);
    else               yfit[i]=log10(Schechter_M_LF(x[i],par[0],par[1],par[2]));
    if(DEBUG) printf(" x %g y %g\n",x[i],yfit[i]);
/*     printf("Valores ajuste x %g y %g\n",x[i],yfit[i]); */
  } 
 
/*   printf(" Schecher de -1. 0.004  -30  -32 %g\n",Schechter_M_LF(-32.,-30.,-1.,0.004)); */
/*   printf(" log Schecher de -1. 0.004  -30  -32 %g\n",log10(Schechter_M_LF(-32.,-30.,-1.,0.004))); */
  
  /* dabreu */
  if(plots) 
  {
    pg2=cpgopen("?");

    ymin=1e38;
    ymax=-1e38;
    for(i=0;i<vvmax.npoint;i++) {
      if(vvmax.lf[i]!=0) {
        if((double)vvmax.loglf[i]>ymax) ymax=(double)vvmax.loglf[i];
        if((double)vvmax.loglf[i]<ymin) ymin=(double)vvmax.loglf[i];
      }
    }
    ymin=ymin-2.5;
    ymax=ymax+2.5;

    printf(" Cuts %f %f \n",ymin,ymax);

    if(vvmax.islum) {
      cpgswin(vvmax.lum_max+0.1,vvmax.lum_min-0.1,ymin,ymax);
      cpgbox("BCTNSL",0,0,"BCTNS",0,0);
    }
    else {
        cpgswin(vvmax.mag_max+1,vvmax.mag_min-1,ymin,ymax);
        cpgbox("BCTNS",0,0,"BCTNS",0,0);
    }
    cpgsci(1);
    cpgsch(1.2); 
    for(i=0;i<nfit;i++) {
      cpgpt1(x[i],y[i],17);
      dum1=x[i];
      dum2=y[i]-sigy[i];
      dum3=y[i]+sigy[i];
      /*     printf(" dum2 %f dum3 %f \n",dum2,dum3); */
        cpgerry_d(1,&dum1,&dum2,&dum3,0.35); 
    }
    cpgsci(4);
    cpgline_d(nfit,x,yfit);
    cpgsci(1);
    cpgsch(1.);
  
    if(vvmax.islum) cpglab("Luminosity (W)","log(Phi\\dlum\\u) (#/Mpc3/W) ","Luminosity function by the V/Vmax method computed with fluxes");
    else            cpglab("Mag","log(Phi\\dmag\\u) (#/Mpc3/Mag) ","Luminosity function by the V/Vmax method computed with magnitudes");
    
  
    cpgclos();
    cpgslct(pg1); 
    cpgclos();
  }  /* final del if del 2º plot */
  free_matrix_d(covarpar,3,3);
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
  static int plots=0; /* para no hacer las gráficas */
  
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

double Schechter_M_LF(double M,double Mstar,double alfa,double phistar)
{
  double tmp;
  double elev;
  
  elev=pow(10.,0.4*(Mstar-M));
  tmp=0.92103404*phistar*pow(elev,alfa+1.);   /* // 0.92103404= 0.4*ln(10) */
/*   //  printf("Mag %f $$ %e $$ %e\n",M,tmp,exp(-elev*(schlf.alfa+1.))); */
  tmp=tmp*exp(-elev);
  return(tmp);

}

/* // Esta no vale (Willmer 97: Estimating Galaxy..) */
/* //double Schechter_M_LF(double M) */
/* //{ */
/* //  double tmp; */
/* //  double elev; */
/* //   */
/* //  elev=pow(10.,0.4*(Mstar-M)*(1.+alfa)); */
/* //  tmp=0.92103404*phistar*elev;   // 0.92103404= 0.4*ln(10) */
/* //  tmp=tmp*exp(-elev); */
/* //  return(tmp); */
/* // */
/* //} */

double Schechter_L_LF(double L,double Lstar,double alfa,double phistar)
{
  double tmp;
/*   double elev; */
  double L_Lstar;

  if(DEBUG2) printf(" L %g Lstar %g\n",L,Lstar);

  L_Lstar=L/Lstar;
  if(DEBUG2) printf(" L_Lstar %g phis %g alfa %f\n",L_Lstar,phistar,alfa);
/*   printf(" L_Lstar %g phis %g alfa %f\n",L_Lstar,phistar,alfa); */
  tmp=phistar/Lstar*pow(L_Lstar,alfa)*exp(-L_Lstar);
/*   printf(" %g %g %g %g\n",log10(phistar),-log10(Lstar),alfa*log10(L_Lstar),-L_Lstar*log10(exp(1))); */
/*   printf(" Todo %g\n",log10(phistar)-log10(Lstar)+alfa*log10(L_Lstar)-L_Lstar*log10(exp(1))); */
  if(DEBUG2) printf(" tmp %g\n",tmp);
  return(log10(tmp));
}



void set_lf_parameters(struct lum_func *lf) {
  
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
/*   //zlow=readf(0.); */
  lf->zlow=(double)readf(lf->zlow); 
  printf(" Input limiting redshift (0=none, magnitude limited sample): ");
  lf->zup=(double)readf(lf->zup); 
/*   //zup=readf(0.); */



/*   printf(" Input minimum redshift  : "); */

/*   printf(" Input maximum redshift  : ");  */


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

void get_sample_sel_cal(struct sample_data_sel_cal *sample) {
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
    collum=readi(collum_sel);
    /*     }else{ */
    printf(" Input column with selection apparent magnitude data: ");
    colmag=readi(colmag_sel);
    printf(" Input column with apparent calculation luminosity (fluxes in W/m2)  data: ");
    collum=readi(collum_cal);
    /*     }else{ */
    printf(" Input column with apparent calculation magnitude data: ");
    colmag=readi(colmag_cal);
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
    ReadDoublecol(datafile,colmag_sel,mag,log2,&ndat);
    ReadDoublecol(datafile,collum_sel,lum,log3,&ndat);
    ReadDoublecol(datafile,colmag_cal,mag,log4,&ndat);
    ReadDoublecol(datafile,collum_cal,lum,log5,&ndat);
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

void get_sample_mag_err(struct sample_data_mag_err *sample)
{
/* dabreu */
/* para leer catálogo con magnitud y error en la magnitud */
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

void set_Schechter_M() {
  printf(" Schechter parameters function: \n");
  printf(" M star: ");
  schlf_M.Mstar=readf(schlf_M.Mstar);
  printf(" Phi star: ");
  schlf_M.phistar=readf(schlf_M.phistar);
  printf(" alfa: ");
  schlf_M.alfa=readf(schlf_M.alfa);

}

void set_Schechter_L() {

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
  static int plots=0; /* para no hacer las gráficas */

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
          "# 4 M_APP (apparent magnitude without observatonal errors)\n"
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
/* toda la función Generate_Cat_L() es una adaptación de Generate_Cat() */

void Generate_Cat_L()
{
  static double fluxlow=0.,fluxup=0.,fluxerr_mean=0.,fluxerr_stddev=0.;
  double Flow,Fup;
  static double zlow=0,zup=0.5;
  static double area=41252.;
  double *Lsample,*zsample,*fsample,*ferror,*fobserved; 
  static int plots=0; /* para no hacer las gráficas */

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

  /* no funcionan los errores, creo que por error en lectura o algo así */
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
  fprintf(fout,"#\n# 1 Z  (redshift)\n# 2 F_APP (apparent flux without observatonal errors)\n# "
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

double Amoe_Funk_d(int n, double *x, double *y, double *p)
{

  int i,j;
  double f,s;
  double prob;              /*  //La probabilidad para cada gal. */
  double Mlow,Mup;
  double elev,tmp;
  double M;
  double sch_int;
  double logL=0;

/*   //printf(" AMOEBA USE %d\n",amoeba_use); */

  if (amoeba_use == 1) {             /* //  Schechter fitting */
    s=0.0;
    for(i=0; i<n; i++) {
      f=Schechter_M_LF(x[i],p[0],p[1],p[2]);
/*       //    f += p[3]+p[4]*x[i]; */
/*       //f += p[4]*x[i]; */
/*       //printf(" f %f y %f\n",f,y[i]); */
      if (f!=0 && y[i]!=0)  s += (log10(f)-y[i])*(log10(f)-y[i]);
    };
    printf(" s = %f  p0 %g p1 %g P2 %g\n",s,p[0],p[1],p[2]); 
    
    return(s);
  }
  
  else if(amoeba_use == 2) {          /* //maximum likelyhood method */
    logL=0;
    for(i=0; i<n; i++) {
      M=Mag(y[i],x[i],cosmo);
      f=Schechter_M_LF(M,p[0],p[1],p[2]); /* //Este es el factor de arriba  */
/*       //Aqui integro la funcion de Schechter por enesima vez */
      sch_int=0;
      Mlow=Mag(y[i],p[3],cosmo);   /* //Truco para pasarle estos parametros */
      Mup=Mag(y[i],p[4],cosmo);   /*  //Ver justo antes de la llamada a Amoeba */
/*       //x[i] es M[i] y y[i] es z[i] */
/*       //printf(" Galaxia %d mag %f z %f\n",i,x[i],y[i]); */


      if(p[0]-Mup> 5 && Mlow-p[0]>0 ) Mup=p[0]-5.;  /* //Para que irse tan lejos!, si es cero */

      for(j=0;j<NSTEP_LF;j++) {
	M=Mlow+j*(Mup-Mlow)/(NSTEP_LF-1.);      
	/*       printf(" M %f Mstar %f \n",M,Mstar); */
	elev=pow(10.,0.4*(p[0]-M));                  /* //p[0] es Mstar!! */
	tmp=0.92103404*p[2]*pow(elev,p[1]+1.);       /* //p[1] es alfa!! */
	tmp=tmp*exp(-elev);                          /* //p[2] es Phistar!! */
	sch_int+=tmp;
	/*       printf(" sch %f\n",sch_int); */
      }
      sch_int*=(Mlow-Mup)/NSTEP_LF;
/*       if(i==100) { */
/* 	printf("  f %g sc %g p %g s %g Mlow %g Mup %g p0 %g p1 %g P2 %g p3 %f\n ",f,sch_int,prob,logL,Mlow,Mup,p[0],p[1],p[2]); */
/*       } */
      prob=f/sch_int;         /* //Calculo la probabilidad para cada gal. */
/*       //printf(" Galaxia %d f %g sc %g p %g s %g\n",i,f,sch_int,prob,logL); */
      logL-=log(prob);/*      //Voy restando el log las probabilidades */
/*       //porque Amoeba lo que sabe es minimizar solo!! */
    }
    printf(" s = %g  p0 %g p1 %g P2 %g p3 %f p4 %f\n",logL,p[0],p[1],p[2],p[3],p[4]); 
    return(logL);
  }
  return(logL);
}

double amofunc_schechterfitmag_d(int n, double *x, double *y, double *p)
{

  int i;
  double f,s;

  if(p[2]<=0) p[2]=-p[2];
  
  s=0.0;
  for(i=0; i<n; i++) {
    f=Schechter_M_LF(x[i],p[0],p[1],p[2]);
    /*       //    f += p[3]+p[4]*x[i]; */
      /*       //f += p[4]*x[i]; */
      /*       //printf(" f %f y %f\n",f,y[i]); */
      if(DEBUG) printf(" f %f log f %f y %f\n",f,log10(f),y[i]); 
      if (f!=0 && y[i]!=0) s += (log10(f)-y[i])*(log10(f)-y[i]);
  };

  if(DEBUG) printf(" s = %f  p0 %g p1 %g p2 %g\n",s,p[0],p[1],p[2]); 
  if(DEBUG) printf(" s = %f  p0 %g p1 %g P2 %g\n",s,p[0],p[1],p[2]); 
  
  return(s);
}


double amofunc_schechterfitlum_d(int n, double *x, double *y, double *p)
{

  int i;
  double f,s;


  if(p[2]<=0) p[2]=-p[2];
 
  s=0.0;
  for(i=0; i<n; i++) {
    f=Schechter_L_LF(pow(10.,x[i]),pow(10.,p[0]),p[1],p[2]);
    /*       //    f += p[3]+p[4]*x[i]; */
      /*       //f += p[4]*x[i]; */
      if(DEBUG2) printf("   f %f y %f\n",f,y[i]); 
      if (f!=0 && y[i]!=0)  s += (f-y[i])*(f-y[i]);
  }
  if(DEBUG2)  printf(" s = %f  p0 %g p1 %g P2 %g\n",s,p[0],p[1],p[2]); 
  
  return(s);
}


void  mrqfunc_schechterfitmag_d(double x,double *p,double *y,double *dyda,int n) {
  double elev;
  /* Valor de la funcion Schechter evaluada en M=x */
  elev=pow(10.,0.4*(p[0]-x));
/*   *y=log10(Schechter_M_LF(x,p[0],p[1],p[2])); */
  *y=log10(0.92103404)+log10(p[2])+(p[1]+1.)*(0.4*(p[0]-x))-elev*log10(exp(1));
  *y=log10(Schechter_M_LF(x,p[0],p[1],p[2]));
  dyda[0]=0.4*(1+p[1]-elev);                    /*Si, tras hacer las derivadas de log10(Phi(M)) queda eso */
  dyda[1]=log10(elev);
  dyda[2]=1./(p[2]*log(10.));
/*   printf("x %f *y %f y2 %f  p0 %f p1 %f p2 %f dyda1 %f dyda2 %f dyda3 %f\n",x,*y,y2,p[0],p[1],p[2],dyda[0],dyda[1],dyda[2]); */

}

void  mrqfunc_schechterfitlum_d(double x,double *p,double *y,double *dyda,int n) {
  double elev;
  double log_L_Lstar;
  log_L_Lstar=x-p[0];

  /* Valor de la funcion Schechter evaluada en L=x */
  elev=pow(10.,0.4*(p[0]-x));
/*   *y=log10(Schechter_M_LF(x,p[0],p[1],p[2])); */
  *y=log10(0.92103404)+log10(p[2])+(p[1]+1.)*(0.4*(p[0]-x))-elev*log10(exp(1));
  *y=Schechter_L_LF(pow(10.,x),pow(10.,p[0]),p[1],p[2]);
  dyda[0]=-1-p[1]+pow(10.,log_L_Lstar);                    /*Si, tras hacer las derivadas de log10(Phi(l)) con respecto a log(Phi_star) queda eso */
  dyda[1]=log_L_Lstar;
  dyda[2]=1./(p[2]*log(10.));
/*   printf("x %f *y %f y2 %f  p0 %f p1 %f p2 %f dyda1 %f dyda2 %f dyda3 %f\n",x,*y,y2,p[0],p[1],p[2],dyda[0],dyda[1],dyda[2]); */

}


void FitSchechter_L(struct lum_func lf, struct Schlf_L *schlf, double *chisq) {

  double *x,*y;
  double *sigy,*yfit;
  double par[3],sigpar[3];
  double **covarpar;
  int nfit;
  int i;
  int iter;
  int ipar[3];


  
  if(!lf.islum) {
    printf(" Internal error. Function was computed in magnitudes. This subrutine fits a luminosity Schechter function\n");
    exit(1);
  }

  x   =vector_d(lf.npoint);
  y   =vector_d(lf.npoint);
  sigy=vector_d(lf.npoint);
  yfit=vector_d(lf.npoint);
  covarpar=matrix_d(3,3);
  ipar[0]=1;ipar[1]=1;ipar[2]=1;
  
  par[2]=0.;
  nfit=0;
  for(i=0;i<lf.npoint;i++) {
    if(lf.lf[i]!=0) {
      x[nfit]=log10((double)lf.lum[i]);
      y[nfit]   =(double)lf.loglf[i];
      sigy[nfit]=(double)lf.errloglf[i];
      par[2]+=lf.lf[i];
      nfit++;
    }
  }
  par[1]=-0.8    ;sigpar[1]=0.1;   /* alfa */
  par[2]=(lf.lum[lf.npoint]-lf.lum[0])*par[2]/lf.npoint;
  sigpar[2]=par[2]/10.;
/*   par[0]=-20.0  ;sigpar[0]=.1; */   /* Mstar */ 
  par[1]=-0.9    ;sigpar[1]=0.01;   /* alfa */ 
  par[2]=0.0033;sigpar[2]=par[2]/10;   /* Phistar */ 
  par[0]=StMedia_d(nfit,x,&(sigpar[0]));

  if(DEBUG) printf(" Llamo amoeba\n");
  iter=Amoeba_d(nfit,x,y,3,par,sigpar,FTOL,10000,amofunc_schechterfitlum_d); 
  iter=Mrq_d(x,y,sigy,nfit,par,ipar,3,covarpar,chisq,mrqfunc_schechterfitlum_d);
  printf("Mrq Fit Schechter_L parameters fit: \n");
  printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",par[0],par[1],par[2],log10(par[2]));
  printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[2][2])/par[2]/log(10)); 
  for(i=0;i<nfit;i++) {
    yfit[i]=Schechter_L_LF(pow(10.,x[i]),pow(10.,par[0]),par[1],par[2]);
    if(DEBUG) printf("Bena x %g y %g y-yteor %f\n",x[i],yfit[i],yfit[i]-Schechter_L_LF(pow(10.,x[i]),pow(10.,36.477),-0.9,0.0033));
  }
  for(i=0;i<nfit;i++) {
    yfit[i]=Schechter_L_LF(pow(10.,x[i]),pow(10.,par[0]),par[1],par[2]);
    if(DEBUG) printf(" x %g y %g\n",x[i],yfit[i]);
  }

  schlf->alfa=par[1];
  schlf->erralfa=sigpar[1];
  schlf->phistar=par[2];
  schlf->errphistar=sigpar[2];
  schlf->Lstar=pow(10.,par[0]);
  schlf->errLstar=pow(10.,sigpar[0]);

  free(x);free(y);
  free(sigy);free(yfit);
  free(covarpar);
}



void FitSchechter_M(struct lum_func lf, struct Schlf_M *schlf, double *chisq) {

  double *x,*y;
  double *sigy,*yfit;
  double par[3],sigpar[3];
  double **covarpar;
  int nfit;
  int i;
  int iter;
  int ipar[3];


  
  if(lf.islum) {
    printf(" Internal error. Function was computed in luminosities. This subrutine fits a magnitude Schechter function\n");
    exit(1);
  }

  x   =vector_d(lf.npoint);
  y   =vector_d(lf.npoint);
  sigy=vector_d(lf.npoint);
  yfit=vector_d(lf.npoint);
  covarpar=matrix_d(3,3);
  ipar[0]=1;ipar[1]=1;ipar[2]=1;
  
  par[2]=0.;
  nfit=0;
  for(i=0;i<lf.npoint;i++) {
    if(lf.lf[i]!=0) {
      x[nfit]=(double)lf.mag[i];
      y[nfit]   =(double)lf.loglf[i];
      sigy[nfit]=(double)lf.errloglf[i];
      par[2]+=lf.lf[i];
      nfit++;
    }
  }
  par[1]=-0.8    ;sigpar[1]=0.1;   /* alfa */
  par[2]=(lf.mag[lf.npoint]-lf.mag[0])*par[2]/lf.npoint;
  sigpar[2]=par[2]/10.;
/*   par[0]=-20.0  ;sigpar[0]=.1; */   /* Mstar */ 
  par[1]=-0.9    ;sigpar[1]=0.01;   /* alfa */ 
  par[2]=0.0033;sigpar[2]=par[2]/10;   /* Phistar */ 
  par[0]=StMedia_d(nfit,x,&(sigpar[0]));

  iter=Amoeba_d(nfit,x,y,3,par,sigpar,FTOL,10000,amofunc_schechterfitmag_d);
  iter=Mrq_d(x,y,sigy,nfit,par,ipar,3,covarpar,chisq,mrqfunc_schechterfitmag_d);
  for(i=0;i<nfit;i++) {
    yfit[i]=log10(Schechter_M_LF(x[i],par[0],par[1],par[2]));
  }
  printf(" Mrq Fit Schechter_M parameters fit: \n");
  printf(" Mstar:    %g   alpha:    %g  Phistar  :  %g (log=%g) \n",par[0],par[1],par[2],log10(par[2]));
  printf(" E_Mstar:  %g   E_alpha:  %g  E_Phistar:  %g (log=%g)  \n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[2][2])/par[2]/log(10));
  
  for(i=0;i<nfit;i++) {
    yfit[i]=log10(Schechter_M_LF(x[i],par[0],par[1],par[2]));
    if(DEBUG) printf(" x %g y %g\n",x[i],yfit[i]);
  }

  schlf->alfa=par[1];
  schlf->erralfa=sigpar[1];
  schlf->phistar=par[2];
  schlf->errphistar=sigpar[2];
  schlf->Mstar=pow(10.,par[0]);
  schlf->errMstar=pow(10.,sigpar[0]);

  free(x);free(y);
  free(sigy);free(yfit);
  free(covarpar);
}


void FitSchechter_Step_L(struct Steplf_L lf, struct Schlf_L *schlf, double *chisq) {

  double *x,*y;
  double *sigy,*yfit;
  double par[3],sigpar[3];
  double **covarpar;
  int nfit;
  int i;
  int iter;
  int ipar[3];

  x   =vector_d(lf.nbin);
  y   =vector_d(lf.nbin);
  sigy=vector_d(lf.nbin);
  yfit=vector_d(lf.nbin);
  covarpar=matrix_d(3,3);
  ipar[0]=1;ipar[1]=1;ipar[2]=1;
  
  par[2]=0.;
  nfit=0;
  for(i=0;i<lf.nbin;i++) {
    if(lf.lf[i]!=0) {
      x[nfit]=log10((double)lf.lumi[i]);
      y[nfit]   =(double)lf.lf[i];
      sigy[nfit]=(double)lf.errlf[i];
      par[2]+=lf.lf[i];
      nfit++;
    }
  }
  par[1]=-0.8    ;sigpar[1]=0.1;   /* alfa */
  par[2]=(lf.lumi[lf.nbin]-lf.lumi[0])*par[2]/lf.nbin;
  sigpar[2]=par[2]/10.;
/*   par[0]=-20.0  ;sigpar[0]=.1; */   /* Mstar */ 
  par[1]=-0.9    ;sigpar[1]=0.01;   /* alfa */ 
  par[2]=0.0033;sigpar[2]=par[2]/10;   /* Phistar */ 
  par[0]=StMedia_d(nfit,x,&(sigpar[0]));

  if(DEBUG) printf(" Llamo amoeba\n");
  iter=Amoeba_d(nfit,x,y,3,par,sigpar,FTOL,10000,amofunc_schechterfitlum_d); 
  iter=Mrq_d(x,y,sigy,nfit,par,ipar,3,covarpar,chisq,mrqfunc_schechterfitlum_d);
  printf("Mrq Fit Schechter_L parameters fit: \n");
  printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",par[0],par[1],par[2],log10(par[2]));
  printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[2][2])/par[2]/log(10)); 
  for(i=0;i<nfit;i++) {
    yfit[i]=Schechter_L_LF(pow(10.,x[i]),pow(10.,par[0]),par[1],par[2]);
    if(DEBUG) printf("Bena x %g y %g y-yteor %f\n",x[i],yfit[i],yfit[i]-Schechter_L_LF(pow(10.,x[i]),pow(10.,36.477),-0.9,0.0033));
  }
  for(i=0;i<nfit;i++) {
    yfit[i]=Schechter_L_LF(pow(10.,x[i]),pow(10.,par[0]),par[1],par[2]);
    if(DEBUG) printf(" x %g y %g\n",x[i],yfit[i]);
  }

  schlf->alfa=par[1];
  schlf->erralfa=sigpar[1];
  schlf->phistar=par[2];
  schlf->errphistar=sigpar[2];
  schlf->Lstar=pow(10.,par[0]);
  schlf->errLstar=pow(10.,sigpar[0]);

  free(x);free(y);
  free(sigy);free(yfit);
  free_matrix_d(covarpar,3,3);
}



void FitSchechter_Step_M(struct Steplf_M lf, struct Schlf_M *schlf, double *chisq) {

  double *x,*y;
  double *sigy,*yfit;
  double par[3],sigpar[3];
  double **covarpar;
  int nfit;
  int i;
  int iter;
  int ipar[3];

  x   =vector_d(lf.nbin);
  y   =vector_d(lf.nbin);
  sigy=vector_d(lf.nbin);
  yfit=vector_d(lf.nbin);
  covarpar=matrix_d(3,3);
  ipar[0]=1;ipar[1]=1;ipar[2]=1;
  
  par[2]=0.;
  nfit=0;
  for(i=0;i<lf.nbin;i++) {
    printf(" mag %f lf %f\n",lf.magni[i],lf.lf[i]);
    if(fabs(lf.errlf[i])<0.5) {
      x[nfit]=(double)lf.magni[i];
      y[nfit]   =(double)lf.lf[i];
      sigy[nfit]=(double)lf.errlf[i];
      par[2]+=lf.lf[i];
      printf(" x %f y %f\n",x[nfit],y[nfit]);
      nfit++;
    }
  }
  par[1]=-0.8    ;sigpar[1]=0.1;   /* alfa */
  par[2]=(lf.magni[lf.nbin]-lf.magni[0])*par[2]/lf.nbin;
  sigpar[2]=par[2]/10.;
/*   par[0]=-20.0  ;sigpar[0]=.1; */   /* Mstar */ 
  par[1]=-0.9    ;sigpar[1]=0.01;   /* alfa */ 
  par[2]=0.0033;sigpar[2]=par[2]/10;   /* Phistar */ 
  par[0]=StMedia_d(nfit,x,&(sigpar[0]));
  for(i=0;i<nfit;i++) printf("x  %d %f y %f sig %f \n",i,x[i],y[i],sigy[i]);

  printf(" INIT\n");
  printf(" Mstar:    %g   alpha:    %g  Phistar  :  %g (log=%g) \n",par[0],par[1],par[2],log10(par[2]));
  iter=Amoeba_d(nfit,x,y,3,par,sigpar,FTOL,10000,amofunc_schechterfitmag_d);
  printf(" First Amoeba\n");
  printf(" Mstar:    %g   alpha:    %g  Phistar  :  %g (log=%g) \n",par[0],par[1],par[2],log10(par[2]));
  iter=Mrq_d(x,y,sigy,nfit,par,ipar,3,covarpar,chisq,mrqfunc_schechterfitmag_d);
  printf(" A partir de aqui\n");   
  for(i=0;i<nfit;i++) {  
    yfit[i]=log10(Schechter_M_LF(x[i],par[0],par[1],par[2]));
  }
  printf(" Mrq Fit Schechter_M parameters fit: \n");
  printf(" Mstar:    %g   alpha:    %g  Phistar  :  %g (log=%g) \n",par[0],par[1],par[2],log10(par[2]));
  printf(" E_Mstar:  %g   E_alpha:  %g  E_Phistar:  %g (log=%g)  \n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[2][2])/par[2]/log(10));
  
  for(i=0;i<nfit;i++) {
    yfit[i]=log10(Schechter_M_LF(x[i],par[0],par[1],par[2]));
    if(DEBUG) printf(" x %g y %g\n",x[i],yfit[i]);
  }

  schlf->alfa=par[1];
  schlf->erralfa=sigpar[1];
  schlf->phistar=par[2];
  schlf->errphistar=sigpar[2];
  schlf->Mstar=par[0];
  schlf->errMstar=sigpar[0]; 

  free(x);free(y);
  free(sigy);free(yfit);
  free(covarpar);
}
