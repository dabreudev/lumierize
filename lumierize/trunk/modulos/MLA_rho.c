#include "modulos.h"
#define FTOL  1e-13
#define MAXITER  600
#define FTOL2  1e-12
#define MAXITER2  100
#define NCONFL 6

#define DEBUG 0



double Amoe_Funk_rho_main(int n, double *x, double *y, double *p);
double Amoe_Funk_rho_conf(int n, double *x, double *y, double *p);
void   EmpiricalCovars_rho(int n,double *ra, double *dec, int *ima, double *par, double *sigpar, struct SurveyDB sdb, double *probdetec, double *dens, double *errdens);

double *pd_rho;
int    *ima_rho;
double *area_rho;
int    *nima_rho;
int     nsur_rho;
double MLmax_rho;
double conflim_rho;

int MLA_rho(int n,double *ra, double *dec, int *ima, struct SurveyDB sdb, double *probdetec, double *dens, double *errdens) {

  int iter_amo=0;
  double par[1];
  double sigpar[1];
  int i,j;
  struct SurveyDB s_item;

  pd_rho=probdetec;
  ima_rho=ima;
  nsur_rho=sdb.nitems;

  area_rho=vector_d(sdb.nitems);
  s_item.nitems=1;
  for(j=0;j<sdb.nitems;j++) {
    s_item.si=&(sdb.si[j]);
    area_rho[j]=s_item.si[0].xdim*s_item.si[0].ydim/3600./3600./180./180.*M_PI*M_PI;
    /* No hay que tener en cuenta el coseno!! Da igual donde este la imagen, ocupara lo mismo */
  }
  nima_rho=vector_i(sdb.nitems);
  for(j=0;j<sdb.nitems;j++) nima_rho[j]=0;

  for(i=0;i<n;i++) { 
    nima_rho[ima_rho[i]]++;
  }
  if(DEBUG) for(j=0;j<sdb.nitems;j++) printf(" Imagen %d. Num obj %d. Area %f\n",j,nima_rho[j],area_rho[j]*180*180/M_PI/M_PI);
  for(j=0;j<sdb.nitems;j++) printf(" Imagen %d. Num obj %d. Area %f\n",j,nima_rho[j],area_rho[j]*180*180/M_PI/M_PI);
 
  while(iter_amo==0) { 
    par[0]=log(1*180*180/M_PI/M_PI);
    sigpar[0]=0.2;
    iter_amo=Amoeba_d(n,ra,dec,1,par,sigpar,FTOL,MAXITER,Amoe_Funk_rho_main);
  }
  MLmax_rho=Amoe_Funk_rho_main(n,ra,dec,par);
  conflim_rho=exp(-.5/10.);
  EmpiricalCovars_rho( n,ra,dec,ima, par,sigpar,sdb,probdetec,dens,errdens);


  *dens=exp(par[0]);
/*   printf(" YA FUERA dens %f\n",*errdens); */
  printf(" Con uno: %f\n",nima_rho[0]/area_rho[0]/pd_rho[0]/180/180*M_PI*M_PI);  
  printf(" Con dos: %f\n",(nima_rho[0]+nima_rho[1])/(area_rho[0]*pd_rho[0]+area_rho[1]*pd_rho[1])/180/180*M_PI*M_PI);   
/*   printf(" Area %f  NIMA %d\n",area_rho[0],nima_rho[0]); */
/*   printf(" Con AMO: %f\n",*dens/180/180*M_PI*M_PI); */
  free(area_rho);
  free(nima_rho);
  return(iter_amo);
}



double Amoe_Funk_rho_main(int n, double *x, double *y, double *p) {
/*   int i;  */
  double logL=0.;
  double rho;
  double nmed;
  int j;

/*   double probexista; */
/*   double probdetec; */
/*   double probnorm; */

  rho=exp(p[0]);
  
/*   for(i=0;i<n;i++) { */
/*     probexista=rho; */
/*     probdetec=pd_rho[ima_rho[i]]; */
/*     probnorm=area_rho[ima_rho[i]]*pd_rho[ima_rho[i]]*rho; */
/*     printf(" ex %f det %f logL %f\n",probexista,probdetec,logL);  */
/*     logL-=log(probexista*probdetec/probnorm); */
    /* Todo esto queda la inversa de la densidad, vamos que no depende
       nada de la densidad. Lo pongo para que se vea. Por eso esta comentado */
  /*   } */
  
  for(j=0;j<nsur_rho;j++) {
    nmed=area_rho[j]*pd_rho[j]*rho;
    logL-=   (nima_rho[j]*log(nmed) - nmed - gammln((double)nima_rho[j]+1.));
  }
  

  
  /* printf(" ITER Rho %f  logL %g\n",rho/180/180*M_PI*M_PI,logL); */

  return(logL);
}

double Amoe_Funk_rho_conf(int n, double *x, double *y, double *p) {
  return(fabs(Amoe_Funk_rho_main(n,x,y,p)-(MLmax_rho-log(conflim_rho)))); 
}



void   EmpiricalCovars_rho(int n,double *ra, double *dec, int *ima, double *par, double *sigpar, struct SurveyDB sdb, double *probdetec, double *dens, double *errdens) {

  int i,j;  
  double *parconf;
  double *sigparconf; 
  double **invcovar;
  double **covar;
  double **parelip;
  double *y;

  double **bb;
  int nconfl,nconflini;
  double first, median, third, *distmax;

  int iter_c;

  nconfl=NCONFL;
  nconflini=NCONFL;
  
  bb=matrix_d(1,1);
  y=vector_d(n);
  parconf=vector_d(1);
  sigparconf=vector_d(1);
  parelip=matrix_d(1,nconflini);
  invcovar=matrix_d(1 ,1 );
  covar=matrix_d(1 ,1 );
  distmax=vector_d(1);

  printf(" Error step  ");
  for(i=0;i<nconfl;i++) printf(".");
  for(i=0;i<nconfl;i++) printf("\b");
  fflush(NULL);

  
  for(i=0;i<nconfl;i++) {
    printf("#");
    fflush(NULL);
    parconf[0]=par[0]+2*sigpar[0]*Gasdev();
    sigparconf[0]=sigpar[0];
    if(i>(int)(nconfl/2.)) {
      parconf[0]=par[0]-((parelip[0])[(int)(i-nconfl/2.)+1]-par[0]);
    }
    iter_c=Amoeba_d(n,ra,dec,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_rho_conf);
    if(DEBUG) printf(" %d SOL iter %d par0 %f \n", i,  iter_c ,parconf[0]);
    (parelip[0])[i]=parconf[0];
    if(iter_c==0 && Amoe_Funk_rho_conf(n,ra,dec,parconf)>FTOL2 ) {
      printf("\b.\b");
      fflush(NULL);
      i--;
    }
  }
  for(i=0;i<nconfl;i++) {
    (parelip[0])[i]-=par[0];
  }
  for(j=0;j<1;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<1;j++) {
      if(fabs((parelip[j])[i])>2*2*distmax[j]/1.35) {
	for(j=0;j<1;j++) {
	  memcpy(parelip[j]+i,parelip[j]+i+1,(nconfl-i-1)*sizeof(double));
	}
	i--;
	nconfl--;
	break;
      }
    }
  }

  MCElipN_d(nconfl,1,parelip,invcovar);
/*   printf(" invcovar %f\n",invcovar[0][0]); */
  for(i=0;i<1;i++) {
    for(j=0;j<1;j++) {
      covar[i][j]=invcovar[i][j];
    }
  }
  gaussj_d(covar,1,bb,1);
/*   printf(" covar %f\n",covar[0][0]); */
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<1;i++) {
    for(j=0;j<1;j++) {
      covar[i][j]/=(-2*log(conflim_rho));
    }
  }

/*   printf(" covar %f\n",covar[0][0]); */

  *errdens=exp(par[0])*sqrt(covar[0][0]);

/*   printf(" errdems %f\n",*errdens); */

  free_matrix_d(bb,1,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,1,nconflini);
  free_matrix_d(invcovar,1  ,1);
  free_matrix_d(   covar,1  ,1);


}
