#include "modulos.h"
      
#define NPOB          500
#define MEANPOB       25.
#define SIGMAPOB       5.
#define SIGMAMED       6.
#define SIGMAMEDDEV    1. 
                                        
#define NBIN           50
     
#define TOL 1.e-7
#define NEXPERIMENTS  100

#define VERBOSE         1
                                 
/* void        ML_g_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma); */
/* void   ML_g_g_corr(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma); */
/* void ML_g_g_corr_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma); */

   
 
int main() { 
 
  int i;
  int n;
  float pi=4*atan(1.); 


  float h[NPOB];
  float x[NPOB];
  float errx[NPOB];
  
  float meanh,sigmah;
  float meanx,sigmax;

  float newmeanh[NEXPERIMENTS],newsigmah[NEXPERIMENTS];
  float errmeanh[NEXPERIMENTS],errsigmah[NEXPERIMENTS],covarmeansigma[NEXPERIMENTS];

  float meannewmeanh,stdnewmeanh;
  float meannewsigmah,stdnewsigmah;
  float meanerrmeanh,stderrmeanh;
  float meanerrsigmah,stderrsigmah; 

  float covar,sigmasigma;

  float dumx,dumy;

  cpgopen("?");
  cpgask(0);


   
  for(n=0;n<NEXPERIMENTS;n++) {
    printf(" EXPER %d\n",n);


    
    for(i=0;i<NPOB;i++) {
      h[i]   =MEANPOB+SIGMAPOB*Gasdev();
/*       errx[i]=Poidev(fabs(SIGMAMED*h[i]/MEANPOB)); */
/*       errx[i]=fabs(SIGMAMED/h[i]*MEANPOB); */
      /*     errx[i]=SIGMAMED*exp((MEANPOB-h[i])/20.); */
      /*     printf(" h %f e %f\n",h[i],errx[i]);  */
      errx[i]=fabs(SIGMAMED+SIGMAMEDDEV*Gasdev());  
/*       errx[i]=SIGMAMED+h[i]+fabs(h[i])*Gasdev(); */
      errx[i]=fabs(SIGMAMED+h[i]*h[i]/MEANPOB/MEANPOB+SIGMAMEDDEV*Gasdev());
      errx[i]=fabs(SIGMAMED+h[i]*h[i]/MEANPOB/MEANPOB+SIGMAMEDDEV*Gasdev());
/*       errx[i]=SIGMAMED;      */
      if(errx[i]==0) errx[i]=1.;
      x[i]   =h[i]+errx[i]*Gasdev();
    }   
                         
    if(VERBOSE) {
      printf(" He generado %d muestras aleatorias de una poblacion\n",NPOB);
      printf(" con dist. gaussiana de media %f y desviacion %f\n",MEANPOB,SIGMAPOB);
      printf(" Histograma rojo. Variable h\n");
      
      printf(" Para cada una de las muestras he realizado una medida\n");
      printf(" con un error medio de %f\n",SIGMAMED);
      printf(" HIstograma azul. Variable x\n");
      
    }
    cpgsvp(0.05,0.65,0.1,0.9);
    cpgeras();
    cpgsci(1);
    
    cpgswin(MEANPOB-(SIGMAPOB+SIGMAMED)*3.,MEANPOB+(SIGMAPOB+SIGMAMED)*3.,0.,(float)NPOB/NBIN*9.);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
    
    cpgsci(2);
    cpghist(NPOB,h,MEANPOB-(SIGMAPOB+SIGMAMED)*3.,MEANPOB+(SIGMAPOB+SIGMAMED)*3.,NBIN,1);
    cpgsci(4);
    cpghist(NPOB,x,MEANPOB-(SIGMAPOB+SIGMAMED)*3.,MEANPOB+(SIGMAPOB+SIGMAMED)*3.,NBIN,1);
    
    meanh=StMedia(NPOB,h,&sigmah); 
    meanx=StMedia(NPOB,x,&sigmax); 
    if(VERBOSE) {
      printf("\n");
      printf(" La media            de h es %f y la desviacion %f  S/N %f\n\n",meanh,sigmah,meanh/sigmah);
      printf(" La media            de x es %f y la desviacion %f  S/N %f\n\n",meanx,sigmax,meanx/sigmax);
      meanx=StWeightMedia(NPOB,x,errx,&sigmax); 
      printf(" La media promediada de x es %f y la desviacion %f  S/N %f\n\n",meanx,sigmax,meanx/sigmax);
      
      covar=StCovar(NPOB,h,errx,&sigmah,&sigmasigma);
      printf(" La covarianza entre h y errx es %f  rho=%f \n",covar,sqrt(covar/(sigmah*sigmah+sigmasigma*sigmasigma)));
      covar=StCovar(NPOB,x,errx,&sigmax,&sigmasigma);
      printf(" La covarianza entre x y errx es %f  rho=%f \n",covar,sqrt(covar/(sigmah*sigmah+sigmasigma*sigmasigma)));
    }

    cpgsci(2);
    cpgmove(MEANPOB-(SIGMAPOB+SIGMAMED)*3.,0.);
    for(i=0;i<1000;i++) {
      dumx=(MEANPOB-(SIGMAPOB+SIGMAMED)*3.)+i/999.*((SIGMAPOB+SIGMAMED)*6.-1.);
      cpgdraw(dumx,NPOB*(SIGMAPOB+SIGMAMED)*6./NBIN*gaussian(dumx,MEANPOB,SIGMAPOB));
      /*     printf(" x %f y %f\n",dumx,NPOB*(SIGMAPOB+SIGMAMED)*6./NBIN*gaussian(dumx,MEANPOB,SIGMAPOB)); */
    }
    cpgsci(4);
    cpgmove(MEANPOB-(SIGMAPOB+SIGMAMED)*3.,0.);
    for(i=0;i<1000;i++) {
      dumx=(MEANPOB-(SIGMAPOB+SIGMAMED*3.))+i/999.*((SIGMAPOB+SIGMAMED)*6.-1.);
      dumy= 1./sqrt(2.*pi)/SIGMAMED/SIGMAPOB/sqrt(1./SIGMAMED/SIGMAMED+1./SIGMAPOB/SIGMAPOB)*
	exp(-dumx*dumx/2./SIGMAMED/SIGMAMED)*exp(-MEANPOB*MEANPOB/2./SIGMAPOB/SIGMAPOB)*
	exp((dumx/SIGMAMED/SIGMAMED+MEANPOB/SIGMAPOB/SIGMAPOB)*(dumx/SIGMAMED/SIGMAMED+MEANPOB/SIGMAPOB/SIGMAPOB)/
	    2./(1./SIGMAMED/SIGMAMED+1./SIGMAPOB/SIGMAPOB));
      cpgdraw(dumx,NPOB*(SIGMAPOB+SIGMAMED)*6./NBIN*dumy);
    }
    
     
    cpgsvp(0.7,0.95,0.1,0.9);
    cpgsci(1);
/*     cpgswin(MEANPOB-(SIGMAPOB+SIGMAMED)*3.,MEANPOB+(SIGMAPOB+SIGMAMED)*3.,SIGMAMED-2.*SIGMAMEDDEV,SIGMAMED+2.*SIGMAMEDDEV); */
    cpgswin(SIGMAMED-10.*SIGMAMEDDEV,SIGMAMED+10.*SIGMAMEDDEV,MEANPOB-(SIGMAPOB+SIGMAMED)*3.,MEANPOB+(SIGMAPOB+SIGMAMED)*3.);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
   
    for(i=0;i<NPOB;i++) {
      cpgpt1(errx[i],h[i],2);  
    }
              
    ML_g_g(NPOB,x,errx,newmeanh+n,newsigmah+n,errmeanh+n,errsigmah+n,covarmeansigma+n); 
    printf("Formulae  Mean %.11f +/- %.11f  Sigma %.11f +/- %.11f\n",newmeanh[n],errmeanh[n],newsigmah[n],errsigmah[n]);
    ML_g_g_corr(NPOB,x,errx,newmeanh+n,newsigmah+n,errmeanh+n,errsigmah+n,covarmeansigma+n);  
    printf("Amoeba    Mean %.11f +/- %.11f  Sigma %.11f +/- %.11f\n",newmeanh[n],errmeanh[n],newsigmah[n],errsigmah[n]);
    ML_g_g_corr_g(NPOB,x,errx,newmeanh+n,newsigmah+n,errmeanh+n,errsigmah+n,covarmeansigma+n);   
    printf("AmoebaCor Mean %.11f +/- %.11f  Sigma %.11f +/- %.11f\n",newmeanh[n],errmeanh[n],newsigmah[n],errsigmah[n]); 

    


     
       
  }
  meanerrmeanh=StMedia(NEXPERIMENTS,errmeanh,&stderrmeanh);       
  meanerrsigmah=StMedia(NEXPERIMENTS,errsigmah,&stderrsigmah); 

 
  meannewmeanh=StMedia(NEXPERIMENTS,newmeanh,&stdnewmeanh);       
  meannewsigmah=StMedia(NEXPERIMENTS,newsigmah,&stdnewsigmah); 

  cpgsvp(0.05,0.95,0.05,0.95);
  cpgpage();
  cpgswin(meannewmeanh-(stdnewmeanh)*3.,meannewmeanh+(stdnewmeanh)*3.,0.,(float)NEXPERIMENTS/NBIN*9.);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  
  cpgsci(2);
  cpghist(NEXPERIMENTS,newmeanh,meannewmeanh-(stdnewmeanh)*2.,meannewmeanh+(stdnewmeanh)*2.,NBIN,1); 


  
  printf(" Al final, valor medio del estimador mean  ML: mean(mean_ML) = %f   std(mean_ML) = %f\n",meannewmeanh,stdnewmeanh);
  printf(" El error medio estimado para mean era %f std = %f\n",meanerrmeanh,stderrmeanh);
  printf(" Al final, valor medio del estimador sigma ML: mean(sigma_ML)= %f   std(sigma_ML)= %f\n",meannewsigmah,stdnewsigmah);
  printf(" El error medio estimado para sigma era %f std = %f\n",meanerrsigmah,stderrsigmah);
  
  return 0;
}  
  
       

void StBayesMedia(int n,float *x,float *errx,float *mean,float *sigma) {
  
  double sold,snew;
  double m;
  int i;

  float smin,smax;

/*   float xmean[NPOB]; */

  float ssum;
 
  double tmp1,tmp2,tmp3;
  double tmp4,tmp5,tmp6;

  MinMax(n,errx,&smin,&smax);
  m=StMedia(n,x,&ssum); 
  snew=ssum;
  sold=snew;
  do {
    sold=snew;
    snew=0;
    tmp1=0;
    tmp2=0;
    tmp4=0;
    tmp5=0;
    for(i=0;i<n;i++) {
      tmp1+=x[i]/errx[i]/errx[i]/(1+sold*sold/errx[i]/errx[i]);
      tmp2+=1./errx[i]/errx[i]/(1+sold*sold/errx[i]/errx[i]);
      tmp4+=x[i]/(1+errx[i]*errx[i]/sold/sold);
      tmp5+=1./(1+sold*sold/errx[i]/errx[i]);

      
    }
/*     snew=sqrt(snew/n); */
    
/*     m=StMedia(n,xmean,&ssum); */
/*     printf(" la sigma y media bien:     %f    %f\n",ssum,m); */

/*     m=0; */
/*     ssum=0; */
    tmp3=0;
    tmp6=0;
    for(i=0;i<n;i++) {
       tmp6+=(x[i]/errx[i]/errx[i]+m/sold/sold)*(x[i]/errx[i]/errx[i]+m/sold/sold)/(1./errx[i]/errx[i]+1./sold/sold)/(1./errx[i]/errx[i]+1./sold/sold); 
       tmp6+=(-2*m*x[i]/errx[i]/errx[i]-2*m*m/sold/sold)/(1./errx[i]/errx[i]+1./sold/sold); 
       tmp3+=(x[i]-m)*(x[i]-m)/sold/sold/(1+errx[i]*errx[i]/sold/sold)/(1+errx[i]*errx[i]/sold/sold);

    }


    snew=sqrt((n*m*m+tmp6)/((float)n-tmp5));   
    m=tmp4/((float)n-tmp5);  

/*     Otra manera de hacerlo: */
/*     m=tmp1/tmp2;  */
/*     snew=sqrt(tmp3/tmp2);  */
 
/*     snew=SIGMAPOB; */
    
/*     printf(" n %d sold %f snew  %f m  %f   s2 %f m2 %f tmp6 %f tmp5 %f tmp3 %f\n",n,sold,snew,m,sqrt((n*m*m+tmp6)/(n-tmp5)),tmp4/((float)n-tmp5),tmp6,tmp5,tmp3); */


  } while (fabs(sold-snew)/fabs(snew)> TOL);

  
  printf("Formulae  Mean %14.10f Sigma %14.10f\n",m,snew);


  *mean=m;
  *sigma=snew;

/*   printf("   Weight  sigma  %f mean  %f\n",ssum,m);  */
  
  return;
}
