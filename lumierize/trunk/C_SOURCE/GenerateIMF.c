#include "modulos.h"

#define NTOT 100000
#define NBIN 25
#define NPOINT 10000

#define XMIN 10
#define XMAX 100

float Scalo(float M,float Mmin, float Mmax,float alfa);

float Scalodev(float Mmin, float Mmax,float alfa);

int main () 
{
  
  int i;
  int ntot=NTOT;
  float Mtot=1e3;
  float media,sigma;
  float x[NTOT];
  float xx[NPOINT],yy[NPOINT];

  float Mmin=0.1;
  float Mmax=100.;
  double M=0;
  float alfa=1.35;
  printf(" Numero total de estrellas (%d maximo): ",NTOT);
  ntot=readi(ntot);
  printf(" Masa total del cumulo (0=generar el numero anterior): ");
  Mtot=readf(Mtot);
  printf(" Limite inferior de la masa: ");
  Mmin=readf(Mmin);
  printf(" Limite superior de la masa: ");
  Mmax=readf(Mmax);
  printf(" Parametro alfa: ");
  alfa=readf(alfa);

  

  
  if(Mtot==0) {
    for(i=0;i<ntot;i++) {
      x[i]=Scalodev(Mmin,Mmax,alfa);
      printf(" x %f\n",x[i]);
    }
  }
  else {
    i=0;
    while (M<Mtot) {
      x[i]=Scalodev(Mmin,Mmax,alfa);
      M+=x[i];
      printf(" x %f   Mtot %f\n",x[i],M);
      i++;
    }
    ntot=i-1;
    printf(" Numero total de estrellas generadas: %d\n",ntot);
  }

  for(i=0;i<NPOINT;i++) {
    xx[i]=Mmin+i/(NPOINT-1.)*(Mmax-Mmin);
    yy[i]=ntot*Scalo(xx[i],Mmin,Mmax,alfa)*(XMAX-XMIN)/NBIN;
    //    printf(" xx %f yy %f\n",xx[i],yy[i]);
  }


  
  cpgopen("?");
  cpgswin(XMIN,XMAX,0.,yy[(int)((XMIN-Mmin)/(Mmax-Mmin)*(NPOINT-1.))]*1.5);
  //  cpgask(0);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgline(NPOINT,xx,yy);
  cpghist(ntot,x,XMIN,XMAX,NBIN,1);
  cpgend(); 
  printf("\n\n");
  printf("Distribucion inicial de x:\n"); 
  media= StMedia(ntot,x,&sigma); 
  printf(" La media es %f y la desviacion %f  S/N %f\n\n",media,sigma,media/sigma); 
  return(0);
}

float Scalo(float M,float Mmin, float Mmax,float alfa) {

  float A;

  A=(1-alfa)/(pow(Mmax,1-alfa)-pow(Mmin,1-alfa));

  return(A*pow(M,-alfa));
  
}

float Scalodev(float Mmin, float Mmax,float alfa) {
  
  float A;
  float tmp;
  float dum;
  static long idum =-1;
  if(idum==-1) idum=-(long)time(NULL)/2;
  do dum=ran2(&idum);
  while (dum == 0.0);

  A=(1-alfa)/(pow(Mmax,1-alfa)-pow(Mmin,1-alfa));

  tmp=dum*(1-alfa)/A+pow(Mmin,1-alfa);
  tmp=pow(tmp,1./(1-alfa));

  return(tmp);
  
}

