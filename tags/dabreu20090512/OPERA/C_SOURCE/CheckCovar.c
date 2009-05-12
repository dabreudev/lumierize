#include "modulos.h"

#define NPOINT 500000
#define NBIN 50

float funeta(float x, float y);
float funpsi(float x, float y);

int main ()
{
  
  float x[NPOINT],y[NPOINT];
  float xmedia,ymedia;
  float sigmax,sigmay,covarxy;
  float eta[NPOINT],psi[NPOINT];
  float etamedia,psimedia;
  float sigmaeta,sigmapsi,covaretapsi;
  float fnul1,fnul2;
  int i;
  static long idum =-1;
  if(idum==-1) idum=-(long)time(NULL)/2;

  sigmax=.1;
  sigmay=.05;
  xmedia=1.;
  ymedia=0.;

  for(i=0;i<NPOINT;i++) {
    x[i]=xmedia+Gasdev()*sigmax;
    y[i]=ymedia+Gasdev()*sigmay;
    y[i]=(x[i]+y[i])/(sqrt(2.));  
    
    eta[i]=funeta(x[i],y[i]);
    psi[i]=funpsi(x[i],y[i]);
/*      printf(" x %f eta %f\n",x[i],eta[i]); */
  }

  xmedia= StMedia(NPOINT,x,&sigmax); 
  ymedia= StMedia(NPOINT,y,&sigmay); 
  covarxy=StCovar(NPOINT,x,y,&fnul1,&fnul2);

  etamedia= StMedia(NPOINT,eta,&sigmaeta); 
  psimedia= StMedia(NPOINT,psi,&sigmapsi); 
  covaretapsi=StCovar(NPOINT,eta,psi,&fnul1,&fnul2);
  
  printf(" Las funciones son:\n eta=x*x\n psi=x+y\n");
  printf(" Desviaciones teoricas:\n sig2(eta)= (d(eta)/d(x) sig(x) )**2 + (d(eta)/d(y) sig(y) )**2 + d(eta) /d(x)d(y) cov(xy)\n");
  printf("  sig2(psi)= (d(psi)/d(x) sig(x) )**2 + (d(psi)/d(y) sig(y) )**2 + d(psi) /d(x)d(y) cov(xy)\n");
  printf(" cov(eta,psi)=d(eta) /d(x)d(y) * d(psi) /d(x) * d(psi) /d(y) * cov(x,y)\n");
  printf(" En nuestro caso:\n sig2(eta)=4*x*x*sig(x)*sig(x)\n");
  printf(" sig2(psi)=(1+y)*(1+y)*sig(x)*sig(x)+(1+x)*(1+x)*sig(y)*sig(y) + (1+x)*(1+y)*cov(x,y)\n");
  printf(" cov(eta,psi)=d(eta)/d(x) * d(psi)/d(x) * sig(x)*sig(x) + d(eta)/d(y) * d(psi)/d(y) * sig(y)*sig(y) + (d(eta)/d(x)*d(psi)/d(y)+d(eta)/d(y)*d(psi)/d(x))*cov(x,y)\n");
  printf(" La media de x es %f y la desviacion %f \n\n",xmedia,sigmax);   
  printf(" La media de y es %f y la desviacion %f \n\n",ymedia,sigmay);   

  printf(" Covarianza entre x e y: %f\n",covarxy);

  printf(" La media de eta es %f y la desviacion %f \n\n",etamedia,sigmaeta);   
  printf(" La media de psi es %f y la desviacion %f \n\n",psimedia,sigmapsi);   

  printf(" Covarianza entre eta y psi: %f\n",covaretapsi);

  printf("Teoricas: \n");
  printf(" La media de eta es %f y la desviacion %f \n\n",xmedia*xmedia+ymedia,sqrt(4*xmedia*xmedia*sigmax*sigmax+sigmay*sigmay+4*xmedia*covarxy));   
  printf(" La media de psi es %f y la desviacion %f \n\n",xmedia+ymedia+xmedia*ymedia,sqrt((1+ymedia)*(1+ymedia)*sigmax*sigmax+(1+xmedia)*(1+xmedia)*sigmay*sigmay+2*(1+xmedia)*(1+ymedia)*covarxy));   
  
  printf(" Covarianza entre eta y psi: %f\n",2*xmedia*(1+ymedia)*sigmax*sigmax+(1+xmedia)*sigmay*sigmay+(2*xmedia*(1+xmedia) + (1+ymedia))*covarxy);
  

/*   cpgopen("?"); */
/*   cpgswin(0.,10.,0.,yy[0]*1.3); */
/*   //  cpgask(0); */
/*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
/*   cpgline(1000,xx,yy); */
/*   cpghist(NTOT,x,0.,10.,NBIN,1); */
/*   cpgend();  */
/*   printf("\n\n"); */
/*   printf("Distribucion inicial de x:\n");  */
/*   media= StMedia(NTOT,x,&sigma);  */

    return 0;
}


float funeta(float x, float y) {
  return(x*x+y);
}

float funpsi(float x, float y) {
  return(x+y+x*y);
}


