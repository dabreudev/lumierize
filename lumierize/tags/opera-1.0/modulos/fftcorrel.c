#include "modulos.h"
#define DEBUG 0
#define DEBUG2 0
#define PLOTDEBUG 0

void fungauss(double x,double *p,double *y,double *dyda,int n);


void fftcorrel(int n, float *data1, float *data2, float *xcorr, float *fcorr) {

  int i;
  int n2;

  float *xr1,*xr2;
  float *xi1,*xi2;

  float *xr,*xi;
  
  if(DEBUG) {
    printf(" DENTRO FFTCORREL\n");
    for(i=0;i<n;i++) 
      printf(" %d  %f %f \n",i,data1[i],data2[i]);
  }

  xr1=vector_f(n);
  xr2=vector_f(n);
  xi1=vector_f(n);
  xi2=vector_f(n);
  xr =vector_f(n);
  xi =vector_f(n);
  
  for(i=0;i<n;i++) {
    xr1[i]=data1[i];
    xr2[i]=data2[i];
    xi1[i]=0.;
    xi2[i]=0.;
  }

  cfft(n,xr1,xi1,1);
  cfft(n,xr2,xi2,1);

  if(DEBUG) {
    printf(" FFT DE LAS DOS\n");
    for(i=0;i<n;i++) 
      printf(" %d  %f %f   %f %f\n",i,xr1[i],xi1[i],xr2[i],xi2[i]);
  }


  for(i=0;i<n;i++)    xi2[i]=-xi2[i];
  
  for(i=0;i<n;i++) {
    xr[i]= xr1[i] * xr2[i] - xi1[i] * xi2[i];
    xi[i]= xr1[i] * xi2[i] + xi1[i] * xr2[i];
  }

  if(DEBUG) {
    printf(" DESPUES MULT\n");
    for(i=0;i<n;i++) 
      printf(" %d  %f %f \n",i,xr[i],xi[i]);
  }


  cfft(n,xr,xi,-1);

  if(DEBUG) {
    printf(" DESPUES CFFT INTEVERSA\n");
    for(i=0;i<n;i++) 
      printf(" %d  %g %g \n",i,xr[i],xi[i]);
  }

  n2=n/2;

  for(i=0;i<n2;i++) {
    xcorr[i]=-(float)(n2-i);
    xcorr[i+n2]=(float)i;
    fcorr[i]=xr[n2+i];
    fcorr[i+n2]=xr[i];
  }

  if(DEBUG) {
    printf(" FINAL FFTCORREL\n");
    for(i=0;i<n;i++) 
      printf(" %d  %f %f \n",i,xcorr[i],fcorr[i]);
  }

  free(xr1);
  free(xi1);
  free(xr2);
  free(xi2);
  free(xr);
  free(xi);
}


float offsetpix(int n, float *spec1, float *spec2) {

  int i;
  float *xcorr;
  float *fcorr;
  float *data1;
  float *data2;
  int nfft;
  float fl=0.1;

  float  fmax;
  int    imax;
  int    nside;
  int    nfit;
  int    ibegin;
  int    iend;
  double  *xfit;
  double  *yfit;
  double  *sigyfit;
  
  double parfit[4];
  int   iparfit[4];
  double **covarpar;
  double chisq;
  
  float first1,median1,third1;
  float first2,median2,third2;
  float stdmed;

  int isfinite;
  
  nfft=fft2power(n);

  if(DEBUG) printf(" Entro nuevo\n");

  if(DEBUG) printf(" nfft %d  n %d\n",nfft,n);

  xcorr=vector_f(nfft);
  fcorr=vector_f(nfft);
  data1=vector_f(nfft);
  data2=vector_f(nfft);
  
  Quartil(n,spec1,&first1,&median1,&third1);
  Quartil(n,spec2,&first2,&median2,&third2);
  if(DEBUG2) printf(" medias %f %f  %f %f %f %f \n",median1,median2,first1,first2,third1,third2);
  
  if(median1<=0)   median1=StMedia(n,spec1,&stdmed);
  if(median1<=0)   median1=1;
  if(median2<=0)   median2=StMedia(n,spec1,&stdmed);
  if(median2<=0)   median2=1;

  for(i=0;i<n;i++) {
    data1[i]=spec1[i]/median1;
    data2[i]=spec2[i]/median2;
    if(DEBUG2) printf(" ENTRADA offset %f %f\n",spec1[i],spec2[i]);
  } 
  for(i=n;i<nfft;i++) {
    data1[i]=0;
    data2[i]=0;
  }
  for(i=0;i<n;i++) {
    data1[i]=data1[i]*cosbell(i,n,fl);
    data2[i]=data2[i]*cosbell(i,n,fl); 
  } 

  fftcorrel(nfft,data1,data2,xcorr,fcorr);

  fmax=-1e30;
  imax=nfft/2;
  for(i=0;i<nfft;i++) {
    if(fabs(fcorr[i])>fmax) {
      fmax=fabs(fcorr[i]);
      imax=i;
    }
  }

  if(DEBUG2) printf(" imax %d  con fmax %f\n",imax,fmax);

  nside=7;
  ibegin=imax-nside;
  iend=imax+nside;
  if(ibegin<0) ibegin=0;
  if(iend>nfft-1) iend=nfft-1;
  if(DEBUG2) printf(" ibegin %d iend %d\n",ibegin,iend);
  nfit=iend-ibegin+1;
  if(DEBUG2) printf(" nfit %d \n",nfit);
  xfit=vector_d(nfit);
  yfit=vector_d(nfit);
  sigyfit=vector_d(nfit);
  covarpar=matrix_d(4,4);

  isfinite=0;
  for(i=ibegin;i<=iend;i++) {
    xfit[i-ibegin]=xcorr[i];
    yfit[i-ibegin]=fabs(fcorr[i]);
    sigyfit[i-ibegin]=fabs(sqrt(fabs(fcorr[i])))/10.;
    if(!(finite(yfit[i-ibegin])) || !(finite(sigyfit[i-ibegin])) ) {
      printf(" ERROR: offsetpix. Any number was not finite. Returning 0.\n");
      return(0.);
    }
  }
  iparfit[0]=1;iparfit[1]=1;iparfit[2]=1;iparfit[3]=1;

  if(DEBUG2) {
    for(i=0;i<nfit;i++) printf(" %d  xfit %f yfit %f sig %f\n",i,xfit[i],yfit[i],sigyfit[i]);
  }

  parfit[0]=0.;
  parfit[2]=xcorr[imax];
  parfit[3]=sqrt((xfit[nfit-1]-xfit[0])*(xfit[nfit-1]-xfit[0])/4./2./log(fmax/((yfit[0]+yfit[nfit-1])/2.)));
  parfit[1]=fmax*sqrt(2*M_PI)*parfit[3];

  
  if(PLOTDEBUG) {
    cpgswin(xcorr[0]-1,xcorr[nfft-1]+1,0.,1.3*fmax);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    cpgline(nfft,xcorr,fcorr);
  }

  if(DEBUG2) printf(" La ini: p0 %f p1 %f p2 %f p3 %f \n",parfit[0],parfit[1],parfit[2],parfit[3]);

  Mrq_d(xfit,yfit,sigyfit,nfit,parfit,iparfit,4,covarpar,&chisq,fungauss);
  
  if(DEBUG) printf(" Ajuste: p0 %f p1 %f p2 %f p3 %f \n",parfit[0],parfit[1],parfit[2],parfit[3]);
  if(DEBUG) printf(" Sigma : p0 %f p1 %f p2 %f p3 %f \n",sqrt(covarpar[0][0]),sqrt(covarpar[1][1]),sqrt(covarpar[2][2]),sqrt(covarpar[3][3]));

  if(PLOTDEBUG) {
    for(i=0;i<nfit;i++) yfit[i]=parfit[0]+parfit[1]*gaussian(xfit[i],parfit[2],parfit[3]);
    for(i=0;i<nfit;i++) printf(" AJUSTA %d  xfit %f yfit %f\n",i,xfit[i],yfit[i]);
    cpgsci(2);
    cpgline_d(nfit,xfit,yfit);
    cpgsci(1);
  }

  free(xcorr);
  free(fcorr);
  free(data1);
  free(data2);
  free(xfit);
  free(yfit);
  free(sigyfit);
  free_matrix_d(covarpar,4,4);
  return(parfit[2]);

}

void fungauss(double x,double *p,double *y,double *dyda,int n) {

  double fac,ex,arg;

  if(p[1]<0)  p[1]=-p[1];
  if(p[3]<0)  p[3]=-p[3];
/*   if(p[3]>2*nspacut) p[3]=2./platescale; */   /* No puede ser una gaussiana muy ancha */
  if(p[3]<0.1) p[3]=1;                     /* Ni mucho menos estrecha que un pixel */
/*   if(p[2]+2*p[3]<0) p[2]=nspacut; */
/*   if(p[2]-2*p[3]>nspacut*2) p[2]=nspacut; */

  arg=(x-p[2])/p[3];
  ex=gaussian(x,p[2],p[3]);
  fac=p[1]*2.0*arg*ex;
  *y=p[0]+p[1]*ex;
  dyda[0]=1;
  dyda[1]=ex;
  dyda[2]=fac/p[3];
  dyda[3]=fac*arg/p[3];

  if(DEBUG) printf(" x0 %g sig %g fac %g const %f \n",p[2],p[3],p[1],p[0]);
}


