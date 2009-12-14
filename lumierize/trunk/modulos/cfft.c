#include "modulos.h"
#define DEBUG 0

int iBitRFFT(int j, int nu);

int fft2power(int n0) {
  int k;
  int n;

  k=0;
  n=1;
  while(n<n0) {
    k++;
    n=n*2;
  }
  return(n);
}
  

void cfft(int n, float *xr, float *xi, int imode) {
  
  int nu,n0;
  int n2,nu1,k,l,i,p,k1,k1n2;
  float xr0,xi0,fn;
  double darg,cc,ss,treal,timag;

  if(DEBUG) {
    printf(" DENTRO CFFT\n");
    for(i=0;i<n;i++) {
      printf(" %d %f %f \n",i,xr[i],xi[i]);
    }
  }


  nu=0;
  n0=n;
  while (n0>2) {
    if((n0-2*floor(n0/2.))!=0) {
      printf(" cfft: Error: N is not a power of 2. Exiting\n");
      exit(1);
    }
    n0=n0/2;
    nu++;
  }
  
  nu++;

  /*   To compute the inverse FFT remember that, given a function f and its Fourier 
       transform F: F = FFT(f) and  f = 1/N [FFT(f*)]*, where the asterisks mean 
       complex conjugate. */
  if(imode==-1) {
    for(i=0;i<n;i++) {
      xi[i]=-xi[i];
    }
  }
  /* Initialization */

  n2=n/2;
  nu1=nu-1;
  k=0;

  for(l=0;l<nu;l++) {
    while(k<n) {
      for(i=1;i<=n2;i++) {
	p=iBitRFFT(k/(int)(pow(2,nu1)),nu);
	darg=2*M_PI*(double)p/((double)n);
	cc=cos(darg);
	ss=sin(darg);
	k1=k+1;
	k1n2=k1+n2;
	treal=(double)(xr[k1n2-1])*cc+(double)(xi[k1n2-1])*ss;
	timag=(double)(xi[k1n2-1])*cc-(double)(xr[k1n2-1])*ss;
	xr[k1n2-1]=xr[k1-1]-treal;
	xi[k1n2-1]=xi[k1-1]-timag;
	xr[k1-1]  =xr[k1-1]+treal;
	xi[k1-1]  =xi[k1-1]+timag;
	k=k+1;
      }
      k=k+n2;
    }
    k=0;
    nu1=nu1-1;
    n2=n2/2;
  }

  if(DEBUG) {
    printf(" PRIMER PASO CFFT\n");
    for(i=0;i<n;i++) {
      printf(" %d %f %f \n",i,xr[i],xi[i]);
    }
  }


  /* unscrambling the FFT */
  for(k=1;k<=n;k++) {
    i=iBitRFFT(k-1,nu)+1;
    if(i>k) {
      xr0=xr[k-1];
      xi0=xi[k-1];
      xr[k-1]=xr[i-1];
      xi[k-1]=xi[i-1];
      xr[i-1]=xr0;
      xi[i-1]=xi0;
    }
  }

  if(DEBUG) {
    printf(" PASO FINAL CFFT\n");
    for(i=0;i<n;i++) {
      printf(" %d %f %f \n",i,xr[i],xi[i]);
    }
  }


  if(imode==-1) {
    for(i=0;i<n;i++)   xi[i]=-xi[i];
    fn=1./n;
    for(i=0;i<n;i++) {
      xr[i]=xr[i]*fn;
      xi[i]=xi[i]*fn;
    }
  }
}


int iBitRFFT(int j, int nu) {
  
  int i,j1,j2;
  int value;

  j1=j;
  value=0;
  for(i=1;i<=nu;i++) {
    j2=j1/2;
    value=value*2+(j1-2*j2);
    j1=j2;
  }
  return(value);
}


float cosbell(int i, int n, float fl) {
  float nl;
  
  if(fl<0 || fl> 0.5 ) {
    printf(" Cosbell Error: fl < 0 or fl > 0.5 (fl = %f). Exiting\n",fl);
    exit(1);
  }
  if(i<0 || i>n ) {
    printf(" Cosbell Error: i < 0 or i > n (i = %d, n = %d). Exiting\n",i,n);
    exit(1);
  }

  nl=(int)(fl*n);

  if(i<=nl-1)       return(0.5*(1.-cos(M_PI*(i+1)/(float)nl)));
  else if(i>=n-nl)  return(0.5*(1.-cos(M_PI*(n-i-1)/(float)nl)));
  else              return(1.); 

}
