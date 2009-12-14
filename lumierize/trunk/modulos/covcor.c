#include "modulos.h"


void covcor(float *a,int k,int l,float *med,float *b,float *c)
{
  /* A: Input matrix of dimensions lxk.
     k: number of data
     l: number of variables
     med: vector with mean values. Dimension k
     B: covariance matrix        . Dimension lxl
     C: correlation matrix       . Dimension lxl
  */
  
  int i,n,m;
/*   //float sig[l],var[l]; */
  float *sig,*var;
  if((sig=malloc(l*sizeof(float)))==NULL) printf("I cannot dimension sig of %d elements \n",l);
  if((var=malloc(l*sizeof(float)))==NULL) printf("I cannot dimension var of %d elements \n",l);
  
  /*   for(i=0;i<30;i++) { */
  /*     printf(" %f %f %f %f %f\n",a[0+5*i],a[1+5*i],a[2+5*i],a[3+5*i],a[4+5*i]); */
  /*     //printf(" %f %f %f %f %f\n",y[0+5*i],y[1+5*i],y[2+5*i],y[3+5*i],y[4+5*i]); */
    /*   } */
    
  for(n=0;n<l;n++) {
    med[n]=0.;
    sig[n]=0.;
    var[n]=0.;
    for(i=0;i<k;i++) {
      med[n]=med[n]+a[n+l*i];
    }
    med[n]=med[n]/(float)k;
    for(i=0;i<k;i++) {
      var[n]=var[n]+(a[n+l*i]-med[n])*(a[n+l*i]-med[n]);
    }
    var[n]=var[n]/(float)k;
    sig[n]=sqrt(var[n]);
/*     //printf(" var %f med %f\n",var[n],med[n]); */
  }
  for(n=0;n<l;n++) {
    for(m=0;m<l;m++) {
      b[n+l*m]=0;
      if(!(m==n)) {
	for(i=0;i<k;i++) {
	  b[n+l*m]=b[n+l*m]+(a[n+l*i]-med[n])*(a[m+l*i]-med[m]);
	}
	b[n+l*m]=b[n+l*m]/(float)k;
	if(sig[n]<1.e-20 || sig[m]<1.e-20) c[n+l*m]=0.;
	else c[n+l*m]=b[n+l*m]/(sig[n]*sig[m]);
      }
      else {
	b[n+l*m]=var[n];
	c[n+l*m]=1.0;
      }
    }
  }
}


/*
  
  SUBROUTINE COVCOR(A,K,L,MED,B,C)
  C A Matriz de K datos con L columnas
  C K Numero de datos (filas)
  C L Numero de variables (columnas)
  C B Matriz LxL de covarianza
  C C Matriz LxL de correlacion
  INTEGER DIM,I,J,K,L,NPOIN,N,M
  REAL A(L,K),B(L,L),C(L,L)
  REAL MED(L),SIG(L),VAR(L)
  
  DO N=1,L
  MED(N)=0.
  SIG(N)=0.
  VAR(N)=0.
  DO I=1,K   
  MED(N)=MED(N)+A(N,I)
  ENDDO
  MED(N)=MED(N)/FLOAT(K)
  DO I=1,K
  VAR(N)=VAR(N)+(A(N,I)-MED(N))**2.
  ENDDO
  VAR(N)=VAR(N)/FLOAT(K)
  SIG(N)=SQRT(VAR(N))
  ENDDO
  
  DO N=1,L
  DO M=1,L
  B(N,M)=0.
  IF (M.NE.N) THEN
  DO I=1,K
  B(N,M)=B(N,M)+(A(N,I)-MED(N))*(A(M,I)-MED(M))
  ENDDO
  B(N,M)=B(N,M)/FLOAT(K)
  IF (SIG(N).LT.1.E-20.OR.SIG(M).LT.1.E-20) THEN
  C(N,M)=0.0
  ELSE
  C(N,M)=B(N,M)/(SIG(N)*SIG(M))
  ENDIF
  ELSE
  B(N,M)=VAR(N)
  C(N,M)=1.0
  ENDIF
  ENDDO
  ENDDO
  
  
  RETURN
  END
  
*/
