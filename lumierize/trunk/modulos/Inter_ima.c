#include "modulos.h"
#define DEBUG 0

float intimapix_lin(float *ima, int nx, int ny, float xp, float yp, float nullval) {

  float cur,cul,cdr,cdl;  /* Contribution from upper/down right/left corners */
  float value;

  if(xp>0.5 && xp<nx-0.5 && yp>0.5 && yp<ny-0.5) {
    cdl=((int)xp+1-xp)*((int)yp+1-yp);
    cdr=(xp - (int)xp)*((int)yp+1-yp);
    cur=(xp - (int)xp)*(yp - (int)yp);
    cul=((int)xp+1-xp)*(yp - (int)yp);
    /* Los -1 son para pasar a pixels C */
    value =cdl*ima[(int)xp  -1 +((int)yp  -1)*nx];
    value+=cdr*ima[(int)xp+1-1 +((int)yp  -1)*nx];
    value+=cur*ima[(int)xp+1-1 +((int)yp+1-1)*nx];
    value+=cul*ima[(int)xp  -1 +((int)yp+1-1)*nx];
    
    return(value);
  }
  else {
    return(nullval);
  }
}


float intspecpix_lin(float *spec, int nx, float xp, float nullval) {

  float cr,cl;  /* Contribution from right/left pixels */
  float value;

  if(xp>0.5 && xp<nx+0.5 ) {
    cl=((int)xp+1-xp);
    cr=(xp - (int)xp);
    /* Los -1 son para pasar a pixels C */
    value =cl*spec[(int)xp  -1 ];
    value+=cr*spec[(int)xp+1-1 ];
    return(value);
  }
  else {
    return(nullval);
  }
}

float sumspecpix_lin(float *spec,int nx,float pix1,float pix2) {

  float sum;
  int ipix;
  int pixent1;
  int pixent2;
 
  pixent1=(int)(pix1+1.5);
  pixent2=(int)(pix2-0.5);
  if(pixent1<0) pixent1=0;
  if(pixent1>nx-1) pixent1=nx-1;
  /* Estos son los pixels que estan enteros. SUmo la contribucion de estos a pelo */
  sum=0;
  for(ipix=pixent1;ipix<=pixent2;ipix++)  sum+=spec[ipix];
  /* Si todo se integra dentro de un pixel: */
  if((pixent1-1)==(pixent2+1)) {
    if(pixent1!=0 && pixent1-0.5-pix1>0 && pixent2!=nx-1 && pix2-pixent2-0.5>0) {
      sum+=(pix2-pix1)*spec[pixent1-1];
    }
    return(sum);
  }
  
  /* Ahora cojo la fraccion izquierda */
  
  if(pixent1!=0 && pixent1-0.5-pix1>0) {
    sum+=(pixent1-0.5-pix1)*spec[pixent1-1];
  }

  /* Ahora la fraccion derecha */
  if(pixent2!=nx-1 && pix2-pixent2-0.5>0) {
    sum+=(pix2-pixent2-0.5)*spec[pixent2+1];
  }

  return(sum);
}

float sumspecpix_vista(float *spec,int nx,float pix1,float pix2) {
  
  float sum;
  int ipix;
  int pixent1;
  int pixent2;
  
  float a0;
  float a1;
  float a2;
  int pixinter;
  float x1,x2;
  /*   pix1=2.1; */
  /*   pix2=3.9; */
  
/*   int j; */
  
  pixent1=(int)(pix1+1.5);
  pixent2=(int)(pix2-0.5);

  if(DEBUG) printf(" pixent1 %d pixent2 %d\n",pixent1,pixent2);
  
  if(pixent1<0) pixent1=0;
  if(pixent1>nx-1) pixent1=nx-1;
  /* Estos son los pixels que estan enteros. SUmo la contribucion de estos a pelo */
  sum=0;
  for(ipix=pixent1;ipix<=pixent2;ipix++)  sum+=spec[ipix];
  if(DEBUG) printf(" Pixeles enteros : %f\n",sum);
  /* Si toda la integracion es en un solo pixel, la hago aqui: */
  if((pixent1-1)==(pixent2+1)) {
    pixinter=pixent1-1;
    if(pixinter>=1 && pixinter<=nx-2) {
      a2=+spec[pixinter-1]/2.+spec[pixinter+1]/2.-spec[pixinter];
      a1=-spec[pixinter-1]/2.+spec[pixinter+1]/2.;
      a0=+spec[pixinter]-a2/12.;
      x1=pix1-pixinter;
      x2=pix2-pixinter;
      sum+=a0*(x2-x1)+a1*(x2*x2-x1*x1)/2.+a2*(x2*x2*x2-x1*x1*x1)/3.;
      if(DEBUG) printf(" Suma en un solo pixel %f\n",sum);
    }
    else if(pixent1!=0 && pixent1-0.5-pix1>0 && pixent2!=nx-1 && pix2-pixent2-0.5>0) {
      /* Entonces interpolacion lineal, que esta en el borde */
      sum+=(pix2-pix1)*spec[pixent1-1];
    }
    return(sum);
  }


  /* Ahora cojo la fraccion izquierda */
  /* Lo hacemos como en VISTA, pero no como en REDUCEME, es un poco mas complicado 
     y no me apetece hacerlo ahora */
  
  pixinter=pixent1-1;
  if(pixinter>=1 && pixinter<=nx-2) {
    a2=+spec[pixinter-1]/2.+spec[pixinter+1]/2.-spec[pixinter];
    a1=-spec[pixinter-1]/2.+spec[pixinter+1]/2.;
    a0=+spec[pixinter]-a2/12.;
    x1=pix1-pixinter;
    x2=+0.5;
    sum+=a0*(x2-x1)+a1*(x2*x2-x1*x1)/2.+a2*(x2*x2*x2-x1*x1*x1)/3.;
    if(DEBUG) printf(" Fraccion izq %f\n",sum);
  }
  else if(pixent1!=0 && pixent1-0.5-pix1>0) {
    /* Entonces interpolacion lineal, que esta en el borde */
    sum+=(pixent1-0.5-pix1)*spec[pixent1-1];
  }
  
  /* Ahora la fraccion derecha */
  pixinter=pixent2+1;
  if(pixinter>=1 && pixinter<=nx-2) {
    a2=+spec[pixinter-1]/2.+spec[pixinter+1]/2.-spec[pixinter];
    a1=-spec[pixinter-1]/2.+spec[pixinter+1]/2.;
    a0=+spec[pixinter]-a2/12.;
    x1=-0.5;
    x2=pix2-pixinter;
    sum+=a0*(x2-x1)+a1*(x2*x2-x1*x1)/2.+a2*(x2*x2*x2-x1*x1*x1)/3.;
    if(DEBUG) printf(" Fraccion izq+der %f\n",sum);
  }
  else  if(pixent2!=nx-1 && pix2-pixent2-0.5>0) {
    sum+=(pix2-pixent2-0.5)*spec[pixent2+1];
  }
  return(sum);
}



