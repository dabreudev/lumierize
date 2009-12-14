#include "modulos.h"


 
void PlateCtes_3(int ndat, float *xdat, float *ydat, float *alfa,float *delta,float alfac,float deltac,struct plt_cte *constantes) {
  /*
    Le meto las coordenadas medidas por mi(xdat,ydat) dimensionado como ndat.DIVIDIDAS ya por la escala de placa
    Le meto las coordenadas ecuatoriales de esas estrellas(alfa,delta) en RADIANES
    Le meto las coordenadas cebtrales de la placa (alfac,deltac) EN RADIANES)
    En la estructura plt_cte salen las constantes de placa:
    psi=bpx+cpy+ap
    eta=bex+cey+ae
  */

  float sx=0,sy=0,sx2=0,sxy=0,sy2=0;
  float sxp=0,syp=0,sp=0,sye=0,sxe=0,se=0;
  
  double psi,eta;
  int i;
  float a[9];
  float c[3];
  float x[3];


/*   float app,dpp; */

  for(i=0;i<ndat;i++) {
    Ecu2Plac((double)alfa[i],(double)delta[i],(double)alfac,(double)deltac,&psi,&eta);
/*         Plac2Ecu(psi,eta,alfac,deltac,&app,&dpp); */
/*     printf("Estrella %d a,d %e %e   psi, erta %f %f  app ,dpp %f %f\n",i,alfa[i],delta[i],psi,eta,app,dpp);  */

/*     //printf("Estrella %d a,d %e %e  ac,dc %e %e psi, erta %f %f \n",i,alfa[i],delta[i],alfac,deltac,psi,eta);  */
    sx+=xdat[i];
    sy+=ydat[i];
    sx2+=xdat[i]*xdat[i];
    sxy+=xdat[i]*ydat[i];
    sy2+=ydat[i]*ydat[i];
    sxp+=xdat[i]*psi;
    syp+=ydat[i]*psi;
    sp+=psi;
    sye+=ydat[i]*eta;
    sxe+=xdat[i]*eta;
    se+=eta;
  }
  a[0]=sx2;
  a[1]=sxy;
  a[2]=sx;
  a[3]=sxy;
  a[4]=sy2;
  a[5]=sy;
  a[6]=sx;
  a[7]=sy;
  a[8]=ndat;
  c[0]=sxp;
  c[1]=syp;
  c[2]=sp;
/*   //printf("a %g %g %g \n",a[0],a[4],a[9]); */
/*   //printf("c %g %g %g \n",c[0],c[1],c[2]); */
/*   //printf("DSAIOJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n"); */
  
  
  SELGauss(3,a,c,x);
/*   //printf("x x x %f %f %f \n",x[0],x[1],x[2]); */
  constantes->bp=x[0];
  constantes->cp=x[1];
  constantes->ap=x[2];
  

  c[0]=sxe;
  c[1]=sye;
  c[2]=se;
  
  
  
  SELGauss(3,a,c,x);
  
  
  constantes->be=x[0];
  constantes->ce=x[1];
  constantes->ae=x[2];

  
  
/*   //  constantes->amdx1=constantes->a*180*3600*1000/pi/constantes->xpixelsz; */
/*   // Esto hay que testearlo sobre todo por el PP03 */
}
