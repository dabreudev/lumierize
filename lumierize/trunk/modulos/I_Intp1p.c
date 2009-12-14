#include "modulos.h"


int I_Intp1p(int nx, int ny, float *ima, float pixx, float pixy, float *value)

{ 

  float f[4],dx,dy;
  int px,py;

/*  Test variados */
/*  ------------- */
/*  Esta era la condicion anterior */
/*   if(pixx < 1.0 || pixx > (float)nx || pixy < 1.0 || pixy > (float)ny) {  */
  if(pixx < 0.0 || pixx > (float)nx +1.0  || pixy < 0.0 || pixy > (float)ny+1.0) { 
/* He cambiado la condicion a que sea menor que 0.0, aunque yo mismo lo habia 
   puesto a menor de 0.5, pero alguna vez petaba */
/*     fprintf(stderr,"I_Intp1p: WARNING. Me pides extrapolar.\n"); */
/*     fprintf(stderr,"          Tamaño de la imagen %dx%d\n",nx,ny); */
/*     fprintf(stderr,"          Pides valor en pixel=%f,%f\n",pixx,pixy); */
    return(0);
    } 


/*   pixx y pixy viene en coordenadas pixel, que no es lo mismo que  */
/*    las componentes de la matriz. */

  px=(int)pixx;
  py=(int)pixy;
/*   printf(" INTE %f %f PIX  %d %d DIM %d %d\n",pixx,pixy,px,py,nx,ny); */

  if(px == nx) px--;
  if(py == ny) py--;
  
  if(px == 0) px=1;
  if(py == 0) py=1;


  f[0]=ima[(px-1)+(py-1)*nx];
  f[1]=ima[(px  )+(py-1)*nx];
  f[2]=ima[(px-1)+(py  )*nx];
  f[3]=ima[(px  )+(py  )*nx];
/*   printf(" VAlues %f %f %f %f\n",ima[(px-1)+(py-1)*nx],ima[(px  )+(py-1)*nx],ima[(px-1)+(py  )*nx],ima[(px  )+(py  )*nx]); */

  dx=pixx-(float)px;
  dy=pixy-(float)py;

  *value = (1.0-dx)*(1.0-dy)*f[0]+dx*(1.0-dy)*f[1]+dy*(1.0-dx)*f[2]+dx*dy*f[3];
  return(1);
}
