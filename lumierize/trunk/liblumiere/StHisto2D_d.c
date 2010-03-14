#include "sthisto.h"

int **StHisto2D_d(int n, double *x, double *y, int nbinx, double *xmin, double *xmax, int nbiny, double *ymin, double *ymax)

{
  double deltx, delty;
  int i,jx,jy,ex,ey,**h;

  /* Creamos h */
  h=matrix_i(nbinx,nbiny);
    
  /* Inicializamos los valores de h a cero */
  for(jx=0; jx<nbinx; jx++) 
  {
     for(jy=0; jy<nbiny; jy++)
     {
         h[jx][jy]=0;
     }
  } 

  /* Si no tenemos definidos los extremos los tomamos de los datos */
  if(*xmin == 0 && *xmax == 0)  MinMax_d(n,x,xmin,xmax);
  if(*ymin == 0 && *ymax == 0)  MinMax_d(n,y,ymin,ymax);
  deltx = (*xmax - *xmin)/(nbinx);
  delty = (*ymax - *ymin)/(nbiny);

  /* Calculamos el histograma */
  for(i=0; i<n; i++) 
  {
    if (x[i] >= *xmin && x[i] < *xmax && y[i] >= *ymin && y[i] < *ymax) 
    {
      ex=(int)((x[i] - *xmin)/deltx);
      ey=(int)((y[i] - *ymin)/delty);
      if(ex>=nbinx || ex<=-1 || ey>=nbiny || ey<=-1) 
      {
	printf(" ERROR: StHisto2, ex=%d ey=%d.\n Exiting\n",ex,ey);
	exit(1);
      }
      h[ex][ey]++;
    }
  }
  return(h);
}

