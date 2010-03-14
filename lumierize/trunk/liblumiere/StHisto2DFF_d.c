#include "sthisto.h"

int **StHisto2DFF_d(int n, double *x, double *y, int nbinx, double *xmin, double *xmax, double yff)

{
  double deltx;
  int i,jx,jy,ex,ey,**h;

  /* Creamos h */
  h=matrix_i(nbinx,2);
    
  /* Inicializamos los valores de h a cero */
  for(jx=0; jx<nbinx; jx++) 
  {
     for(jy=0; jy<2; jy++)
     {
         h[jx][jy]=0;
     }
  } 

  /* Si no tenemos definidos los extremos los tomamos de los datos */
  if(*xmin == 0 && *xmax == 0)  MinMax_d(n,x,xmin,xmax);
  deltx = (*xmax - *xmin)/(nbinx);

  /* Calculamos el histograma */
  for(i=0; i<n; i++) 
  {
    if (x[i] >= *xmin && x[i] < *xmax) 
    {
      ex=(int)((x[i] - *xmin)/deltx);
      if(ex>=nbinx || ex<=-1) 
      {
	printf(" ERROR: StHisto2, ex=%d.\n Exiting\n",ex);
	exit(1);
      }
      if (y[i] >= yff)
           ey = 1;
      else
           ey = 0;
      h[ex][ey]++;
    }
  }
  return(h);
}

