#include "modulos.h"

#define DEBUG 0

int  VVmax_M(int n,double *magn,double *z, double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf) {
  
  int i,j;
  int *ngalbin;
  double M;
  double zmax;
  
  ngalbin=vector_i(lf->nbin);
  
  if(DEBUG) printf(" Entro en vvamx_M\n");

  for(i=0;i<lf->nbin;i++) {
    lf->lf[i]=0.;
    ngalbin[i]=0;
    
    for(j=0;j<n;j++) {
      M=Mag(z[j],magn[j],cosmo);
      if(DEBUG) printf(" L %g   %g - %g\n",M,lf->magni[i],lf->magni[i+1]);
      if(M>lf->magni[i] && M<lf->magni[i+1]) {
	zmax=Z_m(mlim,M,cosmo);      /* Esto no es asi, tengo que comprobarlo */
	if(DEBUG) printf(" zmax %f \n",zmax);
	if(zmax>zup && zup!=0) zmax=zup;
	if(zmax>zlow) {
	  lf->lf[i]+=1./(Vol(zmax,cosmo)-Vol(zlow,cosmo));
	  ngalbin[i]++;
	}
      }
    }
    
    lf->lf[i]=lf->lf[i]/strrad*4*M_PI/(lf->magni[i+1]-lf->magni[i])*1.e18; /* El 1.e18 es para pasar de parsecs3 a megaparsecs3. El dl/vvmax.lum[i] es diferencial de L ( ya que dl es en realidad d(logL)) */
    lf->errlf[i]=lf->lf[i]/sqrt(ngalbin[i]); /* Esto es como suponer que la funcion de lum es proporcional a una variable poissonian que es el numero de de galaxias en ese bin*/

    if(DEBUG) printf(" bin %d L %f - %f lf %g (%f) +/- %g (%f)  ngal %d\n",i,lf->magni[i],lf->magni[i+1],lf->lf[i],log10(lf->lf[i]),lf->errlf[i],lf->errlf[i]/lf->lf[i]/log(10),ngalbin[i]);


    /* Paso a logaritmos */
    lf->errlf[i]=lf->errlf[i]/lf->lf[i];
    lf->lf[i]=log(lf->lf[i]);

    if(ngalbin[i]==0) {
      lf->lf[i]=0;
      lf->errlf[i]=0;
    }
    
  }
  return(0);
}

