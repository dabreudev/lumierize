#include "modulos.h"

#define DEBUG 0

int  VVmax_L(int n,double *flux,double *z, double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf) {
  
  int i,j;
  int *ngalbin;
  double L;
  double zmax;
  
  ngalbin=vector_i(lf->nbin);
  
  if(DEBUG) printf(" Entro en vvamx_l\n");

  for(i=0;i<lf->nbin;i++) {
    lf->lf[i]=0.;
    ngalbin[i]=0;
    
    for(j=0;j<n;j++) {
      L=Lum(z[j],flux[j],cosmo);
      if(DEBUG) printf(" L %g   %g - %g\n",log(L),lf->lumi[i],lf->lumi[i+1]);
      if(log(L)>lf->lumi[i] && log(L)<lf->lumi[i+1]) {
	zmax=Z_l(flim,L,cosmo);      /* Esto no es asi, tengo que comprobarlo */
	if(DEBUG) printf(" zmax %f \n",zmax);
	if(zmax>zup && zup!=0) zmax=zup;
	if(zmax>zlow) {
	  lf->lf[i]+=1./(Vol(zmax,cosmo)-Vol(zlow,cosmo));
	  ngalbin[i]++;
	}
      }
    }
    
    if(DEBUG) printf(" lf %g fac %g\n",lf->lf[i],1./((exp(lf->lumi[i])+exp(lf->lumi[i+1]))/2.));
    lf->lf[i]=lf->lf[i]/strrad*4*M_PI/(exp(lf->lumi[i+1])-exp(lf->lumi[i]))*1.e18; /* El 1.e18 es para pasar de parsecs3 a megaparsecs3. El dl/vvmax.lum[i] es diferencial de L ( ya que dl es en realidad d(logL)) */
    lf->errlf[i]=lf->lf[i]/sqrt(ngalbin[i]); /* Esto es como suponer que la funcion de lum es proporcional a una variable poissonian que es el numero de de galaxias en ese bin*/

    if(DEBUG) printf(" bin %d L %f - %f lf %g (%f) +/- %g (%f)  ngal %d\n",i,lf->lumi[i]/log(10),lf->lumi[i+1]/log(10),lf->lf[i],log10(lf->lf[i]),lf->errlf[i],lf->errlf[i]/lf->lf[i]/log(10),ngalbin[i]);


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

