#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vvmax.h"
#include "alloc.h"
#include "cosmology.h"

#define DEBUG 0

int  VVmax_L(int n,double *flux_sel, double *flux_cal,double *z, double flim_sel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf) {
  
  int i,j;
  int *ngalbin;
  double L_cal;
  double L_sel;
  double zmax;
  
  ngalbin=vector_i(lf->nbin);
  
  if(DEBUG) printf(" Entro en vvamx_l\n");

  for(i=0;i<lf->nbin;i++) {
    lf->lnlf[i]=0.;
    ngalbin[i]=0;
    
    for(j=0;j<n;j++) {
      L_cal=Lum(z[j],flux_cal[j],cosmo);
      L_sel=Lum(z[j],flux_sel[j],cosmo);
      if(DEBUG) printf(" L %g   %g - %g\n",log(L_cal),lf->lumi[i],lf->lumi[i+1]);
      if(log(L_cal)>lf->lumi[i] && log(L_cal)<lf->lumi[i+1]) {
	zmax=Z_l(flim_sel,L_sel,cosmo);    
	if(DEBUG) printf(" zmax %f \n",zmax);
	if(zmax>zup && zup!=0) zmax=zup;
	if(zmax>zlow) {
	  lf->lnlf[i]+=1./(Vol(zmax,cosmo)-Vol(zlow,cosmo));
	  ngalbin[i]++;
	}
      }
    }
    
    if(DEBUG) printf(" lf %g fac %g\n",lf->lnlf[i],1./((exp(lf->lumi[i])+exp(lf->lumi[i+1]))/2.));
    lf->lnlf[i]=lf->lnlf[i]/strrad*4*M_PI/(exp(lf->lumi[i+1])-exp(lf->lumi[i]))*1.e18; /* El 1.e18 es para pasar de parsecs3 a megaparsecs3. El dl/vvmax.lum[i] es diferencial de L ( ya que dl es en realidad d(logL)) */
    lf->errlnlf[i]=lf->lnlf[i]/sqrt(ngalbin[i]); /* Esto es como suponer que la funcion de lum es proporcional a una variable poissonian que es el numero de de galaxias en ese bin*/

    if(DEBUG) printf(" bin %d L %f - %f lf %g (%f) +/- %g (%f)  ngal %d\n",i,lf->lumi[i]/log(10),lf->lumi[i+1]/log(10),lf->lnlf[i],log10(lf->lnlf[i]),lf->errlnlf[i],lf->errlnlf[i]/lf->lnlf[i]/log(10),ngalbin[i]);


    /* Paso a logaritmos */
    lf->errlnlf[i]=lf->errlnlf[i]/lf->lnlf[i];
    lf->lnlf[i]=log(lf->lnlf[i]);

    if(ngalbin[i]==0) 
    {
      lf->lnlf[i]=-1/0.;
      lf->errlnlf[i]=-1/0.;
    }
    
  }
  return(0);
}

