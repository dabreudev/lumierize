#include "modulos.h"

#define DEBUG 0

/**
 * This function will compute the LF using VVmax method in magnitudes.
 * @param mag_sel   A vector with the apparent magnitudes in the band 
 *                  where selection was performed.
 * @param mag_cal   A vector with the apparent magnitudes in the band 
 *                  where we want to compute the LF.
 * @param strrad    The area covered by the survey in radians
 */
int  VVmax_M(int n,double *mag_sel, double *mag_cal,double *z, double mlim_sel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf) 
{
  
  int i,j;
  int *ngalbin;
  double M_cal, M_sel;
  double zmax;
  
  ngalbin=vector_i(lf->nbin);
  
  if(DEBUG) printf(" Entro en vvamx_M\n");

  for(i=0;i<lf->nbin;i++) 
  {
    lf->lf[i]=0.;
    ngalbin[i]=0;
    
    for(j=0;j<n;j++)
    {
      M_cal=Mag(z[j],mag_cal[j],cosmo);
      M_sel=Mag(z[j],mag_sel[j],cosmo);
      if(DEBUG) printf(" z  %g M_cal %g M_sel %g  %g - %g\n",
                       z[j],M_cal,M_sel,lf->magni[i],lf->magni[i+1]);
      if(M_cal>lf->magni[i] && M_cal<lf->magni[i+1]) 
      {
	zmax=Z_m(mlim_sel,M_sel,cosmo);
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


    /* Paso a logaritmos */
    lf->errlnlf[i]=lf->errlf[i]/lf->lf[i];
    lf->lnlf[i]=log(lf->lf[i]);


    if(ngalbin[i]==0) 
    {
      lf->lnlf[i]=-1/0.;
      lf->errlnlf[i]=-1/0.;
    }
    
    if(DEBUG) printf(" bin %d M %g - %g lf %g (%g) +/- %g (%g)  ngal %d\n",i,lf->magni[i],lf->magni[i+1],lf->lf[i],lf->lnlf[i],lf->errlf[i],lf->errlnlf[i],ngalbin[i]);

  }
  return(0);
}

