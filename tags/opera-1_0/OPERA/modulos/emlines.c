#include "modulos.h"


double Flux_ew_mag(double ew, double mag, char photband[51], float gamma, float delta, float Kcoc) {

  double fluxband;
  double fluxcont;
  double fluxline;
  double zeropoint;

  if(!strcmp(photband,"U")) zeropoint=28.36;
  else if(!strcmp(photband,"B")) zeropoint=27.97;
  else if(!strcmp(photband,"V")) zeropoint=28.52;
  else if(!strcmp(photband,"R")) zeropoint=29.39;
  else if(!strcmp(photband,"I")) zeropoint=29.66;
  else if(!strcmp(photband,"J")) zeropoint=31.25;
  else if(!strcmp(photband,"H")) zeropoint=32.35;
  else if(!strcmp(photband,"K")) zeropoint=33.50;
  else {
    printf(" ERROR: Photometric band %-s not defined\n",photband);
    printf(" Photometric bands defined are: U B V R I J H K. \n EXITING\n");
    exit(1);
  }

  fluxband=pow(10.,-0.4*(mag+zeropoint));
  fluxcont=fluxband/(Kcoc+gamma/delta*ew);
  fluxline=ew*fluxcont;
  
  return(fluxline);
  
}



void   Flux_ew_mag_err(double ew, double errew, double mag, double errmag, char photband[51], float gamma, float delta,float Kcoc, double *fluxline, double *errfluxline) {
  
  
  double fluxband,errfluxband;
  double fluxcont,errfluxcont;
  /*   float fluxline,errfluxline; */
  double zeropoint;

  if(!strcmp(photband,"U")) zeropoint=28.36;
  else if(!strcmp(photband,"B")) zeropoint=27.97;
  else if(!strcmp(photband,"V")) zeropoint=28.52;
  else if(!strcmp(photband,"R")) zeropoint=29.39;
  else if(!strcmp(photband,"I")) zeropoint=29.66;
  else if(!strcmp(photband,"J")) zeropoint=31.25;
  else if(!strcmp(photband,"H")) zeropoint=32.35;
  else if(!strcmp(photband,"K")) zeropoint=33.50;
  else {
    printf(" ERROR: Photometric band %-s not defined\n",photband);
    printf(" Photometric bands defined are: U B V R I J H K. \n EXITING\n");
    exit(1);
  }

  fluxband=pow(10.,-0.4*(mag+zeropoint));
  errfluxband=0.4*log(10)*errmag*fluxband;
  fluxcont=fluxband/(Kcoc+gamma/delta*ew);
  errfluxcont=sqrt( (errfluxband/(Kcoc+gamma/delta*ew))*(errfluxband/(Kcoc+gamma/delta*ew)) + (errew*gamma/delta*fluxband/(Kcoc+gamma/delta*ew)/(Kcoc+gamma/delta*ew))*(errew*gamma/delta*fluxband/(Kcoc+gamma/delta*ew)/(Kcoc+gamma/delta*ew)));
  *fluxline=ew*fluxcont;
  *errfluxline=sqrt( (errew*fluxcont)*(errew*fluxcont) + (ew*errfluxcont)*(ew*errfluxcont));

}


double mag_ew_flux(double ew, double flux, char photband[51], float gamma, float delta, float Kcoc) {

  double fluxband;
  double fluxcont;
  double fluxline;
  double zeropoint;
  double magap;
  

  if(!strcmp(photband,"U")) zeropoint=28.36;
  else if(!strcmp(photband,"B")) zeropoint=27.97;
  else if(!strcmp(photband,"V")) zeropoint=28.52;
  else if(!strcmp(photband,"R")) zeropoint=29.39;
  else if(!strcmp(photband,"I")) zeropoint=29.66;
  else if(!strcmp(photband,"J")) zeropoint=31.25;
  else if(!strcmp(photband,"H")) zeropoint=32.35;
  else if(!strcmp(photband,"K")) zeropoint=33.50;
  else {
    printf(" ERROR: Photometric band %-s not defined\n",photband);
    printf(" Photometric bands defined are: U B V R I J H K. \n EXITING\n");
    exit(1);
  }
 
  fluxcont=flux/ew;
  fluxband=fluxcont*(Kcoc+gamma/delta*ew);
  magap=-2.5*log10(fluxband)-zeropoint;

  return(magap);

}

