#include "modulos.h"

double Vol(double z, struct cosmo_param cosmo) {
  /*   Devuelve el resultado en parsecs al cubo */
  double c=299792.46; /*  En km/s */
  double Ok;
  double dm_dh;
  double dH3;
  double Vc=0; /* Volumen comovil */
  Ok=1-2*cosmo.q0;   /* Omega_k, que para Omega_lambda=0 vale 1-2q0 */
  dH3=(c/cosmo.H0*1.e6)*(c/cosmo.H0*1.e6)*(c/cosmo.H0*1.e6); /* La distancia de Hubble al cubo */
  dm_dh=2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z));
  if(cosmo.q0<0.5) {    /* Omega_k > 0 */
    Vc=(4*M_PI*dH3)/2/Ok*(dm_dh*sqrt(1+Ok*(dm_dh)*(dm_dh))-1/sqrt(Ok)*asinh(sqrt(Ok)*dm_dh));
  }
  else if(cosmo.q0>0.5) {   /* Omega_k < 0 */
    Vc=(4*M_PI*dH3)/2/Ok*(dm_dh*sqrt(1+Ok*(dm_dh)*(dm_dh))-1/sqrt(fabs(Ok))*asin(sqrt(fabs(Ok))*dm_dh));
  }
  else if(cosmo.q0==0.5) {   /* Omega_k = 0 */
    /*     printf(" Q0 = 0.5\n"); */
    Vc=(4*M_PI*dH3)/3*dm_dh*dm_dh*dm_dh;
  }
  return(Vc); 
}



double dVdz(double z, struct cosmo_param cosmo) {
  double c=299792.46; /*  En km/s */
  double Dang_DH,DH;  /* Distancia angular segun  Hogg astroph/9905116 */
  double E;        /* Funcion E. */
  /*   Esto es dV/dz, que no es igual a dV/dz=4pi dV/dz/dOmega, Omega: angulo so
       lido */
  Dang_DH=2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0
								      *cosmo.q0*(1+z)*(1+z)); /* La distancia angular */
  /*   Formula solo valida para Omega_lambda=0. Pero Dang=DM/(1+z). Eso siempre.
   */
  DH=(c/cosmo.H0*1.e6); /* La distancia de Hubble */
  E=sqrt(2*cosmo.q0*(1+z)*(1+z)*(1+z)+(1-2*cosmo.q0)*(1+z)*(1+z));
  /*   Esta formula solo es valida para Omega_lambda=0 */
  return(4*M_PI*(1+z)*(1+z)*Dang_DH*Dang_DH/E*DH*DH*DH);
  /*   El resultado en parcsecs al cuadrado */
}


double Lum(double z, double flux, struct cosmo_param cosmo) {
  /* Devuelve la luminosidad dado el flujo aparente */
  /* El flujo debe ir en watt/m2 y 
     la luminosidad la devuelve en Watt */
  double c=299792.46; /*  En km/s */
  double dlum;
  double L;
  double pi=3.1415926535897932384;
  /*   Calculo la distancia de luminosidad relativista. */
  dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z));
  /* Esta dlum esta dada en Megaparsec. 
     Para pasar a metros: 3.08567818589e22 metros  = 1 Mpc */
  L=4*pi*dlum*3.08567818589e22*dlum*3.08567818589e22*flux;
  return(L);
}

double dLumdflux(double z, struct cosmo_param cosmo) {
  /* Devuelve la derivada de la luminosidad con respecto al flujo aparente */
  /* la derivada de la luminosidad la devuelve en Watt / (watt/m2 ) = m2*/
  double c=299792.46; /*  En km/s */
  double dlum;
  double dLdf;
  double pi=3.1415926535897932384;
  /*   Calculo la distancia de luminosidad relativista. */
  dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z
))/(4*cosmo.q0*cosmo.q0*(1+z));
  /* Esta dlum esta dada en Megaparsec. 
     Para pasar a metros: 3.08567818589e22 metros  = 1 Mpc */
  dLdf=4*pi*dlum*3.08567818589e22*dlum*3.08567818589e22;
  return(dLdf);
}

double Flux(double z, double Lum, struct cosmo_param cosmo) {
  /* Devuelve el flujo aparente dada la luminosidad */
  /* La luminosidad debe ir in watt y
     el flujo lo devuelve en watt/m2 */

  double c=299792.46; /*  En km/s */
  double dlum;
  double flux;
  double pi=3.1415926535897932384;
  /*   Calculo la distancia de luminosidad relativista. */
  dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z));
  /* Esta dlum esta dada en Megaparsec. 
     Para pasar a metros: 3.08567818589e22 metros  = 1 Mpc */
  flux=Lum/4/pi/dlum/3.08567818589e22/dlum/3.08567818589e22;
  return(flux);
}


double mag(double z, double M, struct cosmo_param cosmo) {
  double c=299792.46; /*  En km/s */
  double dlum;
  double deltaM;
  dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z));
  deltaM=-5*log10(dlum/10e-6);  /* 10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
  return(M-deltaM);
}

double Mag(double z, double m, struct cosmo_param cosmo) {
  static double c=299792.46; /*  En km/s */
  double dlum;
  /* double deltaM; */
  /* dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z)); */
  dlum=c/cosmo.H0*2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0);
  /*   Esto es la distancia de luminosidad según Hogg. */
  /* deltaM=-5*log10(dlum/10e-6); */  /* 10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
  /* deltaM=-5*log10(dlum)-25.; */  /* solución más simple, pero no aumenta significativamente la vel. */
  /* return(deltaM+m); */
  return(-5.*log10(dlum)-25.+m);
}

void Mag_err(double z, double errz, double m, double errm, struct cosmo_param cosmo, double *Mag, double *errMag) {
  double c=299792.46; /*  En km/s */
  double dlum;
  double deltaM;
  double ddldz;
  dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z));
  ddldz=dlum/(1+z)+c/cosmo.H0*2*(1+z)*(2*cosmo.q0-2*cosmo.q0*(1-1*cosmo.q0)/sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z))+c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z))/(4*cosmo.q0*cosmo.q0*(1+z))*(4*cosmo.q0*cosmo.q0);
  /*   Esto es la distancia de luminosidad según Hogg. */
  deltaM=-5*log10(dlum/10e-6);  /* 10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
  *Mag=deltaM+m;
  *errMag=sqrt(errm*errm+errz*errz*(5/dlum/log(10))*(5/dlum/log(10))*ddldz*ddldz);
}


double Z_l(double flux, double L, struct cosmo_param cosmo) {
  double c=299792.46; /*  En km/s */
  double dlum;
  double dlum_dh; 
  double OM;
  double ztemp;
  double pi=3.1415926535897932384;
  dlum=sqrt(L/4./pi/flux)/3.08567818589e22;
  dlum_dh=dlum/c*cosmo.H0;
  OM=2*cosmo.q0;  /*        //Esto solo es valido para Omega_lambda=0 */
  ztemp=(dlum_dh*OM-2.+OM+(2.-OM)*sqrt(1.+2.*dlum_dh*OM*OM))/2.;
  return(ztemp);
}

double Z_m(double m, double M, struct cosmo_param cosmo) {
  double c=299792.46; /*  En km/s */
  double dlum;
  double dlum_dh; 
  double OM;
  double ztemp;
  double deltam;
  deltam=M-m;
  dlum=pow(10.,(-deltam+5.)/5.-6.);
  dlum_dh=dlum/c*cosmo.H0;
  OM=2*cosmo.q0;  /*        Esto solo es valido para Omega_lambda=0 */
  ztemp=(dlum_dh*OM-2.+OM+(2.-OM)*sqrt(1.+2.*dlum_dh*OM*OM))/2.;
  return(ztemp);
}


double D_co(double z, struct cosmo_param cosmo) {
  double c=299792.46; /*  En km/s */
  return(c/cosmo.H0*2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z)));

}

double D_lum(double z, struct cosmo_param cosmo) {
  return((1+z)*D_co(z,cosmo));
}

double D_ang(double z, struct cosmo_param cosmo) {
  double c=299792.46; /*  En km/s */
  return(c/cosmo.H0*1.e6*(2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z)*(1+z))));
}


