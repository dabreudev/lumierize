#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

#define LSPEED 299792.46 /* LSPEED: light speed in km/sec */

#define ZMAX_CACHED 30.      /* For cached values in cosmo_init */
#define NUMZ_CACHED 1000000 

#define DEBUG1 0
#define DEBUG2 0
#define DEBUG3 0

/* If you want to use cosmoloy functions you must call firts cosmo_init */

void cosmo_init(struct cosmo_param *cosmo, double H0, double OM, double OL)
{
  double zstep;
  double dc_int=0.;
  double DC,DM,DA;
  int i,j;
  double z,E;
  double min=0,max=0;
  double log10dlum=0.,log10dlum_step;

  zstep=ZMAX_CACHED/NUMZ_CACHED;

  /* free(cosmo->cached.z);
  free(cosmo->cached.log10dlum);
  free(cosmo->cached.dVdz);
  free(cosmo->cached.invlog10dlum);
  free(cosmo->cached.invz);
  free(cosmo->cached.DM); */

  cosmo->H0 = H0;
  cosmo->OM = OM;
  cosmo->OL = OL;
  cosmo->Ok = 1. - OM - OL;
  cosmo->DH = LSPEED / H0;
  cosmo->cached.N = NUMZ_CACHED;
  cosmo->cached.z            = malloc(NUMZ_CACHED*sizeof(double));
  cosmo->cached.log10dlum    = malloc(NUMZ_CACHED*sizeof(double));
  cosmo->cached.dVdz         = malloc(NUMZ_CACHED*sizeof(double));
  cosmo->cached.invlog10dlum = malloc(NUMZ_CACHED*sizeof(double));
  cosmo->cached.invz         = malloc(NUMZ_CACHED*sizeof(double));
  cosmo->cached.DM           = malloc(NUMZ_CACHED*sizeof(double));

  /* loop for the integral */
  for(i=0;i<=NUMZ_CACHED;i++)
  { /* first step is zstep */
    z=(i+1)*zstep;
    E= Efunction(z,*cosmo);
    dc_int+= cosmo->DH / E * zstep; /* dz = zstep */
    cosmo->cached.z[i] = z;
    DC = dc_int;
    if(cosmo->Ok == 0)
    {
      DM = DC;
    }
    else if(cosmo->Ok > 0)
    {
      DM = sinh(sqrt(cosmo->Ok)*DC/cosmo->DH)*cosmo->DH / sqrt(cosmo->Ok);
    }
    else
    {
      DM = sin(sqrt(fabs(cosmo->Ok))*DC/cosmo->DH)*cosmo->DH / sqrt(fabs(cosmo->Ok));
    }

    DA = DM / (1.+z);
    cosmo->cached.log10dlum[i] = log10((1+z)*DM);
    cosmo->cached.dVdz[i] = cosmo->DH*1e6*1e6*1e6*(1.+z)*(1+z)*DA*DA*4.*M_PI / E; /* dVdz in pc**3 */
    cosmo->cached.DM[i] = DM;
    if(DEBUG2) printf("dVdz %g z %g\n",cosmo->cached.dVdz[i],cosmo->cached.z[i]);

  } /* end of loop for the integral */
  MinMax_d(NUMZ_CACHED, cosmo->cached.log10dlum, &min, &max);	
  cosmo->cached.maxlog10dlum=max;
  cosmo->cached.minlog10dlum=min;
  log10dlum_step = (max-min) / NUMZ_CACHED;

  j=0;
  log10dlum = min-log10dlum_step;
  /* solo para log10dlum monotonamente creciente */  
  for(i=0;i<NUMZ_CACHED;i++)
  {
    log10dlum+=log10dlum_step;
    while(j<=NUMZ_CACHED)
    {
      if(DEBUG3)
      {
        printf("Inside the while j %d\n",j);
        printf("log10dlum %g cached.log10dlum %g\n",log10dlum,cosmo->cached.log10dlum[j]);
      }
      if(log10dlum > cosmo->cached.log10dlum[j]) j+=1;
      else
      {
        if(DEBUG3) printf("entering else j %d\n",j);
        cosmo->cached.invlog10dlum[i] = log10dlum;
        cosmo->cached.invz[i] = cosmo->cached.z[j];
        //j+=1;
        break;
      }
    } /* end of loop searching for a value */
  }
  if(DEBUG1)
  {
    printf("cosmo: H0 %g OM %g OL %g Ok %g DH %g\n",cosmo->H0,cosmo->OM,cosmo->OL,cosmo->Ok,cosmo->DH);
    printf("cosmo.cached: N %d maxlog10dlum %g minlog10dlum %g\n",cosmo->cached.N,cosmo->cached.maxlog10dlum,cosmo->cached.minlog10dlum);
    printf("log10dlum_step %g\n",log10dlum_step);
  }

  if(DEBUG2)
  {
    for(i=0;i<NUMZ_CACHED;i++)
    {
      printf("cached.z %g\t",cosmo->cached.z[i]);
      printf("cached.log10dlum %g\n",cosmo->cached.log10dlum[i]);
      printf("cached.invlog10dlum %g\t",cosmo->cached.invlog10dlum[i]);
      printf("cached.invz %g\n",cosmo->cached.invz[i]);
    }
  }
  if(DEBUG1) printf("cosmo_init done.\n");
}

void cosmo_free(struct cosmo_param *cosmo)
{
  /* if you want to change your cosmology, first you must do a free */
  free(cosmo->cached.z);
  free(cosmo->cached.log10dlum);
  free(cosmo->cached.dVdz);
  free(cosmo->cached.invlog10dlum);
  free(cosmo->cached.invz);
  free(cosmo->cached.DM);
}

double Efunction(double z, struct cosmo_param cosmo)
{
  return(sqrt(cosmo.OM*(1.+z)*(1.+z)*(1.+z)+cosmo.Ok*(1.+z)*(1.+z)+cosmo.OL)); 
}

double Vol(double z, struct cosmo_param cosmo)
{
  /*   Devuelve el resultado en parsecs al cubo */
  double DM;
  double dH3;
  double Vc=0; /* Volumen comovil */
  double zstep;
  int i;
  double DM_DH;

  zstep=ZMAX_CACHED/NUMZ_CACHED;
  i=rint(z/zstep);
  if (i >= NUMZ_CACHED || i < 0)
  {
    printf("Out of range for cached values. z = %g\n",z);
    exit(-1);
  }
  DM = cosmo.cached.DM[i];
  DM_DH = DM / cosmo.DH;
  dH3=cosmo.DH*1.e6*cosmo.DH*1.e6*cosmo.DH*1.e6; /* La distancia de Hubble al cubo */
  
  if(cosmo.Ok>0) {    /* Omega_k > 0 */
    Vc=(4*M_PI*dH3)/2/cosmo.Ok*(DM_DH*sqrt(1+cosmo.Ok*DM_DH*DM_DH)-1/sqrt(cosmo.Ok)*asinh(sqrt(cosmo.Ok)*DM_DH));
  }
  else if(cosmo.Ok<0) {   /* Omega_k < 0 */
    Vc=(4*M_PI*dH3)/2/cosmo.Ok*(DM_DH*sqrt(1+cosmo.Ok*DM_DH*DM_DH)-1/sqrt(fabs(cosmo.Ok))*asin(sqrt(fabs(cosmo.Ok))*DM_DH));
  }
  else if(cosmo.Ok==0) {   /* Omega_k = 0 */
    if(DEBUG2) printf(" Ok = 0\n");
    Vc=(4*M_PI)/3*DM*DM*DM*1e6*1e6*1e6; /* pc**3 */
  }
  return(Vc); 
}

double Vol_OmegaLambda0(double z, struct cosmo_param cosmo) {
  /*   Devuelve el resultado en parsecs al cubo */
  double c=299792.46; /*  En km/s */
  double Ok;
  double dm_dh;
  double dH3;
  double Vc=0; /* Volumen comovil */
  double q0;
  q0=cosmo.OM/2.;

  if(cosmo.OL != 0)
  {
    printf("WARNING: Results only valid for OmegaLambda = 0\n");
  }

  Ok=1-2*q0;   /* Omega_k, que para Omega_lambda=0 vale 1-2q0 */
  dH3=(c/cosmo.H0*1.e6)*(c/cosmo.H0*1.e6)*(c/cosmo.H0*1.e6); /* La distancia de Hubble al cubo */
  dm_dh=2*(2-2*q0*(1-z)-(2-2*q0)*sqrt(1+2*q0*z))/(4*q0*q0*(1+z));
  if(q0<0.5) {    /* Omega_k > 0 */
    Vc=(4*M_PI*dH3)/2/Ok*(dm_dh*sqrt(1+Ok*(dm_dh)*(dm_dh))-1/sqrt(Ok)*asinh(sqrt(Ok)*dm_dh));
  }
  else if(q0>0.5) {   /* Omega_k < 0 */
    Vc=(4*M_PI*dH3)/2/Ok*(dm_dh*sqrt(1+Ok*(dm_dh)*(dm_dh))-1/sqrt(fabs(Ok))*asin(sqrt(fabs(Ok))*dm_dh));
  }
  else if(q0==0.5) {   /* Omega_k = 0 */
    /*     printf(" Q0 = 0.5\n"); */
    Vc=(4*M_PI*dH3)/3*dm_dh*dm_dh*dm_dh;
  }
  return(Vc); 
}

double dVdz(double z, struct cosmo_param cosmo)
{
//  double c=299792.46; /*  En km/s */
//  double Dang_DH,DH;  /* Distancia angular segun  Hogg astroph/9905116 */
//  double E;        /* Funcion E. */
//  /*   Esto es dV/dz, que no es igual a dV/dz=4pi dV/dz/dOmega, Omega: angulo so
//       lido */
//  Dang_DH=2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0
//								      *cosmo.q0*(1+z)*(1+z)); /* La distancia angular */
//  /*   Formula solo valida para Omega_lambda=0. Pero Dang=DM/(1+z). Eso siempre.
//   */
//  DH=(c/cosmo.H0*1.e6); /* La distancia de Hubble */
//  E=sqrt(2*cosmo.q0*(1+z)*(1+z)*(1+z)+(1-2*cosmo.q0)*(1+z)*(1+z));
//  /*   Esta formula solo es valida para Omega_lambda=0 */
//  return(4*M_PI*(1+z)*(1+z)*Dang_DH*Dang_DH/E*DH*DH*DH);
//  /*   El resultado en parcsecs al cuadrado */
  double zstep;
  int i;

  zstep=ZMAX_CACHED/NUMZ_CACHED;
  i=rint(z/zstep);
  if (i >= NUMZ_CACHED || i < 0)
  {
    printf("Out of range for cached values. z = %g\n",z);
    exit(-1);
  }
  return(cosmo.cached.dVdz[i]);
}

double dVdz_OmegaLambda0(double z, struct cosmo_param cosmo)
{
  double c=299792.46; /*  En km/s */
  double Dang_DH,DH;  /* Distancia angular segun  Hogg astroph/9905116 */
  double E;        /* Funcion E. */
  double q0;

  if(cosmo.OL != 0)
  {
    printf("WARNING: Results only valid for OmegaLambda = 0\n");
  }
  q0 = cosmo.OM / 2.;
  /*   Esto es dV/dz, que no es igual a dV/dz=4pi dV/dz/dOmega, Omega: angulo so
       lido */
  Dang_DH=2*(2-2*q0*(1-z)-(2-2*q0)*sqrt(1+2*q0*z))/(4*q0*q0*(1+z)*(1+z)); /* La distancia angular */
								      
  /*   Formula solo valida para Omega_lambda=0. Pero Dang=DM/(1+z). Eso siempre.
   */
  DH=(c/cosmo.H0*1.e6); /* La distancia de Hubble */
  E=sqrt(2*q0*(1+z)*(1+z)*(1+z)+(1-2*q0)*(1+z)*(1+z));
  /*   Esta formula solo es valida para Omega_lambda=0 */
  return(4*M_PI*(1+z)*(1+z)*Dang_DH*Dang_DH/E*DH*DH*DH);
  /*   El resultado en parcsecs al cuadrado */
}

double Lum(double z, double flux, struct cosmo_param cosmo)
{
  /* Devuelve la luminosidad dado el flujo aparente */
  /* El flujo debe ir en watt/m2 y 
     la luminosidad la devuelve en Watt */
//  double c=299792.46; /*  En km/s */
  double dlum;
  double L;
  double zstep;
  int i;

  zstep=ZMAX_CACHED/NUMZ_CACHED;
  i=rint(z/zstep);
  if (i >= NUMZ_CACHED || i < 0)
  {
    printf("Out of range for cached values. z = %g\n",z);
    exit(-1);
  }
 
//  double pi=3.1415926535897932384;
  /*   Calculo la distancia de luminosidad relativista. */
//  dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z));
  dlum = pow(10.,cosmo.cached.log10dlum[i]);
  /* Esta dlum esta dada en Megaparsec. 
     Para pasar a metros: 3.08567818589e22 metros  = 1 Mpc */
  L=4*M_PI*dlum*3.08567818589e22*dlum*3.08567818589e22*flux;
  return(L);
}

double dLumdflux(double z, struct cosmo_param cosmo)
{
  /* Devuelve la derivada de la luminosidad con respecto al flujo aparente */
  /* la derivada de la luminosidad la devuelve en Watt / (watt/m2 ) = m2*/
  double dlum;
  double dLdf;
  double zstep;
  int i;

  zstep=ZMAX_CACHED/NUMZ_CACHED;
  i=rint(z/zstep);
  if (i >= NUMZ_CACHED || i < 0)
  {
    printf("Out of range for cached values. z = %g\n",z);
    exit(-1);
  }
  /*   Calculo la distancia de luminosidad relativista. */
  //dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z
//))/(4*cosmo.q0*cosmo.q0*(1+z));
  dlum = pow(10.,cosmo.cached.log10dlum[i]);
  /* Esta dlum esta dada en Megaparsec. 
     Para pasar a metros: 3.08567818589e22 metros  = 1 Mpc */
  dLdf=4*M_PI*dlum*3.08567818589e22*dlum*3.08567818589e22;
  return(dLdf);
}

double Flux(double z, double Lum, struct cosmo_param cosmo)
{
  /* Devuelve el flujo aparente dada la luminosidad */
  /* La luminosidad debe ir in watt y
     el flujo lo devuelve en watt/m2 */
  double dlum;
  double flux;
  double zstep;
  int i;

  zstep=ZMAX_CACHED/NUMZ_CACHED;
  i=rint(z/zstep);
  if (i >= NUMZ_CACHED || i < 0)
  {
    printf("Out of range for cached values. z = %g\n",z);
    exit(-1);
  }
  dlum = pow(10.,cosmo.cached.log10dlum[i]);
  /*   Calculo la distancia de luminosidad relativista. */
//  dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z));
  /* Esta dlum esta dada en Megaparsec. 
     Para pasar a metros: 3.08567818589e22 metros  = 1 Mpc */
  flux=Lum/4/M_PI/dlum/3.08567818589e22/dlum/3.08567818589e22;
  return(flux);
}

double mag(double z, double M, struct cosmo_param cosmo)
{
  double log10dlum;
  double zstep;
  int i;

  zstep=ZMAX_CACHED/NUMZ_CACHED;
  i=rint(z/zstep);
  if (i >= NUMZ_CACHED || i < 0)
  {
    printf("Out of range for cached values. z = %g\n",z);
    exit(-1);
  }
  log10dlum = cosmo.cached.log10dlum[i];
  //double deltaM;
  /* dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z)); */
  //dlum=c/cosmo.H0*(1-cosmo.q0*(1-z)-(1-cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(cosmo.q0*cosmo.q0);
  //deltaM=-5*log10(dlum/10e-6);  /* 10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
  //return(M-deltaM);
  return(M+5*log10dlum+25.);
}

double mag_OmegaLambda0(double z, double M, struct cosmo_param cosmo)
{
  double c=299792.46; /*  En km/s */
  double q0;
  double dlum;

  if(cosmo.OL != 0)
  {
    printf("WARNING: Results only valid for OmegaLambda = 0\n");
  }
  q0 = cosmo.OM / 2.;
  dlum=c/cosmo.H0*(1-q0*(1-z)-(1-q0)*sqrt(1+2*q0*z))/(q0*q0);
  return(M+5*log10(dlum)+25.);
}

double Mag(double z, double m, struct cosmo_param cosmo)
{
  double log10dlum;
  double zstep;
  int i;

  zstep=ZMAX_CACHED/NUMZ_CACHED;
  i=rint(z/zstep);
  if (i >= NUMZ_CACHED || i < 0)
  {
    printf("Out of range for cached values. z = %g\n",z);
    exit(-1);
  }
  log10dlum = cosmo.cached.log10dlum[i];
  /* double deltaM; */
  /* dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z)); */
  /* dlum=c/cosmo.H0*2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0); */
  //dlum=c/cosmo.H0*(1-cosmo.q0*(1-z)-(1-cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(cosmo.q0*cosmo.q0);
  /*   Esto es la distancia de luminosidad seg�n Hogg. */
  /* deltaM=-5*log10(dlum/10e-6); */  /* 10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
  /* deltaM=-5*log10(dlum)-25.; */  /* soluci�n m�s simple, pero no aumenta significativamente la vel. */
  /* return(deltaM+m); */
  return(-5.*log10dlum-25.+m);
}

double Mag_OmegaLambda0(double z, double m, struct cosmo_param cosmo)
{
  double c=299792.46; /*  En km/s */
  double dlum;
  double q0;

  if(cosmo.OL != 0)
  {
    printf("WARNING: Results only valid for OmegaLambda = 0\n");
  }
  q0 = cosmo.OM / 2.;

  /* double deltaM; */
  /* dlum=c/cosmo.H0*2*(1+z)*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z)); */
  /* dlum=c/cosmo.H0*2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0); */
  dlum=c/cosmo.H0*(1-q0*(1-z)-(1-q0)*sqrt(1+2*q0*z))/(q0*q0);
  /*   Esto es la distancia de luminosidad seg�n Hogg. */
  /* deltaM=-5*log10(dlum/10e-6); */  /* 10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
  /* deltaM=-5*log10(dlum)-25.; */  /* soluci�n m�s simple, pero no aumenta significativamente la vel. */
  /* return(deltaM+m); */
  return(-5.*log10(dlum)-25.+m);
}

void Mag_err(double z, double errz, double m, double errm, struct cosmo_param cosmo, double *Mag, double *errMag)
{
  double c=299792.46; /*  En km/s */
  double dlum;
  double deltaM;
  double ddldz;

  if(cosmo.OL != 0)
  {
    printf("Not implemented\n");
  }
  else
  {
    dlum=c/cosmo.H0*2*(1+z)*(2-cosmo.OM*(1-z)-(2-cosmo.OM)*sqrt(1+cosmo.OM*z))/(cosmo.OM*cosmo.OM*(1+z));
    ddldz=dlum/(1+z)+c/cosmo.H0*2*(1+z)*(cosmo.OM-cosmo.OM*(1-0.5*cosmo.OM)/sqrt(1+cosmo.OM*z))/(cosmo.OM*cosmo.OM*(1+z))+c/cosmo.H0*2*(1+z)*(2-cosmo.OM*(1-z)-(2-cosmo.OM)*sqrt(1+cosmo.OM*z))/(cosmo.OM*cosmo.OM*(1+z))/(cosmo.OM*cosmo.OM*(1+z))*(cosmo.OM*cosmo.OM);
    /*   Esto es la distancia de luminosidad seg�n Hogg. */
    deltaM=-5*log10(dlum/10e-6);  /* 10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
    *Mag=deltaM+m;
    *errMag=sqrt(errm*errm+errz*errz*(5/dlum/log(10))*(5/dlum/log(10))*ddldz*ddldz);
  }
}

double Z_l(double flux, double L, struct cosmo_param cosmo) {
  double c=299792.46; /*  En km/s */
  double dlum;
  double dlum_dh; 
  double OM;
  double ztemp=0.;
  double pi=3.1415926535897932384;
  if(cosmo.OL != 0)
  {
    printf("Not implemented\n");
  }
  else
  {
    dlum=sqrt(L/4./pi/flux)/3.08567818589e22;
    dlum_dh=dlum/c*cosmo.H0;
    OM=cosmo.OM;  /*        //Esto solo es valido para Omega_lambda=0 */
    ztemp=(dlum_dh*OM-2.+OM+(2.-OM)*sqrt(1.+2.*dlum_dh*OM*OM))/2.;
  }
  return(ztemp);
}

double Z_m(double m, double M, struct cosmo_param cosmo)
{
  double log10dlum;
  double ztemp;
  double deltam;
  double log10dlum_step;
  int i;

  log10dlum_step=(cosmo.cached.maxlog10dlum-cosmo.cached.minlog10dlum)/NUMZ_CACHED;
  deltam=M-m;
  log10dlum=(-deltam+5.)/5.-6.;
  if(DEBUG1) printf("Z_m: log10dlum_step %g\n",log10dlum_step);
  if(DEBUG1) printf("Z_m: log10dlum %g\n",log10dlum);
  i=rint((log10dlum-cosmo.cached.minlog10dlum)/log10dlum_step);
  if (i >= NUMZ_CACHED || i < 0)
  {
    printf("Out of range for cached values. i %d\n",i);
    if(DEBUG1) printf("maxlog10dlum %g\n",cosmo.cached.maxlog10dlum);
    exit(-1);
  }
  ztemp = cosmo.cached.invz[i];
  return(ztemp);
}

double Z_m_OmegaLambda0(double m, double M, struct cosmo_param cosmo) {
  double c=299792.46; /*  En km/s */
  double dlum;
  double dlum_dh; 
  double OM;
  double ztemp;
  double deltam;
  if(cosmo.OL != 0)
  {
    printf("WARNING: Results only valid for OmegaLambda = 0\n");
  }
  deltam=M-m;
  dlum=pow(10.,(-deltam+5.)/5.-6.);
  if(DEBUG2) printf("Z_m_old: log10(dlum) %g\n",log10(dlum));
  dlum_dh=dlum/c*cosmo.H0;
  OM=cosmo.OM;  /*        Esto solo es valido para Omega_lambda=0 */
  ztemp=(dlum_dh*OM-2.+OM+(2.-OM)*sqrt(1.+2.*dlum_dh*OM*OM))/2.;
  return(ztemp);
}


//double D_co(double z, struct cosmo_param cosmo) {
//  double c=299792.46; /*  En km/s */
//  return(c/cosmo.H0*2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z)));

//}

//double D_lum(double z, struct cosmo_param cosmo) {
//  return((1+z)*D_co(z,cosmo));
//}

//double D_ang(double z, struct cosmo_param cosmo) {
//  double c=299792.46; /*  En km/s */
//  return(c/cosmo.H0*1.e6*(2*(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(4*cosmo.q0*cosmo.q0*(1+z)*(1+z))));
//}
