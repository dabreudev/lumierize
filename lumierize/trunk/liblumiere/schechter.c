#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_gamma.h>
#include "cosmology.h"
#include "functions.h"
#include "schechter.h"

//#define NSTEP_LF 200
#define NSTEP_Z  200
#define NSTEP_COLOR  300
//#define NSTEP_Z  80
#define NSTEP_MAG 500
#define NSTEP_MAG_FSEL 400
//#define NSTEP_MAGAP_FSEL 300
#define NSTEP_MAGAP_FSEL 50
#define NSTEP_LUM_FSEL 500
//#define NSTEP_EW 200
#define NSTEP_EW 50
#define ZMIN 0.00001
#define ZMAX 10.
#define DEBUG 0
#define DEBUG2 0
#define DEBUG3 0


double Schechter_M(double M, struct Schlf_M lf) 
{
  double tmp;
  double elev;
  
  elev=pow(10.,0.4*(lf.Mstar-M));
  tmp=0.92103404*lf.phistar*pow(elev,lf.alfa+1.);   /*  0.92103404= 0.4*ln(10) */
  tmp=tmp*exp(-elev);
  /* Lo devuelvo SIN logaritmos */
  return(tmp);
}

double lnSchechter_M(double M, struct Schlf_M lf) 
{
  double tmp;
  double elev;
  
  elev=pow(10.,0.4*(lf.Mstar-M));
  tmp=-0.08225828358 + log(lf.phistar) + (lf.alfa + 1) *
    (lf.Mstar-M) * 0.92103404 - elev;
  /* -0.08225828358 = ln(0.4*ln(10)) */
  /* Este devuelve el logaritmo natural de la func de Schecter */
  return(tmp);
}

double Schechter_L(double L, struct Schlf_L lf) {
  double tmp;
  double L_Lstar;
  
  L_Lstar=L/lf.Lstar; 
/*   printf(" L_Lstar %g phis %g alfa %f\n",L_Lstar,lf.phistar,lf.alfa); */
  tmp=log(lf.phistar)-log(lf.Lstar)+lf.alfa*log(L_Lstar)-L_Lstar;
/*   printf(" %g %g %g %g\n",log10(lf.phistar),-log10(lf.Lstar),lf.alfa*log10(L_Lstar),-L_Lstar*log10(exp(1))); */
/*   printf(" Todo %g\n",tmp); */
  /* Lo devuelvo CON logaritmos neperianos!!! Cuidado con eso*/
  return(tmp);
}


double Sch_rhoz_L(double zint, struct Schlf_L lf, double zlow,double zup,double fluxlim,struct cosmo_param cosmo) {
  
  int nstep_z;
  int i;
  double z,logz;
  double Llow;
  double xmin,xmax;
  double sch_int=0;
  double Dang_DH;  /* Distancia angular segun  Hogg astroph/9905116 */
  double E;        /* Funcion E. */
  double dVdz;      /* Los factores anteriores , para calcular dVdz */
  static double fz[NSTEP_Z],flogz[NSTEP_Z]; /* //fz es la funcion distribucion en z, y fz_sum es la integral */
  static double fz_int;
  double dz;
  static double zplot[NSTEP_Z],logzplot[NSTEP_Z];
  static double Lstarold,alfaold,zupold,zlowold,fluxlimold;
  int oldflag=0;
  double rhoz;
  double q0;

  nstep_z =NSTEP_Z;

  if (cosmo.OL != 0)
  {
    printf("Method only valid for OL = 0\n");
    return(-1);
  }
  q0 = cosmo.OM/2.;
  if(Lstarold==lf.Lstar && alfaold==lf.alfa && zupold==zup && zlowold==zlow && fluxlimold==fluxlim) {
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    oldflag=1;
  }
  else {
    Lstarold=lf.Lstar;alfaold=lf.alfa;zupold=zup;zlowold=zlow;fluxlimold=fluxlim;
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    fz_int=0;
    for(i=0;i<nstep_z;i++) {
      logz=log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.);
      z=exp(logz);
      zplot[i]=z;
      logzplot[i]=logz;
      Llow=Lum(z,fluxlim,cosmo);
      xmin=Llow/lf.Lstar;
      xmax=100;
      sch_int=0;
      if(xmax> 1000 && xmin<1 ) xmax=100.; /*  Para que irse tan lejos!, si es cero */
      sch_int=(incom(1+lf.alfa,xmax)-incom(1+lf.alfa,xmin));
/*       printf(" sch_int %g\n",sch_int); */
      /*       Ahora calculamos dVdz, pero sin los factores, que no influyen */
      Dang_DH=(2-2*q0*(1-z)-(2-2*q0)*sqrt(1+2*q0*z))/(1+z); /* La distancia angular * (1+z) */
      E=sqrt(2*q0*(1+z)*(1+z)*(1+z)+(1-2*q0)*(1+z)*(1+z));
      dVdz=(Dang_DH*Dang_DH/E);
      /*       Aqui iria rho(z) si la densidad comovil variase con el z */
      /*       El z del final es porque estoy integrando en log(z), y dz=z*d(log(z)) */
      flogz[i]=sch_int*dVdz*z;
      dz=exp(log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.))-exp(log(zlow)+(i-1)*(log(zup)-log(zlow))/(nstep_z-1.));
      fz[i]=sch_int*dVdz;
      fz_int+=flogz[i];
/*       printf(" z %f sch_int %g fz %g\n",zplot[i],sch_int,fz[i]); */
    }
  }
  rhoz=Lagr2_d(zplot,fz,nstep_z,zint);
  return(rhoz);
}



double Int_sch_M(struct Schlf_M lf, double zlow,double zup,double mlim,struct cosmo_param cosmo) 
{
  double z;
  int nz,nM;
  int iz;
  double Mlow;
  double Llow,Lstar;
  double N,Ngal;
  double zlow_l;
  double zup_l;

  zlow_l = (zlow < ZMIN ? ZMIN : zlow);
  zup_l = (zup < ZMIN ? ZMAX : zup);

  nz=NSTEP_Z;
  nM=NSTEP_MAG;

  Lstar=pow(10.,-0.4*lf.Mstar);

  N=0;

/*  printf("\n===================\n");
  printf("Antes del bucle de la integral\n");
  printf("------------------------------\n");
  printf(" N %i\n",N);
  printf(" zlow %g\n",zlow_l);
  printf(" zup %g\n",zup_l);
  printf(" milm %g\n",mlim);
  printf(" nz %i\n",nz);
  printf(" nM %i\n",nM);
  printf("===================\n\n");  */

  for(iz=0;iz<nz;iz++) 
  {
    z=zlow_l+iz*(zup_l-zlow_l)/(nz-1.);
    Mlow=Mag(z,mlim,cosmo);

    Llow=pow(10.,-0.4*Mlow);

    /* Y esto es con la funcion gamma incompleta incom */
    /* estamos haciendo la integral desde Llow/Lstar hasta infinito */
    /* debido a un underflow, tuvimos que poner este if */
    if(Llow/Lstar > 0.25 && (lf.alfa*log(Llow/Lstar) - Llow/Lstar) <= GSL_LOG_DBL_MIN)
    {
      Ngal=GSL_DBL_MIN;
    }
    else
    {
      Ngal=lf.phistar*(gsl_sf_gamma_inc(1+lf.alfa,Llow/Lstar))*dVdz(z,cosmo)/1.e18;
    }

    /* printf(" i %d Ngal %g z %g\n",i,Ngal,z); */ 
        
    N=N+Ngal;
  }
  N=N/nz*(zup_l-zlow_l);
  /* printf(" N %g phi %g mst %g\n",N,lf.phistar,lf.Mstar); */
  return(N);
}


double Int_sch_M_wC(struct Schlf_M lf, double zlow,double zup, double color_mean, double color_stddev, double mDetectLim,struct cosmo_param cosmo) 
{
  double z;
  int nz,nM,nColor;
  int i,iColor;
  double Mlow;
  double Llow,Lstar;
  double N,Ncolor,Ngal;  /* WARNING: Ncolor != nColor */
  double zlow_l;
  double zup_l;
  double color, colorLow, colorUp;
  double mDistLim;

  zlow_l = (zlow < ZMIN ? ZMIN : zlow);
  zup_l = (zup < ZMIN ? ZMAX : zup);
  colorLow=-5.*color_stddev+color_mean;
  colorUp=5.*color_stddev+color_mean;

  nz=NSTEP_Z;
  nColor=NSTEP_COLOR;  /* nColor -> nStepColor ... */
  nM=NSTEP_MAG;

  Lstar=pow(10.,-0.4*lf.Mstar);

  N=0;

  if(DEBUG3)
  {
    printf("\n===================\n");
    printf("Antes del bucle de la integral\n");
    printf("------------------------------\n");
    printf(" N %g\n",N);
    printf(" zlow %g\n",zlow_l);
    printf(" zup %g\n",zup_l);
    printf(" mDetectLim %g\n",mDetectLim);
    printf(" nz %i\n",nz);
    printf(" nM %i\n",nM);
    printf("===================\n\n");
  }

  for(i=0;i<nz;i++) 
  {
    z=zlow_l+i*(zup_l-zlow_l)/(nz-1.);
    Ncolor=0;
    for(iColor=0; iColor<nColor;iColor++)
    {

      /* color=distrib - detect -> J-K */

      color=colorLow+iColor*(colorUp-colorLow)/(nColor-1.);

      mDistLim=mDetectLim+color;
      Mlow=Mag(z,mDistLim,cosmo);

      Llow=pow(10.,-0.4*Mlow);

      /* Y esto es con la funcion gamma incompleta incom */
      /* parece que hab�a un error con las gamma incompleta*/
      /* Npar=lf.phistar*(gsl_sf_gamma(lf.alfa+1.)-incom(1+lf.alfa,Llow/Lstar))*dVdz(z,cosmo)/1.e18; */
      /* estamos haciendo la integral desde Llow/Lstar hasta infinito */
      /* debido a un underflow, tuvimos que poner este if */
      /* 0.25 es por la implementaci�n de gamma_inc.c de gsl */
      if(Llow/Lstar > 0.25 && (lf.alfa*log(Llow/Lstar) - Llow/Lstar) <= GSL_LOG_DBL_MIN)
      {
        Ngal=GSL_DBL_MIN;
      }
      else
      {
        if(DEBUG2) printf("alfa %g Llow/Lstar %g dVdz(z,cosmo) %g gsl_sf_gamma %g gsl_sf_gamma %g\n", lf.alfa, Llow/Lstar, dVdz(z,cosmo), gsl_sf_gamma_inc(1.+lf.alfa,Llow/Lstar), gsl_sf_gamma_inc(1.+lf.alfa,1.03162e-7));
        Ngal=lf.phistar*(gsl_sf_gamma_inc(1.+lf.alfa,Llow/Lstar))*dVdz(z,cosmo)/1.e18;
      }
      Ncolor += gaussian(color, color_mean, color_stddev)*Ngal;
      if(DEBUG2) printf(" color %g color_mean %g color_stddev %g\n",color,color_mean,color_stddev);
      if(DEBUG2) printf(" icolor %d Ngal %g Ncolor %g z %g\n",iColor,Ngal,Ncolor,z); 
    }
    Ncolor = Ncolor/nColor*(colorUp - colorLow);
    if(DEBUG2) printf(" iz %d Ncolor %g z %g\n",i,Ncolor,z);  

    N=N+Ncolor;
  }
  N=N/nz*(zup_l-zlow_l);
  if(DEBUG2) printf(" N %g phi %g mst %g\n",N,lf.phistar,lf.Mstar);
  return(N);
}


double Int_sch_L(struct Schlf_L lf, double zlow,double zup,double fluxlim, struct cosmo_param cosmo) {
  double z;
  int nz,nM;
  int i;
  double Llow;
  double xmin;
  double N,Npar;
  double zlow_l;
  double zup_l;

  zlow_l = (zlow < ZMIN ? ZMIN : zlow);
  zup_l =  (zup  < ZMIN ? ZMAX : zup);

  nz=NSTEP_Z;
  nM=NSTEP_MAG;

  N=0;
  
  for(i=0;i<nz;i++) {
    z=zlow+i*(zup_l-zlow_l)/(nz-1.);
    /* En luminosidades: */
/*     printf(" Ho %f q0 %f\n",cosmo.H0,cosmo.q0); */
    Llow=Lum(z,fluxlim,cosmo);
/*     printf("Ante z %f Llow/Lstar %f Lup/Lstar %f\n",z,Llow/lf.Lstar,Lup/lf.Lstar);  */
    xmin=Llow/lf.Lstar;
    Npar=lf.phistar*(gsl_sf_gamma(lf.alfa+1.)-incom(1+lf.alfa,xmin))*dVdz(z,cosmo)/1.e18; 
    /*if(DEBUG)   printf("N %f Npar %f  in1 %f in2 %f xmin %g xmax %g z %g\n",N,Npar,incom(1+lf.alfa,xmax),incom(1+lf.alfa,xmin),xmin,xmax,z);*/
    N=N+Npar;
  }
  N=N/nz*(zup_l-zlow_l);
/*  printf(" Devuelvo %f\n",N); */
  return(N);
}

/* void PlotSchLF_L( struct Schlf_L lfsch) { */

/*   float x, y, sigy; */
/*   double dum1,dum2,dum3; */
/*   double ymin,ymax; */
/*   float xmin,xmax; */
/*   int i; */

/*   int nline=100; */

/*   cpgpage(); */
  
/*   ymin=1e38; */
/*   ymax=-1e38; */
/*   for(i=0;i<lfstep.nbin;i++) { */
/*     if(!(lfstep.lf[i]==0 && lfstep.errlf[i]==0)) {  */
/*       if((double)log10(exp(lfstep.lf[i]))>ymax) ymax=(double)log10(exp(lfstep.lf[i])); */
/*       if((double)log10(exp(lfstep.lf[i]))<ymin) ymin=(double)log10(exp(lfstep.lf[i])); */
/*     } */
/*   } */
  
/*   ymin=ymin-2.5; */
/*   ymax=ymax+2.5; */
/*   pgLimits_d(lfstep.nbin,lfstep.lumi,&xmin,&xmax);  */
/*   xmin=xmin-1.5; */
/*   xmax=xmax+1.5; */
/*   cpgsci(1); */
/*   cpgsch(1.2); */
/*   cpgswin((float)xmax/log(10),(float)xmin/log(10),(float)ymin,(float)ymax); */
/*   cpgbox("BCTNSL",0,0,"BCTNS",0,0); */
  
  
/*   for(i=0;i<lfstep.nbin;i++) { */
/*     if(!(lfstep.lf[i]==0 && lfstep.errlf[i]==0)) { */
/*       x=log10(exp((lfstep.lumi[i]+lfstep.lumi[i+1])/2.)); */
/*       y=log10(exp(lfstep.lf[i])); */
/*       sigy=lfstep.errlf[i]/log(10.); */
/*       cpgpt1(x,y,17); */
/*       dum1=x; */
/*       dum2=y-sigy; */
/*       dum3=y+sigy; */
/*       cpgerry_d(1,&dum1,&dum2,&dum3,0.35);  */
/*     } */
/*   } */
  
/*   cpgsci(4); */
/*   cpgmove(xmin,Schechter_L(pow(10.,xmin)/log(10),lfsch)); */
  
/*   for(i=0;i<nline;i++) { */
/*     x=xmin+(xmax-xmin)*i/(nline-1); */
/*     y=Schechter_L(exp(x),lfsch); */
/*     x=x/log(10); */
/*     y=y/log(10); */

/*     cpgdraw(x,y); */
/*   } */

/*   cpgsci(1); */
/*   cpglab("Luminosity (W)","log\\dlum\\u(Phi) (#/Mpc3/W) ",""); */
  
/* } */

double Int_sch_f_M(struct Schlf_M lf, double zlow, double zup, struct fermifsel_M fsel,struct cosmo_param cosmo) {
  double z;
  double zlow_l, zup_l;
  int nz,nM_fs;
  int iz,iMag;
  double Mleft,Mright;
  double M,m;
  double Lleft,Lright,Lstar;
  double N,Ngal;

  if(DEBUG) printf(" Entro en Int_sch_f_M\n");
  zlow_l = (zlow < ZMIN ? ZMIN : zlow);
  zup_l = (zup < ZMIN ? ZMAX : zup);

  nz=NSTEP_Z;
  nM_fs=NSTEP_MAG_FSEL;

  Lstar=pow(10.,-0.4*lf.Mstar);

  N=0;

  for(iz=0;iz<nz;iz++) 
  {
    z=zlow_l+iz*(zup_l-zlow_l)/(nz-1.);
    /* Integro primero hasta 6 veces m�s a la izquierda del corte en Fermi */
    Mleft =Mag(z,fsel.magcut-6*fsel.deltamag,cosmo);
    Mright=Mag(z,fsel.magcut+6*fsel.deltamag,cosmo);
    Lleft=pow(10.,-0.4*Mleft);
    Lright=pow(10.,-0.4*Mright);
    if(DEBUG) printf(" mag cu %f delt %f\n",fsel.magcut,fsel.deltamag);
    if(DEBUG) printf(" z %f Mleft %f Mrith %f\n",z,Mleft,Mright);
    /* Esto que viene es haciendo la integral a pelo: de Mleft hasta Mright*/
    Ngal=0;
    for(iMag=0;iMag<nM_fs;iMag++) {
      M=Mleft+iMag*(Mright-Mleft)/(nM_fs-1.);
      m=mag(z,M,cosmo);
      Ngal=Ngal+Schechter_M(M,lf)*Fermi(m,fsel.magcut,fsel.deltamag);
      if(DEBUG2) printf(" M %f Npar %g Sch %g Fer %f \n",M,Ngal,Schechter_M(M,lf),Fermi(m,fsel.magcut,fsel.deltamag));
    }
    Ngal=Ngal/nM_fs*(Mright-Mleft)*dVdz(z,cosmo)/1.e18;  
    if(DEBUG) printf(" Primer Npar %g \n",Ngal);
      
    /* Y esto es con la funcion gamma incompleta incom desde Mleft a -inf*/
    /* debido a un underflow, tuvimos que poner este if */
    /* 0.25 es por la implementaci�n de gamma_inc.c de gsl */
    if(Lleft/Lstar > 0.25 && 
       (lf.alfa*log(Lleft/Lstar) - Lleft/Lstar) <= GSL_LOG_DBL_MIN)
    {
      Ngal+=GSL_DBL_MIN;
    }
    else
    {
      Ngal+=lf.phistar*(gsl_sf_gamma_inc(1.+lf.alfa,Lleft/Lstar))*
            dVdz(z,cosmo)/1.e18;
      //Ngal+=lf.phistar*(incom(1+lf.alfa,Lup/Lstar)-incom(1+lf.alfa,Lleft/Lstar))*dVdz(z,cosmo)/1.e18; 
      if(DEBUG) printf(" Segundo Npar %g \n",Ngal);
    }
    N=N+Ngal;
  }
  N=N/nz*(zup_l-zlow_l);
  return(N);
}

double Int_sch_f_M_wC(struct Schlf_M lf, double zlow,double zup, double color_mean, double color_stddev, struct fermifsel_M fsel, struct cosmo_param cosmo) 
{
  double z;
  double zlow_l, zup_l;
  int nz,nColor,nM_fs;
  int iz, iColor,iMag;
  double Mleft,Mright;
  double M,m;
  double Lleft,Lright,Lstar;
  double Ngal,Ngal_color,Ngal_mag;
  double color, colorLow, colorUp;
  double mDistMagCut;

  if(DEBUG) printf(" Entro en Int_sch_f_M\n");
  zlow_l = (zlow < ZMIN ? ZMIN : zlow);
  zup_l = (zup < ZMIN ? ZMAX : zup);
  colorLow = color_mean-6.*color_stddev;
  colorUp  = color_mean+6.*color_stddev;
  if(DEBUG) printf(" zlow_l %g zup_l %g zlow %g zup %g\n",zlow_l, zup_l, zlow, zup);

  nz=NSTEP_Z;
  nColor=NSTEP_COLOR;
  nM_fs=NSTEP_MAG_FSEL;

  Lstar=pow(10.,-0.4*lf.Mstar);

  Ngal=0;

  for(iz=0;iz<nz;iz++) 
  {
    z=zlow_l+iz*(zup_l-zlow_l)/(nz-1.);
    Ngal_color=0;
    for(iColor=0; iColor<nColor;iColor++)
    {
      /* color =distrib - detect */
      color=colorLow+iColor*(colorUp-colorLow)/(nColor -1.);
      mDistMagCut = fsel.magcut + color;
      /* Integro primero hasta 5 veces m�s a la izquierda 
         del corte en Fermi */
      Mleft =Mag(z,mDistMagCut-6*fsel.deltamag,cosmo);
      Mright=Mag(z,mDistMagCut+6*fsel.deltamag,cosmo);
      Lleft=pow(10.,-0.4*Mleft);
      Lright=pow(10.,-0.4*Mright);
      if(DEBUG) printf(" mag cu %f del %f\n",fsel.magcut,fsel.deltamag);
      if(DEBUG) printf(" z %f color %f Mleft %f Mrith %f\n",z, color, Mleft,Mright);
      /* Esto que viene es haciendo la integral a pelo: 
         de Mleft hasta Mright*/
      Ngal_mag=0;
      for(iMag=0;iMag<nM_fs;iMag++) {
        M=Mleft+iMag*(Mright-Mleft)/(nM_fs-1.);
        m=mag(z,M,cosmo);
        Ngal_mag=Ngal_mag+
          Schechter_M(M,lf)*Fermi(m-color,fsel.magcut,fsel.deltamag);
        if(DEBUG2) printf(" M %f Npar %g Sch %g Fer %f \n",M,Ngal_mag,Schechter_M(M,lf),Fermi(m,fsel.magcut,fsel.deltamag));
      }
      Ngal_mag=Ngal_mag/nM_fs*(Mright-Mleft);  
      if(DEBUG) printf(" Primer Npar %g \n",Ngal_mag);
      
      /* Y esto es con la funcion gamma incompleta incom desde Mleft a -inf*/
      /* debido a un underflow, tuvimos que poner este if */
      /* 0.25 es por la implementaci�n de gamma_inc.c de gsl */
      if(Lleft/Lstar > 0.25 && 
         (lf.alfa*log(Lleft/Lstar) - Lleft/Lstar) <= GSL_LOG_DBL_MIN)
      {
        Ngal_mag+=GSL_DBL_MIN;
      }
      else
      {
        Ngal_mag+=lf.phistar*(gsl_sf_gamma_inc(1.+lf.alfa,Lleft/Lstar));
        if(DEBUG) printf(" Segundo Npar %g \n",Ngal_mag);
      }
      Ngal_color+=Ngal_mag * gaussian(color, color_mean, color_stddev);
    }
    Ngal_color = Ngal_color/nColor*(colorUp - colorLow);
    if(DEBUG) printf(" Ngal_color %g\n",Ngal_color);
    Ngal=Ngal+Ngal_color*dVdz(z,cosmo)/1.e18;
  }
  Ngal=Ngal/nz*(zup_l-zlow_l);
  return(Ngal);
}

double Int_sch_f_L(struct Schlf_L lf, double zlow,double zup,struct fermifsel_L fsel, struct cosmo_param cosmo) {
  double z;
  int nz,nL_fs;
  int i,j;
  double Lleft,Lright;
  double Lup=1e60;
  double L,flux,fmin;
  double xmax;
  double xright;
  double N,Npar;

  zlow = (zlow < ZMIN ? ZMIN : zlow);

  nz=NSTEP_Z;
  nL_fs=NSTEP_LUM_FSEL;

  N=0;
  
  if(Lup/lf.Lstar>100) Lup=100*lf.Lstar; /* Para que irse tan lejos!, si es cero */
  for(i=0;i<nz;i++) {
    z=zlow+i*(zup-zlow)/(nz-1.);
    /* En luminosidades: */
    fmin=fsel.fluxcut-5*fsel.deltaflux;
    if(fmin<=0) fmin=fsel.fluxcut/20;
    Lleft=Lum(z,fmin,cosmo);
    Lright=Lum(z,fsel.fluxcut+5*fsel.deltaflux,cosmo);
    xright=Lright/lf.Lstar;
    xmax=200;
    if(xright> 150) xright=150;
    if(Lleft/lf.Lstar> 150) Lleft =150*lf.Lstar;
    /* Esto que viene es haciendo la integral a pelo: de Lleft hasta Lright*/
    Npar=0;
    if(DEBUG) printf(" Lleft %g Lright %g Lstar %g \n",Lleft/lf.Lstar, Lright/lf.Lstar, lf.Lstar);
    for(j=0;j<nL_fs;j++) {
      L=Lleft+j*(Lright-Lleft)/(nL_fs-1.);
      flux=Flux(z,L,cosmo);
      Npar=Npar+exp(Schechter_L(L,lf))*Fermi(flux,fsel.fluxcut,-fsel.deltaflux);
      /* 	printf("j %d L %g flux %g  Npar %g Sch %f Fer %f \n",j,L,flux,Npar,Schechter_L(L,lf),Fermi(flux,fsel.fluxcut,-fsel.deltaflux)); */
    }
    Npar=Npar/nL_fs*(Lright-Lleft)*dVdz(z,cosmo)/1.e18;  
    if(DEBUG) printf(" Primera parte Npar %g \n",Npar);

    if(xright> 190) xright=190.; /*  Para que irse tan lejos!, si es cero */
    Npar+=lf.phistar*(incom(1+lf.alfa,xmax)-incom(1+lf.alfa,xright))*dVdz(z,cosmo)/1.e18; 
    if(DEBUG) printf(" Al final Npar %g \n",Npar);
    N=N+Npar;
  }
  N=N/nz*(zup-zlow);
/*  printf(" Devuelvo %f\n",N); */
  return(N);
}



void PlotSchLF_L_ov( struct Schlf_L lfsch) {

  float x, y;
  float xmin,xmax;
  float ymin,ymax;
  static int color=2;
  static int linestyle=2;
  int i;

  int nline=300;

  
  cpgsci(color);
  cpgsls(linestyle);
  cpgsch(1.2);
  cpgqwin(&xmax,&xmin,&ymin,&ymax);

  xmin*=log(10);
  xmax*=log(10);

 printf(" xmin %f xmax %f ymin %f ymax %f\n",xmin,xmax,ymin,ymax);

  
  for(i=0;i<nline;i++) {
    x=xmin+(xmax-xmin)*i/(nline-1);
    y=Schechter_L(exp(x),lfsch);
    x=x/log(10);
    y=y/log(10);

    printf(" x %f y %f\n",x,y); 

    cpgdraw(x,y);
  }

  cpgsci(1);
  color++;
  linestyle++;
  
}
