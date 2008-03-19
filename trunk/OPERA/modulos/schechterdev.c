#include "modulos.h"

#define NSTEP_LF 200
#define NSTEP_Z  500
#define NSTEP_MAG 500
#define NSTEP_MAG_FSEL 500
#define NSTEP_MAGAP_FSEL 200
#define NSTEP_LUM_FSEL 500
#define NSTEP_EW 200
#define ZMIN 0.00001
#define DEBUG 0
#define DEBUG2 0



double Schechterdev_M(struct Schlf_M lf, double Mlow,double Mup) {
  double xmin,xmax; /* Corresponden a los valores limites de x=L/L* */
  /*   xmax puede ser infinito. */
  double ymin,ymax; /* ymin=ln(xmin) */
  double e;
  double xe;
  double mag;
  static long idum =-1;
  struct timeval tv;
  if(idum==-1) 
  {
    gettimeofday(&tv,NULL);
    idum=tv.tv_usec;
  }
/*  if(idum==-1) idum=-(long)time(NULL)/2; */
  xmin=pow(10.,-0.4*(Mlow-lf.Mstar));
  xmax=pow(10.,-0.4*(Mup -lf.Mstar));

  if(xmax> 1000 && xmin<1 ) xmax=100.; /*  Para que irse tan lejos!, si es cero */
  ymin=log(xmin);ymax=log(xmax);
  if(lf.alfa<0) {
    do {
      do  {
        xe=Expdev();
      } while (xe < xmin || xe > xmax);
      e=pow(xe,lf.alfa)/(pow(xmin,lf.alfa));
    } while (ran2(&idum) > e);
    mag=lf.Mstar-2.5*log10(xe); /* Esta es la buena */
    return mag;
  }
  else return(0);
}

double zSchdev_M(struct Schlf_M lf, double zlow,double zup,double mlow,double mup, struct cosmo_param cosmo) {
  
  int nstep_lf=1000;
  int nstep_z=1000;
  int i,j;
  double ri;
  double z,logz,M;
  double rr;
  double Mlow,Mup;
  double sch_int=0;
  double Dang_DH;  /* Distancia angular segun  Hogg astroph/9905116 */
  double E;        /* Funcion E. */
  double dVdz;      /* Los factores anteriores , para calcular dVdz */
  static double fz[NSTEP_Z],fz_sum[NSTEP_Z],flogz[NSTEP_Z]; /* //fz es la funcion distribucion en z, y fz_sum es la integral */
  static double fz_int;
  double dz;
  double elev,tmp;
  static double zplot[NSTEP_Z],logzplot[NSTEP_Z];
  static double Mstarold,alfaold,zupold,zlowold,mlowold,mupold;
  int oldflag=0;
  static long idum =-1;
  struct timeval tv;
  if(idum==-1) 
  {
    gettimeofday(&tv,NULL);
    idum=tv.tv_usec;
  }
  //if(idum==-1) idum=-(long)time(NULL)/2;
  nstep_lf=NSTEP_LF;
  nstep_z =NSTEP_Z;

  
  if(Mstarold==lf.Mstar && alfaold==lf.alfa && zupold==zup && zlowold==zlow && mlowold==mlow && mupold==mup) {
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    oldflag=1;
  }
  else {
    Mstarold=lf.Mstar;alfaold=lf.alfa;zupold=zup;zlowold=zlow;mlowold=mlow;mupold=mup;
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    fz_int=0;
    for(i=0;i<nstep_z;i++) {
      logz=log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.);
      z=exp(logz);
      zplot[i]=z;
      logzplot[i]=logz;
      Mlow=Mag(z,mlow,cosmo);
      Mup =Mag(z,mup,cosmo);
      /*       Aqui integro la funcion de Schechter por enesima vez */
      sch_int=0;
      if(lf.Mstar-Mup> 5 && Mlow-lf.Mstar>0 ) Mup=lf.Mstar-5.; /*  Para que irse tan lejos!, si es cero */
      for(j=0;j<nstep_lf;j++) {
        M=Mlow+j*(Mup-Mlow)/(nstep_lf-1.);      
        elev=pow(10.,0.4*(lf.Mstar-M));
        tmp=pow(elev,lf.alfa+1.);  
        tmp=tmp*exp(-elev);
        sch_int+=tmp;
      }
      sch_int*=(Mlow-Mup)/nstep_lf;
      /*       Ahora calculamos dVdz, pero sin los factores, que no influyen */
      Dang_DH=(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(1+z); /* La distancia angular * (1+z) */
      E=sqrt(2*cosmo.q0*(1+z)*(1+z)*(1+z)+(1-2*cosmo.q0)*(1+z)*(1+z));
      dVdz=(Dang_DH*Dang_DH/E);
      /*       Aqui iria rho(z) si la densidad comovil variase con el z */
      /*       El z del final es porque estoy integrando en log(z), y dz=z*d(log(z)) */
      flogz[i]=sch_int*dVdz*z;
      dz=exp(log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.))-exp(log(zlow)+(i-1)*(log(zup)-log(zlow))/(nstep_z-1.));
      fz[i]=sch_int*dVdz;
      fz_int+=flogz[i];
    }
  }

  rr=fz_int*ran2(&idum);
  fz_sum[0]=0.;
  i=0;
  do {
    i++;
    fz_sum[i]=fz_sum[i-1]+flogz[i];
  } while (rr>fz_sum[i] && i < nstep_z-1);
  /*   Interpolamos linealmente la i: */
  ri=((i-1.)*(fz_sum[i]-rr)+i*(rr-fz_sum[i-1]))/(fz_sum[i]-fz_sum[i-1]);
  logz=log(zlow)+(log(zup)-log(zlow))/(nstep_z-1)*ri;
  z=exp(logz);
  return(z);
}

double Schechterdev_L(struct Schlf_L lf, double Llow,double Lup) {
  double xmin,xmax; /* Corresponden a los valores limites de x=L/L* */
  /*   xmax puede ser infinito. */
  double ymin,ymax; /* ymin=ln(xmin) */
  double e;
  double xe;
  double lum;
  static long idum =-1;
  if(idum==-1) idum=-(long)time(NULL)/2;
  xmin=Llow/lf.Lstar;
  xmax=Lup/lf.Lstar;

  if(xmax> 1000 && xmin<1 ) xmax=100.; /*  Para que irse tan lejos!, si es cero */
  ymin=log(xmin);ymax=log(xmax);
  if(lf.alfa<0) {
    do {
      do  {
        xe=Expdev();
      } while (xe < xmin || xe > xmax);
      e=pow(xe,lf.alfa)/(pow(xmin,lf.alfa));
    } while (ran2(&idum) > e);
    lum=xe*lf.Lstar;
    return(lum);
  }
  else return(0);
}

double zSchdev_L(struct Schlf_L lf, double zlow,double zup,double ffaint,double fbright, struct cosmo_param cosmo) {
  
  int nstep_lf=1000;
  int nstep_z=1000;
  int i;
  double ri;
  double z,logz;
  double rr;
  double Llow,Lup;
  double xmin,xmax;
  double sch_int=0;
  double Dang_DH;  /* Distancia angular segun  Hogg astroph/9905116 */
  double E;        /* Funcion E. */
  double dVdz;      /* Los factores anteriores , para calcular dVdz */
  static double fz[NSTEP_Z],fz_sum[NSTEP_Z],flogz[NSTEP_Z]; /* //fz es la funcion distribucion en z, y fz_sum es la integral */
  static double fz_int;
  double dz;
  static double zplot[NSTEP_Z],logzplot[NSTEP_Z];
  static double Lstarold,alfaold,zupold,zlowold,ffaintold,fbrightold;
  int oldflag=0;
  static long idum =-1;
  if(idum==-1) idum=-(long)time(NULL)/2;
  nstep_lf=NSTEP_LF;
  nstep_z =NSTEP_Z;

  
  if(Lstarold==lf.Lstar && alfaold==lf.alfa && zupold==zup && zlowold==zlow && ffaintold==ffaint && fbrightold==fbright) {
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    oldflag=1;
  }
  else {
    if(DEBUG) printf(" ESTE ES NUEVO\n"); 
    Lstarold=lf.Lstar;alfaold=lf.alfa;zupold=zup;zlowold=zlow;ffaintold=ffaint;fbrightold=fbright;
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    fz_int=0;
    for(i=0;i<nstep_z;i++) {
      logz=log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.);
      z=exp(logz);
      zplot[i]=z;
      logzplot[i]=logz;

      Llow=Lum(z,ffaint,cosmo);
      Lup =Lum(z,fbright,cosmo);
      xmin=Llow/lf.Lstar;
      xmax=Lup /lf.Lstar;
      if(DEBUG) printf(" iz %d Llow %g Lup %g\n",i,Llow,Lup); 
      /*       Aqui integro la funcion de Schechter por enesima vez */
      sch_int=0;
      if(xmax> 200) xmax=200.; /*  Para que irse tan lejos!, si es cero */
      if(xmin>xmax) xmin=xmax-50.;
      if(DEBUG)printf(" xmin %g xmax %g\n",xmin,xmax); 
      
      sch_int=(incom(1+lf.alfa,xmax)-incom(1+lf.alfa,xmin));
/*       printf(" sch_int %g\n",sch_int); */
      /*       Ahora calculamos dVdz, pero sin los factores, que no influyen */
      Dang_DH=(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(1+z); /* La distancia angular * (1+z) */
      E=sqrt(2*cosmo.q0*(1+z)*(1+z)*(1+z)+(1-2*cosmo.q0)*(1+z)*(1+z));
      dVdz=(Dang_DH*Dang_DH/E);
      /*       Aqui iria rho(z) si la densidad comovil variase con el z */
      /*       El z del final es porque estoy integrando en log(z), y dz=z*d(log(z)) */
      if(DEBUG) printf(" z %f sch_int %f dVdz %f \n",z,sch_int,dVdz);
      flogz[i]=sch_int*dVdz*z;
      dz=exp(log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.))-exp(log(zlow)+(i-1)*(log(zup)-log(zlow))/(nstep_z-1.));
      fz[i]=sch_int*dVdz;
      fz_int+=flogz[i];
    }
  }

  rr=fz_int*ran2(&idum);
  fz_sum[0]=0.;
  i=0;
  do {
    i++;
    fz_sum[i]=fz_sum[i-1]+flogz[i];
  } while (rr>fz_sum[i] && i < nstep_z-1);
  /*   Interpolamos linealmente la i: */
  ri=((i-1.)*(fz_sum[i]-rr)+i*(rr-fz_sum[i-1]))/(fz_sum[i]-fz_sum[i-1]);
  logz=log(zlow)+(log(zup)-log(zlow))/(nstep_z-1)*ri;
  z=exp(logz);
  return(z);
}

