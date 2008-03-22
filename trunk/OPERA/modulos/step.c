#include "modulos.h"

#define NSTEP_LF 200
#define NSTEP_Z  500
#define NSTEP_MAG 500
#define ZMIN 0.00001
#define DEBUG 0
#define DEBUG2 0

void mrqfunc_fitsch2steplf_L(double x,double *p,double *y,double *dyda,int n);
void mrqfunc_fitsch2steplf_M(double x,double *p,double *y,double *dyda,int n);
/* Estas no se están usando de momento, ya que se usa mrqfunc */
double amofunc_schechterfitmag_d(int n, double *x, double *y, double *p);
double amofunc_schechterfitlum_d(int n, double *x, double *y, double *p);


double Step_M(double M, struct Steplf_M lf) 
{
  
  int j;

  if(M<lf.magni[0] || M>lf.magni[lf.nbin]) 
    return(0);
  else 
    for(j=0;j<lf.nbin;j++) 
      if(M<=lf.magni[j+1]) return(exp(lf.lnlf[j]));
  return(0);
}

double Step_L(double L, struct Steplf_L lf) 
{
  
  int j;

  if(log(L)<lf.lumi[0] || log(L)>lf.lumi[lf.nbin]) 
    return(0);
  else 
    for(j=0;j<lf.nbin;j++) 
      if(log(L)<=lf.lumi[j+1]) return(exp(lf.lnlf[j]));
  return(0);
}


double Steplfdev_M(struct Steplf_M lf, double Mlow,double Mup)
{

  double r;
  double M=0;
  static long idum =-1;
  static int nbin=-1;
  static double *PLFcum;
  int i;
  /* Cuidado!! El bin j va desde magni[j] a magni[j+1], luego magni esta dimensionado
     de nbin+1 y lf de nbin. De esta manera, el intervalo es magni[j+1]-magni[j] */
  if(idum==-1) idum=-(long)time(NULL)/2;

  nbin=lf.nbin;
  if(nbin!=-1) free(PLFcum);
  PLFcum=vector_d(nbin);

  PLFcum[0]=exp(lf.lnlf[0])*(lf.magni[1]-lf.magni[0]);
  for(i=1;i<nbin;i++) {
    PLFcum[i]=exp(lf.lnlf[i])*(lf.magni[i+1]-lf.magni[i])+PLFcum[i-1];
  }
  /* Normalizo */
  for(i=0;i<nbin;i++)  PLFcum[i]=PLFcum[i]/PLFcum[nbin-1];

  if(Mup >lf.magni[0])       
  {
    printf(" Mup >lf.lumi[0] shouldn't happen Mup %f magni[0] %f. Exiting\n",Mup,lf.magni[0]);
    exit(1);
   return(lf.magni[0]);
  }
  if(Mlow<lf.magni[lf.nbin]) 
  {
    printf(" Mlow <lf.lumi[nbin] shouldn't happen Mlow %f magni[0] %f. Exiting\n",Mlow,lf.magni[nbin]);
    exit(1);
    return(lf.magni[lf.nbin]);
  }

   
  do {
    r=ran2(&idum);
    for(i=0;i<nbin;i++) {
      if(r<=PLFcum[i]) {
	if(i==0) M=(r-0.)/(PLFcum[0]-0)*(lf.magni[1]-lf.magni[0])+lf.magni[0];
	else     M=(r-PLFcum[i-1])/(PLFcum[i]-PLFcum[i-1])*(lf.magni[i+1]-lf.magni[i])+lf.magni[i];
	/* Estas a continuacion es si fuera como lo hacia antes: magni[j] es el centro del intervalo j.
	   Tal y como es ahora, el intervalo va de magni[j] a magni[j+1]. */
	/*       if(i==0) M=(r-0.)/(PLFcum[0]-0)*(lf.magni[1]-lf.magni[0])+lf.magni[0]-(lf.magni[1]-lf.magni[0])/2.; */
	/*       else     M=(r-PLFcum[i-1])/(PLFcum[i]-PLFcum[i-1])*((lf.magni[i+1]-lf.magni[i-1])/2.)+(lf.magni[i]+lf.magni[i-1])/2.; */
	if(DEBUG) printf(" De esta ya sale con %g  Mlow %g Mup %g\n",M,Mlow,Mup);
	break;
      }
    }
  } while (!(M>Mup && M<Mlow));

  return(M);

  printf(" Steplfdev: Never reaches this point\n");
  exit(1);

  return(lf.magni[nbin]);
}

double zSteplfdev_M(struct Steplf_M lf, double zlow, double zup, double mlow, double mup, struct cosmo_param cosmo) {
  
  int nstep_lf=1000;
  int nstep_z=1000;
  int i,j;
  int jstart,jend;
  double ri;
  double z,logz;
  double rr; 
  double Mlow,Mup;
  double step_int=0;
  double Dang_DH;  /* Distancia angular segun  Hogg astroph/9905116 */
  double E;        /* Funcion E. */
  double dVdz;      /* Los factores anteriores , para calcular dVdz */
  static double fz[NSTEP_Z],fz_sum[NSTEP_Z],flogz[NSTEP_Z]; /* //fz es la funcion distribucion en z, y fz_sum es la integral */
  static double fz_int;
  double dz;
  static double zplot[NSTEP_Z],logzplot[NSTEP_Z];
  static double zupold,zlowold,mlowold,mupold;
  static struct Steplf_M lfold;
  int oldflag=0;
  static long idum =-1;
  if(idum==-1) idum=-(long)time(NULL)/2;
  nstep_lf=NSTEP_LF;
  nstep_z =NSTEP_Z;

  
  if(!memcmp(&lf,&lfold,sizeof(struct Steplf_M)) && zlow==zlowold && zup==zupold && mlow==mlowold && mup==mupold) {
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    oldflag=1;
  }
  else {
    memcpy(&lfold,&lf,sizeof(struct Steplf_M));
    zupold=zup;zlowold=zlow;mlowold=mlow;mupold=mup;
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    if(DEBUG) printf(" NUEVO\n");
    fz_int=0;
    for(i=0;i<nstep_z;i++) {
      logz=log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.);
      z=exp(logz);
      zplot[i]=z;
      logzplot[i]=logz;
      Mlow=Mag(z,mlow,cosmo);
      Mup =Mag(z,mup,cosmo);
      if(DEBUG) printf(" %d z %f Mlow %f Mup %f mlow %f mup %f zlow %f zup %f\n",i,z,Mlow,Mup,mlow,mup,zlow,zup); 
      /*       Aqui integro la funcion de lumi de tipo Step */
      step_int=0;
      jstart=0;
      jend=lf.nbin;
      if(Mup<lf.magni[0])  jstart=-1;
      if(Mup>lf.magni[lf.nbin])	jstart=lf.nbin;
      if(Mlow<lf.magni[0])  jend=-1;
      if(Mlow>lf.magni[lf.nbin])   jend=lf.nbin;
      for(j=0;j<lf.nbin;j++) {
	if(Mup>lf.magni[j] && Mup<lf.magni[j+1]) jstart=j;
	if(Mlow>lf.magni[j] && Mlow<lf.magni[j+1]) jend=j;
      }
      if(jstart!=-1 && jstart!=lf.nbin) {
	step_int+=exp(lf.lnlf[jstart])*(lf.magni[jstart+1]-Mup);
 	if(DEBUG) printf(" %d Sumando en mag %f-%f Mup %f   %g  : %g\n",jstart,lf.magni[jstart],lf.magni[jstart+1],Mup,lf.lnlf[jstart],step_int); 
      }
      if(jend!=-1   && jend!=lf.nbin)   {
	step_int+=exp(lf.lnlf[jend])*(Mlow-lf.magni[jend]);
       	if(DEBUG) printf(" %d Sumando en mag %f-%f Mlow %f   %g  : %g\n",jend,lf.magni[jend],lf.magni[jend+1],Mlow,lf.lnlf[jend],step_int); 
      }
      for(j=jstart+1;j<jend;j++) {
	step_int+=exp(lf.lnlf[j])*(lf.magni[j+1]-lf.magni[j]);
 	if(DEBUG) printf(" %d Sumando en mag %f-%f   %g  : %g\n",j,lf.magni[j],lf.magni[j+1],lf.lnlf[j],step_int); 
      }
      if(DEBUG) printf(" Al final step %g\n",step_int); 
      /* Esto que viene a continuacion es como lo hacia antes, con magni[j] el centro del intervalo j */
      /*       if(Mup<(lf.magni[0]-(lf.magni[1]-lf.magni[0])/2.)) { */
      /* 	jstart=-1; */
      /* 	step_int+=exp(lf.lf[0])*(lf.magni[1]-lf.magni[0]); */
      /*       } */
      /*       else if(Mup<(lf.magni[0]+(lf.magni[1]-lf.magni[0])/2.)) { */
      /* 	jstart=0; */
      /* 	step_int+=exp(lf.lf[0])*((lf.magni[0]+(lf.magni[1]-lf.magni[0])/2.)-Mup); */
      /*       } */
      /*       if(Mup>(lf.magni[lf.nbin-1]+(lf.magni[lf.nbin-1]-lf.magni[lf.nbin-2])/2.)) { */
      /* 	jstart=lf.nbin; */
      /*       } */
      /*       else if(Mup>(lf.magni[lf.nbin-1]-(lf.magni[lf.nbin-1]-lf.magni[lf.nbin-2])/2.)) { */
      /* 	jstart=lf.nbin-1; */
      /* 	step_int+=exp(lf.lf[lf.nbin-1])*((lf.magni[lf.nbin-1]+(lf.magni[lf.nbin-1]-lf.magni[lf.nbin-2])/2.)-Mup); */
      /*       } */
      /*       if(Mlow<(lf.magni[0]-(lf.magni[1]-lf.magni[0])/2.)) jend=-1; */
      /*       else if(Mlow<(lf.magni[0]+(lf.magni[1]-lf.magni[0])/2.)) { */
      /* 	jend=0; */
      /* 	step_int+=exp(lf.lf[0])*(Mlow-(lf.magni[0]-(lf.magni[1]-lf.magni[0])/2.)); */
      /*       }	       */
      /*       if(Mlow>(lf.magni[lf.nbin-1]+(lf.magni[lf.nbin-1]-lf.magni[lf.nbin-2])/2.)) { */
      /* 	jend=lf.nbin; */
      /* 	step_int+=exp(lf.lf[lf.nbin-1])*(lf.magni[lf.nbin-1]-lf.magni[lf.nbin-2]); */
      /*       } */
      /*       else if(Mlow>(lf.magni[lf.nbin-1]-(lf.magni[lf.nbin-1]-lf.magni[lf.nbin-2])/2.)) { */
      /* 	jend=lf.nbin-1; */
      /* 	step_int+=exp(lf.lf[lf.nbin-1])*(Mlow-(lf.magni[lf.nbin-1]-(lf.magni[lf.nbin-1]-lf.magni[lf.nbin-2])/2.)); */
      /*       } */
      /*       for(j=1;j<lf.nbin-1;j++) { */
      /* 	if((Mup>(lf.magni[j]+lf.magni[j-1])/2.) && (Mup<(lf.magni[j]+lf.magni[j+1])/2.)) jstart=j; */
      /* 	if((Mlow>(lf.magni[j]+lf.magni[j-1])/2.) && (Mlow<(lf.magni[j]+lf.magni[j+1])/2.)) jend=j; */
      /*       } */
      /*       for(j=jstart+1;j<jend-1;j++) { */
      /* 	step_int+=exp(lf.lf[j])*((lf.magni[j+1]-lf.magni[j-1])/2.); */
      /*       } */
      
      /*       Ahora calculamos dVdz, pero sin los factores, que no influyen */
      Dang_DH=(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(1+z); /* La distancia angular * (1+z) */
      E=sqrt(2*cosmo.q0*(1+z)*(1+z)*(1+z)+(1-2*cosmo.q0)*(1+z)*(1+z));
      dVdz=(Dang_DH*Dang_DH/E);
      /*       Aqui iria rho(z) si la densidad comovil variase con el z */
      /*       El z del final es porque estoy integrando en log(z), y dz=z*d(log(z)) */
      flogz[i]=step_int*dVdz*z;
      dz=exp(log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.))-exp(log(zlow)+(i-1)*(log(zup)-log(zlow))/(nstep_z-1.));
      fz[i]=step_int*dVdz;
      fz_int+=flogz[i];
      if(DEBUG2) printf(" z %g fz %g \n",z,fz[i]);
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

double Steplfdev_L(struct Steplf_L lf, double Llow,double Lup) {

  double r;
  double L=0;
  static long idum =-1;
  static int nbin=-1;
  static double *PLFcum;
  int i;
  /* Cuidado!! El bin j va desde lumi[j] a lumi[j+1], luego lumi esta dimensionado
     de nbin+1 y lf de nbin. De esta manera, el intervalo es lumi[j+1]-lumi[j] */
  if(idum==-1) idum=-(long)time(NULL)/2;
  nbin=lf.nbin;
  if(nbin!=-1) free(PLFcum);
  PLFcum=vector_d(nbin);

  PLFcum[0]=exp(lf.lnlf[0]+lf.lumi[1])-exp(lf.lnlf[0]+lf.lumi[0]);
  for(i=1;i<nbin;i++) 
  {
    PLFcum[i]=(exp(lf.lnlf[i]+lf.lumi[i+1])-exp(lf.lnlf[i]+lf.lumi[i]))+PLFcum[i-1];
  }
  /* Normalizo */
  for(i=0;i<nbin;i++)  {
    PLFcum[i]=PLFcum[i]/PLFcum[nbin-1];
    if(DEBUG) printf(" %d PLF %f lum %f-%f\n",i,PLFcum[i],lf.lumi[i]/log(10),lf.lumi[i+1]/log(10));
  }

  if(log(Lup) <lf.lumi[0])      {
    printf(" Lup <lf.lumi[0] shouldn't happen Lup %f lumi[0] %f. Exiting\n",log(Lup),lf.lumi[0]);
    exit(1);
    return(lf.lumi[0]);
  }
  if(log(Llow)>lf.lumi[lf.nbin]) {
    printf(" Llow >lf.lumi[nbin] shouldn't happen Llow %f lumi[0] %f. Exiting\n",log(Llow),lf.lumi[lf.nbin]);
    exit(1);
    return(lf.lumi[lf.nbin]);
  }
   
  do {
    r=ran2(&idum);
    for(i=0;i<nbin;i++) {
      if(r<=PLFcum[i]) {
	if(i==0) L=(r-0.)/(PLFcum[0]-0)*(exp(lf.lumi[1])-exp(lf.lumi[0]))+exp(lf.lumi[0]);
	else     L=(r-PLFcum[i-1])/(PLFcum[i]-PLFcum[i-1])*(exp(lf.lumi[i+1])-exp(lf.lumi[i]))+exp(lf.lumi[i]);
	if(DEBUG) printf(" De esta ya sale con %g  Llow %g Lup %g i %d r %f\n",L,Llow,Lup,i,r);
	break;
      }
    }
    if(DEBUG) printf(" L %g Lup %g Llow %g\n",L,Lup,Llow);
  } while (!(L>Llow && L<Lup));

  return(L);

  printf(" Steplf_Ldev: Never reaches this point\n");
  exit(1);

  return(exp(lf.lumi[nbin]));
}

double zSteplfdev_L(struct Steplf_L lf, double zlow, double zup, double ffaint, double fbright, struct cosmo_param cosmo) {
  
  int nstep_lf=1000;
  int nstep_z=1000;
  int i,j;
  int jstart,jend;
  double ri;
  double z,logz;
  double rr; 
  double Llow,Lup;
  double step_int=0;
  double Dang_DH;  /* Distancia angular segun  Hogg astroph/9905116 */
  double E;        /* Funcion E. */
  double dVdz;      /* Los factores anteriores , para calcular dVdz */
  static double fz[NSTEP_Z],fz_sum[NSTEP_Z],flogz[NSTEP_Z]; /* //fz es la funcion distribucion en z, y fz_sum es la integral */
  static double fz_int;
  double dz;
  static double zplot[NSTEP_Z],logzplot[NSTEP_Z];
  static double zupold,zlowold,ffaintold,fbrightold;
  static struct Steplf_L lfold;
  int oldflag=0;
  static long idum =-1;
  if(idum==-1) idum=-(long)time(NULL)/2;
  nstep_lf=NSTEP_LF;
  nstep_z =NSTEP_Z;

  
  if(!memcmp(&lf,&lfold,sizeof(struct Steplf_L)) && zlow==zlowold && zup==zupold && ffaint==ffaintold && fbright==fbrightold) {
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    oldflag=1;
  }
  else {
    memcpy(&lfold,&lf,sizeof(struct Steplf_L));
    zupold=zup;zlowold=zlow;ffaintold=ffaint;fbrightold=fbright;
    zlow = (zlow < ZMIN ? ZMIN : zlow);
    if(DEBUG) printf(" NUEVO\n");
    fz_int=0;
    for(i=0;i<nstep_z;i++) {
      logz=log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.);
      z=exp(logz);
      zplot[i]=z;
      logzplot[i]=logz;
      Llow=Lum(z,ffaint,cosmo);
      Lup =Lum(z,fbright,cosmo);
      if(DEBUG) printf(" %d z %f Llow %g Lup %g ffaint %g fbright %g zlow %f zup %f\n",i,z,Llow,Lup,ffaint,fbright,zlow,zup); 
      /*       Aqui integro la funcion de lumi de tipo Step */
      step_int=0;
      jstart=0;
      jend=lf.nbin;
      if(log(Lup)<lf.lumi[0])  jend=-1;
      if(log(Lup)>lf.lumi[lf.nbin])	jend=lf.nbin;
      if(log(Llow)<lf.lumi[0])  jstart=-1;
      if(log(Llow)>lf.lumi[lf.nbin])   jstart=lf.nbin;
      for(j=0;j<lf.nbin;j++) {
	if(log(Lup) >lf.lumi[j] && log(Lup) <lf.lumi[j+1]) jend  =j;
	if(log(Llow)>lf.lumi[j] && log(Llow)<lf.lumi[j+1]) jstart=j;
      }
      if(jstart!=-1 && jstart!=lf.nbin) {
	step_int+=exp(lf.lnlf[jstart]+lf.lumi[jstart+1])-exp(lf.lnlf[jstart]+log(Llow));
 	if(DEBUG) printf(" %d Sumando en mag %f-%f Lup %g   %g  : %g\n",jstart,lf.lumi[jstart],lf.lumi[jstart+1],Lup,lf.lnlf[jstart],step_int); 
      }
      if(jend!=-1   && jend!=lf.nbin)   {
	step_int+=exp(lf.lnlf[jend]+log(Lup))-exp(lf.lnlf[jend]+lf.lumi[jend]);
       	if(DEBUG) printf(" %d Sumando en mag %f-%f Llow %g   %g  : %g\n",jend,lf.lumi[jend],lf.lumi[jend+1],Llow,lf.lnlf[jend],step_int); 
      }
      for(j=jstart+1;j<jend;j++) {
	step_int+=exp(lf.lnlf[j]+lf.lumi[j+1])-exp(lf.lnlf[j]+lf.lumi[j]);
 	if(DEBUG) printf(" %d Sumando en lumi %f-%f   %g  : %g\n",j,log10(exp(lf.lumi[j])),lf.lumi[j+1]/log(10),lf.lnlf[j],step_int); 
      }
      if(DEBUG) printf(" Al final step %g\n",step_int);
      Dang_DH=(2-2*cosmo.q0*(1-z)-(2-2*cosmo.q0)*sqrt(1+2*cosmo.q0*z))/(1+z); /* La distancia angular * (1+z) */
      E=sqrt(2*cosmo.q0*(1+z)*(1+z)*(1+z)+(1-2*cosmo.q0)*(1+z)*(1+z));
      dVdz=(Dang_DH*Dang_DH/E);
      /*       Aqui iria rho(z) si la densidad comovil variase con el z */
      /*       El z del final es porque estoy integrando en log(z), y dz=z*d(log(z)) */
      flogz[i]=step_int*dVdz*z;
      dz=exp(log(zlow)+i*(log(zup)-log(zlow))/(nstep_z-1.))-exp(log(zlow)+(i-1)*(log(zup)-log(zlow))/(nstep_z-1.));
      fz[i]=step_int*dVdz;
      fz_int+=flogz[i];
      if(DEBUG2) printf(" z %g fz %g \n",z,fz[i]);
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



double Int_step_M(struct Steplf_M lf, double zlow,double zup,double mlim,struct cosmo_param cosmo) {
  double z;
  int nz,nM;
  int i,j;
  int jend;
  double Mlow;
  double N,Npar;

  zlow = (zlow < ZMIN ? ZMIN : zlow);

  nz=NSTEP_Z;
  nM=NSTEP_MAG;
  
  
  N=0;
  for(i=0;i<nz;i++) {
    z=zlow+i*(zup-zlow)/(nz-1.);
    Mlow=Mag(z,mlim,cosmo);
    /*         Esto que viene es haciendo la integral a pelo: */
    Npar=0;
    jend=lf.nbin;
    if(Mlow<lf.magni[0])  jend=-1;
    if(Mlow>lf.magni[lf.nbin])   jend=lf.nbin;
    for(j=0;j<lf.nbin;j++) if(Mlow>lf.magni[j] && Mlow<lf.magni[j+1]) jend=j;
    if(jend!=-1   && jend!=lf.nbin)   Npar+=exp(lf.lnlf[jend])*(Mlow-lf.magni[jend]);
    for(j=0;j<jend;j++) 	Npar+=exp(lf.lnlf[j])*(lf.magni[j+1]-lf.magni[j]);

    /*     printf("i  %d  Npar %f\n",i,Npar); */
    Npar=Npar*dVdz(z,cosmo)/1.e18;
    N=N+Npar;
  }
  N=N/nz*(zup-zlow);
  return(N);
}



double Int_step_L
(struct Steplf_L lf, double zlow,double zup,double fluxlim,
 struct cosmo_param cosmo) 
{
  double z;
  int nz,nM;
  int i,j;
  int jini;
  double Llow;
/*   double Lup=1e60; */
  double N,Npar;

  nz=NSTEP_Z;
  nM=NSTEP_MAG;

  N=0;
  
/*   if(Lup/lf.Lstar>100) Lup=100*lf.Lstar; */ /* Para que irse tan lejos!, si es cero */
  for(i=0;i<nz;i++) 
  {
    z=zlow+i*(zup-zlow)/(nz-1.);
    /* En luminosidades: */
/*     printf(" Ho %f q0 %f\n",cosmo.H0,cosmo.q0); */
    Llow=Lum(z,fluxlim,cosmo);
    Npar=0;
    jini=lf.nbin;
    if(log(Llow)<lf.lumi[0])  jini=-1;
    if(log(Llow)>lf.lumi[lf.nbin])   jini=lf.nbin;
    for(j=0;j<lf.nbin;j++) 
      if(log(Llow)>lf.lumi[j] && log(Llow)<lf.lumi[j+1]) jini=j;
    if(jini!=-1   && jini!=lf.nbin)   
    {
      Npar+=exp(lf.lnlf[jini]+lf.lumi[jini+1])-exp(lf.lnlf[jini]+log(Llow));
    }
    for(j=jini+1;j<lf.nbin;j++) 	
    {
      Npar+=exp(lf.lnlf[j]+lf.lumi[j+1])-exp(lf.lnlf[j]+lf.lumi[j]);
    }

    Npar=Npar*dVdz(z,cosmo)/1.e18;
    N=N+Npar;
  }
  N=N/nz*(zup-zlow);
  return(N);
}

void PrintStepLF_L(struct Steplf_L lf) 
{
  unsigned int i;
  for(i=0;i<lf.nbin;i++) 
  {
    printf(" log(Lum) %11g  LF %11g (log=%9g) Err_LF %11g (log=%9g)\n",
           lf.lumi[i]/log(10),exp(lf.lnlf[i]),lf.lnlf[i]/log(10),
           exp(lf.lnlf[i])*lf.errlnlf[i],lf.errlnlf[i]);
  }
}

void PrintStepLF_M(struct Steplf_M lf) 
{
  unsigned int i;
  for(i=0;i<lf.nbin;i++) 
  {
    printf(" Mag %11g    LF %11g (log=%9g) Err_LF %11g (log=%9g)\n",
           lf.magni[i],exp(lf.lnlf[i]),lf.lnlf[i]/log(10),
           exp(lf.lnlf[i])*lf.errlnlf[i],lf.errlnlf[i]);
  }
}

void PlotStepLF_L(struct Steplf_L lf) {

  float x, y, sigy;
  double dum1,dum2,dum3;
  double ymin,ymax;
  float xmin,xmax;
  int i;
  
  cpgpage();

  ymin=1e38;
  ymax=-1e38;
  for(i=0;i<lf.nbin;i++)
  {
    if(!(lf.lnlf[i]==-1/0.))
    {
      if((double)log10(exp(lf.lnlf[i]))>ymax) 
        ymax=(double)log10(exp(lf.lnlf[i]));
      if((double)log10(exp(lf.lnlf[i]))<ymin)
        ymin=(double)log10(exp(lf.lnlf[i]));
    }
  }
  
  ymin=ymin-2.5;
  ymax=ymax+2.5;
  pgLimits_d(lf.nbin,lf.lumi,&xmin,&xmax); 
  xmin=xmin-1.5;
  xmax=xmax+1.5;
  cpgsci(1);
  cpgsch(1.2);
  cpgswin((float)xmin/log(10),(float)xmax/log(10),(float)ymin,(float)ymax);
  cpgbox("BCTNSL",0,0,"BCTNS",0,0);
  
  
  for(i=0;i<lf.nbin;i++) 
  {
    if(!(lf.lnlf[i]==-1/0.))
    {
      x=log10(exp((lf.lumi[i]+lf.lumi[i+1])/2.));
      y=lf.lnlf[i]/log(10);
      sigy=lf.errlnlf[i]/log(10.);
      cpgpt1(x,y,17);
      dum1=x;
      dum2=y-sigy;
      dum3=y+sigy;
      /*     printf(" dum2 %f dum3 %f \n",dum2,dum3); */
      cpgerry_d(1,&dum1,&dum2,&dum3,0.35); 
    }
  }
  cpglab("Luminosity (W)","log(Phi\\dlum\\u) (#/Mpc\\u3\\d/W) ","");

}

void PlotStepLF_L_ov(struct Steplf_L lf) {

  float x, y, sigy;
  double dum1,dum2,dum3;
/*   double ymin,ymax; */
/*   float xmin,xmax; */
  static int color=2;
  static int symbol=18;
  int i;
  
  cpgpage();


  cpgsci(color);
  cpgsch(1.2);
  
  for(i=0;i<lf.nbin;i++) 
  {
    if(!(lf.lnlf[i]==-1/0.)) 
    {
      x=log10(exp((lf.lumi[i]+lf.lumi[i+1])/2.));
      y=lf.lnlf[i]/log(10);
      sigy=lf.errlnlf[i]/log(10.);
      cpgpt1(x,y,symbol);
      dum1=x;
      dum2=y-sigy;
      dum3=y+sigy;
      /*     printf(" dum2 %f dum3 %f \n",dum2,dum3); */
      cpgerry_d(1,&dum1,&dum2,&dum3,0.35); 
    }
  }

}


void PlotStepLF_M(struct Steplf_M lf) 
{

  float x, y, sigy;
  double dum1,dum2,dum3;
  double ymin,ymax;
  float xmin,xmax;
  int i;

  cpgpage();

  ymin=1e38;
  ymax=-1e38;
  printf(" El bin %d\n",lf.nbin);
  for(i=0;i<lf.nbin;i++)
  {
    printf(" dentro %g \n",lf.lnlf[i]);
    if(!(lf.lnlf[i]==-1/0.))
    {
      if((double)(log10(exp(lf.lnlf[i])))>ymax)
         ymax=(double)log10(exp((lf.lnlf[i])));
      if((double)(log10(exp(lf.lnlf[i])))<ymin) 
         ymin=(double)log10(exp((lf.lnlf[i])));
    }
  }

  ymin=ymin-2.5;
  ymax=ymax+2.5;
  
  pgLimits_d(lf.nbin,lf.magni,&xmin,&xmax);
  cpgsci(1);
  cpgsch(1.2);
  
  printf("xmin %g xmax %g\n",xmin,xmax);
    
  cpgswin((float)xmax,(float)xmin,(float)ymin,(float)ymax);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  
  
  for(i=0;i<lf.nbin;i++) 
  {
    if(!(lf.lnlf[i]==-1/0.))
    {
      x=(lf.magni[i]+lf.magni[i+1])/2.;
      y=log10(exp(lf.lnlf[i]));
      sigy=lf.errlnlf[i] / log(10.);
      cpgpt1(x,y,17);
      dum1=x;
      dum2=y-sigy;
      dum3=y+sigy;
      printf(" dum1 %f dum2 %f dum3 %f \n",dum1,dum2,dum3); 
      cpgerry_d(1,&dum1,&dum2,&dum3,0.35); 
    }
  }
  cpglab("Mag","log(Phi\\dmag\\u) (#/Mpc\\u3\\d/Mag) ","");

}



void PlotStepSchLF_L(struct Steplf_L lfstep, struct Schlf_L lfsch) 
{

  float x, y, sigy;
  double dum1,dum2,dum3;
  double ymin,ymax;
  float xmin,xmax;
  int i;

  int nline=100;

  cpgpage();
  
  ymin=1e38;
  ymax=-1e38;
  for(i=0;i<lfstep.nbin;i++) 
  {
    if(!(lfstep.lnlf[i]==-1/0.)) 
    { 
      if((double)log10(exp(lfstep.lnlf[i]))>ymax) 
        ymax=(double)log10(exp(lfstep.lnlf[i]));
      if((double)log10(exp(lfstep.lnlf[i]))<ymin) 
        ymin=(double)log10(exp(lfstep.lnlf[i]));
    }
  }
  
  ymin=ymin-2.5;
  ymax=ymax+2.5;
  pgLimits_d(lfstep.nbin,lfstep.lumi,&xmin,&xmax); 
  xmin=xmin-1.5;
  xmax=xmax+1.5;
  cpgsci(1);
  cpgsch(1.2);
  cpgswin((float)xmin/log(10),(float)xmax/log(10),(float)ymin,(float)ymax);
  cpgbox("BCTNSL",0,0,"BCTNS",0,0);
  
  
  for(i=0;i<lfstep.nbin;i++) 
  {
    if(!(lfstep.lnlf[i]==-1/0.)) 
    {
      x=log10(exp((lfstep.lumi[i]+lfstep.lumi[i+1])/2.));
      y=log10(exp(lfstep.lnlf[i]));
      sigy=lfstep.errlnlf[i]/log(10.);
      cpgpt1(x,y,17);
      dum1=x;
      dum2=y-sigy;
      dum3=y+sigy;
      /*     printf(" dum2 %f dum3 %f \n",dum2,dum3); */
      cpgerry_d(1,&dum1,&dum2,&dum3,0.35); 
    }
  }
  
  cpgsci(4);
  cpgmove(xmin,Schechter_L(pow(10.,xmin)/log(10),lfsch));
  
  for(i=0;i<nline;i++) {
    x=xmin+(xmax-xmin)*i/(nline-1);
    y=Schechter_L(exp(x),lfsch);
    x=x/log(10);
    y=y/log(10);
/*     printf("2  x %g y %g\n",x,y); */

    cpgdraw(x,y);
  }

  cpgsci(1);
  cpglab("Luminosity (W)","log\\dlum\\u(Phi) (#/Mpc\\u3\\d/W) ","");
  
}


void PlotStepSchLF_M(struct Steplf_M lfstep, struct Schlf_M lfsch)
{

  float x, y, sigy;
  double dum1,dum2,dum3;
  double ymin,ymax;
  float xmin,xmax;
  int i;
  int nline=100;

  cpgpage();

  ymin=1e38;
  ymax=-1e38;
  for(i=0;i<lfstep.nbin;i++) 
  {
    if(!(lfstep.lnlf[i]==-1/0.)) 
    {
      if((double)(log10(exp(lfstep.lnlf[i])))>ymax) 
        ymax=(double)log10(exp(lfstep.lnlf[i]));
      if((double)(log10(exp(lfstep.lnlf[i])))<ymin)
        ymin=(double)log10(exp(lfstep.lnlf[i]));
    }
  }
  
  ymin=ymin-2.5;
  ymax=ymax+2.5;
  /* dabreu: hay que mirar pgLimits_d.c y MinMax_d.c */
  pgLimits_d(lfstep.nbin,lfstep.magni,&xmin,&xmax);
  xmin=xmin-1.;
  xmax=xmax+1.;
  cpgsci(1);
  cpgsch(1.2);
  cpgswin((float)xmax,(float)xmin,(float)ymin,(float)ymax);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  
  
  for(i=0;i<lfstep.nbin;i++) 
  {
    if(!(lfstep.lnlf[i]==-1/0.)) 
    {
      x=(lfstep.magni[i]+lfstep.magni[i+1]) / 2.;
      y=lfstep.lnlf[i] / log(10);
      sigy=lfstep.errlnlf[i] / log(10.);
      cpgpt1(x,y,17);
      dum1=x;
      dum2=y-sigy;
      dum3=y+sigy;
      printf(" x %g y %g dum2 %g dum3 %g \n",x,y,dum2,dum3); 
      cpgerry_d(1,&dum1,&dum2,&dum3,0.35); 
    }
  }

  cpgsci(4);
  cpgmove(xmin,Schechter_M(xmin,lfsch));
  for(i=0;i<nline;i++) 
  {
    x=xmin+(xmax-xmin)*i/(nline-1);
    y=log10(Schechter_M(x,lfsch));
    cpgdraw(x,y);
    /* printf(" Pinto %f %g \n",x,y); */
  }
  printf(" La sch %g %g %g\n",lfsch.Mstar,lfsch.alfa,lfsch.phistar);
  
  cpgsci(1);
  cpglab("Mag","log\\d10\\u \\gF (#/Mpc\\u3\\d/Mag) ","LF and Sch fit");  

}

int  FitSch2StepLF_L
(struct Steplf_L lfstep, struct Schlf_L *lfsch, double *chisq) 
{
  int iter;
  double par[3],sigpar[3];
  double **covarpar;
  int ipar[3];
  
  double *lfx,*lfy,*lfsigy;
  int i;
  int nfit;
  double first,median,third;

  lfx   =vector_d(lfstep.nbin);
  lfy   =vector_d(lfstep.nbin);
  lfsigy=vector_d(lfstep.nbin);
  covarpar=matrix_d(3,3);
  ipar[0]=1;ipar[1]=1;ipar[2]=1;


  par[2]=0;
  nfit=0;
  for(i=0;i<lfstep.nbin;i++) 
  {
    if(!(lfstep.lnlf[i]==-1/0.)) 
    {
      lfx[nfit]=(lfstep.lumi[i+1]+lfstep.lumi[i])/2.;
      lfy[nfit]=lfstep.lnlf[i]/log(10);
      lfsigy[nfit]=lfstep.errlnlf[i]/log(10);
/*       par[2]+=exp(lfstep.lnlf[i]); */
      nfit++;
    }
  }

  par[1]=-0.8; sigpar[1]=0.3;
/*   printf(" par2 %g\n",par[2]); */
/*   par[2]=(exp(lfstep.lumi[lfstep.nbin])-exp(lfstep.lumi[0]))*par[2]/lfstep.nbin; */
/*   printf(" %g %g \n",lfstep.lumi[lfstep.nbin+1],lfstep.lumi[0]); */
  
  sigpar[2]=par[2]/10.;
  par[0]=StMedia_d(nfit,lfx,&(sigpar[0]));
  Quartil_d(nfit,lfy,&first,&median,&third);
  par[2]=exp(first*log(10)+par[0]+1);

/*   printf(" first %g par[0] %g \n",first,par[0],par[2]); */

  
    
/*   printf(" Ats par %g %g %g\n",par[0],par[1],par[2]); */

/*   par[0]=34*log(10);  */
/*   par[1]=-0.7;  */
/*   par[2]=0.0033; */
/*   printf(" nfit %d lfstep.nbin %d\n",nfit,lfstep.nbin); */
  iter=Mrq_d(lfx,lfy,lfsigy,nfit,par,ipar,3,covarpar,chisq,mrqfunc_fitsch2steplf_L);

/*   printf(" Despue s\n"); */

  lfsch->alfa=par[1];
  lfsch->erralfa=sqrt(covarpar[1][1]);
  lfsch->Lstar=exp(par[0]);
  lfsch->errLstar=sqrt(covarpar[0][0])*lfsch->Lstar;
  lfsch->phistar=par[2];
  lfsch->errphistar=sqrt(covarpar[2][2]);
  lfsch->covaralfaphistar=covarpar[1][2];
  lfsch->covaralfaLstar=covarpar[1][0]*lfsch->Lstar;
  lfsch->covarphistarLstar=covarpar[2][0]*lfsch->Lstar;

  free(lfx);free(lfy);free(lfsigy);
  free_matrix_d(covarpar,3,3);
  return(iter);
}

int  FitSch2StepLF_M
(struct Steplf_M lfstep, struct Schlf_M *lfsch, double *chisq) 
{
  int iter;
  double par[3],sigpar[3];
  double **covarpar;
  int ipar[3];
  
  double *lfx,*lfy,*lfsigy;
  int i;
  int nfit;
  double first,median,third;

  lfx   =vector_d(lfstep.nbin);
  lfy   =vector_d(lfstep.nbin);
  lfsigy=vector_d(lfstep.nbin);
  covarpar=matrix_d(3,3);
  ipar[0]=1;ipar[1]=1;ipar[2]=1;


  par[2]=0;
  nfit=0;
  for(i=0;i<lfstep.nbin;i++) 
  {
    if(!(lfstep.lnlf[i]==-1/0.)) 
    {
      lfx[nfit]=(lfstep.magni[i+1]+lfstep.magni[i])/2.;
      lfy[nfit]=lfstep.lnlf[i]/log(10);
      lfsigy[nfit]=lfstep.errlnlf[i]/log(10);
/*       par[2]+=exp(lfstep.lnlf[i]); */
      nfit++;
    }
  }

  par[1]=-0.8; sigpar[1]=0.3;
/*   printf(" par2 %g\n",par[2]); */
/*   par[2]=(exp(lfstep.lumi[lfstep.nbin])-exp(lfstep.lumi[0]))*par[2]/lfstep.nbin; */
/*   printf(" %g %g \n",lfstep.lumi[lfstep.nbin+1],lfstep.lumi[0]); */
  
  sigpar[2]=par[2]/10.;
  par[0]=StMedia_d(nfit,lfx,&(sigpar[0]));
  Quartil_d(nfit,lfy,&first,&median,&third);
  par[2]=exp(first*log(10)+1);
  
  if(DEBUG) printf(" first %f par[2] %g\n",first,par[2]);


  iter=Mrq_d(lfx,lfy,lfsigy,nfit,par,ipar,3,covarpar,chisq,mrqfunc_fitsch2steplf_M);

/*   printf(" Despue s\n"); */

  lfsch->Mstar=par[0];
  lfsch->errMstar=sqrt(covarpar[0][0]);
  lfsch->alfa=par[1];
  lfsch->erralfa=sqrt(covarpar[1][1]);
  lfsch->phistar=par[2];
  lfsch->errphistar=sqrt(covarpar[2][2]);
  lfsch->covaralfaphistar=covarpar[1][2];
  lfsch->covaralfaMstar=covarpar[1][0];
  lfsch->covarphistarMstar=covarpar[2][0];

  free(lfx);free(lfy);free(lfsigy);
  free_matrix_d(covarpar,3,3);
  return(iter);
}

void mrqfunc_fitsch2steplf_L(double x,double *p,double *y,double *dyda,int n) {
  double log_L_Lstar;

  struct Schlf_L lf;

/*   printf(" AA\n"); */

  if(p[2]<0) p[2]=-p[2];
  if(p[2]==0) p[2]=1;

  log_L_Lstar=x-p[0];
  lf.alfa=p[1];
  lf.Lstar=exp(p[0]);
  lf.phistar=p[2];  
/*   printf(" p[0] %g p[1] %g p[2] %g\n",p[0],p[1],p[2]); */

  /* Valor de la funcion Schechter evaluada en L=x */
  *y=Schechter_L(exp(x),lf)/log(10);
  dyda[0]=-(1+p[1])*log(exp(1))+exp(log_L_Lstar)*log(exp(1));
  dyda[1]=log_L_Lstar/log(10);
  dyda[2]=1./(p[2]*log(10.));
/*   printf("x %f *y %f y2 %f  p0 %f p1 %f p2 %f dyda1 %f dyda2 %f dyda3 %f\n",x,*y,y2,p[0],p[1],p[2],dyda[0],dyda[1],dyda[2]); */

/*   printf(" y %g dyda0 %g dyda1 %g dyda2 %g\n",*y,dyda[0],dyda[1],dyda[2]); */

/*   printf(" BB\n"); */

}


void mrqfunc_fitsch2steplf_M(double x,double *p,double *y,double *dyda,int n) 
{
  double log_L_Lstar;

  struct Schlf_M lf;

/*   printf(" AA\n"); */

  if(p[2]<0) p[2]=-p[2];
  if(p[2]==0) p[2]=1;

  log_L_Lstar=x-p[0];
  lf.alfa=p[1];
  lf.Mstar=p[0];
  lf.phistar=p[2];  


  /* Valor de la funcion Schechter evaluada en M=x */
  *y=log(Schechter_M(x,lf))/log(10);
  dyda[0]=0.4*(1+p[1])-0.4*pow(10.,0.4*(lf.Mstar-x));
  dyda[1]=0.4*(lf.Mstar-x);
  dyda[2]=1./(p[2]*log(10.));
}

 
double amofunc_schechterfitmag_d(int n, double *x, double *y, double *p)
{

  int i;
  double f,s;

  if(p[2]<=0) p[2]=-p[2];
  
  s=0.0;
  for(i=0; i<n; i++)
  {
    /* dabreu: Creo que ya no se usa */
    /* f=Schechter_M(x[i],p[0],p[1],p[2]); */
    /*       //    f += p[3]+p[4]*x[i]; */
      /*       //f += p[4]*x[i]; */
      /*       //printf(" f %f y %f\n",f,y[i]); */
      if(DEBUG) printf(" f %f log f %f y %f\n",f,log10(f),y[i]); 
      if (f!=0 && y[i]!=0) s += (log10(f)-y[i])*(log10(f)-y[i]);
  };

  if(DEBUG) printf(" s = %f  p0 %g p1 %g p2 %g\n",s,p[0],p[1],p[2]); 
  if(DEBUG) printf(" s = %f  p0 %g p1 %g P2 %g\n",s,p[0],p[1],p[2]); 
  
  return(s);
}


double amofunc_schechterfitlum_d(int n, double *x, double *y, double *p)
{

  int i;
  double f,s;


  if(p[2]<=0) p[2]=-p[2];
 
  s=0.0;
  for(i=0; i<n; i++)
  {
    /* dabreu: creo que ya no se usa */
    /* f=Schechter_L_LF(pow(10.,x[i]),pow(10.,p[0]),p[1],p[2]); */
    /*       //    f += p[3]+p[4]*x[i]; */
      /*       //f += p[4]*x[i]; */
      if(DEBUG2) printf("   f %f y %f\n",f,y[i]); 
      if (f!=0 && y[i]!=0)  s += (f-y[i])*(f-y[i]);
  }
  if(DEBUG2)  printf(" s = %f  p0 %g p1 %g P2 %g\n",s,p[0],p[1],p[2]); 
  
  return(s);
}
