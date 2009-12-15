#include "modulos.h"
#define YMIN ymin
#define YMAX ymax
/* //#include "cpgplot.h" */
/* Esta version 8 lo que hace es calcular los pixelss del CCD que se ven contribuidos por el obejto. Solo integra en esa zona 
 */


struct th {
/*   //Esto es para ajustar un polinomio */
  int ngrad;
  int ngradsig;
  float *coef;
  float *coefsig;
/*   //Esto es para splines ajustados */
  int nbin;
  float xmin,xmax;
  float *x;
  float *ym;
  float *sig;
  float nsig;
};


struct instr {
  float    diameter;         /* //Telescope diameter ( in meters) */
  float    focal;            /* //Telescope focal    (in meters) */
  char     filter_file[50];  /* // File with response function of filter range:0-1 */
  char     qe_file[50];      /* // File with quantum eff in the range:0-1(0%-100%) */
  char     prism_file[50];   /* //File with refraction index: Lambda  Ref_index */
  float    alfa_prism;       /* //Angle of  the prism (degrees) */
  float    texp;             /* //Exposure (seconds) */
  float    pixsize;          /* //Pixel size (microns) */
  float    gain;             /* // CCD gain (e-/ADU) */
  float    noise_e;          /* //CCD Read_out_noise RON in electrons */
  float    bias;
  float    dark;
  float    eff;              /* //Total efficiency of the system. Just a factor */
  int      xnpix;
  int      ynpix;
  int      satura;        /* //Saturation level of CCD (o if no saturation) */
};

struct field {
  float ewmin;
  float ewmax;
  int new;
  float mrmin;
  float mrmax;
  int nmr;
  float zmin;
  float zmax;
  int nz;
  float fwhmmin;
  float fwhmmax;
  int nfwhm;
  float sizemin;
  float sizemax;
  int nsize;
  float excenmin;
  float excenmax;
  int nexcen;
  float APmin;
  float APmax;
  int nAP;
  int ngrid;
  float spadist;
  int profile;
};

struct gal {
  float ew;
  float mr;
  float fwhm;
  float reff;
  float z;
  float AP;
  float excen;
  float x;
  float y;
};

struct star {
  float mr;
  float teff;
};

struct other {
  float seeing;
  float ldomin;
  float ldomax;
  int nldo;
  char  sky_file[50];
  float ldoline;
  char photband[51];
};
struct qe
{
  float *ldo;
  float *y;
  int n;
};
struct filter
{
  float *ldo;
  float *y;
  int n;
};
struct skytab
{
  float *ldo;
  float *y;
  int n;
};
struct prism
{
  float *ldo;
  float *y;
  int n;
};
struct qe Q;
struct filter F;
struct skytab Sky;
struct instr I;
struct gal G;
struct star E;
struct other O;
struct prism P;
struct field FD;
struct th TH;


char catfile[51];
char inputcat[51];
float zeropoint;
void Band2ZP();
float addgal();
float addsky();
float magsky();
void SaveOptions();
void ReadOptions(FILE *filepar);
void InputOptions();
float QE(float ldo);
float T(float ldo);
float Skyesp(float ldo);
float n(float ldo);
void Readfilter( char file[100]);
void Readqe(char file[100]);
void Readsky(char file[100]);
void Readprism(char file[100]);
int min(int x1,int x2);
int max(int x1,int x2);
float maxf(float x1,float x2);
int Disper2( float cosx, float cosy, float cocin, float alfa, float *cosxs, float *cosys );
void Numbernet2param(int number, int *jew, int *jmr, int *jz, int *jfwhm, int *jsize,int  *jexcen,int *jAP,int *jgrid);
void Numbercat2param(int number, FILE *fileicat);
void ThressholdMean(int n, float *x,float *y);
void ThressholdTheorcoc(int n, float *x,float K,float nsg);
void ThressholdTheorlogcoc(int n, float *x,float K,float nsg);
void PlotThress();
float DeltaLamb();


int main(int argc, char **argv)
{
  float hc=1.979e-25;  /*  J·m Julio x metro */

/*     float *imgp; */
/*     char fileimg[50]="salida.fits"; */
  float pi=4*atan(1.);
  float skyfot;
  FILE *filepar;
  FILE *fileocat;
  FILE *fileicat ;
  int number=0;   /* //Contador del numero de galaxia por el que vamos */
  int ntot;       /* // Numero total de galaxias simuladas. */
  int jew,jmr,jz,jfwhm,jsize,jexcen,jAP,jgrid;

  int fromcat=0;
/*   // Para saber si los parametros de las galaxias vendran del catalogo. */


  float *galcounts,*maggal,*ewgal,*cocflux,*logcocflux;

  
  float areaobj,pixelobj;

  float xmin,xmax,ymin,ymax;

  float nsg;

  float Kinst;

  float cocew20,cocew50,cocew100,cocew1000;

  float deltalambda;

  float totfotflux,totfot,totelec,totelecnoise,totcounts,fotgal;

  float fg;

  char xlabel[100];

/*   printf(" erf(1.) %f %f %f\n",erf(0.5*sqrt(2.)),erf(0.5*2.*sqrt(2.)),erf(0.5*3.*sqrt(2.))); */
/*   exit(1); */


  /********************************************
  **********************************************
  LECTURA DE PARAMETROS
  *********************
  **********************/
  srandom((unsigned int)time(NULL)/2); 
  printf("\n Welcome to SimulaOP_field Version 5\n");
  if(argc < 2) InputOptions();
  else {
    printf("This program has been called with parameter file %s\n",argv[1]);
    if((filepar=fopen(argv[1],"r"))==NULL) {
      printf("ERROR: Can't open options file %s\n",argv[1]);
      exit(1);
    }
    ReadOptions(filepar);
  }
/*   printf(" Bein %s\n",catfile); */
  if((fileocat=fopen(catfile,"w"))==NULL) {
    printf("ERROR: Can't open log file %s\n",catfile);
    exit(1);
  }
  /*     printf(" AMDOJO\n"); */
  

					    
  Band2ZP();




  fprintf(fileocat,"#Catalogue file generated by  SimulaOP_narrow with parameter file %s\n",argv[1]);
  fprintf(fileocat,"#Number       x         y        mr     EW   fwhm    z    reff     excen    PA\n");
  Readfilter(I.filter_file);
  Readqe(I.qe_file);
  Readsky(O.sky_file);
  Readprism(I.prism_file);

/*   //printf(" efeicnen   %f\n",I.eff); */
  
  /* Contribucion del cielo
  ***********************
  ********************/

  deltalambda=DeltaLamb();
/*   printf(" deltalambda %f\n",deltalambda); */


  skyfot=addsky();

  printf("Sky contribution: %e phot/s/pixel ---> %e ADUs in %f s.\n",skyfot,skyfot*I.texp*I.eff/I.gain,I.texp);



  printf(" Plate scale in the output image %f arcsec/pix\n",I.pixsize/I.focal*206265); 

  setvbuf(stdin,"",_IOLBF,0);




/*   //  printf("INPU %s \n",inputcat); */

  if(strcmp(inputcat,"NONE")) {
    fromcat=1;    
    if((fileicat=fopen(inputcat,"r"))==NULL) {
      printf("ERROR: Can't open input catalogue %s\n",inputcat);
      exit(1);
    }
    printf(" Galaxy parameters from catalogue file %s\n",inputcat);
  }
  
  if (fromcat) ntot=FileNLin(inputcat)-2;
  else ntot=(FD.ngrid*FD.new*FD.nmr*FD.nz*FD.nfwhm*FD.nsize*FD.nexcen*FD.nAP);
  
  printf("\n Total number of galaxies simulated %d\n",ntot);
/*   printf(" Using galaxy profile %d\n",FD.profile); */

  printf("********************************************\n");
  printf("******** Beginning of computation **********\n");
  
  
/*   printf(" %5d Net     %5d %5d %5d %5d %5d %5d %5d %5d \n",ntot,FD.ngrid,FD.new,FD.nmr,FD.nsize,FD.nz,FD.nfwhm,FD.nexcen,FD.nAP); */



  if((galcounts=malloc(ntot*sizeof(float)))==NULL )  { printf("I cannot dimension  of %d bytes",ntot);exit(1);} 
  if((maggal=malloc(ntot*sizeof(float)))==NULL )  { printf("I cannot dimension  of %d bytes",ntot);exit(1);} 
  if((ewgal=malloc(ntot*sizeof(float)))==NULL )  { printf("I cannot dimension  of %d bytes",ntot);exit(1);} 
  if((cocflux=malloc(ntot*sizeof(float)))==NULL )  { printf("I cannot dimension  of %d bytes",ntot);exit(1);} 
  if((logcocflux=malloc(ntot*sizeof(float)))==NULL )  { printf("I cannot dimension  of %d bytes",ntot);exit(1);} 


/*   printf(" SWDS %f\n",magsky()); */
  //exit(1);

  for(number=0;number<ntot;number++) {  
    
    printf(" %06d / %06d",number,ntot);
    
    if(fromcat) Numbercat2param(number,fileicat );
    else Numbernet2param(number,&jew,&jmr,&jz,&jfwhm, &jsize,&jexcen,&jAP,&jgrid);
    
    /*        printf(" Hasta aqui  %d at X=%f Y=%f with EW=%f m_r=%f z=%f\n",number,G.x,G.y,G.ew,G.mr,G.z); */
    /*     printf(" %5d Numeros %5d %5d %5d %5d %5d %5d %5d %5d \n",number,jgrid,jew,jmr,jsize,jz,jfwhm,jexcen,jAP); */
    
    /*         fprintf(fileocat," De aqui %6d  %6.2f  %6.2f  %5.2f  %4.1f  %4.2f %6.5f  %f  %3.1f  %5.1f\n",number,G.x,G.y,G.mr,G.ew,G.fwhm,G.z,G.reff,G.excen,G.AP); */
    /* Posicion del objeto */
    
/*     printf(" Galaxy  %d at X=%f Y=%f with EW=%f m_r=%f z=%f\n",number+1,G.x,G.y,G.ew,G.mr,G.z); */
    
    /*     exit(1); */
    fprintf(fileocat," %6d  %8.2f  %8.2f  %8.4f  %7.2f  %7.2f %6.5f  %f  %3.1f  %5.1f\n",number,G.x,G.y,G.mr,G.ew,G.fwhm,G.z,G.reff,G.excen,G.AP);



    if(ntot==1 ) printf(" Theorical number of counts from object (Efficiency = 1, Effective wavelength: %f ): %e ADU\n",(O.ldomax+O.ldomin)/2,pow(10.,-0.4*(G.mr+29.39))*(pi*I.diameter*I.diameter/4)*I.texp*I.eff*(O.ldomax-O.ldomin)*(O.ldomax+O.ldomin)/1.e10/2/hc/I.gain);
    
    
    
    
    
    /* Aqui vienen ahora las contribuciones finales, Poisson, ADU,...,bias,..*/
    
    nsg=2.;
    
    areaobj=3.1415*((O.seeing/2.35)*nsg)*((O.seeing/2.35)*nsg); /*el 2 es porque suponemos que cojemos dos sigmas del perfil gaussiano*/
    
    pixelobj=areaobj/I.pixsize/I.pixsize*I.focal*I.focal/206265/206265;
    
    /*     printf("are %f  npix %f galph %f sky %f \n",areaobj,pixelobj,addgal(),skyfot); */
    
    fotgal=addgal();

    fg=1.-exp(-nsg*nsg/2.);

    totfotflux=fg*fotgal + skyfot * pixelobj;  /* Sumamos el numero de fotones por s del cielo; */
    /*     printf(" galc 1 %f %f \n",galcounts[number],fnul); */
    totfot = totfotflux * I.texp * I.eff; /* Esta ya son fotones; */
    /*     printf(" galc 2 %f\n",galcounts[number]); */
    totelec = totfot + I.texp*I.dark*I.gain*pixelobj;
/*     printf(" galc 3 %f\n",galcounts[number]); */
    
    if(number==0) {
      Kinst=-(2.5*log10(fotgal / I.gain )+G.mr);
/*       printf(" Kinst= %f\n",Kinst); */

    }

/*     printf(" fot*fg*I.texp %f \n",fotgal*I.texp*I.eff*fg); */
/*     printf(" totele %f\n",totelec); */
/*     printf(" ruido %f\n",sqrt(totelec)); */

    if(totelec<1.e7)     totelec=Poidev((int)totelec);/*  Ruido poissoniano */
    else totelec=Gasdev()*sqrt(totelec)+totelec; /*  Ruido gaussiano; */

/*     printf(" totele %f\n",totelec);  */
    
/*      printf(" galc 4 %f\n",galcounts[number]); */
   
    totelecnoise = totelec + I.noise_e*Gasdev();     /* Anadimos ruido de lectura; */
/*     printf(" totelecnoise %f\n",totelecnoise); */
    totcounts = totelecnoise / I.gain + I.bias;                 /*   Pasamos a ADU; */
/*     printf(" totcounts %f\n",totcounts);  */
    /* Restamos el cielo */
    if(I.satura != 0)  if(totcounts>I.satura) totcounts=I.satura; /* Saturacion; */

    galcounts[number] = totcounts - I.bias - (( skyfot * I.eff + I.dark * I.gain) * pixelobj * I.texp ) / I.gain;

/*     printf(" gal[number]  %f\n",galcounts[number]);  */
    
    maggal[number]=G.mr;
    ewgal[number]=G.ew;
    //    if(galcounts[number]>=0 ) {
    logcocflux[number]=Kinst+2.5*log10(fabs(galcounts[number]/I.texp/I.eff/fg))+maggal[number];
    cocflux[number]=(galcounts[number]/I.texp/I.eff/fg)*pow(10.,0.4*(maggal[number]+Kinst));
/*     printf(" gal[number]/texp/fg  %f\n",galcounts[number]/I.texp/fg);  */

    //cocflux[number]=Gasdev()*.2;
/*     printf(" mag %f  counts %f  coc %f\n",maggal[number],galcounts[number],cocflux[number]); */
      //    }
      //else {
      //number--;
      //ntot--;
      // }

    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
  }
  printf("\n");

  deltalambda=DeltaLamb();
  TH.nbin=100;
  TH.nsig=3.;


  pgLimits(ntot,cocflux,&ymin,&ymax);
  pgLimits(ntot,maggal,&xmin,&xmax);
  cpgopen("?");
  cpgpage();


  sprintf(xlabel,"%-s magnitude",O.photband);  

  pgLimits(ntot,logcocflux,&ymin,&ymax);
  cpgswin(xmax,xmin,ymin,ymax);

  cpglab(xlabel,"-2.5 * log(Flux ratio)","Selection diagram");
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgpt(ntot,maggal,logcocflux,1);
  cocew20=log10(1+20./deltalambda);
  cocew50=log10(1+50./deltalambda);
  cocew100=log10(1+100./deltalambda);
  cocew1000=log10(1+1000./deltalambda);

  cpgsci(5);
  cpgmove(xmin,cocew20);
  cpgdraw(xmax,cocew20);
  cpgmove(xmin,cocew50);
  cpgdraw(xmax,cocew50);
  cpgmove(xmin,cocew100);
  cpgdraw(xmax,cocew100);
  cpgmove(xmin,cocew1000);
  cpgdraw(xmax,cocew1000);
  ThressholdMean(ntot,maggal,logcocflux);
  ThressholdTheorlogcoc(ntot,maggal,Kinst,nsg);


  
/*   printf(" N = %d \n",ntot); */
  //exit(1);
  cpgask(1);
  cpgpage();

  cpgsci(1);

  pgLimits(ntot,cocflux,&ymin,&ymax);
  cpgswin(xmax,xmin,deltalambda*(ymin-1),deltalambda*(ymax-1));
  cpgbox("BCTNS",0,0,"CTMS",0,0);
  cpgmtxt("R",2.5,0.5,0.5,"Anchura equivalente");

  cpgsci(5);
  cpgmove(xmin,20.);
  cpgdraw(xmax,20.);
  cpgmove(xmin,50.);
  cpgdraw(xmax,50.);
  cpgmove(xmin,100.);
  cpgdraw(xmax,100.);
  cpgmove(xmin,1000.);
  cpgdraw(xmax,1000.);
  cpgsci(1);

  cpgswin(xmax,xmin,ymin,ymax);
  cpglab(xlabel,"Flux ratio (narrow/broad)","Selection diagram");
  cpgbox("BCTNS",0,0,"BTNS",0,0);
  cpgpt(ntot,maggal,cocflux,1);

  cpgswin(xmax,xmin,ymin,ymax);
  ThressholdMean(ntot,maggal,cocflux);

  pgLimits(ntot,cocflux,&ymin,&ymax);
  cpgswin(xmax,xmin,ymin,ymax);
  cpgswin(xmax,xmin,YMIN,YMAX);
  ThressholdTheorcoc(ntot,maggal,Kinst,nsg);

  cpgend();

  return  0;
}


float addgal()
{ 

  int l; 
  float ldo; 
  float pi=4*atan(1.);
  float fcont=0;
  float halfa=O.ldoline*(1+G.z);
  double contri=0;
  float galphot=0;
  float hc=1.979e-25;  /*  J·m Julio x metro */
  float *espec;
  if((espec=malloc(O.nldo*sizeof(float)))== NULL) {printf("sin memoria\n");exit(1);}
  fcont=pow(10.,-0.4*(G.mr+31.25));
/*   printf("Continuum flux  %e W/m2/A %e erg/cm2/A  log(erg/cm2/A) %f\n",pow(10.,-0.4*(G.mr+39.39))*1.e4,pow(10.,-0.4*(G.mr+39.39))*1.e7,-0.4*(G.mr+39.39)+7); */
  
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    espec[l]=(fcont*(1+G.ew/sqrt(2.*pi)/(G.fwhm/2.35)*exp(-(ldo-halfa)*(ldo-halfa)/2./(G.fwhm/2.35)/(G.fwhm/2.35))));
    contri+=espec[l]*ldo*T(ldo)*QE(ldo);
  } 
  galphot=contri*(O.ldomax-O.ldomin)/O.nldo*1e-10; /*  1e-10 From angstroms to m */
  galphot=galphot/hc;

  galphot*=(pi*I.diameter*I.diameter/4.);
  
  return(galphot);
  
  
}





float addsky()
{
  float hc=1.979e-25;  /* // J·m Julio x metro */
  float pi=4*atan(1.);
  float skycont;
  int l;
  double ldo;
  double contri;
/*   //  printf(" Entra\n"); */
  contri=0;
  

  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
/*     //printf("Aqui tam %e\n",Skyesp(ldo)); */
/*     //printf(" Que pasa ldo %e contri %e sky %e t %f qe %f\n",ldo,contri,Skyesp(ldo),T(ldo),QE(ldo)); */

    contri=contri+Skyesp(ldo)*ldo*T(ldo)*QE(ldo);
/*     //printf(" Que pasa ldo %f contri %e sky %e t %f qe %f\n",ldo,contri,Skyesp(ldo),T(ldo),QE(ldo)); */
/*     //    printf("  %f  %e %e %f %f\n",ldo,contri,Skyesp(ldo),T(ldo),QE(ldo)); */
  }
  /* printf(" Contri %g \n",contri); */


  skycont=contri*(O.ldomax-O.ldomin)/O.nldo*1e-10; /*  //1e-10 From angstroms to m */
  skycont=skycont/hc;
/*   //skycont:  fot/m/m/s/arcsec */
/*   // Skyesp must have  W/m/m/s/Angstrom/arcsec/arcsec  unities */
/*   // ldo:  Angstroms */
/*   // T(ldo), Q(ldo) : dimensionless */

/*   //I.pixsize:  meters */
/*   //I.focal:    meters */
  skycont=skycont*(206165*I.pixsize/I.focal)*(206265*I.pixsize/I.focal); 
/*   // skycont: fot/m/m/s/pix */

  skycont=skycont*(pi*I.diameter*I.diameter/4);
/*   // skycont: fot/s/pix */

/*   printf(" Sky flux: %e  fot/s/pix\n",skycont); */
  return(skycont);
}

float magsky()
{
  int l;
  double ldo;
  double contri=0.;
  float filtro=0.;
  contri=0;
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    contri=contri+Skyesp(ldo)*T(ldo)*QE(ldo);
    filtro+=T(ldo)*QE(ldo);
  }

  contri*=(206165*I.pixsize/I.focal)*(206265*I.pixsize/I.focal);
/*   printf(" coc %g brillo %g\n",contri/filtro,-31.25-2.5*log10(contri)); */


  return(-31.25-2.5*log10(contri/filtro));
}

float DeltaLamb()
{ 
  int   l; 
  float ldo; 

  double filt=0.,lfilt=0.;

  float lc;

  

  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    filt+=T(ldo)*QE(ldo);
  } 
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    lfilt+=ldo*T(ldo)*QE(ldo);
  } 
  
  lc=lfilt/filt;

  filt*=(O.ldomax-O.ldomin)/O.nldo;

  filt/=(T(lc)*QE(lc));
/*   printf(" xfil %f %f\n",filt,lc); */
  return(filt);
  
  
}


void SaveOptions(struct instr I)
{
  FILE *fp;
  char file[100];
  char ch21[51];
  printf("Name of file: ");
  scanf("%s",file);
  if((fp=fopen(file,"w")) ==NULL) {
    printf("ERROR: Can't open file\n");
    return;
  }
  fprintf(fp,"COMMENT  Parameter file for SimulaOP                                           \n");
  fprintf(fp,"COMMENT  Please do not insert any extra character in this file. Just overwrite \n");
  fprintf(fp,"COMMENT  Instrumentation parameters                                            \n");
  fprintf(fp,"DIAMETER=%21f / Telescope diameter ( in meters)               \n",I.diameter);
  fprintf(fp,"FOCAL   =%21f / Telescope focal    ( in meters)               \n",I.focal   );
  fprintf(fp,"PRISMANG=%21f / Prism angle        ( in degrees)              \n",-I.alfa_prism*206265/3600);
  fprintf(fp,"EXPTIME =%21f / Exposure time      ( in seconds)              \n",I.texp    );
  sprintf(ch21,"'%s'",O.sky_file);
  fprintf(fp,"SKYFILE = %-20s / Sky spectrum file W/(m·m·s·Angstrom·arc·arc)  \n",ch21);
  fprintf(fp,"PIXSIZE =%21f / Pixel size of the CCD (in microns)            \n",I.pixsize*1e6 );
  fprintf(fp,"GAIN    =%21f / Gain of the CCD    ( in e/ADU  )              \n",I.gain    );
  fprintf(fp,"NOISE_E =%21f / Electron noise     ( in e/s    )              \n",I.noise_e );
  fprintf(fp,"BIAS    =%21f / Bias of the CCD    ( in   ADU  )              \n",I.bias    );
  fprintf(fp,"DARK    =%21f / Darkness current of the CCD ( in ADU/s)       \n",I.dark    );
  fprintf(fp,"SATURA  =%21d / Saturation level of CCD (0 if no saturation)  \n",I.satura  );
  fprintf(fp,"TOTEFF  =%21f / Total efficiency of the system                \n",I.eff  );
  fprintf(fp,"NAXIS1  =%21d / number of CCD pixels in X direction           \n",I.xnpix);
  fprintf(fp,"NAXIS2  =%21d / number of CCD pixels in X direction           \n",I.ynpix);
  sprintf(ch21,"'%s'",I.filter_file);
  fprintf(fp,"FILTERFI= %-20s / File with filter response                     \n",ch21);
  sprintf(ch21,"'%s'",I.qe_file);
  fprintf(fp,"QEFFFILE= %-20s / Quantum efficiency file                       \n",ch21); 
  sprintf(ch21,"'%s'",I.prism_file);
  fprintf(fp,"NFILE   = %-20s / Prism refraction index                        \n",ch21); 
  fprintf(fp,"SEEING  =%21f / Seeing of the observations ( arcseconds)      \n",O.seeing  );
  fprintf(fp,"COMMENT  Galaxy parameters                                                     \n");
  fprintf(fp,"LDOLINE =%21f / Wavelenght of emission line to simulate       \n",FD.spadist);
  fprintf(fp,"MRMIN   =%21f / Apparent  magnitude ( integrated)  Min        \n",FD.mrmin  );
  fprintf(fp,"MRMAX   =%21f / Apparent  magnitude ( integrated)  Max        \n",FD.mrmax  );
  fprintf(fp,"EWMIN   =%21f / EW of Halpha (Angstroms)            Min       \n",FD.ewmin  );
  fprintf(fp,"EWMAX   =%21f / EW of Halpha (Angstroms)            Max       \n",FD.ewmax  );
  fprintf(fp,"REFFMIN =%21f / Efective radius of the galaxy (arcsec)  Min   \n",FD.sizemin);
  fprintf(fp,"REFFMAX =%21f / Efective radius of the galaxy (arcsec)  Max   \n",FD.sizemax);
  fprintf(fp,"FWHMMIN =%21f / FWHM of the Halpha line (Angstroms)   Min     \n",FD.fwhmmin);
  fprintf(fp,"FWHMMAX =%21f / FWHM of the Halpha line (Angstroms)   Max     \n",FD.fwhmmax);
  fprintf(fp,"ZMIN    =%21f / Redshift of the galaxy                 Min    \n",FD.zmin);
  fprintf(fp,"ZMAX    =%21f / Redshift of the galaxy                 Max    \n",FD.zmax);
  fprintf(fp,"ECCMIN  =%21f / Eccentricity of the galaxy             Min    \n",FD.excenmin);
  fprintf(fp,"ECCMAX  =%21f / Eccentricity of the galaxy             Max    \n",FD.excenmax);
  fprintf(fp,"PAMAX   =%21f / Position angle (degrees)               Max    \n",FD.APmax);
  fprintf(fp,"PAMIN   =%21f / Position angle (degrees)               Min    \n",FD.APmin);
  fprintf(fp,"SPADIST =%21f / Spatial Distribution (0=random, -1=net)       \n",FD.spadist);
  fprintf(fp,"PROFILE =%21d / Spatial Profile (0=Step 1=Exp 2=r1/4 3=Gauss) \n",FD.profile);
  fprintf(fp,"PHOTBAND= %-20s / Photometric band for magnitude                \n",O.photband); 
  fprintf(fp,"COMMENT  U,B,V,R,I,J,H,K parameters                                            \n");
  fprintf(fp,"INPUTCAT= %-20s / Catalogue file. It overrides previous settings\n",ch21);
  sprintf(ch21,"'%s'",inputcat );
  fprintf(fp,"COMMENT  This catalogue file must obey the same structure than CATFILE.        \n");
  fprintf(fp,"COMMENT  Set INPUTCAT to 'NONE' if you want to use MRMIN,...,PAMIN instead.    \n");

  
  fprintf(fp,"COMMENT  Algorithm parameters. Do not change if you are not sure what you do   \n");
  fprintf(fp,"LDOMIN  =%21f / Beginning wavelenght of the algorithm         \n",O.ldomin);
  fprintf(fp,"LDOMAX  =%21f / End wavelentgh of the algorithm               \n",O.ldomax);
  fprintf(fp,"NLDO    =%21d / Number of steps in wavelentgh calculation     \n",O.nldo);
  fprintf(fp,"NEW     =%21d / Number of steps in EW grid                    \n",FD.new);
  fprintf(fp,"NMR     =%21d / Number of steps in apparente magnitude grid   \n",FD.nmr);
  fprintf(fp,"NZ      =%21d / Number of steps in redshift grid              \n",FD.nz);
  fprintf(fp,"NFWHM   =%21d / Number of steps in FWHM grid                  \n",FD.nfwhm);
  fprintf(fp,"NREFF   =%21d / Number of steps in effective radius grid      \n",FD.nsize);
  fprintf(fp,"NECCEN  =%21d / Number of steps in eccentrity grid            \n",FD.nexcen);
  fprintf(fp,"NPA     =%21d / Number of steps in position angle grid        \n",FD.nAP);
  fprintf(fp,"NGRID   =%21d / Number of galaxies in every dot of the grid   \n",FD.ngrid);
  fprintf(fp,"COMMENT  Output parameters                                                     \n");
  sprintf(ch21,"'%s'",catfile);
  fprintf(fp,"CATFILE = %-20s / Log file with information about each galaxy   \n",ch21); 
  fprintf(fp,"END                                                           ");
  fclose(fp);



  
}
  
void ReadOptions(FILE *filepar)
{
  int c=0;
/*   //printf("Estoy dentro\n"); */
  c+=f_kff(filepar,"DIAMETER",&(I.diameter));
  c+=f_kff(filepar,"FOCAL   ",&(I.focal));
  c+=f_kff(filepar,"EXPTIME",&(I.texp));
  c+=f_kff(filepar,"PRISMANG",&(I.alfa_prism));
  c+=f_kff(filepar,"PIXSIZE",&(I.pixsize));
  I.pixsize=I.pixsize/1e6;

  c+=f_kff(filepar,"GAIN",&(I.gain));
  c+=f_kff(filepar,"NOISE_E",&(I.noise_e));
  c+=f_kff(filepar,"BIAS",&(I.bias));
  c+=f_kff(filepar,"DARK",&(I.dark));
  c+=f_kff(filepar,"TOTEFF",&(I.eff));
/*   //  printf(" efeicnen   %f gainign %f\n",I.eff,I.gain); */

  c+=f_kfls(filepar,"FILTERFI",I.filter_file);
  c+=f_kfls(filepar,"QEFFFILE",I.qe_file);
  c+=f_kfls(filepar,"SKYFILE",O.sky_file);
  c+=f_kfls(filepar,"NFILE",I.prism_file);
  c+=f_kff(filepar,"LDOLINE",&(O.ldoline));
  c+=f_kff(filepar,"EWMIN",&(FD.ewmin));
  c+=f_kff(filepar,"MRMIN",&(FD.mrmin));
  c+=f_kff(filepar,"FWHMMIN",&(FD.fwhmmin));
  c+=f_kff(filepar,"REFFMIN",&(FD.sizemin));
  c+=f_kff(filepar,"ZMIN  ",&(FD.zmin));
  c+=f_kff(filepar,"ECCMIN",&(FD.excenmin));
  c+=f_kff(filepar,"PAMIN",&(FD.APmin));
  c+=f_kff(filepar,"EWMAX",&(FD.ewmax));
  c+=f_kff(filepar,"MRMAX",&(FD.mrmax));
  c+=f_kff(filepar,"FWHMMAX",&(FD.fwhmmax));
  c+=f_kff(filepar,"REFFMAX",&(FD.sizemax));
  c+=f_kff(filepar,"ZMAX  ",&(FD.zmax));
  c+=f_kff(filepar,"ECCMAX",&(FD.excenmax));
  c+=f_kff(filepar,"PAMAX",&(FD.APmax));
  c+=f_kfls(filepar,"PHOTBAND",O.photband);

  c+=f_kff(filepar,"SEEING",&(O.seeing));
  c+=f_kfi(filepar,"NAXIS1",&(I.xnpix));
  c+=f_kfi(filepar,"NAXIS2",&(I.ynpix));
  c+=f_kfi(filepar,"SATURA",&(I.satura));
  c+=f_kff(filepar,"LDOMAX",&(O.ldomax));
  c+=f_kff(filepar,"LDOMIN",&(O.ldomin));
  c+=f_kfi(filepar,"NLDO",&(O.nldo));

  c+=f_kfi(filepar,"NEW ",&(FD.new));
  c+=f_kfi(filepar,"NMR ",&(FD.nmr));
  c+=f_kfi(filepar,"NFWHM",&(FD.nfwhm));
  c+=f_kfi(filepar,"NREFF",&(FD.nsize));
  c+=f_kfi(filepar,"NZ",&(FD.nz));
  c+=f_kfi(filepar,"NECCEN",&(FD.nexcen));
  c+=f_kfi(filepar,"NPA",&(FD.nAP));
  c+=f_kfi(filepar,"NGRID",&(FD.ngrid));

/*   //FD.ngrid=1; */

  c+=f_kfls(filepar,"CATFILE",catfile);
  strcpy(inputcat,"NONE");
  c+=f_kfls(filepar,"INPUTCAT",inputcat);

  
/*   //printf("El tm  %f \n",G.reff); */
/*   //printf("El tm  %f \n",G.mr); */
/*   //printf("El   %f \n",G.mr); */
  I.alfa_prism=I.alfa_prism/206265*3600;
  
  if(c!=47) {
    printf("Not enough parameters in parameters file\n EXITING\n");
/*     printf("La prxima vez........\n"); */
    exit(1); 
  }

}

void InputOptions()
{
  printf("*****************\nParameters Input\n****************\n"); 
  printf("Telescope diameter (in meters): ");
  scanf("%f",&(I.diameter)); 
  printf("Telescope focal (in meters): ");
  scanf("%f",&(I.focal)); 
  printf("Prism angle (in degrees): ");
  scanf("%f",&(I.alfa_prism)); 
  I.alfa_prism=I.alfa_prism/206265*3600;
  printf("Exposure time (in seconds): ");
  scanf("%f",&(I.texp)); 
  printf("pixel size (in microns): ");
  scanf("%f",&(I.pixsize)); 
  I.pixsize=I.pixsize/1e6;


  printf("CCD gain:");
  scanf("%f",&(I.gain)); 
  printf("CCD read-out noise (electrons / second): ");
  scanf("%f",&(I.noise_e)); 
  printf("CCD bias: ");
  scanf("%f",&(I.bias)); 
  printf("CCD dark: ");
  scanf("%f",&(I.dark)); 
  printf("CCD pixels in X direction: ");
  scanf("%d",&(I.xnpix)); 
  printf("CCD pixels in Y direction: ");
  scanf("%d",&(I.ynpix)); 
  printf("CCD saturation level (counts): ");
  scanf("%d",&(I.satura)); 

  printf("File with quantum efficiency: ");
  scanf("%s",I.qe_file); 
  

  printf("File with filter response: ");
  scanf("%s",I.filter_file); 
  printf("File with prism refraction index: ");
  scanf("%s",I.prism_file); 
  printf("Total efficiency of the system: ");
  scanf("%f",&(I.eff)); 
  printf("Seeing of the observations: ");
  scanf("%f",&(O.seeing)); 
  printf("File with sky brightness in W/(m·m·s·Angstrom·arc·arc): ");
  scanf("%s",O.sky_file); 
  


  printf("Minimum Apparent R magnitude of the galaxy integrated: ");
  scanf("%f",&(FD.mrmin)); 
  printf("Maximum Apparent R magnitude of the galaxy integrated: ");
  scanf("%f",&(FD.mrmax)); 
  printf("Minimum Equivalent width of the Halpha line ( Angstroms): ");
  scanf("%f",&(FD.ewmin)); 
  printf("Maxmum Equivalent width of the Halpha line ( Angstroms): ");
  scanf("%f",&(FD.ewmin)); 
  printf("Efective radius of the galaxy (arcseconds): ");
  scanf("%f",&(G.reff)); 
  printf("Redshift of the galaxy: ");
  scanf("%f",&(G.z)); 


  printf("Beginning wavelenght of the algorithm: ");
  scanf("%f",&(O.ldomin)); 
  printf("End  wavelenght of the algorithm: ");
  scanf("%f",&(O.ldomax)); 
  printf("Number of steps in wavelentgh calculation: ");
  scanf("%d",&(O.nldo)); 

  printf("Output catalogue file: ");
  scanf("%s",catfile); 


  
  
  /*
    printf("EW (Å): "); 
    scanf("%f",&(G.ew)); 
    printf("%f\n",G.ew);
    printf("G.Fwhmma (Å): "); 
    scanf("%f",&(G.fwhm)); 
    printf("Continuum flux (W/arcsec2/m2/Å): "); 
    scanf("%f",&(G.mr)); 
    printf("Sky flux (W/arcsec2/m2/Å): "); 
    //  scanf("%f",&(s); 
    printf("Galaxy diameter (arcsec): "); 
    scanf("%f",&(G.reff)); 
    printf("X size of the image (arcsec): "); 
    scanf("%f",&(I.xnpix)); 
    printf("Y size of the image (arcsec): "); 
    scanf("%f",&(I.ynpix)); 
    printf("Minimum wavelenght (Å) :");
    scanf("%f",&ldomin);
    printf("Minimum wavelenght (Å) :");
    scanf("%f",&ldomax);
    printf("Number of steps in wavelenght: "); 
    scanf("%d",&nldo); 
    printf("Number of steps in X dimension: "); 
    scanf("%d",&nx); 
    printf("Number of steps in Y dimension: "); 
    scanf("%d",&ny); 
    */
}


void Readfilter( char file[100])
{
  FILE *fp;
  int i;       
/*   //  float *lx,*ly; */
  if((fp=fopen(file,"r"))==NULL) {
    printf("Cannot open file %s\n",file);
    exit(1);
  }
  F.n=FileNLin(file);
/*   //printf("El numero fe s%d   \n",F.n); */
  if((F.ldo=malloc(F.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector F.ldo of %d elements",F.n);
    exit(1);
  }
  if((F.y=malloc(F.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector F.y of %d elements",F.n);
    exit(1);
  }

  for (i=0;i<F.n;i++) {
    fscanf(fp," %f %f",F.ldo+i,F.y+i);
/*     //    F.y[i]=F.y[i]/100; */
  }  
}


void Readqe( char file[100])
{
  FILE *fp;
  int i;       
  Q.n=FileNLin(file);
  if((fp=fopen(file,"r"))==NULL) {
    printf("Cannot open file %s\n",file);
    exit(1);
  }
  if((Q.ldo=malloc(Q.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector Q.ldo of %d elements",Q.n);
    exit(1);
  }
  if((Q.y=malloc(Q.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector Q.y of %d elements",Q.n);
    exit(1);
  }
  for (i=0;i<Q.n;i++) {
    fscanf(fp," %f %f",Q.ldo+i,Q.y+i);
/*     //    Q.y[i]=Q.y[i]/100; */
  }
}
void Readsky( char file[100])
{
  FILE *fp;
  int i;       
  Sky.n=FileNLin(file);
  if((fp=fopen(file,"r"))==NULL) {
    printf("Cannot open file %s\n",file);
    exit(1);
  }
  if((Sky.ldo=malloc(Sky.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector Sky.ldo of %d elements",Sky.n);
    exit(1);
  }
  if((Sky.y=malloc(Sky.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector Sky.y of %d elements",Sky.n);
    exit(1);
  }
  for (i=0;i<Sky.n;i++) {
    fscanf(fp," %f %f",Sky.ldo+i,Sky.y+i);
    Sky.y[i]=Sky.y[i];
  }
}
void Readprism( char file[100])
{
  FILE *fp;
  int i;       
  P.n=FileNLin(file);
  if((fp=fopen(file,"r"))==NULL) {
    printf("Cannot open file %s\n",file);
    exit(1);
  }
  if((P.ldo=malloc(P.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector P.ldo of %d elements",P.n);
    exit(1);
  }
  if((P.y=malloc(P.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector P.y of %d elements",P.n);
    exit(1);
  }
  for (i=0;i<P.n;i++) {
    fscanf(fp," %f %f",P.ldo+i,P.y+i);
    P.y[i]=P.y[i];
  }
}
float T(float ldo)
{
  return(Lagr2(F.ldo,F.y,F.n,ldo));
}

float QE(float ldo)
{
  return(Lagr2(Q.ldo,Q.y,Q.n,ldo));
}

float Skyesp(float ldo)
{
  return(Lagr2(Sky.ldo,Sky.y,Sky.n,ldo));
}

float n(float ldo)
{
  return(Lagr2(P.ldo,P.y,P.n,ldo));
}







int min(int x1,int x2) {
  if(x1>x2) return(x2);
  else return(x1);
}
int max(int x1,int x2) {
  if(x1>x2) return(x1);
  else return(x2);
}

float maxf(float x1,float x2) {
  if(x1>x2) return(x1);
  else return(x2);
} 
int Disper2( float cosx, float cosy, float cocin, float alfa, float *cosxs, float *cosys ) {
  *cosxs=cosx+I.alfa_prism*(cocin-1);
  *cosys=cosy;
  return 1;
}




void Numbernet2param(int number,int *jew, int *jmr, int *jz, int *jfwhm, int *jsize,int  *jexcen,int *jAP,int *jgrid) {

  
  div_t divi;
  int nrest;
  int ntot;
  int nx,ny,i,j;
  ntot=(FD.ngrid*FD.new*FD.nmr*FD.nz*FD.nfwhm*FD.nsize*FD.nexcen*FD.nAP);
  
  divi=div(number,FD.ngrid*FD.new*FD.nmr*FD.nz*FD.nfwhm*FD.nsize*FD.nexcen);
  *jAP=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new*FD.nmr*FD.nz*FD.nfwhm*FD.nsize);
  *jexcen=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new*FD.nmr*FD.nz*FD.nsize);
  *jfwhm=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new*FD.nmr*FD.nsize);
  *jz=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new*FD.nmr);
  *jsize=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new);
  *jmr=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid);
  *jew=divi.quot;
  *jgrid=divi.rem;


  if(FD.nAP!=1) G.AP=FD.APmin+(FD.APmax-FD.APmin)**jAP/(FD.nAP-1);
  else G.AP=FD.APmin;
  if(FD.nexcen!=1) G.excen=FD.excenmin+(FD.excenmax-FD.excenmin)**jexcen/(FD.nexcen-1);
  else G.excen=FD.excenmin;
  if(FD.nfwhm!=1) G.fwhm=FD.fwhmmin+(FD.fwhmmax-FD.fwhmmin)**jfwhm/(FD.nfwhm-1);
  else G.fwhm=FD.fwhmmin;
  if(FD.nz!=1) G.z=FD.zmin+(FD.zmax-FD.zmin)**jz/(FD.nz-1);
  else G.z=FD.zmin;
  if(FD.nsize!=1) G.reff=FD.sizemin+(FD.sizemax-FD.sizemin)**jsize/(FD.nsize-1);
  else G.reff=FD.sizemin;
  if(FD.nmr!=1) G.mr=FD.mrmin+(FD.mrmax-FD.mrmin)**jmr/(FD.nmr-1);
  else G.mr=FD.mrmin;
  if(FD.new!=1) G.ew=FD.ewmin+(FD.ewmax-FD.ewmin)**jew/(FD.new-1);
  else G.ew=FD.ewmin;


  if(FD.spadist==0) {
    G.x=(float)I.xnpix*random()/2147483647.;
    G.y=(float)I.ynpix*random()/2147483647.;
  }
  else if(FD.spadist==-1) {
/*     // Determino nx,ny, las dimensiones de la malla como las mas cercanas  */
/*     // a un cuadrado. */
    nx=(int)(sqrt(ntot));
    ny=(int)(ntot/nx)+1;
/*     // Calculo los indices de esta  galaxia dentro de la red. Van desde i=1,,nx */
    j=(int)(number/nx)+1;
    i=number-(j-1)*nx+1;
/*     //Calculo las posiciones para que esten centradas y le pongo sumo */
/*     //una unidad aleatoria para que el centro caiga */
/*     // en una fraccion de pixel diferente */
    G.x=(float)(I.xnpix*(i-0.5)/nx)+random()/2147483647-.5;
    G.y=(float)(I.ynpix*(j-0.5)/ny)+random()/2147483647-.5;

/*     printf(" Position %d %d    %d  %d %f %f \n",nx,ny,i,j,G.x,G.y); */
  }
  else {
    printf(" Spatial distribution value %f not coherent\n",FD.spadist);
    exit(1);
  }





}


void Numbercat2param(int number, FILE *fileicat) {
  char cnul='#';
  char snul[2000];
  int inul;
  while(cnul=='#') {
    fgets(snul,2000,fileicat);
/*     //    printf("Leyendo dentro:  %s\n",snul); */
    cnul=snul[0];
/*     //printf(" cara %c\n",cnul); */
  }
/*   //printf("Leyendo:  %s\n",snul); */
  sscanf(snul," %d %f %f %f %f %f %f %f %f %f ",&inul,&G.x,&G.y,&G.mr,&G.ew,&G.fwhm,&G.z,&G.reff,&G.excen,&G.AP);

}



void ThressholdMean(int n,float *x, float *y)
{
  float yobjth,sigobjth;
  int j,ib; 
  float *xbuf,*ybuf,*nsig;
  int nbuf;
  float xdelt;
  float nrej=10;
  int inul;


  if((xbuf=malloc(n*sizeof(float)))==NULL) printf("I cannot dimension xbuf   of %d elements \n",n);
  if((ybuf=malloc(n*sizeof(float)))==NULL) printf("I cannot dimension ybuf   of %d elements \n",n);
  if((nsig=malloc(n*sizeof(float)))==NULL) printf("I cannot dimension nsig   of %d elements \n",n);
  for(j=0;j<n;j++) {
    nsig[j]=0;
  }


  if((TH.x    =malloc(TH.nbin*sizeof(float)))==NULL) printf("I cannot dimension x        of %d elements \n",TH.nbin);
  if((TH.ym   =malloc(TH.nbin*sizeof(float)))==NULL) printf("I cannot dimension ym       of %d elements \n",TH.nbin);
  if((TH.sig  =malloc(TH.nbin*sizeof(float)))==NULL) printf("I cannot dimension sig      of %d elements \n",TH.nbin);
  MinMax(n,x,&(TH.xmin),&(TH.xmax));
  
  /*     printf(" cmin %f xmax %f\n",TH.xmin,TH.xmax); */
  for(ib=0;ib<TH.nbin;ib++) TH.x[ib]=(TH.xmax-TH.xmin)*(ib+.5)/TH.nbin+TH.xmin;
  xdelt=(TH.xmax-TH.xmin)/TH.nbin;
  
  for(ib=0;ib<TH.nbin;ib++) {
    nbuf=0;
    for(j=0;j<n;j++) {
      
      if(x[j]>=TH.x[ib]-xdelt/2 && x[j]<TH.x[ib]+xdelt/2 && nsig[j]<nrej) {
	xbuf[nbuf]=x[j];
	ybuf[nbuf]=y[j];
	nbuf++;
      }
      
    }

    /* for(j=0;j<nbuf;j++) printf(" %d y %f\n",j,ybuf[j]); */
    
    /*       printf(" nbuf= %d\n",nbuf); */
    TH.ym[ib]=StMedia(nbuf,ybuf,&(TH.sig[ib]));
    /* TH.ym[ib]=1.; */
/*     printf(" x %f y %f sig %f\n",TH.x[ib],TH.ym[ib],TH.sig[ib]);  */
    
  }
  
  /*     Calculo nsig para todos lo objetos */
  inul=0;
  for(j=0;j<n;j++) {
    
    yobjth=Lagr2(TH.x,TH.ym,TH.nbin,x[j]);
    sigobjth=Lagr2(TH.x,TH.sig,TH.nbin,x[j]);
    nsig[j]=(y[j]-yobjth)/sigobjth;
    /* printf(" %d y %f nsig %f\n",n,y[j],nsig[j]); */
    
    if(fabs(nsig[j])<1.) inul++;
    
  }
  
/*   printf(" inul %d tanto %f\n",inul,((float)inul)/n); */
  PlotThress();
  
  
  
  
  free(xbuf);free(ybuf);
}


void PlotThress(){
  int ib;
  float *sigd,*sigu;
  
  if((sigd=malloc(100*sizeof(float)))==NULL) printf("I cannot dimension sigd   of %d elements \n",100);
  if((sigu=malloc(100*sizeof(float)))==NULL) printf("I cannot dimension sigu   of %d elements \n",100);
  
  for(ib=0;ib<TH.nbin;ib++) {
    sigu[ib]=TH.ym[ib]+TH.nsig*TH.sig[ib];
    sigd[ib]=TH.ym[ib]-TH.nsig*TH.sig[ib];  
/*     printf(" x %f y %f sig %f\n",TH.x[ib],TH.ym[ib],TH.sig[ib]);  */
  }
/*   //printf(" nbin %d\n",thress.nbin); */
  
  cpgsci(2);
  cpgline(TH.nbin,TH.x,TH.ym);
  cpgsci(4);
  cpgline(TH.nbin,TH.x,sigu);
  cpgline(TH.nbin,TH.x,sigd);
  cpgsci(1);
  free(sigu);free(sigd);
  
}



void ThressholdTheorcoc(int n, float *x,float K,float nsg)
{
  int j; 
  float *xteor,*yteor,*nteorsig,*ysig;
  float xmin,xmax;
  int npt=100;


  float SJ;
  float J;

  float fg, po;

  float areaobj;

/*   printf(" Entroa en toe\n"); */


  SJ=magsky();
  areaobj=3.1415*((O.seeing/2.35)*nsg)*((O.seeing/2.35)*nsg); /*el 2 es porque suponemos que cojemos dos sigmas del perfil gaussiano*/
  po=areaobj/I.pixsize/I.pixsize*I.focal*I.focal/206265/206265;

  fg=erf(nsg/sqrt(2.));
  
  if((xteor=malloc(npt*sizeof(float)))==NULL) printf("I cannot dimension xteor   of %d elements \n",n);
  if((yteor=malloc(npt*sizeof(float)))==NULL) printf("I cannot dimension yteor   of %d elements \n",n);
  if((ysig=malloc(npt*sizeof(float)))==NULL) printf("I cannot dimension ysig   of %d elements \n",n);
  if((nteorsig=malloc(npt*sizeof(float)))==NULL) printf("I cannot dimension nteorsig   of %d elements \n",n);

  MinMax(n,x,&(xmin),&(xmax));
  
  
  cpgsci(3);

/*   printf(" K %f SJ %f fg %f po %f\n",K,SJ,fg,po); */


  for(j=0;j<npt;j++) {
    
    xteor[j]=xmin+j*(xmax-xmin)/npt;
    J=xteor[j];
    yteor[j]=1.;
    nteorsig[j]=3.*sqrt(I.noise_e*I.noise_e/I.gain+(I.texp*I.eff*(I.dark*po+pow(10.,-0.4*K)*(pow(10.,-0.4*J)*fg+pow(10.,-0.4*SJ)*po))));
/*     printf(" cac %f  %f\n",I.texp*(pow(10.,-0.4*K)*(pow(10.,-0.4*J)*fg)),I.texp*pow(10.,-0.4*K)*pow(10.,-0.4*SJ)*po); */
    //nteorsig[j]/=(I.texp*(pow(10.,-0.4*K)*(pow(10.,-0.4*J)*fg)));
    //nteorsig[j]/=(I.texp*pow(10.,-0.4*(K+J)));
    ysig[j]=yteor[j]+nteorsig[j]/(sqrt(I.gain)*fg*I.texp*I.eff*pow(10.,-0.4*(K+J)));
/*     printf(" xteor[j] %f sig %g y %g\n",J,nteorsig[j],ysig[j]); */
  }
  cpgline(npt,xteor,ysig);
  for(j=0;j<npt;j++) {
    
    xteor[j]=xmin+j*(xmax-xmin)/npt;
    J=xteor[j];
    yteor[j]=1.;
    nteorsig[j]=-3.*sqrt(I.noise_e*I.noise_e/I.gain+(I.texp*I.eff*(I.dark*po+pow(10.,-0.4*K)*(pow(10.,-0.4*J)*fg+pow(10.,-0.4*SJ)*po))));
/*     printf(" cac %f  %f\n",I.texp*(pow(10.,-0.4*K)*(pow(10.,-0.4*J)*fg)),I.texp*pow(10.,-0.4*K)*pow(10.,-0.4*SJ)*po); */
    ysig[j]=yteor[j]+nteorsig[j]/(sqrt(I.gain)*fg*I.texp*I.eff*pow(10.,-0.4*(K+J)));
  }
  cpgline(npt,xteor,ysig);

}


void ThressholdTheorlogcoc(int n, float *x,float K,float nsg)
{
  int j; 
  float *xteor,*yteor,*nteorsig;
  float xmin,xmax;
  int npt=100;


  float SJ;
  float J;

  float fg, po;

  float areaobj;

/*   printf(" Entroa en toe\n"); */


  SJ=magsky();
  areaobj=3.1415*((O.seeing/2.35)*nsg)*((O.seeing/2.35)*nsg); /*el 2 es porque suponemos que cojemos dos sigmas del perfil gaussiano*/
  po=areaobj/I.pixsize/I.pixsize*I.focal*I.focal/206265/206265;

  fg=erf(nsg/sqrt(2.));
  
  if((xteor=malloc(npt*sizeof(float)))==NULL) printf("I cannot dimension xbuf   of %d elements \n",n);
  if((yteor=malloc(npt*sizeof(float)))==NULL) printf("I cannot dimension ybuf   of %d elements \n",n);
  if((nteorsig=malloc(npt*sizeof(float)))==NULL) printf("I cannot dimension nsig   of %d elements \n",n);

  MinMax(n,x,&(xmin),&(xmax));
  
  
  cpgsci(3);

  for(j=0;j<npt;j++) {
    
    xteor[j]=xmin+j*(xmax-xmin)/npt;
    J=xteor[j];
    yteor[j]=0.;
    nteorsig[j]=3.*sqrt(1./K*(pow(10.,0.4*J)/fg+pow(10.,0.8*J-0.4*SJ)*po/fg/fg));
  }
  cpgline(npt,yteor,nteorsig);
  for(j=0;j<npt;j++) {
    
    xteor[j]=xmin+j*(xmax-xmin)/npt;
    J=xteor[j];
    yteor[j]=0.;
    nteorsig[j]=-3.*sqrt(1./K*(pow(10.,0.4*J)/fg+pow(10.,0.8*J-0.4*SJ)*po/fg/fg));
  }
  cpgline(npt,yteor,nteorsig);

}



void Band2ZP() {

  if(!strcmp(O.photband,"U")) zeropoint=28.36;
  else if(!strcmp(O.photband,"B")) zeropoint=27.97;
  else if(!strcmp(O.photband,"V")) zeropoint=28.52;
  else if(!strcmp(O.photband,"R")) zeropoint=29.39;
  else if(!strcmp(O.photband,"I")) zeropoint=29.66;
  else if(!strcmp(O.photband,"J")) zeropoint=31.25;
  else if(!strcmp(O.photband,"H")) zeropoint=32.35;
  else if(!strcmp(O.photband,"K")) zeropoint=33.50;
  else {
    printf(" ERROR: Photometric band %-s not defined\n",O.photband);
    printf(" Photometric bands defined are: U B V R I J H K. \n EXITING\n");
    exit(1);
  }
}
