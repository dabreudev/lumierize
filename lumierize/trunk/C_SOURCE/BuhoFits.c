#include "modulos.h"

#define MAXPIX 10000
#define MAXOBJECTS 50
#define nvarsig 15

#define DEBUG 0

static int icompare(x1,x2)
     float *x1,*x2; 
{  
  if(*x1  < *x2) return(-1);
  if(*x1 == *x2) return(0);
  if(*x1  > *x2) return(+1);
  return(0);
}


struct obj {
  int area;
  int nb;
  float xobj[MAXPIX];
  float yobj[MAXPIX];
  int ipix[MAXPIX];
  int jpix[MAXPIX];
  float f[MAXPIX];
  float xb[MAXPIX];
  float yb[MAXPIX];
  float sky,sig;

  float xcpix,ycpix;
  float xpv,ypv;
 
  float xc,yc,ft,apfot1,apfot2,apfot25,fixapfot;
  float xcelip,ycelip,elipa,elipb,elipt;
  float mx2,my2,mxy,rkron;

  int heritage;
  int nheritaged;
  int herflag;
  struct obj *objher;
  int iobj,isig;
  int ideb;

  int deblenflag;
};


struct objerr {
  float errsky,errsig;

  float errf[MAXPIX];
  float errxcpix,errycpix;
  float errxc,erryc,errft,errarea,errapfot1,errapfot2,errapfot25,errfixapfot;
  float errelipa,errelipb,errelipt;
  float errmx2,errmy2,errmxy,errrkron;
  
  struct objerr *errobjher;
};




int mierdaflag=0;

int errflag;
int deblenflag=0;
int minpix;

int Detecta(float *v, struct headfits *hv, int i, int j, struct obj *object);

int errDetecta(float *v, float *errv,struct headfits *hv, int i, int j, struct obj *object, struct objerr *errobject);

void CFM(struct obj *object);

void errCFM(struct obj *object, struct objerr *errobject);

void AjustaElipse(struct obj *object);

void ApertureFot(float *ima, int nx, int ny,float rfixap, struct obj *object);

void errApertureFot(float *ima, float *errima,int nx, int ny, float rfixap, struct obj *object,struct objerr *errobject);

int Deblending(float *ima, struct headfits hv, float *imasky, struct headfits hsky,float *imasig, struct headfits hsig,struct obj objprimal, float nsig,struct obj **objsdeb, int *ndeb);

int errDeblending(float *ima, float *errima, struct headfits hv, float *imasky, struct headfits hsky,float *imasig, struct headfits hsig, float *imaerrsky, struct headfits herrsky, struct obj objprimal, float nsig,struct obj **objsdeb, struct objerr **errobjsdeb, int *ndeb); 

int LookParent(struct obj *obj,struct obj **objparent);

int errLookParent(struct obj *obj, struct objerr *objerr, struct obj **objparent,struct objerr **objerrparent);
    
int main(int argc, char **argv)

{
/*   float xc,yc,ft,apfot1,apfot2,apfot25,fixapfot; */
/*   float xcelip,ycelip,elipa,elipb,elipt; */
/*   float mx2,my2,mxy,rkron; */

/*   float errxc,erryc,errft,errarea,errapfot1,errapfot2,errapfot25,errfixapfot; */
/*   float errelipa,errelipb,errelipt; */
/*   float errmx2,errmy2,errmxy,errrkron; */
/*   int area,nb; */
/*   float xcpix,ycpix; */
/*   float errxcpix,errycpix; */
/*   float xpv,ypv; */ /* //El centro del objeto en cood pixel de windetec */
  /*  Todo esto esta ahora metido en la estructura obj */
  struct obj object;
  struct objerr errobject;

  struct obj *objsdeb;
  struct objerr *errobjsdeb;
  int ndeb;
  int ideb;


  FILE *fpout,*fpar;
/*   FILE *fperrout; */
  float *imasky,*imasig,*imaerrsky=NULL,*linea,*errlinea=NULL,x,y,sky,sig,errsky;
/*   float *imaerrin,*errlinea; */
  float *windetec,*winerrdetec=NULL,*winfot,*winerrfot=NULL;
  char fitsin[19],fitssky[19],fitssig[19],fitserrin[19],fitserrsky[19],device[19];
  char fileout[19],fileerrout[19];
  struct headfits hin, hsky, hsig, herrin, herrsky, hv;
  int vy,nobj=0,idx,liny;
  int c=0,x1,x2,y1,y2,i,j;
  int jdetect;
  long blc[2],trc[2],incr[2];

  float xmin,xmax,ymin,ymax;
  float tr[6]={0,1,0,0,0,1};
  float nsig;

/*   float xb[MAXPIX],yb[MAXPIX]; */

  int pgbegx,pgbegy,pgask;
  float gvpx1,gvpy1,gvpx2,gvpy2,pxmin,pxmax,pymin,pymax;
  char graphmod,chlab[30];
  float sn;
  float rfixap;

  /* Variables para el FITS  */
  fitsfile *infits,*errinfits,*skyfits,*sigfits,*errskyfits;
  long naxes_in[2],naxes_errin[2];
  int status=0;
  int nfound, anynull;

  long  fpixel=1;
  float  nullval;




/*   char cnul; */
/* // ------------------- */
/* // Inicio del programa */
/* // ------------------- */
/*   //system("clear"); */
  printf("This is BuhoFits  ...\n\n");
/*   //system("banner BuhoFits"); */
/*   //SlpTime(stdout); */
  if(argc > 2 ) {
    printf("Use: BuhoFits [parameter_file]\n\n");
    exit(1);
  }

  /*   printf(" Nuero par %d \n",argc); */
  /*  --------------------- */
  /*  Entrada de parametros */
  /*  --------------------- */
  
		  
  if(argc == 2 ) {
    

    if((fpar=fopen(argv[1],"r")) == NULL) {
      fprintf(stderr,"BuhoFits: ERROR. No such parameter file %s\n", argv[1]);
      exit(1);
    }


    c += f_kfs(fpar,"IMAGE   =",fitsin);
    c += f_kfs(fpar,"SKY     =",fitssky);
    c += f_kfs(fpar,"SIG     =",fitssig);
    c += f_kfs(fpar,"ERRIMAGE=",fitserrin);
    c += f_kfs(fpar,"ERRSKY  =",fitserrsky);
    c += f_kfi(fpar,"VY      =",&vy);
    c += f_kfs(fpar,"DEVICE  =",device);
    c += f_kff(fpar,"NSIG    =",&nsig);
    c += f_kff(fpar,"SNWHOLE =",&sn);
    c += f_kfi(fpar,"MINPIX  =",&minpix);
    c += f_kfs(fpar,"OUTCAT  =",fileout);
/*     c += f_kfs(fpar,"ERRCAT  =",fileerrout); */
    c += f_kfc(fpar,"GRAPHMOD=",&graphmod);
    c += f_kfi(fpar,"PGBEGX  =",&pgbegx);
    c += f_kfi(fpar,"PGBEGY  =",&pgbegy);
    c += f_kfi(fpar,"PGASK   =",&pgask);
    c += f_kff(fpar,"RFIXAP  =",&rfixap);

    if(strcmp(fitserrin,"NONE"))    errflag=1;
    else                            errflag=0;
    
    
    fclose(fpar);
    if(c != 16) {
      fprintf(stderr,"BuhoFits: ERROR. No he encontrado todos los parametros ");
      fprintf(stderr,"en BuhoFits.par (%d/16)\n",c);
      exit(1);
    }
  }
  else if(argc == 1 ) {
    printf(" Input FITS file: ");
    scanf("%s",fitsin);
    printf(" Input sky fit FITS file: ");
    scanf("%s",fitssky);
    printf(" Input sigma fit FITS file: ");
    scanf("%s",fitssig);
    printf(" Input error image for image (NONE=no errors will be computed): ");
    scanf("%s",fitserrin);
    if(strcmp(fitserrin,"NONE")) {
      errflag=1;
      printf(" Input error image for sky image: ");
      scanf("%s",fitserrsky);
    }
    else {
      errflag=0;
    }
    
    printf(" Input window size in Y direction ( big enough to hold an object: ");
    vy=readi(150);
    printf(" Number of sigmas above the sigma of the sky to detect objects: ");
    nsig=readf(3.);
    printf(" Minimum signal-to-noise ratio of the whole object: ");
    sn=readf(sn);
    printf(" Minimum number of pixels of an object: ");
    minpix=readi(10);
    printf(" Radius in pixels to compute fixd aperture photometry: ");
    rfixap=readf(3.);
    printf(" Output file with detected objects: ");
    scanf("%s",fileout);
    if(errflag) {
      printf(" Output file with errors in object parameters: ");
      scanf("%s",fileerrout);
    }
    strcpy(device,"?");
    graphmod='F';
    pgask=0;
  }



/* // ------------------------- */
/* // Apertura del modo grafico */
/* // ------------------------- */
  if(graphmod == 'C') {
    cpgbeg(0,device,pgbegx,pgbegy);
    cpgask(pgask);
    }
  if(graphmod == 'F') {
    cpgbeg(0,device,1,1);
    cpgask(pgask);
    }

/* // ----------------------------------------------------- */
/* // Inicializacion de variables y lectura de las imagenes */
/* // ----------------------------------------------------- */
  printf("Leyendo la imagen de fondo de cielo ...\n");
  if(fitsrh(fitssky,&hsky)) {   
    fprintf(stderr,"BuhoFits: ERROR. Error leyendo la imagen %s\n",fitssky);
    exit(1);
  }
  printf(" Salio satisfactorio\n"); 

  if( ffopen(&skyfits, fitssky, READONLY, &status)) fits_report_error(stderr,status);
  if((imasky=malloc(hsky.naxis1*hsky.naxis2*sizeof(float)))==NULL )  { printf("I cannot dimension imaksy of %d bytes",hsky.naxis1*hsky.naxis2);exit(1);} 
  if(fits_read_img(skyfits, TFLOAT, fpixel, hsky.naxis1*hsky.naxis2, &nullval, imasky, &anynull, &status )) fits_report_error(stderr,status);

  fitsinfo(fitssky,stdout);
  printf("-----------------------------------------\n\n");

  printf("Leyendo la imagen de fondo de sigma ...\n");
  if(fitsrh(fitssig,&hsig)) {
    fprintf(stderr,"BuhoFits: ERROR. Error leyendo la imagen %s\n",fitssig);
    exit(1);
  }
  if( ffopen(&sigfits, fitssig, READONLY, &status)) fits_report_error(stderr,status);
  if((imasig=malloc(hsig.naxis1*hsig.naxis2*sizeof(float)))==NULL )  { printf("I cannot dimension imaksy of %d bytes",hsig.naxis1*hsig.naxis2);exit(1);} 
  if(fits_read_img(sigfits, TFLOAT, fpixel, hsig.naxis1*hsig.naxis2, &nullval, imasig, &anynull, &status )) fits_report_error(stderr,status);
  
  fitsinfo(fitssig,stdout);
  printf("-----------------------------------------\n\n");
  
  if(errflag) {
    printf("Leyendo la imagen de error de fondo de cielo ...\n");
    if(fitsrh(fitserrsky,&herrsky)) {
      fprintf(stderr,"BuhoFits: ERROR. Error leyendo la imagen %s\n",fitserrsky);
      exit(1);
    }
    if( ffopen(&errskyfits, fitserrsky, READONLY, &status)) fits_report_error(stderr,status);
    if((imaerrsky=malloc(herrsky.naxis1*herrsky.naxis2*sizeof(float)))==NULL )  { printf("I cannot dimension imaksy of %d bytes",herrsky.naxis1*herrsky.naxis2);exit(1);} 
    if(fits_read_img(errskyfits, TFLOAT, fpixel, herrsky.naxis1*herrsky.naxis2, &nullval, imaerrsky, &anynull, &status )) fits_report_error(stderr,status);
    
    if (hsky.naxis1 != herrsky.naxis1 || hsky.naxis2 != herrsky.naxis2 ||
	hsky.crval1 != herrsky.crval1 || hsky.crval2 != herrsky.crval2 ||
	hsky.cdelt1 != herrsky.cdelt1 || hsky.cdelt2 != herrsky.cdelt2 ) {
      fprintf(stderr,"Parece que los descriptores de las imagenes de cielo ");
      fprintf(stderr,"y de error de fondo de cielo no coinciden \n");
      exit(1);
    }    
  }
  
  if (hsky.naxis1 != hsig.naxis1 || hsky.naxis2 != hsig.naxis2 ||
      hsky.crval1 != hsig.crval1 || hsky.crval2 != hsig.crval2 ||
      hsky.cdelt1 != hsig.cdelt1 || hsky.cdelt2 != hsig.cdelt2 ) {
    fprintf(stderr,"Parece que los descriptores de las imagenes de cielo ");
    fprintf(stderr,"y sigma no coinciden \n");
    exit(1);
    }


/* // Limites para la interpolacion del fondo y sigma */
/* // ----------------------------------------------- */

/* 	// KIKE Estos primeros son los que quiero. */
  xmin=hsky.crval1-hsky.cdelt1/2.;
  xmax=hsky.crval1+hsky.cdelt1*(hsky.naxis1-.5);
  ymin=hsky.crval2-hsky.cdelt2/2.;
  ymax=hsky.crval2+hsky.cdelt2*(hsky.naxis2-.5);
  
/*   //  xmin=hsky.crval1; */
/*   //xmax=hsky.crval1+hsky.cdelt1*(hsky.naxis1-1.); */
/*   //ymin=hsky.crval2; */
/*   //ymax=hsky.crval2+hsky.cdelt2*(hsky.naxis2-1.); */
  
  printf("Limits by sky and sigma images:\n");
  printf("X: %f - %f\n",xmin,xmax);
  printf("Y: %f - %f\n\n",ymin,ymax);

  
/* // ----------------------------------------------------- */
/* // Lectura de la primera windetec de la imagen de entrada */
/* // ----------------------------------------------------- */
  printf("Leyendo la imagen de entrada ...\n");
  fitsrh(fitsin,&hin);
  fitsinfo(fitsin,stdout);
  

  if(errflag) {  
    fitsrh(fitserrin,&herrin);
    if (hin.naxis1 != herrin.naxis1 || hin.naxis2 != herrin.naxis2 ||
	hin.crval1 != herrin.crval1 || hin.crval2 != herrin.crval2 ||
	hin.cdelt1 != herrin.cdelt1 || hin.cdelt2 != herrin.cdelt2 ) {
      fprintf(stderr,"Parece que los descriptores de la imagen de entrada ");
      fprintf(stderr,"y de error no coinciden \n");
      exit(1);
    }
  }
  

/*   //Esto lo cambio por la cagada de la WFC: Asi comentado esta bien.  */
/*   hin.bscale=1.0; */
/*   hin.bzero=0; */
/*   hin.cdelt1=1.0; */
/*   hin.cdelt2=1.0; */
/*   hin.crval1=0.0; */
/*   hin.crval2=0.0; */

/*   //printf(" SKY: crval1 %f cdelt1 %f crval2 %f cdelt2 %f\n",hsky.crval1,hsky.cdelt1,hsky.crval2,hsky.cdelt2); */
/*   //printf(" SKY: crval1 %f cdelt1 %f crval2 %f cdelt2 %f\n",hsig.crval1,hsig.cdelt1,hsig.crval2,hsig.cdelt2); */
  printf(" IMAGE: crval1 %f cdelt1 %f crval2 %f cdelt2 %f\n",hin.crval1,hin.cdelt1,hin.crval2,hin.cdelt2);

/*   // exit(1); */



  printf("Leyendo la primera windetec ....\n");
  x1=1;
  x2=hin.naxis1;
  y1=1;
  y2=vy;
  windetec=fitsrfb(fitsin,&hin,&x1,&x2,&y1,&y2);
  if(errflag) winerrdetec=fitsrfb(fitserrin,&herrin,&x1,&x2,&y1,&y2);
/*   //Esto esta comentado por lo de la WFC. Asi esta bien para todo. */
/*  hin.bscale=1.0;
  hin.bzero=0;
  hin.cdelt1=1.0;
  hin.cdelt2=1.0;
  hin.crval1=0.0;
  hin.crval2=0.0; */
/*   //Dimensiono otra ventana identica para poder hacer fotometria sobre ella */
/*   // y copio una sobre la otra antes de que a windetect se le hagan cambios */
  if( (winfot=malloc((x2-x1+1)*(y2-y1+1)*sizeof(float)))== NULL) {
    printf(" ERROR: Cannot dimension winfot of %d bytes\n",(x2-x1+1)*(y2-y1-1)*sizeof(float));
    exit(1);
  }
  memcpy(winfot,windetec,(x2-x1+1)*(y2-y1+1)*sizeof(float));

  if(errflag) {
    if( (winerrfot=malloc((x2-x1+1)*(y2-y1+1)*sizeof(float)))== NULL) {
      printf(" ERROR: Cannot dimension winerrfot of %d bytes\n",(x2-x1+1)*(y2-y1-1)*sizeof(float));
      exit(1);
    }
    memcpy(winerrfot,winerrdetec,(x2-x1+1)*(y2-y1+1)*sizeof(float));
  }

/* // Unos Test */
/* // --------- */
  if(x1 != 1 || y1 != 1 || x2 != hin.naxis1) {
    fprintf(stderr,"BuhoFits: ERROR. Algo raro falla en la primera lectura\n");
    fprintf(stderr,"          Se ha redimensionado la windetec\n");
    fprintf(stderr,"          x1=%d, x2=%d, y1=%d, y2=%d\n",x1,x2,y1,y2);
    exit(1);
  }
  if(y2 != vy) {
    fprintf(stderr,"BuhoFits: WARNING. Imagen mas pequena que VY\n");
    fprintf(stderr,"          NAXIS2=%d. VY=%d. Nuevo VY=%d\n",
               hin.naxis2,vy,y2);
    vy=y2;
    }
 
/* // -------------------------------------------- */
/* // Abro el fichero de salida y escribo cabecera */
/* // -------------------------------------------- */
  if((fpout=fopen(fileout,"w")) == NULL) {
    printf("BuhoFits: ERROR. No puedo abrir el fichero de salida %s",
		   fileout);
    exit(1);
  }
/*   Esto lo comento porque los errores en los parametros los voy a poner en el mismo fichero  */
/*   if(errflag) { */
/*     if((fperrout=fopen(fileerrout,"w")) == NULL) { */
/*       printf("BuhoFits: ERROR. No puedo abrir el fichero de salida %s", */
/* 	     fileerrout); */
/*       exit(1); */
/*     } */
/*   } */



/*   //t=time(NULL); */
/*   //tm=localtime(&t); */
/*   //fprintf(fpout,"#%19s (%19s %19s) %2d/%2d/%2d %2d:%2d:%2d  %4d %5.2f %4d\n", */
/*   //fitsin,fitssky,fitssig,tm->tm_mday,tm->tm_mon,tm->tm_year, */
/*   //  tm->tm_hour,tm->tm_min,tm->tm_sec,minpix,nsig,vy); */
/*   //printf("Aqui si \n%19s (%19s %19s) %2d/%2d/%2d %2d:%2d:%2d  %4d %5.2f %4d\n",fitsin,fitssky,fitssig,tm->tm_mday,tm->tm_mon,tm->tm_year, */
/*   //tm->tm_hour,tm->tm_min,tm->tm_sec,minpix,nsig,vy); */
/*   //free(&tm); */
  fprintf(fpout,"#Objects found by BuhoFits from image %s, with sky image %s and sigma image %s. Number of sigma detection: %f, rfixap= %f. Last parameters are errors in parameters\n",fitsin,fitssky,fitssig,nsig,rfixap);
  fprintf(fpout,"#Nobj        Xc          Yc       Flux    Area     Elip_A    Elip_B    Elip_T        mx2            my2          mxy       Sky         Sig      Xcpix    Ycpix      rKron     Flux_1rk     Flux2rk     Flux2.5rk    Fixapfot  ndeb");  
  if(errflag) {
    fprintf(fpout,"    err_Xc     err_Yc  err_Flux err_Area err_Elip_A err_Elip_B err_Elip_T err_mx2         err_my2      err_mxy    err_Sky    err_Sig  err_Xcpix  err_Ycpix err_rKron err_Flux_1rk err_Flux2rk err_Flux2.5rk err_Fixapfot");
/*     fprintf(fperrout,"#Errors in parameters for objects in %s from image %s, with sky image %s and sigma image %s. Number of sigma detection: %f, rfixap= %f\n",fileout,fitsin,fitssky,fitssig,nsig,rfixap); */
/*     fprintf(fperrout,"#Nobj     err_Xc      err_Yc  err_Flux err_Area err_Elip_A err_Elip_B err_Elip_T err_mx2         err_my2      err_mxy    err_Sky    err_Sig  err_Xcpix  err_Ycpix err_rKron err_Flux_1rk err_Flux2rk err_Flux2.5rk err_Fixapfot\n"); */
/*     Esto lo comento porque los errores los meto en el mismo fichero */
  }
  fprintf(fpout,"\n");

/* // Creamos una cabecera para la imagen de ventana */
/* // ---------------------------------------------- */
  memcpy(&hv,&hin,sizeof(struct headfits));
  hv.naxis2=vy;

/* // Restamos el fondo de cielo a la ventana (-1: no cielo; -0.9: < nsig sigmas) */
/* // --------------------------------------------------------------------------- */
  printf("Procesando  la primera windetec ....\n");
  printf("Procesando linea %3d/%3d",0,hv.naxis2);
  fflush(stdout);
  idx=0;
  
  for(j=0; j < hv.naxis2; j++) {
    y=hv.crval2+j*hv.cdelt2;
    printf("\b\b\b\b\b\b\b%3d/%3d",j+1,hv.naxis2); 
    fflush(stdout);
    if(y < ymin || y > ymax) {
      for(i=0; i < hv.naxis1; i++) {
        windetec[idx]=-1.0;
        /* if(errflag) winerrdetec[idx]=-1.0; */ /* Creo que no hace falta */
        idx++;
      }
    }
    else {
      for(i=0; i < hv.naxis1; i++) {
	x=hv.crval1+i*hv.cdelt1;
	if(x < xmin || x > xmax) {
	  windetec[idx]=-1.0;
	  /* if(errflag) winerrdetec[idx]=-1.0; */  /* Ni aqui tampoco  */
	}
	else {
	  fits_pv(imasky,&hsky,x,y,&sky);
	  windetec[idx] -= sky;
	  if(errflag) {
	    fits_pv(imaerrsky,&herrsky,x,y,&errsky);
	    winerrdetec[idx] = sqrt(winerrdetec[idx]*winerrdetec[idx]+errsky*errsky);
	  }
	  fits_pv(imasig,&hsig,x,y,&sig);
	  if(windetec[idx] < nsig*sig) {
	    windetec[idx]=-0.9; /* KIKE */
	    /* 	    if(errflag) winerrdetec[idx]=-0.9;  */
	  }
	  /*   windetec[idx]/=nsig*sig*3;    Antiguo */
	  /* 	  printf("VENNN 333 %f\n",windetec[idx]); */
	}
        idx++;
      }
    }
  }

  printf("\n\n");
  if(ffopen(&infits, fitsin, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(infits, "NAXIS", 1, 2, naxes_in, &nfound, &status)) fits_report_error(stderr,status);

  if(errflag) {
    if( ffopen(&errinfits, fitserrin, READONLY, &status)) fits_report_error(stderr,status);
    if(fits_read_keys_lng(errinfits, "NAXIS", 1, 2, naxes_errin, &nfound, &status)) fits_report_error(stderr,status);
  }

  

/*   fpin=fopen(fitsin,"rb"); */
/*   fperrin=fopen(fitserrin,"rb"); */

/* // ------------------------------- */
/* // Dimensionado de la matriz linea */
/* // ------------------------------- */
  if((linea=malloc(hv.naxis1*sizeof(float))) == NULL) {
    fprintf(stderr,"BuhoFits: ERROR. Error dimensionando la metriz linea\n");
    exit(1);
    }
  if(errflag) {
    if((errlinea=malloc(hv.naxis1*sizeof(float))) == NULL) {
      fprintf(stderr,"BuhoFits: ERROR. Error dimensionando la metriz errlinea\n");
      exit(1);
    }
  }

/* // ---------------------------------- */
/* // Programilla Principal: Preparacion */
/* // ---------------------------------- */
  if(graphmod == 'F') {
    cpgqvpXY(0.05,0.95,0.1,0.4,0.,hv.naxis1+1.0,0.,hv.naxis2+1.0,
	     &gvpx1,&gvpx2,&gvpy1,&gvpy2);
    cpgsvp(gvpx1,gvpx2,gvpy1,gvpy2);
    cpgswin(0.,hv.naxis1+1.0,0.,hv.naxis2+1.0);
    cpgbox("bctns",0.0,0,"bctns",0.0,0);
    }

  printf("Busqueda de objetos ...\n\n");
  printf(" Linea     No. de objetos\n");
  printf("-------------------------------------\n");
  printf(" %6d/%6d     %7d",0,hin.naxis2,nobj);
  fflush(stdout);
  liny=hv.naxis2;

  /*  ----------------------------------------------- */
  /*  Bucle principal a todas las lineas de la imagen */
  /*  =============================================== */
  
  jdetect=-1;  /* Es la linea dentro de la ventana en la que esta detectando objetos */

  for(j=0; j<hin.naxis2; j++) {
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"); 
    printf("%6d/%6d     %7d",j+1,hin.naxis2,nobj); 
    fflush(stdout);
    
    jdetect++;  /*Aumenta uno la linea de detección */
    if(j<=hin.naxis2-(int)(hv.naxis2/2))
      if(jdetect>(int)(hv.naxis2/2.)) jdetect=(int)(hv.naxis2/2.);


    /* Lo que he hecho ahora es que la linea de detección esté 
       siempre a la mitad de la ventana uan vez que ya ha pasado por todas 
       las líneas de abajo */
/*     jdetect=0; */
/*     printf("\n Linea real: %d %f %f jdet %d liny %d hinacis2 %d\n",j,jdetect+hv.crval2,hv.crval2,jdetect,liny,hin.naxis2); */


    if(graphmod == 'F') {

      cpgsvp(gvpx1,gvpx2,gvpy1,gvpy2);
      cpgswin(0.,hv.naxis1+1.0,0.,hv.naxis2+1.0);
      cpggray(windetec,hv.naxis1,hv.naxis2,1,hv.naxis1,1,hv.naxis2,-1.0,0.0,tr); /* // KIKE  */
/*       //cpggray(windetec,hv.naxis1,hv.naxis2,1,hv.naxis1,1,hv.naxis2,-0.1,1.0,tr); */
      }
    for(i=0; i<hv.naxis1; i++) {

      object.area=0;
      object.nb=0;

	if(errflag) errDetecta(windetec,winerrdetec,&hv,i,jdetect,&object,&errobject);
	else       Detecta(windetec,&hv,i,jdetect,&object);

/* // Sera un objeto si tiene mas de minpix pixels */
/* // -------------------------------------------- */
      if(object.area > minpix) {
	CFM(&object);
	object.xcpix=((object.xc-hin.crval1)/hin.cdelt1)+1.0;
	object.ycpix=((object.yc-hin.crval2)/hin.cdelt2)+1.0;
	fits_pv(imasky,&hsky,object.xc,object.yc,&sky);
	fits_pv(imasig,&hsig,object.xc,object.yc,&sig);
	if(errflag) fits_pv(imaerrsky,&herrsky,object.xc,object.yc,&errsky);
	object.sky=sky;
	object.sig=sig;
	if(errflag) errobject.errsky=errsky;
	object.xpv=((object.xc-hv.crval1)/hv.cdelt1)+1.0;
	object.ypv=((object.yc-hv.crval2)/hv.cdelt2)+1.0;

	ndeb=0;
 	if(DEBUG) printf(" 666 ENRO \n"); 
/* 	ndeb=readi(ndeb); */
	if(deblenflag) {
	  if(errflag)  object.deblenflag=errDeblending(winfot,winerrfot,hv,imasky,hsky,imasig,hsig,imaerrsky,herrsky,object,nsig,&objsdeb,&errobjsdeb,&ndeb);
	  else         object.deblenflag=   Deblending(winfot,hv,imasky,hsky,imasig,hsig,object,nsig,&objsdeb,&ndeb);
	}
	else {
	  ndeb=1;
	  objsdeb=malloc(1*sizeof(struct obj));
	  memcpy(objsdeb,&(object),sizeof(struct obj));
	  if(errflag) {
	    errobjsdeb=malloc(1*sizeof(struct objerr));
	    memcpy(errobjsdeb,&(errobject),sizeof(struct objerr));
	  }
	}
/*  	if(DEBUG) printf(" 666 PUNTERO %d\n",objsdeb);  */

	for(ideb=0;ideb<ndeb;ideb++) {
/* 	  printf("  FINAL OBJECT: %d  xc %f yc %f sig %d  \n",ideb,objsdeb[ideb].xc,objsdeb[ideb].yc,objsdeb[ideb].isig); */
	}	    
	
        nobj=nobj+ndeb;
        printf("\b\b\b\b\b\b\b");
        printf("%7d",nobj);
        fflush(stdout);
	
	/*        Calculo de los parametros */
	/*        ------------------------- */
	/* 	printf(" empieza calculo \n"); */

	for(ideb=0;ideb<ndeb;ideb++) {
	  
	  
	  if(errflag) errCFM(objsdeb+ideb,errobjsdeb+ideb);
	  else           CFM(objsdeb+ideb);
	  objsdeb[ideb].xcpix=((objsdeb[ideb].xc-hin.crval1)/hin.cdelt1)+1.0;
	  objsdeb[ideb].ycpix=((objsdeb[ideb].yc-hin.crval2)/hin.cdelt2)+1.0;
	  fits_pv(imasky,&hsky,objsdeb[ideb].xc,objsdeb[ideb].yc,&sky);
	  fits_pv(imasig,&hsig,objsdeb[ideb].xc,objsdeb[ideb].yc,&sig);
	  objsdeb[ideb].sky=sky;
	  objsdeb[ideb].sig=sig;
	  objsdeb[ideb].xpv=((objsdeb[ideb].xc-hv.crval1)/hv.cdelt1)+1.0;
	  objsdeb[ideb].ypv=((objsdeb[ideb].yc-hv.crval2)/hv.cdelt2)+1.0;

	  if(errflag) {
	    errobjsdeb[ideb].errxcpix=errobjsdeb[ideb].errxc/hin.cdelt1;
	    errobjsdeb[ideb].errycpix=errobjsdeb[ideb].erryc/hin.cdelt2;
	    fits_pv(imaerrsky,&herrsky,objsdeb[ideb].xc,objsdeb[ideb].yc,&errsky);
	    errobjsdeb[ideb].errsky=errsky;
	  }

	  AjustaElipse(objsdeb+ideb);
	 

	  if(errflag) errApertureFot(winfot,winerrfot,hv.naxis1,hv.naxis2,rfixap,objsdeb+ideb,errobjsdeb+ideb); 
	  else           ApertureFot(winfot,hv.naxis1,hv.naxis2,rfixap, objsdeb+ideb);


	  /* 	xcpix, ycpix son las coordenadas en PIXELES. NO? */
	  /* 	mientras que las coordenadas en unidades de matriz son estas -1 */
	  
	  /* 	 Escritura de los parametros */
	  /* 	 --------------------------- */
	  fprintf(fpout," %6d %10.4f %10.4f %10.3g %6d ",  nobj-ndeb+ideb+1,objsdeb[ideb].xc,objsdeb[ideb].yc,objsdeb[ideb].ft,objsdeb[ideb].area);
	  fprintf(fpout,"%9.3f %9.3f %9.3f ",              objsdeb[ideb].elipa,objsdeb[ideb].elipb,objsdeb[ideb].elipt);
	  fprintf(fpout,"%13.5g %13.5g %13.5g ",           objsdeb[ideb].mx2,objsdeb[ideb].my2,objsdeb[ideb].mxy);
	  fprintf(fpout,"%10.4f %10.4f %8.2f %8.2f %10.4f",objsdeb[ideb].sky,objsdeb[ideb].sig,objsdeb[ideb].xcpix,objsdeb[ideb].ycpix,objsdeb[ideb].rkron);
	  fprintf(fpout," %12.4g %12.4g %12.4g %12.4g %4d",    objsdeb[ideb].apfot1,objsdeb[ideb].apfot2,objsdeb[ideb].apfot25,objsdeb[ideb].fixapfot,objsdeb[ideb].ideb);
	  if(errflag) {
	    errobjsdeb[ideb].errarea=0.; errobjsdeb[ideb].errelipa=0.; errobjsdeb[ideb].errelipb=0.; errobjsdeb[ideb].errelipt=0.; errobjsdeb[ideb].errsig=0.;
	    fprintf(fpout," %10.4f %10.4f %10.3g %6f ",      errobjsdeb[ideb].errxc,errobjsdeb[ideb].erryc,errobjsdeb[ideb].errft,errobjsdeb[ideb].errarea);
	    fprintf(fpout,"%9.3f %9.3f %9.3f ",              errobjsdeb[ideb].errelipa,errobjsdeb[ideb].errelipb,errobjsdeb[ideb].errelipt);
	    fprintf(fpout,"%13.5g %13.5g %13.5g ",           errobjsdeb[ideb].errmx2,errobjsdeb[ideb].errmy2,errobjsdeb[ideb].errmxy);
	    fprintf(fpout,"%10.4f %10.4f %8.2f %8.2f %10.4f",errobjsdeb[ideb].errsky,errobjsdeb[ideb].errsig,errobjsdeb[ideb].errxcpix,errobjsdeb[ideb].errycpix,errobjsdeb[ideb].errrkron);
	    fprintf(fpout," %12.4g %12.4g %12.4g %12.4g"    ,errobjsdeb[ideb].errapfot1,errobjsdeb[ideb].errapfot2,errobjsdeb[ideb].errapfot25,errobjsdeb[ideb].errfixapfot);
	  }
	  fprintf(fpout,"\n");
	  
	}
	
	/*   Dibujo */
	/*   ------ */
	if(graphmod == 'F'|| graphmod == 'C') {
	  pgLimits(object.area,object.xobj,&pxmin,&pxmax);
	  pgLimits(object.area,object.yobj,&pymin,&pymax);
	  if(graphmod == 'F') {
	    cpgsvp(0.0,1.0,0.41,1.0);
    	    cpgswin(0.,1.,0.,1.);
    	    cpgsci(0);
	    cpgrect(0.,1.,0.,1.);
	    cpgsci(1);
	    cpgsvpXY(0.2,0.9,0.5,0.9,pxmin,pxmax,pymin,pymax);
	  }
	  else {
	    cpgenv(pxmin,pxmax,pymin,pymax,1,-2);
	  }
	  cpgswin(pxmin,pxmax,pymin,pymax);
	  cpgbox("bctns",0.0,0,"bctns",0.0,0);
	  sprintf(chlab,"%s#%05d",fitsin,nobj);
	  cpglab("X\\dWORLD\\u","Y\\dWORLD\\u",chlab);

	  for(ideb=0;ideb<ndeb;ideb++) {
	    
	    cpgpt(objsdeb[ideb].area,objsdeb[ideb].xobj,objsdeb[ideb].yobj,1);
	    cpgsci(3); 
	    cpgpt(objsdeb[ideb].nb,objsdeb[ideb].xb,objsdeb[ideb].yb,2);
	    
	    cpgsci(2);
	    cpgpt(1,&(objsdeb[ideb].xc),&(objsdeb[ideb].yc),6);
	    cpgpt(1,&(objsdeb[ideb].xcelip),&(objsdeb[ideb].ycelip),7);
	    cpgelip(objsdeb[ideb].xcelip,objsdeb[ideb].ycelip,objsdeb[ideb].elipa,objsdeb[ideb].elipb,objsdeb[ideb].elipt);
	  } 
	  cpgsci(1);
	}
	
	free(objsdeb);
	if(errflag) free(errobjsdeb);

      }
    }
    
    

    /*  Actualizacion de la ventana: Scroll, Lectura y Resta del cielo */
    /*  -------------------------------------------------------------- */
    /*     Scroll de la ventana de deteccion y de la de fotometria de forma identica */
    /*     printf(" Antes del scroll\n"); */
		
    if(j>=(int)(hv.naxis2/2.) && liny < hin.naxis2) {
      memcpy(windetec,&(windetec[hv.naxis1]),hv.naxis1*(hv.naxis2-1)*sizeof(float));
      if (errflag)  memcpy(winerrdetec,&(winerrdetec[hv.naxis1]),hv.naxis1*(hv.naxis2-1)*sizeof(float));
      memcpy(winfot,&(winfot[hv.naxis1]),hv.naxis1*(hv.naxis2-1)*sizeof(float));
      if (errflag)  memcpy(winerrfot,&(winerrfot[hv.naxis1]),hv.naxis1*(hv.naxis2-1)*sizeof(float));
      hv.crval2 = hin.crval2 + (liny+1-hv.naxis2) * hin.cdelt2;
/*       printf(" mas real %f %d %f %f %d\n",hin.crval2,j,hin.cdelt2,hv.crval2,liny); */
      idx=0;
      liny++;
      
      /*       f_rfb(fpin,&hv,linea,1,hv.naxis1,liny,liny); */
      blc[0]=1;blc[1]=liny;trc[0]=hv.naxis1;trc[1]=liny;incr[0]=1;incr[1]=1;
      if(fits_read_subset_flt(infits,1,2,naxes_in,blc,trc,incr,0.,linea,&anynull,&status))fits_report_error(stderr,status) ;
      if(errflag)  if(fits_read_subset_flt(errinfits,1,2,naxes_in,blc,trc,incr,0.,errlinea,&anynull,&status))fits_report_error(stderr,status) ;
      
      /*       Hay que copiar este trozo en winfot antes de que sea modificado */
      idx = (hv.naxis2-1)*hv.naxis1;
      for(i=0; i<hv.naxis1; i++) winfot[idx+i]=linea[i];
      if(errflag) for(i=0; i<hv.naxis1; i++) winerrfot[idx+i]=errlinea[i];
      /*       Seguimos modificando linea que posteriormente se copia en windetec */
      idx=0;
      y=hv.crval2+(hv.naxis2-1)*hv.cdelt2;
      if(y < ymin  || y > ymax) {
        for(i=0; i < hv.naxis1; i++) {
          linea[idx]=-1.0;
          idx++;
	}
      }
      else {
        x=hv.crval1;
        for(i=0; i < hv.naxis1; i++) {
          if(x < xmin || x > xmax) {
            linea[idx]=-1.0;
	  }
          else {
            fits_pv(imasky,&hsky,x,y,&sky);
	    /* 	    printf(" SSKY %f   x  %f y  %f\n",sky,x,y); */
	    linea[idx] -= sky;
	    if(errflag) {
	      fits_pv(imaerrsky,&herrsky,x,y,&errsky);
	      errlinea[idx] = sqrt(errlinea[idx]*errlinea[idx]+errsky*errsky);
	    }
            fits_pv(imasig,&hsig,x,y,&sig);
	    /* 	    printf(" SIG %f   x  %f y  %f\n",sig,x,y); */
	    /* 	    printf(" sky %f sig %f\n",sky,sig); */
            if(linea[idx] < nsig*sig) linea[idx]=-0.9;
	  }
          idx++;
          x += hv.cdelt1;
	}
      }
      idx = (hv.naxis2-1)*hv.naxis1;
      for(i=0; i<hv.naxis1; i++) windetec[idx+i]=linea[i];
      if(errflag) for(i=0; i<hv.naxis1; i++) winerrdetec[idx+i]=errlinea[i];
      
      
    }
    else if(liny>=hin.naxis2) {
      /*       for(i=0; i < hv.naxis1; i++) { */
      /* 	Si se ha acabado la imagen, rellenamos con ceros la ventana fotometrica */
      /* 	winfot[(hv.naxis2-1)*hv.naxis1+i]=0; */
      /* 	if(errflag) winerrfot[(hv.naxis2-1)*hv.naxis1+i]=0; */
      /* 	Y ademas la linea de siempre */
      /* 	linea[idx]=-1.0; */
      /*         idx++; */
      /*       } */
      /*       idx = (hv.naxis2-1)*hv.naxis1; */
      /*       for(i=0; i<hv.naxis1; i++) windetec[idx+i]=linea[i]; */
      /*       if(errflag) for(i=0; i<hv.naxis1; i++) winerrdetec[idx+i]=errlinea[i]; */
      /*     Creo que esto hay que ponerlo para que la fotometria no haga tonterias */
      /*     Y solo en este caso, en el que liny>hv.naxis2 */
      /*       hv.naxis2--;       */
/*       printf("real pasa menos todavia\n"); */
    }
    else {
      /*Estoy todavia en las primeras vy/2 lineas */
      /* Luego no hago nada, todas ls ventanas deben estar estaticas */
/*       printf(" real pasa nada\n"); */
    }





    
    
  }
  
  
  cpgend();
  free(linea);
  if(errflag) free(errlinea);
  free(windetec);
  free(winfot);
  if(errflag) free(winerrfot);
  free(imasky);
  free(imasig);
/*   fclose(fpin); */
  fclose(fpout);
  
  printf("\n\n");
  SlpTime(stdout);
  printf("Fin del programa\n");
  return(0);
} 



/* //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int Detecta(float *v, struct headfits *hv, int i, int j, struct obj *O) 

{
/*   int k; */

  int p=0;
/*   float tr[6]; */
/*   tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;   */

  
/*   cpgwnad(0,hv->naxis1,0,hv->naxis2); */
/*   cpggray(v,hv->naxis1,hv->naxis2,1,hv->naxis1,1,hv->naxis2,0.,-1.0,tr);  */
/*   cpgbox("bctns",0.0,0,"bctns",0.0,0);  */
/*   cpgswin((-1)*hv->cdelt1+hv->crval1,(hv->naxis1-1)*hv->cdelt1+hv->crval1,(-1)*hv->cdelt2+hv->crval2,(hv->naxis2-1)*hv->cdelt2+hv->crval2); */

  if(mierdaflag) {

  cpgsci(6);
  cpgpt1(hv->crval1+i*hv->cdelt1,hv->crval2+j*hv->cdelt2,6);
  p=readi(p);
  cpgsci(5);
  cpgpt1(hv->crval1+i*hv->cdelt1,hv->crval2+j*hv->cdelt2,6);
  cpgsci(1);

  }

/*   printf(" 666 0\n"); */

  if((*O).area > MAXPIX - 1) return(0);
  if(i < 0 || i > hv->naxis1-1 || j < 0 || j > hv->naxis2-1 ) return(0);
  if(v[i+j*hv->naxis1] < -0.8) return(0);
  if(v[i+j*hv->naxis1] < 0) {

    if(mierdaflag) {
      cpgsci(4);
      cpgpt1(hv->crval1+i*hv->cdelt1,hv->crval2+j*hv->cdelt2,6);
/*       printf(" Ya incluido %d %d\n",i,j); */
    }
    return(1);
  }
  
  /*  Llego aqui si hay pixel por encima del fondo */
  /*  -------------------------------------------- */
  /*   for(k=0;k<*area;k++)   if(f[k]=v[i+j*hv->naxis1]) return(1); */
  

/*   printf(" 666 1 area %d  punt %d\n",(*O).area,O); */

  (*O).f[(*O).area]=v[i+j*hv->naxis1];
  (*O).xobj[(*O).area]=hv->crval1+i*hv->cdelt1;
  (*O).yobj[(*O).area]=hv->crval2+j*hv->cdelt2;
  (*O).ipix[(*O).area]=i;
  (*O).jpix[(*O).area]=j;

/*   printf(" 666 2\n"); */

  v[i+j*hv->naxis1] = -0.7;
  ((*O).area)++;
  if(mierdaflag) {
    cpgsci(3);
    cpgpt1(hv->crval1+i*hv->cdelt1,hv->crval2+j*hv->cdelt2,6);
    
/*     printf(" %d  Nuevo %d %d\n",(*O).area,i,j); */
  }
  /* // Recursivo a los pixels de alrededor */
/* // ----------------------------------- */
  p += Detecta(v,hv,i  ,j+1,O);
  p += Detecta(v,hv,i+1,j+1,O);
  p += Detecta(v,hv,i+1,j  ,O);
  p += Detecta(v,hv,i+1,j-1,O);
  p += Detecta(v,hv,i  ,j-1,O);
  p += Detecta(v,hv,i-1,j-1,O);
  p += Detecta(v,hv,i-1,j  ,O);
  p += Detecta(v,hv,i-1,j+1,O);

  if(p != 8) {
    (*O).xb[(*O).nb]=hv->crval1+i*hv->cdelt1;
    (*O).yb[(*O).nb]=hv->crval2+j*hv->cdelt2;
    if(mierdaflag) { 
      cpgsci(4);
      cpgpt1(hv->crval1+i*hv->cdelt1,hv->crval2+j*hv->cdelt2,6);
/*       printf(" %d BORDE i %d j %d\n",(*O).nb,i,j);  */
    }
    ((*O).nb)++;
    v[i+j*hv->naxis1] = -0.1;
    }

  return(1);

}

/* //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void CFM(struct obj *O)


{
  int i;
  float rf,r;

  float *f,*x,*y;

  f=(*O).f;
  x=(*O).xobj;
  y=(*O).yobj;

  (*O).ft = 0.0;
  (*O).xc = 0.0;
  (*O).yc = 0.0;
  (*O).mx2 = 0.0;
  (*O).my2 = 0.0;
  (*O).mxy = 0.0;
  rf  = 0.0;
  for(i=0; i<(*O).area; i++) {
    (*O).ft += f[i];
    (*O).xc += x[i]*f[i];
    (*O).yc += y[i]*f[i];
    (*O).mx2 += x[i]*x[i]*f[i];
    (*O).my2 += y[i]*y[i]*f[i];
    (*O).mxy += x[i]*y[i]*f[i];
  }
  (*O).xc /= (*O).ft;
  (*O).yc /= (*O).ft;
  (*O).mx2 /= (*O).ft;
  (*O).my2 /= (*O).ft;
  (*O).mxy /= (*O).ft;

  for(i=0;i<(*O).area;i++) {
    r=sqrt((x[i]-(*O).xc)*(x[i]-(*O).xc)+(y[i]-(*O).yc)*(y[i]-(*O).yc));
    rf+=r*f[i];
  }
  (*O).rkron=rf / (*O).ft;
}

/* //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void AjustaElipse(struct obj *O)


{
  int i;
  float a,b,c,d,f,e;
/*   float dum1; */
  int ctrl1,ctrl2;

  float *xb,*yb;

  if(DEBUG) printf(" ntes nada ajus\n");
  xb=(*O).xb;
  yb=(*O).yb;

  
  
  /*  Restamos las coordenadas del centro para que no se disparen las sumas */
  /*  --------------------------------------------------------------------- */
  for(i=0; i<(*O).nb; i++) {
    xb[i] -= (*O).xc;
    yb[i] -= (*O).yc;
    if(DEBUG) printf(" Pix %d  xb %f yb %f\n",i,xb[i],yb[i]);
  }
  
  /*   if(MCElip(nb,xb,yb,&a,&b,&c,&d,&f,&e)) { */
  /*     if(ElipPar(a,b,c,d,f,e,xcelip,ycelip,elipa,elipb,elipt)) { */
  if(DEBUG) printf(" Antes elip\n");
  ctrl1=MCElip((*O).nb,xb,yb,&a,&b,&c,&d,&f,&e);
  if(DEBUG) printf(" despues elip\n");
  
  /* Volvemos a restablecer los valores de los bordes */
  for(i=0; i<(*O).nb; i++) {
    xb[i] += (*O).xc;
    yb[i] += (*O).yc;
  }
   
  if(ctrl1) {
    ctrl2=ElipPar(a,b,c,d,f,e,&(*O).xcelip,&(*O).ycelip,&(*O).elipa,&(*O).elipb,&(*O).elipt);
    if(ctrl2) {
      if(DEBUG) printf("Ajuste a %g  b %g  c %g d %g e  %g f %g\n",a,b,c,d,e,f); 
      if(DEBUG) printf("Ajuste xcelip %g ycelip %g elipa %g elipb %g elipt %g\n",(*O).xcelip,(*O).ycelip,(*O).elipa,(*O).elipb,(*O).elipt); 
      (*O).xcelip += (*O).xc;
      (*O).ycelip += (*O).yc;
      return;
      }
    }



  (*O).xcelip=0;
  (*O).ycelip=0;
  (*O).elipa=1;
  (*O).elipb=1;
  (*O).elipt=0;
  return;
}



void ApertureFot(float *ima, int nx, int ny,float rfixap, struct obj *O) 
{

  /*   //Subrutina para calcular fotometria de apertura */
  /*   // De momento lo hago copiando lo de ventana en otra matriz de */
  /*   // modo que se necesita el doble de memoria. */
  /*   float tr[6]={0.,1.,0.,0.,0.,1.}; */
  /*   int pg; */
  /*   float r; */
  int npix;
  /*   int i; */
  /*   float fnul; */
  /*   char cnul; */
  /*   cpgqid(&pg); */
  /*   cpgopen("/xserve"); */
  /*   cpgsvp(0.1,0.4,0.1,0.9); */
  /*   cpgwnad((int)xc-30,(int)xc+30,1.,ny-1); */
  /*   cpggray(ima,nx,ny,(int)xc-30,(int)xc+30,1,ny-1,sky-1.5*sqrt(sky),sky+3.5*sqrt(sky),tr); */
  /*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
  /*   cpgsci(2); */
  /*   cpgslw(3); */
  /*   cpgsls(2); */
  /*   cpgelip(xc,yc,elipa,elipb,elipt); */
  /*   cpgslw(1); */
  /*   cpgsls(1); */
  /*   cpgsci(1); */
  /*   cpgsvp(0.5,0.9,0.1,0.9); */
  /*   cpgswin(0.,rkron*5,0.,sky*npixobj*5); */
  /*   printf(" Aqui no se quen\n"); */
  /*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
  /*     printf(" Seque da\n"); */
  /*   cpgsch(2.); */
  /*   printf("Entra\n"); */
  
  /*   cpgsfs(2); */
  /*   for(i=0;i< 15;i++) { */
  /*     cpgsci(i); */
  /*     r=rkron*5/15*i; */
  /*     printf(" IMAAAA 590 21 %f\n",ima[590+nx*21]); */
  /*     *apfot1=circ_aper(ima,nx,ny,xc,yc,r,&npix); */
  /*     cpgsvp(0.1,0.4,0.1,0.9); */
  /*     cpgwnad((int)xc-30,(int)xc+30,1.,ny-1); */
  /*     cpgcirc(xc,yc,r); */
  /*     cpgsvp(0.5,0.9,0.1,0.9); */
  /*     cpgswin(0.,rkron*5,0.,ft*2.5); */
  /*     printf(" radio %f fot %f fot-sky %f npix %d\n",r,*apfot1,*apfot1-npix*sky,npix); */
  /*     cpgpt1(r,*apfot1-npix*sky,3); */
  /*   } */

/*   printf(" RKRON %f rfixap %f\n",(*O).rkron,rfixap); */


  (*O).apfot1=circ_aper(ima,nx,ny,(*O).xpv,(*O).ypv,(*O).rkron,&npix);
  /*  *apfot1-=npix*sky; */
  (*O).apfot1-=M_PI*(*O).rkron*(*O).rkron*(*O).sky; 
  (*O).apfot2=circ_aper(ima,nx,ny,(*O).xpv,(*O).ypv,2*(*O).rkron,&npix);
  (*O).apfot2-=M_PI*2*(*O).rkron*2*(*O).rkron*(*O).sky;
  (*O).apfot25=circ_aper(ima,nx,ny,(*O).xpv,(*O).ypv,2.5*(*O).rkron,&npix);
  (*O).apfot25-=M_PI*2.5*(*O).rkron*2.5*(*O).rkron*(*O).sky;

/* circ_aper_bak(ima,nx,ny,xc,yc,2.5*rkron,&npix); */
/*   printf(" F1 %g F2 %g FS1 %g FS2 %g NPIX %d NNEW %f SKY %f\n",circ_aper(ima,nx,ny,(*O).xpv,(*O).ypv,2.5*(*O).rkron,&npix),circ_aper_bak(ima,nx,ny,(*O).xc,(*O).yc,2.5*(*O).rkron,&npix),(*O).apfot25,circ_aper_bak(ima,nx,ny,(*O).xpv,(*O).ypv,2.5*(*O).rkron,&npix)-npix*(*O).sky,npix,M_PI*2.5*(*O).rkron*2.5*(*O).rkron,(*O).sky); */

  (*O).fixapfot=circ_aper(ima,nx,ny,(*O).xpv,(*O).ypv,rfixap,&npix);
  (*O).fixapfot-=M_PI*rfixap*rfixap*(*O).sky;

/*   cpgsci(1); */
/*   cpgmove(rkron,0); */
/*   cpgdraw(rkron,sky*npixobj); */
/*   cpgsci(2); */
/*   cpgmove(2.5*rkron,0); */
/*   cpgdraw(2.5*rkron,sky*npixobj); */
/*   cpgsci(5); */
/*   cpgmove(sqrt(elipa*elipb),0); */
/*   cpgdraw(sqrt(elipa*elipb),ft*2.5); */
/*   cpgmove(0,ft); */
/*   cpgdraw(5*rkron,ft); */
/*   cpgslw(1); */

/*   printf(" Los pixeles cnetrle s %f %f SKY %f\n",xc,yc,sky); */
/*   printf(" Rkron %f\n",rkron); */
/*   //cpgcurs(&fnul,&fnul,&cnul); */
/*   cpgclos(); */
/*   cpgslct(pg); */
}


int errDetecta(float *v, float *errv,struct headfits *hv, int i, int j, struct obj *O, struct objerr *EO)

{
/*   int k; */

  int p=0;

  if((*O).area > MAXPIX - 1) return(0);
  if(i < 0 || i > hv->naxis1-1 || j < 0 || j > hv->naxis2-1 ) return(0);
  if(v[i+j*hv->naxis1] < -0.8) return(0);
  if(v[i+j*hv->naxis1] < 0) return(1);

/* // Llego aqui si hay pixel por encima del fondo */
/* // -------------------------------------------- */
/*   //for(k=0;k<(*O).area;k++)   if(f[k]=v[i+j*hv->naxis1]) return(1); */
  
  (*O).f[(*O).area]=v[i+j*hv->naxis1];
  (*O).xobj[(*O).area]=hv->crval1+i*hv->cdelt1;
  (*O).yobj[(*O).area]=hv->crval2+j*hv->cdelt2;
  (*O).ipix[(*O).area]=i;
  (*O).jpix[(*O).area]=j;
  v[i+j*hv->naxis1] = -0.7;
  (*EO).errf[(*O).area]=errv[i+j*hv->naxis1];

/*   printf(" errf[%d] = %f x %f y %f\n",(*O).area,errf[(*O).area],x[(*O).area],y[(*O).area]); */

  ((*O).area)++;


/* // Recursivo a los pixels de alrededor */
/* // ----------------------------------- */
  p += errDetecta(v,errv,hv,i  ,j+1,O,EO);
  p += errDetecta(v,errv,hv,i+1,j+1,O,EO);
  p += errDetecta(v,errv,hv,i+1,j  ,O,EO);
  p += errDetecta(v,errv,hv,i+1,j-1,O,EO);
  p += errDetecta(v,errv,hv,i  ,j-1,O,EO);
  p += errDetecta(v,errv,hv,i-1,j-1,O,EO);
  p += errDetecta(v,errv,hv,i-1,j  ,O,EO);
  p += errDetecta(v,errv,hv,i-1,j+1,O,EO);

  if(p != 8) {
    (*O).xb[(*O).nb]=hv->crval1+i*hv->cdelt1;
    (*O).yb[(*O).nb]=hv->crval2+j*hv->cdelt2;
    ((*O).nb)++;
    v[i+j*hv->naxis1] = -0.1;
    }

  return(1);
}



void errCFM(struct obj *O,struct objerr *EO)


{
  int i,j;
  float rf,r;
  float tmp1=0.0,tmp2=0.0,tmp3=0.0,tmp4=0.0,tmp5=0.0,tmp6=0.0;
  float *f,*errf,*x,*y;

  f=(*O).f;
  errf=(*EO).errf;
  x=(*O).xobj;
  y=(*O).yobj;

  (*O).ft = 0.0;
  (*O).xc = 0.0;
  (*O).yc = 0.0;
  (*O).mx2 = 0.0;
  (*O).my2 = 0.0;
  (*O).mxy = 0.0;
  rf  = 0.0;
  (*EO).errft=0.0;
  for(i=0; i<(*O).area; i++) {
    (*O).ft    += f[i];
    (*O).xc    += x[i]*f[i];
    (*O).yc    += y[i]*f[i];
    (*O).mx2   += x[i]*x[i]*f[i];
    (*O).my2   += y[i]*y[i]*f[i];
    (*O).mxy   += x[i]*y[i]*f[i];
/*     printf(" errf %f errft  %f \n",errf[i],*errft); */
    (*EO).errft += errf[i]*errf[i];
  }
  for(j=0; j<(*O).area; j++) {
    tmp1 += ((*O).ft*x[j]-(*O).xc)*((*O).ft*x[j]-(*O).xc)*errf[j]*errf[j];
    tmp2 += ((*O).ft*y[j]-(*O).yc)*((*O).ft*y[j]-(*O).yc)*errf[j]*errf[j];
    tmp3 += ((*O).ft*x[j]*x[j]-(*O).mx2)*((*O).ft*x[j]*x[j]-(*O).mx2)*errf[j]*errf[j];
    tmp4 += ((*O).ft*y[j]*y[j]-(*O).my2)*((*O).ft*y[j]*y[j]-(*O).my2)*errf[j]*errf[j];
    tmp5 += ((*O).ft*x[j]*y[j]-(*O).mxy)*((*O).ft*x[j]*y[j]-(*O).mxy)*errf[j]*errf[j];
  }

  (*EO).errxc=sqrt(tmp1/((*O).ft*(*O).ft*(*O).ft*(*O).ft));
  (*EO).erryc=sqrt(tmp2/((*O).ft*(*O).ft*(*O).ft*(*O).ft));
  (*EO).errmx2=sqrt(tmp3/((*O).ft*(*O).ft*(*O).ft*(*O).ft));
  (*EO).errmy2=sqrt(tmp4/((*O).ft*(*O).ft*(*O).ft*(*O).ft));
  (*EO).errmxy=sqrt(tmp5/((*O).ft*(*O).ft*(*O).ft*(*O).ft));
  (*EO).errft=sqrt((*EO).errft);


  (*O).xc /= (*O).ft;
  (*O).yc /= (*O).ft;
  (*O).mx2 /= (*O).ft;
  (*O).my2 /= (*O).ft;
  (*O).mxy /= (*O).ft;

  for(i=0;i<(*O).area;i++) {
    r=sqrt((x[i]-(*O).xc)*(x[i]-(*O).xc)+(y[i]-(*O).yc)*(y[i]-(*O).yc));
    rf+=r*f[i];
  }
  for(j=0;j<(*O).area;j++) {
    r=sqrt((x[j]-(*O).xc)*(x[j]-(*O).xc)+(y[j]-(*O).yc)*(y[j]-(*O).yc));
    tmp6 += ((*O).ft*r-rf)*((*O).ft*r-rf)*errf[j]*errf[j];
  }
  (*EO).errrkron=sqrt(tmp6/((*O).ft*(*O).ft*(*O).ft*(*O).ft));
  (*O).rkron=rf / (*O).ft;
}





void errApertureFot(float *ima, float *errima,int nx, int ny,float rfixap, struct obj *O, struct objerr *EO)

{

/*   //Subrutina para calcular fotometria de apertura */
  /* Incluyendo tratamiento de errores */
  int npix;
  (*O).apfot1=circ_aper_err(ima,errima,nx,ny,(*O).xpv,(*O).ypv,(*O).rkron,&npix,&(*EO).errapfot1);
  (*O).apfot1-=M_PI*(*O).rkron*(*O).rkron*(*O).sky; 
  (*O).apfot2=circ_aper_err(ima,errima,nx,ny,(*O).xpv,(*O).ypv,2*(*O).rkron,&npix,&(*EO).errapfot2);
  (*O).apfot2-=M_PI*2*(*O).rkron*2*(*O).rkron*(*O).sky;
  (*O).apfot25=circ_aper_err(ima,errima,nx,ny,(*O).xpv,(*O).ypv,2.5*(*O).rkron,&npix,&(*EO).errapfot25);
  (*O).apfot25-=M_PI*2.5*(*O).rkron*2.5*(*O).rkron*(*O).sky;
  (*O).fixapfot=circ_aper_err(ima,errima,nx,ny,(*O).xpv,(*O).ypv,rfixap,&npix,&(*EO).errfixapfot);
  (*O).fixapfot-=M_PI*rfixap*rfixap*(*O).sky;


  (*EO).errapfot1=sqrt(((*EO).errapfot1)*((*EO).errapfot1)+(M_PI*(*O).rkron*(*O).rkron*(*EO).errsky)*(M_PI*(*O).rkron*(*O).rkron*(*EO).errsky));
  (*EO).errapfot2=sqrt(((*EO).errapfot2)*((*EO).errapfot2)+(M_PI*2*(*O).rkron*2*(*O).rkron*(*EO).errsky)*(M_PI*2*(*O).rkron*2*(*O).rkron*(*EO).errsky));
  (*EO).errapfot25=sqrt(((*EO).errapfot25)*((*EO).errapfot25)+(M_PI*2.5*(*O).rkron*2.5*(*O).rkron*(*EO).errsky)*(M_PI*2.5*(*O).rkron*2.5*(*O).rkron*(*EO).errsky));
  (*EO).errfixapfot=sqrt(((*EO).errfixapfot)*((*EO).errfixapfot)+(M_PI*rfixap*rfixap*(*EO).errsky)*(M_PI*rfixap*rfixap*(*EO).errsky));


}


int Deblending(float *ima, struct headfits hv, float *imasky, struct headfits hsky,float *imasig, struct headfits hsig,struct obj objprimal, float nsig, struct obj **objsdeb, int *ndeb) {



  float maxnsig;


  int   i,j,k;
  int   iobj,iobj2;
  float newsig;
  int isig=0;
  float *flux;
  float *windebl;
  float xmin,xmax,ymin,ymax;
  int imin,imax,jmin,jmax;
  float x,y;
  float sky,sig;
  struct headfits hwdb;
  int   nx,ny;
  int nobj=0;
  float nobj_vs_sig[nvarsig],sig_vs_nobj[nvarsig];
  struct obj *objects[nvarsig];
  static struct obj *objdeblen;
  struct obj *objparent;
  float apcut,xcut,ycut,mindist;
 
  float relip,robj;
  int iherit;

  int ntotobj=0;

  float significancy;
  int minpix_this_sig;

  float a,b,t,re;

  float tr[6];

  float graphflag=0;
  int debugflag=0;
  int interactflag=0;

  int pgid;

/*   if(objprimal.xc>820 && objprimal.xc<840 && objprimal.yc>510 && objprimal.yc<550)   debugflag=graphflag=1; */


/*   printf(" 666 Antes tamaño obj %d objerr %d float %d\n",sizeof(struct obj),sizeof(struct objerr),sizeof(float)); */

  
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;  
  if(graphflag) cpgqid(&pgid);
  /*   significancy=sqrt(minpix)*nsig; */ /* No esta muy claro */
  significancy=minpix*nsig;

/*   printf(" SI\n"); */

  memcpy(&hwdb,&hv,sizeof(struct headfits));

/*   printf(" NO666\n"); */
  if((flux=malloc(objprimal.area*sizeof(float))) == NULL) {
    printf("BuhoFits: ERROR. Cannot dimension flux of %d bytes\n ",objprimal.area*sizeof(float)) ;
    exit(1);
  }
  memcpy(flux,objprimal.f,objprimal.area*sizeof(float));
  
  qsort(flux,objprimal.area,sizeof(float),icompare);
  fits_pv(imasig,&hsig,objprimal.xc,objprimal.yc,&sig);
  
  maxnsig=(flux[objprimal.area-minpix])/objprimal.sig;
  if(debugflag) printf("\n %f %f %f %f\n",flux[objprimal.area-minpix],objprimal.sig,objprimal.sky,maxnsig); 
  
  MinMax(objprimal.area,objprimal.xobj,&xmin,&xmax);
  MinMax(objprimal.area,objprimal.yobj,&ymin,&ymax);
  imin=(int)((xmin-hv.crval1)/hv.cdelt1)-1;
  imax=(int)((xmax-hv.crval1)/hv.cdelt1)+2;
  jmin=(int)((ymin-hv.crval2)/hv.cdelt2)-1;
  jmax=(int)((ymax-hv.crval2)/hv.cdelt2)+2;
  if(debugflag)  printf("X  %f-%f Y %f-%f\n",xmin,xmax,ymin,ymax); 
  if(debugflag)  printf("I  %d-%d J %d-%d\n",imin,imax,jmin,jmax); 
  if(imin<0) imin=0;
  if(jmin<0) jmin=0;
  if(imax>hv.naxis1-1) imax=hv.naxis1-1;
  if(jmax>hv.naxis2-1) jmax=hv.naxis2-1;
  

  nx=imax-imin+1;
  ny=jmax-jmin+1;

  hwdb.naxis1=nx;
  hwdb.naxis2=ny;
  hwdb.crval1=hv.crval1+imin*hv.cdelt1;
  hwdb.crval2=hv.crval2+jmin*hv.cdelt2;
  if(debugflag)   printf(" %d-%d %d-%d\n",imin,imax,jmin,jmax); 
  
  if((windebl=malloc(nx*ny*sizeof(float)))==NULL) {
    printf("I cannot dimension windebl of %d elements \n",nx*ny);
    exit(1);
  }


  
  for(isig=0;isig<nvarsig;isig++) {
    if((objects[isig]=malloc(1*sizeof(struct obj)))==NULL) {
      printf("I cannot dimension objects[%d] of %d elements \n",isig,sizeof(struct obj));
      exit(1);
    }
    newsig=exp(isig*(log(maxnsig)-log(nsig))/nvarsig+log(nsig));
    /*     minpix_this_sig=(int)((significancy/newsig)*(significancy/newsig));  */  /* No lo tengo muy claro */
    minpix_this_sig=(int)(significancy/newsig);
    if( minpix_this_sig<2) minpix_this_sig=2;
    if(debugflag)    printf(" SIG %f NPIX %d\n",newsig,minpix_this_sig);
    
    for(i=0;i<nx;i++)  for(j=0;j<ny;j++) windebl[i+j*nx]=-1.0;
    for(k=0;k<objprimal.area;k++) {
      i=-imin+objprimal.ipix[k];
      j=-jmin+objprimal.jpix[k];
      x=(imin+i)*hv.cdelt1+hv.crval1;
      y=(jmin+j)*hv.cdelt2+hv.crval2;
      fits_pv(imasky,&hsky,x,y,&sky);
      fits_pv(imasig,&hsig,x,y,&sig);
      
      if(ima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1] < sky+ newsig*sig)  windebl[i+j*nx]=-0.9;
      else windebl[i+j*nx]=ima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1]-sky;
    }

    nobj=0;    

    if(graphflag) {
      cpgopen("/xserve"); 
      cpgwnad(0,nx+1,0,ny+1);
    }


/*     printf(" ESTA ES LA %d\n",isig); */
/* 	  cpgwnad(0,nx,0,ny); */
/* 	  cpggray(windebl,nx,ny,1,nx,1,ny,0.,-1.0,tr);  */
/* 	  cpgbox("bctns",0.0,0,"bctns",0.0,0);  */
/* 	  cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin-1)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin-1)*hv.cdelt2+hv.crval2); */



    /*!!!!!!!!!!!!!!!!!!!!PASARLO SOLO SOBRE LOS PIXELES DEL OBJETO
      PRIMORDIAL!!!!!! LOS DEMAS YA SE SABEN QUE VAN A DAR Detecta CERO*/
    for(j=0;j<ny;j++) {
      for(i=0;i<nx;i++) {
/* 	printf(" 666 Asigno isig %d  nobj %d\n",isig,nobj); */
	objects[isig][nobj].area=0;
	objects[isig][nobj].nb=0;
/* 	printf(" 666 se mete en detecta y apunta %d  %d\n",&(objects[isig][nobj]),objects[isig]); */
	Detecta(windebl,&hwdb,i,j,(&(objects[isig][nobj]))); 
/* 	printf(" 666 ya uqi\n"); */
	if(objects[isig][nobj].area > minpix_this_sig)  {
/* 	  printf(" 666 calculos\n"); */
	  CFM(&(objects[isig][nobj]));
 	  AjustaElipse(&(objects[isig][nobj])); 
	  nobj++;
	  if(nobj>MAXOBJECTS) nobj--;
/* 	  printf(" 666 realocateo\n"); */
	  if((objects[isig]=realloc(objects[isig],(nobj+1)*sizeof(struct obj)))==NULL) {
	    printf("I cannot dimension objects[%d] of %d elements \n",isig,nobj*sizeof(struct obj));
	    exit(1);
	  }

/* 	  printf(" 666 Realocateado isig %d: %d  nobj %d\n",isig,objects[isig],nobj); */

	}
      }
    }

/*     printf(" 666 Salio detecta\n"); */


    nobj_vs_sig[isig]=nobj;
    
    if(debugflag)      printf("\nISIG %d NSIG %f/%f NOBJ  %f\n",isig,newsig,maxnsig,nobj_vs_sig[isig]); 
    
    if(graphflag) {
      cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin-1)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin-1)*hv.cdelt2+hv.crval2);
      cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin-1)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin-1)*hv.cdelt2+hv.crval2);
      cpgbox("bctns",0.0,0,"bctns",0.0,0); 
      cpglab("X\\dWORLD\\u","Y\\dWORLD\\u","DEBLENDING");
      for(iobj=0;iobj<nobj;iobj++) {
	cpgsci(iobj+3);
	cpgpt(objects[isig][iobj].area,objects[isig][iobj].xobj,objects[isig][iobj].yobj,1);
	cpgpt(objects[isig][iobj].nb,objects[isig][iobj].xb,objects[isig][iobj].yb,2);
	cpgpt(1,&(objects[isig][iobj].xc),&(objects[isig][iobj].yc),6);
	cpgpt(1,&(objects[isig][iobj].xcelip),&(objects[isig][iobj].ycelip),7);
	cpgelip(objects[isig][iobj].xcelip,objects[isig][iobj].ycelip,objects[isig][iobj].elipa,objects[isig][iobj].elipb,objects[isig][iobj].elipt);
	cpgsci(1);
      }
      
      cpgclos();     
    }  
  }
  
  if(debugflag) printf(" Asigno las sigmas\n");

  for(isig=0;isig<nvarsig;isig++) {
    newsig=exp(isig*(log(maxnsig)-log(nsig))/nvarsig+log(nsig));
    sig_vs_nobj[isig]=newsig;
  }
  
  if(graphflag) {
    cpgopen("/xserve"); 
/*     cpgask(0); */
    cpgwnad(0,nx+1,0,ny+1); 
    
    /*     cpggray(windebl,nx,ny,1,nx,1,ny,00.,-2.,tr);   */
    cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin-1)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin-1)*hv.cdelt2+hv.crval2);
    cpgbox("bctns",0.0,0,"bctns",0.0,0);   
  }
/*   cpgline(nvarsig,sig_vs_nobj,nobj_vs_sig); */

  for(iobj=0;iobj<nobj_vs_sig[0];iobj++) {
    objects[0][iobj].heritage=-1;
    objects[0][iobj].nheritaged=0; 
  }

  for(isig=1;isig<nvarsig;isig++) {
    for(iobj=0;iobj<nobj_vs_sig[isig];iobj++) {
      if(debugflag) printf("####ISIG %d IOBJ %d\n",isig,iobj);
      if(graphflag) {
	cpgpage();
	cpgsch(3);
	cpgsci(1);
	cpgpt1(objects[isig][iobj].xcelip,objects[isig][iobj].ycelip,7);
      }
      mindist=1.e38;
      iherit=-1;
      objects[isig][iobj].nheritaged=0; 
      for(iobj2=0;iobj2<nobj_vs_sig[isig-1];iobj2++) {
	if(graphflag) {
	  cpgsci(iobj2+3);
	  cpgelip(objects[isig-1][iobj2].xcelip,objects[isig-1][iobj2].ycelip,objects[isig-1][iobj2].elipa,objects[isig-1][iobj2].elipb,objects[isig-1][iobj2].elipt);
	  cpgpt1(objects[isig-1][iobj2].xcelip,objects[isig-1][iobj2].ycelip,9);
	}
 	apcut=atan2(objects[isig][iobj].ycelip-objects[isig-1][iobj2].ycelip,objects[isig][iobj].xcelip-objects[isig-1][iobj2].xcelip);
	a=objects[isig-1][iobj2].elipa;
	b=objects[isig-1][iobj2].elipb;
	t=objects[isig-1][iobj2].elipt/180*M_PI;
	re=a*b/(sqrt(a*a*sin(apcut-t)*sin(apcut-t)+b*b*cos(apcut-t)*cos(apcut-t)));
	xcut=+re*cos(apcut-t)*cos(t)-re*sin(apcut-t)*sin(t); 
	ycut=+re*cos(apcut-t)*sin(t)+re*sin(apcut-t)*cos(t); 
/* 	if(debugflag) printf("comp %f\n",(xcut*cos(t)+y*sin(t))*(xcut*cos(t)+y*sin(t))/a/a+(xcut*sin(t)-y*cos(t))*(xcut*sin(t)-y*cos(t))/b/b); */
	xcut=xcut+objects[isig-1][iobj2].xcelip; 
	ycut=ycut+objects[isig-1][iobj2].ycelip; 
	if(graphflag) 	cpgpt1(xcut,ycut,8);
	if(debugflag) printf("isig %d i1 %d i2 %d xc1 %f xc2 %f yc1 %f yc2 %f  apcut %f xcut %f ycut %f \n",isig,iobj,iobj2,objects[isig][iobj].xcelip,objects[isig-1][iobj2].xcelip,objects[isig][iobj].ycelip,objects[isig-1][iobj2].ycelip,apcut*180/M_PI,xcut,ycut);
	relip=sqrt((objects[isig-1][iobj2].xcelip-xcut)*(objects[isig-1][iobj2].xcelip-xcut)+(objects[isig-1][iobj2].ycelip-ycut)*(objects[isig-1][iobj2].ycelip-ycut));
	robj =sqrt((objects[isig-1][iobj2].xcelip-objects[isig][iobj].xcelip)*(objects[isig-1][iobj2].xcelip-xcut)+(objects[isig-1][iobj2].ycelip-ycut)*(objects[isig-1][iobj2].ycelip-objects[isig][iobj].ycelip));
	if(debugflag) printf(" relip %f robj %f\n",relip,robj);
	if(mindist> fabs(robj/relip) && fabs(robj/relip)< 1.5 ) {
	  mindist=fabs(robj/relip);
	  iherit=iobj2;
	  
	  if(debugflag) printf("MIN  %d d %f\n",iobj2,mindist);
	}
	if(iherit>-1) {
	  if(debugflag) printf(" Voy a pintar i %d\n",iherit);
	  if(graphflag) {
	    cpgsci(1);
	    cpgpt1(objects[isig-1][iherit].xcelip,objects[isig-1][iherit].ycelip,2);
	  }
	}
	
      }
      
      if(debugflag) printf(" Asigno iheri %d \n",iherit);
      
      if(iherit!=-1) {
	objects[isig][iobj].heritage=iherit;
	objects[isig][iobj].objher=&(objects[isig-1][iherit]);
      }
      else {
	objects[isig][iobj].heritage=-1;
	objects[isig][iobj].objher=NULL;
      }

/*       isig=readi(isig); */

    }
  }

  if(debugflag) printf(" Salio de heritage\n");
  
  if(graphflag) cpgclos();



  
  if(graphflag) {
    cpgopen("/xserve");
    cpgswin(-1.,5.,-1.,nvarsig+1.);
    cpgbox("bctns",0.0,0,"bctns",0.0,0);
  }

  


  for(isig=1;isig<nvarsig;isig++) {
    for(iobj=0;iobj<nobj_vs_sig[isig];iobj++) {
      if(graphflag) cpgsci(isig);
      if(objects[isig][iobj].heritage!=-1) {
	objects[isig-1][objects[isig][iobj].heritage].nheritaged++;
	if(graphflag) {
	  cpgpt1((float)iobj,(float)isig,5);
	  cpgdraw((float)(objects[isig][iobj].heritage),(float)(isig-1));
	}
      }
      objects[isig][iobj].herflag=1;
      objects[isig][iobj].iobj=iobj;
      objects[isig][iobj].isig=isig;

      }
  }

 for(isig=0;isig<nvarsig;isig++) {
    for(iobj=0;iobj<nobj_vs_sig[isig];iobj++) {
      if(debugflag)  printf(" ISIG %d IOBJ %d NHERIT %d HERITAGE %d\n",isig,iobj,objects[isig][iobj].nheritaged,objects[isig][iobj].heritage);
    }
  }


 if(debugflag)  printf(" Mira parientes\n");
 
 
/*  isig=readi(isig); */
 if(graphflag) {
   cpgsci(1);
   cpgsch(1.);
 }
 
 
 if((objdeblen=malloc(1*sizeof(struct obj)))==NULL) {
   printf("I cannot dimension objdeblen of %d elements \n",sizeof(struct obj));
   exit(1);
 }
 for(isig=nvarsig-1;isig>=0;isig--) { 
   for(iobj=0;iobj<nobj_vs_sig[isig];iobj++) { 
     if(debugflag)   printf(" MIRANDO ARBOL De %d %d\n",iobj,isig);
     if(objects[isig][iobj].herflag) { 
       if(debugflag)    printf(" Mem eto en Look\n");
       if(LookParent(&(objects[isig][iobj]),&(objparent))) {
	 if(objects[isig][iobj].nheritaged<=1) {
	   if(debugflag)  printf(" salido  el padre era isig %d iobj %d\n",(*objparent).isig,(*objparent).iobj);
	   ntotobj++;
	   if((objdeblen=realloc(objdeblen,ntotobj*sizeof(struct obj)))==NULL) {
	     printf("I cannot dimension objdeblen of %d elements \n",ntotobj*sizeof(struct obj));
	     exit(1);
	   }
	   memcpy(objdeblen+ntotobj-1,objparent,sizeof(struct obj));
	   if(interactflag) ntotobj=readi(ntotobj); 
	   objdeblen[ntotobj-1].ideb=ntotobj;
	 }
       }
     }
   }
 }
 
 if(debugflag)  printf(" TERMINO LOS PARIENTES\n");
 
 if(graphflag) cpgclos(); 
 
 if(graphflag) {
   cpgopen("/xserve"); 
   for(k=0;k<objprimal.area;k++) {
     i=-imin+objprimal.ipix[k];
     j=-jmin+objprimal.jpix[k];
     x=(imin+i)*hv.cdelt1+hv.crval1;
     y=(jmin+j)*hv.cdelt2+hv.crval2;
     fits_pv(imasky,&hsky,x,y,&sky);
     fits_pv(imasig,&hsig,x,y,&sig);
     if(ima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1] < sky+ nsig*sig)  windebl[i+j*nx]=-0.9;
     else windebl[i+j*nx]=ima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1]-sky;
   }
   cpgwnad(0,nx+1,0,ny+1);
   cpggray(windebl,nx,ny,1,nx,1,ny,sig*maxnsig,-3,tr);  
   cpgbox("bctns",0.0,0,"bctns",0.0,0);  
   cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin)*hv.cdelt2+hv.crval2);
   
 }
 
 
 for(iobj=0;iobj<ntotobj;iobj++) {
   
   if(debugflag)  printf("  FINAL OBJECT: %d  xc %f yc %f sig %d  \n",iobj,objdeblen[iobj].xc,objdeblen[iobj].yc,objdeblen[iobj].isig);
   if(graphflag) {
     cpgsci(iobj+2);
     cpgpt(objdeblen[iobj].area,objdeblen[iobj].xobj,objdeblen[iobj].yobj,1);
     cpgpt(objdeblen[iobj].nb,objdeblen[iobj].xb,objdeblen[iobj].yb,2);
     cpgpt(1,&(objdeblen[iobj].xcelip),&(objdeblen[iobj].ycelip),7);
     cpgelip(objdeblen[iobj].xcelip,objdeblen[iobj].ycelip,objdeblen[iobj].elipa,objdeblen[iobj].elipb,objdeblen[iobj].elipt);
   }
 }
 if(graphflag){
   
   cpgclos();
 }
 
 if(interactflag) ntotobj=readi(ntotobj); 

 
 *ndeb=ntotobj;
 *objsdeb=objdeblen;
 
 
 free(windebl);
 free(flux);
 for(isig=0;isig<nvarsig;isig++)    free(objects[isig]);
 

 return(ntotobj);
 
}


int errDeblending(float *ima, float *errima, struct headfits hv, float *imasky, struct headfits hsky,float *imasig, struct headfits hsig, float *imaerrsky, struct headfits herrsky, struct obj objprimal, float nsig,struct obj **objsdeb, struct objerr **errobjsdeb, int *ndeb) 

 { 
  
  
  
  float maxnsig;


  int   i,j,k;
  int   iobj,iobj2;
  float newsig;
  int isig=0;
  float *flux;
  float *windebl;
  float *winerrdebl;
  float xmin,xmax,ymin,ymax;
  int imin,imax,jmin,jmax;
  float x,y;
  float sky,sig,errsky;
  struct headfits hwdb;
  int   nx,ny;
  int nobj=0;
  float nobj_vs_sig[nvarsig],sig_vs_nobj[nvarsig];
  struct obj *objects[nvarsig];
  struct objerr *errobjects[nvarsig]; 
  struct obj *objdeblen; 
  static struct objerr *errobjdeblen; 
  static struct obj *objparent; 
  struct objerr *errobjparent; 
  float apcut,xcut,ycut,mindist;
 
  float relip,robj;
  int iherit;

  int ntotobj=0;

  float significancy;
  int minpix_this_sig;

  float a,b,t,re;

  float tr[6];

  float graphflag=0;
  int debugflag=1;
  int interactflag=0;

  int pgid;




  
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;  
  if(graphflag) cpgqid(&pgid);
  significancy=sqrt(minpix)*nsig;
  memcpy(&hwdb,&hv,sizeof(struct headfits));
  if((flux=malloc(objprimal.area*sizeof(float))) == NULL) {
    printf("BuhoFits: ERROR. Cannot dimension flux of %d bytes\n ",objprimal.area*sizeof(float)) ;
    exit(1);
  }
  memcpy(flux,objprimal.f,objprimal.area*sizeof(float));
  
  qsort(flux,objprimal.area,sizeof(float),icompare);
  fits_pv(imasig,&hsig,objprimal.xc,objprimal.yc,&sig);
  
  maxnsig=(flux[objprimal.area-minpix])/objprimal.sig;
  if(debugflag) printf("\n %f %f %f %f\n",flux[objprimal.area-minpix],objprimal.sig,objprimal.sky,maxnsig); 
  
  pgLimits(objprimal.area,objprimal.xobj,&xmin,&xmax);
  pgLimits(objprimal.area,objprimal.yobj,&ymin,&ymax);
  imin=(int)((xmin-hv.crval1)/hv.cdelt1)-1;
  imax=(int)((xmax-hv.crval1)/hv.cdelt1)+2;
  jmin=(int)((ymin-hv.crval2)/hv.cdelt2)-1;
  jmax=(int)((ymax-hv.crval2)/hv.cdelt2)+2;
  if(debugflag)  printf("X  %f-%f Y %f-%f\n",xmin,xmax,ymin,ymax); 
  if(debugflag)  printf("I  %d-%d J %d-%d\n",imin,imax,jmin,jmax); 
  if(imin<0) imin=0;
  if(jmin<0) jmin=0;
  if(imax>hv.naxis1-1) imax=hv.naxis1-1;
  if(jmax>hv.naxis2-1) jmax=hv.naxis2-1;
  

  nx=imax-imin+1;
  ny=jmax-jmin+1;

  hwdb.naxis1=nx;
  hwdb.naxis2=ny;
  hwdb.crval1=hv.crval1+imin*hv.cdelt1;
  hwdb.crval2=hv.crval2+jmin*hv.cdelt2;


  if(debugflag)   printf(" %d-%d %d-%d\n",imin,imax,jmin,jmax); 
  
  if((windebl=   malloc(nx*ny*sizeof(float)))==NULL) {
    printf("I cannot dimension windebl of %d elements \n",nx*ny);
    exit(1);
  }
  if((winerrdebl=malloc(nx*ny*sizeof(float)))==NULL) {
    printf("I cannot dimension winerrdebl of %d elements \n",nx*ny);
    exit(1);
  }
  
  for(isig=0;isig<nvarsig;isig++) {
    printf(" isig %d ANTES ALLOC OBJECTS\n",isig);
    if((objects[isig]=malloc(1*sizeof(struct obj)))==NULL) {
      printf("I cannot dimension objects[%d] of %d elements \n",isig,sizeof(struct obj));
      exit(1);
    }
    if((errobjects[isig]=malloc(1*sizeof(struct objerr)))==NULL) {
      printf("I cannot dimension errobjects[%d] of %d elements \n",isig,sizeof(struct objerr));
      exit(1);
    }

    newsig=exp(isig*(log(maxnsig)-log(nsig))/nvarsig+log(nsig));
    minpix_this_sig=(int)((significancy/newsig)*(significancy/newsig));
    if( minpix_this_sig<2) minpix_this_sig=2;
    if(debugflag)    printf(" SIG %f NPIX %d\n",newsig,minpix_this_sig);
    
    for(i=0;i<nx;i++)  for(j=0;j<ny;j++) windebl[i+j*nx]=-1.0;
    for(i=0;i<nx;i++)  for(j=0;j<ny;j++) winerrdebl[i+j*nx]=-1.0;
    for(k=0;k<objprimal.area;k++) {
      i=-imin+objprimal.ipix[k];
      j=-jmin+objprimal.jpix[k];
      x=(imin+i)*hv.cdelt1+hv.crval1;
      y=(jmin+j)*hv.cdelt2+hv.crval2;
      fits_pv(imasky,&hsky,x,y,&sky);
      fits_pv(imasig,&hsig,x,y,&sig);
      
      if(ima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1] < sky+ newsig*sig) {
	windebl[i+j*nx]=-0.9;
	winerrdebl[i+j*nx]=0.;
      }
      else {
	windebl[i+j*nx]=ima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1]-sky;
	fits_pv(imaerrsky,&herrsky,x,y,&errsky);
	winerrdebl[i+j*nx] = sqrt(errima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1]*errima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1]+errsky*errsky);
      }
    }

    nobj=0;    


    if(graphflag) {
      cpgopen("/xserve"); 
      cpgwnad(0,nx+1,0,ny+1);
    }


    if(debugflag) printf(" Antes realloc\n");


    /*!!!!!!!!!!!!!!!!!!!!PASARLO SOLO SOBRE LOS PIXELES DEL OBJETO
      PRIMORDIAL!!!!!! LOS DEMAS YA SE SABEN QUE VAN A DAR Detecta CERO*/

    for(j=0;j<ny;j++) {
      for(i=0;i<nx;i++) {
	objects[isig][nobj].area=0;
	objects[isig][nobj].nb=0;
	errDetecta(windebl,winerrdebl,&hwdb,i,j,(&(objects[isig][nobj])),(&(errobjects[isig][nobj]))); 
	if(objects[isig][nobj].area > minpix_this_sig)  {
	  if(debugflag) printf(" Antes CFM\n");
	  CFM(&(objects[isig][nobj]));
	  if(debugflag) printf(" despues cfm nobj %d isig %d\n",nobj,isig);
 	  AjustaElipse(&(objects[isig][nobj])); 
	  if(debugflag) printf(" despues ajusta\n");
	  nobj++;
	  if(nobj>MAXOBJECTS) nobj--;
	  else {
	    if(debugflag) printf(" inmedia realloc nobj %d\n",nobj);
	    if((objects[isig]=realloc(objects[isig],(nobj+1)*sizeof(struct obj)))==NULL) {
	      printf("I cannot dimension objects[%d] of %d elements \n",isig,nobj*sizeof(struct obj));
	      exit(1);
	    }
	    if((errobjects[isig]=realloc(errobjects[isig],(nobj+1)*sizeof(struct objerr)))==NULL) {
	      printf("I cannot dimension errobjects[%d] of %d elements \n",isig,nobj*sizeof(struct objerr));
	      exit(1);
	    }
	  }
	}
      }
    }

    if(debugflag) printf(" Despeus realloc\n");

    nobj_vs_sig[isig]=nobj;
    
    if(debugflag)      printf("\nISIG %d NSIG %f/%f NOBJ  %f\n",isig,newsig,maxnsig,nobj_vs_sig[isig]); 
    
    if(graphflag) {
      cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin-1)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin-1)*hv.cdelt2+hv.crval2);
      cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin-1)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin-1)*hv.cdelt2+hv.crval2);
      cpgbox("bctns",0.0,0,"bctns",0.0,0); 
      cpglab("X\\dWORLD\\u","Y\\dWORLD\\u","DEBLENDING");
      for(iobj=0;iobj<nobj;iobj++) {
	cpgsci(iobj+3);
	cpgpt(objects[isig][iobj].area,objects[isig][iobj].xobj,objects[isig][iobj].yobj,1);
	cpgpt(objects[isig][iobj].nb,objects[isig][iobj].xb,objects[isig][iobj].yb,2);
	cpgpt(1,&(objects[isig][iobj].xc),&(objects[isig][iobj].yc),6);
	cpgpt(1,&(objects[isig][iobj].xcelip),&(objects[isig][iobj].ycelip),7);
	cpgelip(objects[isig][iobj].xcelip,objects[isig][iobj].ycelip,objects[isig][iobj].elipa,objects[isig][iobj].elipb,objects[isig][iobj].elipt);
	cpgsci(1);
      }
      
      cpgclos();     
    }  
  }

  if(debugflag) printf(" Asigno las sigmas\n");

  for(isig=0;isig<nvarsig;isig++) {
    newsig=exp(isig*(log(maxnsig)-log(nsig))/nvarsig+log(nsig));
    sig_vs_nobj[isig]=newsig;
  }
  
  if(graphflag) {
    cpgopen("/xserve"); 
    cpgask(0);
    cpgwnad(0,nx+1,0,ny+1); 
    
    cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin-1)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin-1)*hv.cdelt2+hv.crval2);
    cpgbox("bctns",0.0,0,"bctns",0.0,0);   
  }

  for(iobj=0;iobj<nobj_vs_sig[0];iobj++) {
    objects[0][iobj].heritage=-1;
    objects[0][iobj].nheritaged=0; 
  }

  for(isig=1;isig<nvarsig;isig++) {
    for(iobj=0;iobj<nobj_vs_sig[isig];iobj++) {
      if(debugflag) printf("####ISIG %d IOBJ %d\n",isig,iobj);
      if(graphflag) {
	cpgpage();
	cpgsch(3);
	cpgsci(1);
	cpgpt1(objects[isig][iobj].xcelip,objects[isig][iobj].ycelip,7);
      }
      mindist=1.e38;
      iherit=-1;
      objects[isig][iobj].nheritaged=0; 
      for(iobj2=0;iobj2<nobj_vs_sig[isig-1];iobj2++) {
	if(graphflag) {
	  cpgsci(iobj2+3);
	  cpgelip(objects[isig-1][iobj2].xcelip,objects[isig-1][iobj2].ycelip,objects[isig-1][iobj2].elipa,objects[isig-1][iobj2].elipb,objects[isig-1][iobj2].elipt);
	  cpgpt1(objects[isig-1][iobj2].xcelip,objects[isig-1][iobj2].ycelip,9);
	}
 	apcut=atan2(objects[isig][iobj].ycelip-objects[isig-1][iobj2].ycelip,objects[isig][iobj].xcelip-objects[isig-1][iobj2].xcelip);
	a=objects[isig-1][iobj2].elipa;
	b=objects[isig-1][iobj2].elipb;
	t=objects[isig-1][iobj2].elipt/180*M_PI;
	re=a*b/(sqrt(a*a*sin(apcut-t)*sin(apcut-t)+b*b*cos(apcut-t)*cos(apcut-t)));
	xcut=+re*cos(apcut-t)*cos(t)-re*sin(apcut-t)*sin(t); 
	ycut=+re*cos(apcut-t)*sin(t)+re*sin(apcut-t)*cos(t); 
/* 	if(debugflag) printf("comp %f\n",(xcut*cos(t)+y*sin(t))*(xcut*cos(t)+y*sin(t))/a/a+(xcut*sin(t)-y*cos(t))*(xcut*sin(t)-y*cos(t))/b/b); */
	xcut=xcut+objects[isig-1][iobj2].xcelip; 
	ycut=ycut+objects[isig-1][iobj2].ycelip; 
	if(graphflag) 	cpgpt1(xcut,ycut,8);
	if(debugflag) printf("isig %d i1 %d i2 %d xc1 %f xc2 %f yc1 %f yc2 %f  apcut %f xcut %f ycut %f \n",isig,iobj,iobj2,objects[isig][iobj].xcelip,objects[isig-1][iobj2].xcelip,objects[isig][iobj].ycelip,objects[isig-1][iobj2].ycelip,apcut*180/M_PI,xcut,ycut);
	relip=sqrt((objects[isig-1][iobj2].xcelip-xcut)*(objects[isig-1][iobj2].xcelip-xcut)+(objects[isig-1][iobj2].ycelip-ycut)*(objects[isig-1][iobj2].ycelip-ycut));
	robj =sqrt((objects[isig-1][iobj2].xcelip-objects[isig][iobj].xcelip)*(objects[isig-1][iobj2].xcelip-xcut)+(objects[isig-1][iobj2].ycelip-ycut)*(objects[isig-1][iobj2].ycelip-objects[isig][iobj].ycelip));
	if(debugflag) printf(" relip %f robj %f\n",relip,robj);
	if(mindist> fabs(robj/relip) && fabs(robj/relip)< 1.5 ) {
	  mindist=fabs(robj/relip);
	  iherit=iobj2;
	  
	  if(debugflag) printf("MIN  %d d %f\n",iobj2,mindist);
	}
	if(iherit>-1) {
	  if(debugflag) printf(" Voy a pintar i %d\n",iherit);
	  if(graphflag) {
	    cpgsci(1);
	    cpgpt1(objects[isig-1][iherit].xcelip,objects[isig-1][iherit].ycelip,2);
	  }
	}
	
      }
      
      if(debugflag) printf(" Asigno iheri %d \n",iherit);
      
      if(iherit!=-1) {
	objects[isig][iobj].heritage=iherit;
	objects[isig][iobj].objher=&(objects[isig-1][iherit]);
	errobjects[isig][iobj].errobjher=&(errobjects[isig-1][iherit]);
      }
      else {
	objects[isig][iobj].heritage=-1;
	objects[isig][iobj].objher=NULL;
	errobjects[isig][iobj].errobjher=NULL;
      }

        if(interactflag) isig=readi(isig); 

    }
  }

  if(debugflag) printf(" Salio de heritage\n");
  
  if(graphflag) cpgclos();



  
  if(graphflag) {
    cpgopen("/xserve");
    cpgswin(-1.,5.,-1.,nvarsig+1.);
    cpgbox("bctns",0.0,0,"bctns",0.0,0);
  }

  


  for(isig=1;isig<nvarsig;isig++) {
    for(iobj=0;iobj<nobj_vs_sig[isig];iobj++) {
      if(graphflag) cpgsci(isig);
      if(objects[isig][iobj].heritage!=-1) {
	objects[isig-1][objects[isig][iobj].heritage].nheritaged++;
	if(graphflag) {
	  cpgpt1((float)iobj,(float)isig,5);
	  cpgdraw((float)(objects[isig][iobj].heritage),(float)(isig-1));
	}
      }
      objects[isig][iobj].herflag=1;
      objects[isig][iobj].iobj=iobj;
      objects[isig][iobj].isig=isig;

      }
  }

 for(isig=0;isig<nvarsig;isig++) {
    for(iobj=0;iobj<nobj_vs_sig[isig];iobj++) {
      if(debugflag)  printf(" ISIG %d IOBJ %d NHERIT %d HERITAGE %d\n",isig,iobj,objects[isig][iobj].nheritaged,objects[isig][iobj].heritage);
    }
  }


 if(debugflag)  printf(" Mira parientes\n");
 
 
 if(graphflag) {
   cpgsci(1);
   cpgsch(1.);
 }
 
 
 if((objdeblen=malloc(1*sizeof(struct obj)))==NULL) {
   printf("I cannot dimension objdeblen of %d elements \n",sizeof(struct obj));
   exit(1);
 }
 if((errobjdeblen=malloc(1*sizeof(struct objerr)))==NULL) {
   printf("I cannot dimension errobjdeblen of %d elements \n",sizeof(struct objerr));
   exit(1);
 }

 for(isig=nvarsig-1;isig>=0;isig--) { 
   for(iobj=0;iobj<nobj_vs_sig[isig];iobj++) { 
     if(debugflag)   printf(" MIRANDO ARBOL De %d %d\n",iobj,isig);
     if(objects[isig][iobj].herflag) { 
       if(debugflag)    printf(" Mem eto en Look\n");
       if(errLookParent(&(objects[isig][iobj]),&(errobjects[isig][iobj]),&(objparent),&(errobjparent))) {
	 if(objects[isig][iobj].nheritaged<=1) {
	   if(debugflag)  printf(" salido  el padre era isig %d iobj %d\n",(*objparent).isig,(*objparent).iobj);
	   ntotobj++;
	   if((objdeblen=realloc(objdeblen,ntotobj*sizeof(struct obj)))==NULL) {
	     printf("I cannot dimension objdeblen of %d elements \n",ntotobj*sizeof(struct obj));
	     exit(1);
	   }
	   memcpy(objdeblen+ntotobj-1,objparent,sizeof(struct obj));
	   if((errobjdeblen=realloc(errobjdeblen,ntotobj*sizeof(struct objerr)))==NULL) {
	     printf("I cannot dimension errobjdeblen of %d elements \n",ntotobj*sizeof(struct objerr));
	     exit(1);
	   }
	   memcpy(errobjdeblen+ntotobj-1,errobjparent,sizeof(struct objerr));
	   if(interactflag) ntotobj=readi(ntotobj); 
	   objdeblen[ntotobj-1].ideb=ntotobj;
	 }
       }
     }
   }
 }
 
 if(debugflag)  printf(" TERMINO LOS PARIENTES\n");
 
 if(graphflag) cpgclos(); 
 
 if(graphflag) {
   cpgopen("/xserve"); 
   for(k=0;k<objprimal.area;k++) {
     i=-imin+objprimal.ipix[k];
     j=-jmin+objprimal.jpix[k];
     x=(imin+i)*hv.cdelt1+hv.crval1;
     y=(jmin+j)*hv.cdelt2+hv.crval2;
     fits_pv(imasky,&hsky,x,y,&sky);
     fits_pv(imasig,&hsig,x,y,&sig);
     if(ima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1] < sky+ nsig*sig)  windebl[i+j*nx]=-0.9;
     else windebl[i+j*nx]=ima[objprimal.ipix[k]+objprimal.jpix[k]*hv.naxis1]-sky;
   }
   cpgwnad(0,nx+1,0,ny+1);
   cpggray(windebl,nx,ny,1,nx,1,ny,sig*maxnsig,-3,tr);  
   cpgbox("bctns",0.0,0,"bctns",0.0,0);  
   cpgswin((imin-1)*hv.cdelt1+hv.crval1,(nx+imin)*hv.cdelt1+hv.crval1,(jmin-1)*hv.cdelt2+hv.crval2,(ny+jmin)*hv.cdelt2+hv.crval2);
   
 }
 
 
 for(iobj=0;iobj<ntotobj;iobj++) {
   
   if(debugflag)  printf("  FINAL OBJECT: %d  xc %f yc %f sig %d  \n",iobj,objdeblen[iobj].xc,objdeblen[iobj].yc,objdeblen[iobj].isig);
   if(graphflag) {
     cpgsci(iobj+2);
     cpgpt(objdeblen[iobj].area,objdeblen[iobj].xobj,objdeblen[iobj].yobj,1);
     cpgpt(objdeblen[iobj].nb,objdeblen[iobj].xb,objdeblen[iobj].yb,2);
     cpgpt(1,&(objdeblen[iobj].xcelip),&(objdeblen[iobj].ycelip),7);
     cpgelip(objdeblen[iobj].xcelip,objdeblen[iobj].ycelip,objdeblen[iobj].elipa,objdeblen[iobj].elipb,objdeblen[iobj].elipt);
   }
 }
 if(graphflag){
   
   cpgclos();
 }
 
 if(interactflag) ntotobj=readi(ntotobj); 

 
 *ndeb=ntotobj;
 *objsdeb=objdeblen;
 *errobjsdeb=errobjdeblen;
 
 
 if(graphflag) cpgslct(pgid); 
 
 if(debugflag) printf(" antes free\n");

 free(windebl);
 free(winerrdebl);
 free(flux);
 for(isig=0;isig<nvarsig;isig++)    free(objects[isig]);
 for(isig=0;isig<nvarsig;isig++)    free(errobjects[isig]);

 if(debugflag) printf(" Despues free\n");

 return(ntotobj);


} 









int  LookParent(struct obj *O, struct obj **OP) {

  int graphflag=0;
  
/*   printf(" SLook con objeto  isig %d iobj %d puntero %d\n",(*obj).isig,(*obj).iobj,obj); */
/*   printf(" SLook con objeto anterior isig %d iobj %d\n",(*(*obj).objher).isig,(*(*obj).objher).iobj); */

  /*   printf(" Estoy que heritage %d nherit %d nherit anterior %d \n",(*obj).heritage,(*obj).nheritaged,(*(*obj).objher).nheritaged);  */

  if((*O).heritage==-1) {
/*     printf(" LookParent: PRINCIPIO DE OBJETO\n"); */
    if(graphflag) cpgpt1((*O).iobj,(*O).isig,8);
    (*O).herflag=0;
    *OP=O;
    return(0);
  }
  else if((((*(*O).objher)).nheritaged)==1) {
/*     printf(" entra Look parent isig %d iobj %d\n",(*O).isig,(*O).iobj); */
    if(graphflag)     cpgpt1((*O).iobj,(*O).isig,8);
    (*O).herflag=0;
    LookParent((*O).objher,OP); 
/*     printf(" Salio\n"); */
/*     printf(" antes de salir el padre era isig %d iobj %d\n",(*(*OP)).isig,(*(*OP)).iobj); */
    return 1;
  }
  else if((((*(*O).objher)).nheritaged)>1) {
/*     printf(" ACABA Look isig %d iobj %d\n",(*O).isig,(*O).iobj); */
   if(graphflag)    cpgpt1((*O).iobj,(*O).isig,8);
    (*O).herflag=0;
    *OP=O;
/*     printf(" esta saliendo y asigno padre a isig %d iobj %d\n",(*(*OP)).isig,(*(*OP)).iO); */
/*     printf(" OP %d\n",OP); */

    return 1;
  }
  else if((((*(*O).objher)).nheritaged)==0) {
/*     printf(" Se corta la cadena Look isig %d iobj %d\n",(*O).isig,(*O).iobj); */
    (*O).herflag=0;
    *OP=O;
    return(0);
  }
  else if(((*(*O).objher).nheritaged)==-1) {
/*     printf(" LookParent: No se que pasa en este caso\n"); */
  }
  (*O).herflag=0;
  return(0);
}


int  errLookParent(struct obj *O, struct objerr *EO, struct obj **OP,struct objerr **EOP) {

  int graphflag=0;
  if((*O).heritage==-1) {
    if(graphflag) cpgpt1((*O).iobj,(*O).isig,8);
    (*O).herflag=0;
    *OP=O;
    *EOP=EO;
    return(0);
  }
  else if((((*(*O).objher)).nheritaged)==1) {
    if(graphflag)     cpgpt1((*O).iobj,(*O).isig,8);
    (*O).herflag=0;
    errLookParent((*O).objher,(*EO).errobjher,OP,EOP); 
    return 1;
  }
  else if((((*(*O).objher)).nheritaged)>1) {
   if(graphflag)    cpgpt1((*O).iobj,(*O).isig,8);
    (*O).herflag=0;
    *OP=O;
    *EOP=EO;
    return 1;
  }
  else if((((*(*O).objher)).nheritaged)==0) {
    (*O).herflag=0;
    *OP=O;
    *EOP=EO;
    return(0);
  }
  else if(((*(*O).objher).nheritaged)==-1) {
  }
  (*O).herflag=0;
  return(0);
}
