#include "modulos.h"
#define DISP Disper2
/* //#include "cpgplot.h" */
/* Esta version 8 lo que hace es calcular los pixelss del CCD que se ven contribuidos por el obejto. Solo integra en esa zona 
 */
#define DEBUG 0

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
  float magmin;
  float magmax;
  int nmag;
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
  float mag;
  float fwhm;
  float reff;
  float z;
  float AP;
  float excen;
  float x;
  float y;
};

struct star {
  float mag;
  float teff;
};

struct other {
  float seeing;
  float ldomin;
  float ldomax;
  int nldo;
  char  sky_file[50];
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

char parfilename[100];
char imgfile[51];
char catfile[51];
char inputcat[51];
float zeropoint;
void addstepgal(float *spa,float *espec,int nx, int ny,float tamx,float tamy);
void addexpgal(float *spa,float *espec,int nx, int ny,float tamx,float tamy);
void addgaussgal(float *spa,float *espec,int nx, int ny,float tamx,float tamy);
void addgalaxy(float *img, struct gal G,int nx, int ny);
void addstar(float *spa,float *espec,int nx, int ny);
float addsky(void);
void SaveOptions(struct instr);
void LoadParam_file(void);
void LoadParam_kbd(void);
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
void Numbernet2param(int number, int *jew, int *jmag, int *jz, int *jfwhm, int *jsize,int  *jexcen,int *jAP,int *jgrid);
void Numbercat2param(int number, FILE *fileicat);
void Band2ZP(void);

int main(int argc, char **argv)
{
  float fnul,fnul2; 
/*   float alfa; */
  float xmin,xmax,xdesp,xs,ys,x,y; 
  float xs1,xs2,ys1,ys2;
/*   int *contador;  */
/*   char cnul;  */
/*   char snul[2000]; */
  float *imgccd;
  float *imgobjccd;
  float ldo,tamx,tamy;
  float psi,eta;
  int ipix,jpix;
  int ipix1,ipix2,jpix1,jpix2;
/*   int ic1,ic2,jc1,jc2; */
  int icmin,icmax,jcmin,jcmax;
  float hc=1.979e-25;  /* // J·m Julio x metro */
  float *img;
  float *direct;
  float *espec;
/*   //  float *imgp; */
  int nx,ny;
/*   FILE *fimg; */
/*   //  char fileimg[50]="salida.fits"; */
  int i,j,l;
  int ii,jj;
  int i1,j1,i2,j2;
  /*   float coef1,coef2,coef3; */  /* Para la interpolacion de Lagrange */ 
  float contri;
  float tr[6];
/*   float b=3.5e6; */
  float imin,imax=0;
  float pi=4*atan(1.);
  float skyfot;
  FILE *fileocat;
  FILE *fileicat=NULL ;
/*   float *xesp; */
/*   float *yesp; */
  float *esp1d;
  int number=0;   /* Contador del numero de galaxia por el que vamos */
  int ntot;       /*  Numero total de galaxias simuladas. */
  int jew,jmag,jz,jfwhm,jsize,jexcen,jAP,jgrid;
  
  int fromcat=0;
  /*    Para saber si los parametros de las galaxias vendran del catalogo. */
  
  int status=0;
  fitsfile *fits;
  long naxes[2];
  /********************************************
 **********************************************
 LECTURA DE PARAMETROS
 *********************
 **********************/
  srandom((unsigned int)time(NULL)/2); 
  printf("\n Welcome to SimulaOP_field Version 5\n");
  if(argc < 2) LoadParam_kbd();
  else {
    strcpy(parfilename,argv[1]);
    LoadParam_file(); 
  }
/*   printf(" Bein %s\n",catfile); */
  if((fileocat=fopen(catfile,"w"))==NULL) {
    printf("ERROR: Can't open log file %s\n",catfile);
    exit(1);
  }
  /*     printf(" AMDOJO\n"); */
  
  
  
  
  
  fprintf(fileocat,"#Catalogue file generated by  SimulaOP_field in image %s with parameter file %s\n",imgfile,argv[1]);
  fprintf(fileocat,"#Number       x         y        m      EW   fwhm    z    reff     excen    PA\n");
  Readfilter(I.filter_file);
  Readqe(I.qe_file);
  Readsky(O.sky_file);
  Readprism(I.prism_file);
  Band2ZP();
  
  
  /*   printf(" efeicnen   %f\n",I.eff); */
  
  /* Contribucion del cielo
***********************
********************/
  I.alfa_prism=-I.alfa_prism; /*  Tiene su razon de ser ; */
  skyfot=addsky();
  
/*   printf(" Elk cielo es %f\n",skyfot); */
/*   printf(" gananxia %f  eff %f  s %f\n",I.gain,I.eff, skyfot*I.texp*I.eff); */
  printf("Contribucion del cielo %e fotones/s/pixel %e ADUS\n",skyfot,skyfot*I.texp*I.eff/I.gain);
  
  /**********************************************
						 BUCLE
						 *********************************************
						 *********************************************
						 ***********************************************/
  DISP(I.xnpix*I.pixsize/I.focal/2,0.,n(O.ldomax),I.alfa_prism,&xmax,&fnul);
/*   printf(" UNO %g %g\n",xmax,fnul); */
  Disper2(I.xnpix*I.pixsize/I.focal/2,0.,n(O.ldomax),I.alfa_prism,&xmax,&fnul);
/*   printf(" DOS %g  %g\n",xmax,fnul); */
  
  
  DISP(0.,0.,n(O.ldomin),I.alfa_prism,&xmin,&fnul);  
  DISP(0.,0.,n(O.ldomax),I.alfa_prism,&xmax,&fnul);  
/*   printf("1st xmax %f xmin %f\n",xmax,xmin); */
  if(xmin>xmax) {
    printf(" Changing alfa to -alfa in order to increase wavelenght rightwards\n");
    I.alfa_prism=-I.alfa_prism;
  }
  DISP(0.,0.,n(O.ldomin),I.alfa_prism,&xmin,&fnul);  
  DISP(0.,0.,n(O.ldomax),I.alfa_prism,&xmax,&fnul);  
/*   printf("2nd xmax %f xmin %f\n",xmax,xmin); */
  
  DISP(-I.xnpix*I.pixsize/I.focal/2,0.,n(O.ldomin),I.alfa_prism,&xmin,&fnul);  
  DISP(I.xnpix*I.pixsize/I.focal/2,0.,n(O.ldomax),I.alfa_prism,&xmax,&fnul);
  
/*   printf("3rd xmax %f xmin %f\n",xmax,xmin); */
  
  xdesp=(xmax+xmin)/2.; 
  
/*   printf("xdesp %f\n",xdesp); */
  
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  
  printf(" Total amount of galaxies simulated %d\n",(FD.ngrid*FD.new*FD.nmag*FD.nz*FD.nfwhm*FD.nsize*FD.nexcen*FD.nAP));
  
  
  
  printf(" Field covered by the output image %f x %f arcsec\n",I.pixsize/I.focal*I.xnpix*206265,I.pixsize/I.focal*I.ynpix*206265); 
  printf(" Plate scale in the output image %f arcsec/pix\n",I.pixsize/I.focal*206265); 
  
  setvbuf(stdin,"",_IOLBF,0);
  
  cpgbeg(0,"?",2,2);
  if((imgccd=malloc(I.xnpix*I.ynpix*sizeof(float)))==NULL) { printf("I cannot dimension imgccd of %d bytes",I.xnpix*I.ynpix*sizeof(float));exit(1);} 
  
  /*     printf("INPU %s \n",inputcat); */

  if(strcmp(inputcat,"NONE")) {
    fromcat=1;    
    if((fileicat=fopen(inputcat,"r"))==NULL) {
      printf("ERROR: Can't open input catalogue %s\n",inputcat);
      exit(1);
    }
    printf(" Galaxy parameters from catalogue file %s\n",inputcat);
  }
  
  if (fromcat) ntot=FileNLin(inputcat)-2;
  else ntot=(FD.ngrid*FD.new*FD.nmag*FD.nz*FD.nfwhm*FD.nsize*FD.nexcen*FD.nAP);
  
  printf(" Total number of galaxies simulated %d\n",ntot);
  printf(" Using galaxy profile %d\n",FD.profile);
  
  printf("********************************************\n");
  printf("******** Beginning of computation **********\n");
  
  for(i=0;i<I.xnpix*I.ynpix;i++) imgccd[i]=0;
  
  printf(" %5d Net     %5d %5d %5d %5d %5d %5d %5d %5d \n",ntot,FD.ngrid,FD.new,FD.nmag,FD.nsize,FD.nz,FD.nfwhm,FD.nexcen,FD.nAP);
  
  
  for(number=0;number<ntot;number++) {  
    
    
    if(fromcat) Numbercat2param(number,fileicat );
    else Numbernet2param(number,&jew,&jmag,&jz,&jfwhm, &jsize,&jexcen,&jAP,&jgrid);
    
    /*        printf(" Hasta aqui  %d at X=%f Y=%f with EW=%f m_r=%f z=%f\n",number,G.x,G.y,G.ew,G.mag,G.z); */
    /*     printf(" %5d Numeros %5d %5d %5d %5d %5d %5d %5d %5d \n",number,jgrid,jew,jmag,jsize,jz,jfwhm,jexcen,jAP); */
    
    /*         fprintf(fileocat," De aqui %6d  %6.2f  %6.2f  %5.2f  %4.1f  %4.2f %6.5f  %f  %3.1f  %5.1f\n",number,G.x,G.y,G.mag,G.ew,G.fwhm,G.z,G.reff,G.excen,G.AP); */
    /* Posicion del objeto */
    
    if(!fromcat && (ntot==1)) {
      G.x=I.xnpix/2.;
      G.y=I.ynpix/2.;
    }
    
    
    
    printf(" Galaxy  %d at X=%f Y=%f with EW=%f m_r=%f z=%f\n",number+1,G.x,G.y,G.ew,G.mag,G.z);
    /*     		  printf("Mac %e %e %e\n",I.pixsize/I.focal,4*(O.seeing+G.reff)/206265,tamx*206265); */
    
    /*      He puesto el 6 que significa que cojo hasta 3 radios efectivos */
    /*     EL tamaño en radianes es el maximo entre el tamaño de un pixel  */
    /*      en el cielo y tres veces la distancia caracteristica:  */
    /*      Reff+Seeing */
    tamx=maxf(I.pixsize/I.focal,6*(G.reff+O.seeing)/206265);
    tamy=maxf(I.pixsize/I.focal,6*(G.reff+O.seeing)/206265);
    
    /*     Para la discretizacion, primero comparamos la discretizacion  */
    /*     correspondiente a un quinto del Reff y del seeing. Me quedo */
    /*     con la mas discretizada. */
    nx=maxf(tamx*206265/G.reff*5,tamx*206265/O.seeing*5);
/*     printf(" Los estos son: %f %f \n", tamx*206265/G.reff*5,tamx*206265/O.seeing*5); */
    ny=maxf(tamx*206265/G.reff*5,tamx*206265/O.seeing*5);
    
    /*      Luego comparo esto con el tamaño del pixel en el cielo y me quedo con */
    /*     el mas discretizado. */
    /*     Es posible que este ultimo paso na haga falta, simplemente puedo  */
    /*     coger una discretizacion mayor, ya me aseguro que esa discretizacion sea */
    /*     suficiente segun las escalas tipicas de la imagen (paso anterior) */
/*     printf(" Al final 111 %d %d \n",nx,ny); */
    nx=maxf(nx,6*(G.reff+O.seeing)/I.pixsize/206265*I.focal);
    ny=maxf(ny,6*(G.reff+O.seeing)/I.pixsize/206265*I.focal);
/*     printf(" Al final 222 %d %d \n",nx,ny); */
/*     printf(" Tamaño en segundos de arco de la imagen directa %f  Reff %f seeing %f\n",6*(G.reff+O.seeing),G.reff,O.seeing); */
/*     printf(" Tamaño del pixel en segundos de arco %f Numero de pixels en imagen directa %f\n",1/I.focal*I.pixsize*206265,6*(G.reff+O.seeing)/I.pixsize/206265*I.focal); */
/*     printf(" Al final %d %d \n",nx,ny); */
    
    /*     exit(1); */
    fprintf(fileocat," %6d  %8.2f  %8.2f  %8.4f  %7.2f  %7.2f %6.5f  %f  %3.1f  %5.1f\n",number,G.x,G.y,G.mag,G.ew,G.fwhm,G.z,G.reff,G.excen,G.AP);
    /*     printf("TAMAASOAOSAOS %f \n",G.reff); */
    if((img=malloc(1*sizeof(float)))== NULL) {
      printf("I can't dimension the matrix img of %d elements",nx*ny*O.nldo);
      exit(1);
    }   
    if((direct=malloc(nx*ny*sizeof(float)))== NULL) {
      printf("I can't dimension the matrix direct of %d elements",nx*ny*O.nldo);
      exit(1);
    }   
    if((espec=malloc(O.nldo*sizeof(float)))== NULL) {
      printf("I can't dimension the vector espec of %d elements",O.nldo);
      exit(1);
    }  
    /*     printf("********************\n"); */
    /*     printf("********************\n"); */
    /*     printf("Generating direct image\n"); */
    /*     printf("********************\n"); */
    /*     printf("********************\n"); */
    /*     printf("Size of direct image: \n   %e x %e pixels \n   %d x %d steps\n   %e x %e arcsecs",tamx/I.pixsize*I.focal,tamy/I.pixsize*I.focal,nx,ny,tamx*206265,tamy*206265); */
    /*     		  printf("Graphics interface in use\n"); */
    /*     addgalaxy(img,G,nx,ny); */
    switch(FD.profile) {
    case 0:
      tamx=maxf(I.pixsize/I.focal,5*(G.reff+O.seeing)/206265);
      tamy=maxf(I.pixsize/I.focal,5*(G.reff+O.seeing)/206265);
      nx=maxf(tamx*206265/G.reff*5,tamx*206265/O.seeing*5);
      ny=maxf(tamx*206265/G.reff*5,tamx*206265/O.seeing*5);
      nx=maxf(nx,5*(G.reff+O.seeing)/I.pixsize/206265*I.focal);
      ny=maxf(ny,5*(G.reff+O.seeing)/I.pixsize/206265*I.focal);
      addstepgal(direct,espec,nx,ny,tamx,tamy);
      break;
    case 1:
      tamx=maxf(I.pixsize/I.focal,6*(G.reff+O.seeing)/206265);
      tamy=maxf(I.pixsize/I.focal,6*(G.reff+O.seeing)/206265);
      nx=maxf(tamx*206265/G.reff*5,tamx*206265/O.seeing*5);
      ny=maxf(tamx*206265/G.reff*5,tamx*206265/O.seeing*5);
      nx=maxf(nx,6*(G.reff+O.seeing)/I.pixsize/206265*I.focal);
      ny=maxf(ny,6*(G.reff+O.seeing)/I.pixsize/206265*I.focal);
      addexpgal(direct,espec,nx,ny,tamx,tamy);
      break;
    case 2:
      /*       addr14gal(direct,espec,nx,ny,tamx,tamy); */
      break;
    case 3:
      tamx=maxf(I.pixsize/I.focal,6*(G.reff+O.seeing)/206265);
      tamy=maxf(I.pixsize/I.focal,6*(G.reff+O.seeing)/206265);
      nx=maxf(tamx*206265/(G.reff+O.seeing)*5,5*(G.reff+O.seeing)/I.pixsize/206265*I.focal);
      ny=maxf(tamx*206265/(G.reff+O.seeing)*5,5*(G.reff+O.seeing)/I.pixsize/206265*I.focal);
      addgaussgal(direct,espec,nx,ny,tamx,tamy);
      break;
    }
    /*     addstepgal(direct,espec,nx,ny,tamx,tamy); */
    /*     addexpgal(direct,espec,nx,ny,tamx,tamy); */
    /*     addgaussgal(direct,espec,nx,ny,tamx,tamy); */
    
/*     printf("End of generating direct image\n"); */
    
    
    E.teff=4500;
    E.mag=G.mag;
    
    /*     addstar(direct,espec,nx,ny); */
      /*************DIBUJO IMAGEN DIRECTA
***********************************
***************************************/
    /*
      cpgpanl(1,2);
      cpgswin(0,nx+1,0,ny+1);
      cpgbox("BCTNS",0,0,"BCTNS",0,0);
      cpggray(direct,nx,ny,1,nx,1,ny,4.,0,tr);
      
      cpgpanl(1,2);
      cpgswin(O.ldomin, O.ldomax,0.,2* espec[O.nldo/2]);
      cpgbox("BCTNS",0,0,"BCTNS",0,0);
      cpgmove(O.ldomin, 0);
      for(l=0;l<O.nldo;l++) {
      ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1);
      cpgdraw(ldo,espec[l]);
      }
      cpgswin(O.ldomin,O.ldomax,0.,1);
      cpgmove(4000., 0);
      for(l=0;l<O.nldo;l++) {
      ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1);
      cpgdraw(ldo,T(ldo)*QE(ldo));
      }
    */
        
    if((tamx/nx > I.pixsize /I.focal) || (tamy/ny > I.pixsize/I.focal)) { 
      printf(" WARNING: Input image not enough discretized\n Do you want to continue anyway?  "); 
      fflush(NULL);
    } 

    imin=1e38;
    imax=0;
    xs=-tamx/2;
    ys=-tamy/2;
    DISP(xs,ys,n(O.ldomin),I.alfa_prism,&x,&y);
    psi=y*I.focal;
    eta=(x-xdesp)*I.focal;
    icmin=(int)(eta/I.pixsize+G.x-.5)-1;
    jcmin=(int)(psi/I.pixsize+G.y-.5)-1;
/*     printf(" meto x %f y %f saco x %f y %f\n",xs,ys,x,y); */
    
/*     printf(" desde %d %d \n",icmin,jcmin);  */
    xs=tamx/2;
    ys=tamy/2;
    DISP(xs,ys,n(O.ldomax),I.alfa_prism,&x,&y);
/*     printf(" meto x %f y %f saco x %f y %f\n",xs,ys,x,y); */
    psi=y*I.focal;
    eta=(x-xdesp)*I.focal;
    icmax=(int)(eta/I.pixsize+G.x+.5)+1;
    jcmax=(int)(psi/I.pixsize+G.y+.5)+1;
    /*     printf(" hasta %d %d \n",icmax,jcmax);		   */
    if(icmin>icmax) {
      printf(" Error: icmin>icmax. Debug program. Exiting\n");
      exit(1);
    }
    if(jcmin>jcmax) {
      printf(" Error: jcmin>jcmax. Debug program. Exiting\n");
      exit(1);
    }
    icmin=max(0,icmin);
    icmax=min(I.xnpix,icmax);
    jcmin=max(0,jcmin);
    jcmax=min(I.ynpix,jcmax);
/*     printf(" Los limites %d %d %d %d :: %d %d\n",icmin,icmax,jcmin,jcmax,(icmax-icmin),(jcmax-jcmin)); */
    /*        printf(" Empieza lo bueno\n"); */
    
    if((imgobjccd=malloc((icmax-icmin)*(jcmax-jcmin)*sizeof(float)))==NULL) { printf("I cannot dimension imgccd of %d bytes",(icmax-icmin)*(jcmax-jcmin)*sizeof(float));exit(1);} 
    
    for(i=0;i<(icmax-icmin)*(jcmax-jcmin);i++) { 
      imgobjccd[i]=0;
    }
    printf("  Prism image size  %d pix x %d pix\n",(icmax-icmin),(jcmax-jcmin)); 
    printf("  Beginning to add object to CCD array\n");
    if(DEBUG) printf("                                              ");
    for(i=icmin;i<icmax;i++) { 
      if(DEBUG) printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05d/%05d                                    ",i,icmax);
/*       //      printf(" x %e  x-des %e  xs %e   i %d\n",x,x-xdesp,xs,i); */
      
      for(j=jcmin;j<jcmax;j++) {
 	if(DEBUG) printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05d/%05d                        ",j,jcmax); 
	for(l=0;l<O.nldo;l++) { 
	  ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
 	  if(DEBUG) printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%11f/%11f",ldo,O.ldomax); 
	  eta=(i-.5-G.x)*I.pixsize;
	  psi=(j-.5-G.y)*I.pixsize;
/* 	  //x=eta/I.focal+xdesp; */
	  x=eta/I.focal+xdesp;
	  y=psi/I.focal;     
	  /*Aqui es necesario poner -I.alfa, porque quiero sacarlo al reves.  
	    Me explico: si pongo -I.alfa saco x,y a partir de xs, ys */
	  DISP(x,y,n(ldo),-I.alfa_prism,&xs1,&ys1); /* He cambiando -I.al */
	  if(DEBUG) printf("i %d j %d eta %f psi %f ",i,j,eta,psi);
 	  if(DEBUG) printf("x %e y %e xs %e ys %e,   ldo  %f desp %e ",x,y,xs1,ys1,ldo,xdesp);  
	  
/* 	  //printf(" x %e  x-des %e",x,x-xdesp); */
	  ipix1=(int)((xs1+tamx/2)*nx/tamx);
	  jpix1=(int)((ys1+tamy/2)*ny/tamy);
	  eta=(i+.5-G.x)*I.pixsize;
	  psi=(j+.5-G.y)*I.pixsize;
	  x=eta/I.focal+xdesp;
	  y=psi/I.focal;
	  DISP(x,y,n(ldo),-I.alfa_prism,&xs2,&ys2);
	  ipix2=(int)((xs2+tamx/2)*nx/tamx);
	  jpix2=(int)((ys2+tamy/2)*ny/tamy);
	  
	  if(FD.profile==3) {
	    /* Con esto me quito todo el rollo de la imagen directa. Hago la integral directamente */
	    if(DEBUG) printf(" xs1 %f xs2 %f   rr %f\n",xs1,xs2,((G.reff)/0.832+O.seeing/2.35)/206265.);
	    contri=int2dgaussian(xs1,xs2,0.,ys1,ys2,0.,((G.reff)/0.832+O.seeing/2.35)/206265.)*espec[l];
	  }
	  else {
	    contri=0;
	    i1=min(ipix1,ipix2);
	    i1=max(0,i1);
	    i2=max(ipix1,ipix2);
	    i2=min(nx,i2);
	    j1=min(jpix1,jpix2);
	    j1=max(0,j1);
	    j2=max(jpix1,jpix2);
	    j2=min(ny,j2);
	    /*  	  printf("Los iniced %d %d %d %d ",ipix1,jpix1,ipix2,jpix2);  */
	    /*  	  printf("Los finea %d %d %d %d de %d %d\n",i1,j1,i2,j2,nx,ny);  */
	    
	    /* 	    if(i1==i2) printf(" Este caso hay que verlo con detenimiento\n"); */
	    /* 	  if(j1==j2) printf(" Este caso hay que verlo con detenimiento\n"); */
	    
	    
	    for(ii=i1;ii<i2+1;ii++) {
	      for(jj=j1;jj<j2+1;jj++) {	    
		/*if(i1==i2) {
		  //Entonces lo que hago es una interpolacion con los
		  // tres puntos mas cercanos. Es el polinomio de Lagrange
		  eta=(i-G.x)*I.pixsize;
		  psi=(j-G.y)*I.pixsize;
		  x=eta/I.focal+xdesp;
		  y=psi/I.focal;
		  DISP(x,y,n(ldo),-I.alfa_prism,&xs,&ys);
		  coef1=(((xs+tamx/2)*nx/tamx)-ii)*(((xs+tamx/2)*nx/tamx)-ii-1);
		  coef2=(((xs+tamx/2)*nx/tamx)-ii+1)*(((xs+tamx/2)*nx/tamx)-ii);
		  coef3=(((xs+tamx/2)*nx/tamx)-ii+1)*(((xs+tamx/2)*nx/tamx)-ii-1);
		  contri+=(coef1*direct[ii-1+jj*nx]+coef2*direct[ii+jj*nx]+coef3*direct[ii+1+jj*nx])*espec[l];
		  }
		*/
		contri=contri+direct[ii+jj*nx]*espec[l];
		
		/* 	      printf("img %e imgold %e\n",direct[ii+jj*nx]*espec[l],img[ii+jj*nx+l*nx*ny]); */
		/* 	      printf("Estoy integrandao %d %d %d %d %d %d\n",ii,jj,i1,i2,j1,j2); */
		/* 	      printf("en i %d j %d\n",i,j); */
		/* 	      if(ii==20 && jj==20 && l==40) printf("En el 20, 20,50: %e %e\n",contri,img[ii+nx*jj+nx*ny*l]); */
		
	      }
	    }
	    
	    /* 	  printf("el producot %d",(i2-i1+1)*(j2-j1+1)); */
	    if((i2-i1)*(j2-j1)!=0) contri=contri/(i2-i1)/(j2-j1);
	    /* 	  printf("Los adas %e %e %d %d %d %d \n",contri,img[ipix1+jpix1*nx+l*nx*ny],ipix1,jpix1,ipix2,jpix2); */
	    /* 	  	  printf(" calculo %d %d\n",i-icmin,(j-jcmin)*(jcmax-jcmin)); */
	    if(i-icmin+(j-jcmin)*(icmax-icmin)>(icmax-icmin)*(jcmax-jcmin)-1) printf(" CATSAASAS %d %d %d %d\n", i-icmin+(j-jcmin)*(icmax-icmin),i-icmin,j-jcmin,(icmax-icmin)*(jcmax-jcmin));
	  }
	  imgobjccd[i-icmin+(j-jcmin)*(icmax-icmin)]+=contri*ldo*T(ldo)*QE(ldo);
	  
	  
	}
	/* 	printf("Esto va ahora despues de todo\n"); */
	/* 	imgccd[i+j*I.xnpix]=imgccd[i+j*I.xnpix]/O.nldo*(O.ldomax-O.ldomin)/hc/1.e10*(pi*I.diameter*I.diameter/4)*I.pixsize*I.pixsize/I.focal/I.focal*206265*206265; */
	
      }
      /*       //exit(1); */
      
    }

    printf(" Finishing to add object\n");

/*     //    printf(" Hasta luego\n"); */
    for(i=0;i<(icmax-icmin);i++) {
      for(j=0;j<(jcmax-jcmin);j++) {
	imgccd[i+icmin+(j+jcmin)*I.xnpix]+=imgobjccd[i+j*(icmax-icmin)];
	if(imgccd[i+icmin+(j+jcmin)*I.xnpix]<imin) imin=imgccd[i+icmin+(j+jcmin)*I.xnpix];
	if(imgccd[i+icmin+(j+jcmin)*I.xnpix]>imax) imax=imgccd[i+icmin+(j+jcmin)*I.xnpix];
      }
    }



/*     //    printf(" Acabao \n"); */
    free(img);
    free(espec);
    free(direct);
    free(imgobjccd);
    cpgpanl(1,1);
    cpgswin(0.,(float)I.xnpix,0.,(float)I.ynpix);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
/*     //		  printf("dibuho ccd\n"); */
    /*    for(i=icmin;i<icmax;i++) { 
	  for(j=jcmin;j<jcmax;j++) {
	  imgccd[i+j*I.xnpix]=imgccd[i+j*I.xnpix]+imgp[i+j*I.xnpix];
	  }
	  }
    */  
    
    cpggray(imgccd,I.xnpix,I.ynpix,1,I.xnpix,1,I.ynpix,imax ,0,tr);
    
  }
  
  
  fclose(fileocat);
  
  
  printf("*****************************************\n");
  printf("******* End of primary computation ******\n");
  
/*   //  printf("sal\n"); */
/*   //printf("imin %e imax %e\n",imin,imax); */
  cpgpanl(1,1);
  cpgswin(0.,(float)I.xnpix,0.,(float)I.ynpix);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
/*   //  printf("dibuho ccd\n"); */
/*   //  printf("SSSSSSSSS AQUI NP \n"); */

  cpggray(imgccd,I.xnpix,I.ynpix,1,I.xnpix,1,I.ynpix,imax ,0,tr);
  cpgpanl(1,1);
  esp1d=malloc(I.xnpix*sizeof(float));
  
  for(i=0;i<I.xnpix;i++) {
    esp1d[i]=0;
    for(j=5-3-1;j<5+3;j++) {
      esp1d[i]=esp1d[i]+imgccd[i+I.xnpix*j];
    }
  }
/*   //  printf(" AQUI NP \n"); */
  cpgswin(0.,(float)I.xnpix,esp1d[1]/1.2,esp1d[I.xnpix/2]*1.2);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgmove(0.,0.);
  for(i=0;i<I.xnpix;i++) {
    cpgdraw(i,esp1d[i]);
  } 
  free(esp1d);
  printf("********************\n");
  printf("********************\n");
  printf("********************\n");
  printf("Ahora vienen las contribuciones finales\n");
  printf("********************\n");
  printf("********************\n");
  printf("********************\n");
  printf("Contribucion del cielo %e fotones/s/pixel %e ADUS\n",skyfot,skyfot*I.texp*I.eff/I.gain);
  if(ntot==1 ) printf(" Theorical number of counts from object (Efficiency = 1, Effective wavelength: %f ): %e ADU\n",(O.ldomax+O.ldomin)/2,pow(10.,-0.4*(G.mag+zeropoint))*(pi*I.diameter*I.diameter/4)*I.texp*(O.ldomax-O.ldomin)*(O.ldomax+O.ldomin)/1.e10/2/hc/I.gain);



  

  /* Aqui vienen ahora las contribuciones finales, Poisson, ADU,...,bias,..*/
  for(i=0;i<I.xnpix*I.ynpix;i++) {   
    imgccd[i]=imgccd[i]/O.nldo*(O.ldomax-O.ldomin)/hc/1.e10*(pi*I.diameter*I.diameter/4)*I.pixsize*I.pixsize/I.focal/I.focal*206265*206265;
    imgccd[i] +=skyfot;  /* Sumamos el numero de fotones por s del cielo; */
    imgccd[i] *=I.texp*I.eff; /* Esta ya son fotones por pixel; */
    imgccd[i] +=I.texp*I.dark*I.gain;
    if(imgccd[i]<1.e9)     imgccd[i]=Poidev((int)imgccd[i]);/* // Ruido poissoniano */
    else imgccd[i]=Gasdev()*sqrt(imgccd[i])+imgccd[i]; /*  //Ruido gaussiano; */
    
    imgccd[i] +=I.noise_e*Gasdev();     /* //Anadimos ruido de lectura; */
    imgccd[i] /=I.gain;                 /*  // Pasamos a ADU; */
    imgccd[i] +=I.bias;                 /*  //Anadimos BIAS; */
    if(I.satura != 0)  if(imgccd[i]>I.satura) imgccd[i]=I.satura; /* //Saturacion; */
  }
  imin=1e38;
  imax=0;
  for(i=0;i<I.xnpix*I.ynpix;i++) {   
    if(imgccd[i]<imin) imin=imgccd[i];
    if(imgccd[i]>imax) imax=imgccd[i];
  }
  printf("Cuts: imin %f imax %f\n",imin,imax);
  fnul=fnul2=0;
  esp1d=malloc(I.xnpix*sizeof(float));
  for(i=0;i<I.xnpix;i++) {
    esp1d[i]=0;
    for(j=5-3-1;j<5+3;j++) {
      esp1d[i]=esp1d[i]+imgccd[i+I.xnpix*j];
    }
  }
  cpgpanl(2,2);
  cpgeras();
  cpgswin(0.,(float)I.xnpix,esp1d[1]/1.2,esp1d[I.xnpix/2]*1.2);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgmove(0.,0.);
  for(i=0;i<I.xnpix;i++) {
    cpgdraw(i,esp1d[i]);
  }
  cpgpanl(2,1);
  cpgswin(0.,(float)I.xnpix+1.,0.,(float)I.ynpix+1.);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpggray(imgccd,I.xnpix,I.ynpix,1,I.xnpix,1,I.ynpix,imax ,imin,tr);  
  DISP(0,0.,n(O.ldomin),I.alfa_prism,&xmax,&fnul);
  DISP(0,0.,n(O.ldomax),I.alfa_prism,&xmin,&fnul);
/*   //  printf(" xmin %e xmax %e\n",xmax,xmin); */
  eta=(xmax-xdesp)*I.focal; 
  ipix=(int)(eta/I.pixsize+I.xnpix/2+.5); 
  eta=(xmin-xdesp)*I.focal; 
  jpix=(int)(eta/I.pixsize+I.xnpix/2+.5); 
/*   //  printf("pixel ldomin %d ldomax %d\n",ipix,jpix); */
  cpgmove(ipix,0);
  cpgdraw(ipix,I.ynpix);
  cpgmove(jpix,0);
  cpgdraw(jpix,I.ynpix); 
  DISP(0,0.,n(6562*(1+G.z)),I.alfa_prism,&xmax,&fnul);
/*   //  printf("desvicaino 0,0:  %e\n",xmax); */
  eta=(xmax-xdesp)*I.focal; 
  ipix=(int)(eta/I.pixsize+I.xnpix/2+.5); 
/*   //  printf("pixel Halfa %d \n",ipix); */
  cpgmove(ipix,0);
  cpgdraw(ipix,I.ynpix);
  /*  while(cnul!='D') {
    printf("Pulsa");
    printf(" %c %f %f %e\n",cnul,hc,fnul,imgccd[(int)hc+I.xnpix*(int)fnul]);
    cpgcurs(&hc,&fnul,&cnul);
    cpggray(imgccd,I.xnpix,I.ynpix,1,I.xnpix,1,I.ynpix,imax/fnul*I.ynpix/2 ,imin,tr);
  }
  */
    
  cpgswin(0.,(float)I.xnpix,0.,(float)I.ynpix);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpggray(imgccd,I.xnpix,I.ynpix,1,I.xnpix,1,I.ynpix,imax/fnul*I.ynpix/2 ,imin,tr);
  cpgend();
  /* Intentando salvar en FITS */
  ffinit(&fits,imgfile,&status);
  if(status) fits_report_error(stderr,status);
  naxes[0]=I.xnpix;
  naxes[1]=I.ynpix;
  fits_create_img(fits,-32,2,naxes,&status);
  fits_write_img(fits,TFLOAT,1,I.xnpix*I.ynpix,imgccd,&status);
  fits_close_file(fits,&status);
  if(status) fits_report_error(stderr,status);
  /*   SaveOptions(I); */
  
  return(0);
}

void addgalaxy(float *img,struct gal G,int nx, int ny)
{ 
  float fcont;
  float tamx,tamy;
  int i,j,l; 
  float ldo,x,y; 
  float halfa;
  float pi=4*atan(1.);
  float tr[6];
  float *imgdirect;
  float *imgseeing;
  long naxes[2];
  int status=0;
  fitsfile  *direct;
/*   fitsfile *seeing; */
  imgdirect=malloc(nx*ny*sizeof(float));
  imgseeing=malloc(nx*ny*sizeof(float));
  printf("Galaxia %f %f %f\n",G.z,G.mag,G.ew);
  halfa=6562*(1+G.z);
  printf("Halfa %f  z %f\n",halfa,G.z);
/*   // Galaxia escalon	   ; */
  fcont=pow(10.,-0.4*(G.mag+zeropoint))/pi/G.reff/G.reff; 
/*   // fcont: W/(m·m·Angstrom·arcsec·arcsec) */
/*   // Galaxia exponencial	    */
/*   //  fcont=pow(10.,-0.4*(G.mag+39.39))*1.e4/pi/G.reff/G.reff;//no esta cambiado; */
/*   // Galaxia gausiana ; */
/*   //fcont=pow(10.,-0.4*(G.mag+39.39))*1.e4/2/pi/G.reff/G.reff; */
  tamx=60;tamy=60;
  tamx=tamx/206265; 
  tamy=tamy/206265; 
  cpgbeg(0,"?",2,2);
  cpgsch(3);
  for(j=0;j<ny;j++) { 
    y=tamy*(j+0.5)/ny-tamy/2; 
    for(i=0;i<nx;i++) { 
      x=tamx*(i+0.5)/nx-tamx/2;
/*       // Galaxia escalon; */
      imgdirect[i+nx*j]=((x*x+y*y)<( sin(G.reff/206265.)*sin(G.reff/206265.)));
/*       // Perfil gaussiano; */
/*       //imgdirect[i+nx*j]=(exp(-(x*x+y*y)/G.reff*206265./G.reff*206265./2)); */
/*       // Perfil exponencial; */
/*       //imgdirect[i+nx*j]=(exp(-sqrt(x*x+y*y)/G.reff*206265.));     */
    }
  }
  printf("la imagen entrada %f %f\n",imgdirect[4+nx*4],imgdirect[5+nx*5]);
  cpgswin(0,nx,0,ny);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  cpggray(imgdirect,nx,ny,1,nx,1,ny,1.5,0,tr);
  printf("Seeing en pixels: %f %f\n",O.seeing/(tamx/nx*206265),O.seeing );
  
  Gausfilter2D(imgdirect,imgseeing,nx,ny,O.seeing/(tamx/nx*206265));
  printf("Es aqui\n");
  cpgpanl(2,1);
  cpgswin(0,nx,0,ny);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  printf("aqui no\n");
  cpggray(imgseeing,nx,ny,1,nx,1,ny,1.5,0,tr);
  
  printf("aqui yes\n");
  ffinit(&direct,"direct.fits",&status);
  naxes[0]=nx;
  naxes[1]=ny;
  /*  fits_create_img(direct, -32,2,naxes,&status);
  fits_report_error(stderr,status);
  fits_write_img(direct,TFLOAT,1,nx*ny,imgdirect,&status);
  fits_report_error(stderr,status);
  fits_close_file(direct,&status);
  ffinit(&seeing,"seeing.fits",&status);
  fits_create_img(seeing, -32,2,naxes,&status);
  fits_write_img(seeing,TFLOAT,1,nx*ny,imgseeing,&status);
  fits_close_file(seeing,&status);
  */
  printf("El continuo es %e\n",fcont);
  printf("\nWavelenght range: %f  --  %f with step %f",O.ldomin,O.ldomax,(O.ldomax-O.ldomin)/O.nldo); 
  
  for(l=0;l<O.nldo;l++) { 
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1);
    for(j=0;j<ny;j++) { 
      y=tamy*(j+0.5)/ny-tamy/2; 
      for(i=0;i<nx;i++) { 
	x=tamx*(i+0.5)/nx-tamx/2;
/* 	//	if(l==0) printf("dentro %f %e %f\n",1.*((x*x+y*y)<( sin(G.reff)*sin(G.reff))),(x*x,y*y), sin(G.reff)*sin(G.reff)); */
/* 	//img[i+nx*j+nx*ny*l]=(fcont*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-halfa)*(ldo-halfa)/2./G.fwhm/G.fwhm))*((x*x+y*y)<( sin(G.reff/206265.)*sin(G.reff/206265.))))*sqrt(1.-x*x-y*y); */
/* 	//img[i+nx*j+nx*ny*l]=(fcont*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-halfa)*(ldo-halfa)/2./G.fwhm/G.fwhm))*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-6790)*(ldo-6790)/2./G.fwhm/G.fwhm))*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-6310)*(ldo-6310)/2./G.fwhm/G.fwhm))*((x*x+y*y)<( sin(G.reff/206265.)*sin(G.reff/206265.))))*sqrt(1.-x*x-y*y); */
	img[i+nx*j+nx*ny*l]=(fcont*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-halfa)*(ldo-halfa)/2./G.fwhm/G.fwhm))*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-5900)*(ldo-5900)/2./G.fwhm/G.fwhm))*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-7000)*(ldo-7000)/2./G.fwhm/G.fwhm))*imgseeing[i+nx*j])*sqrt(1.-x*x-y*y);
/* 	//	img[i+nx*j+nx*ny*l]=(fcont*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-halfa)*(ldo-halfa)/2./G.fwhm/G.fwhm))*imgseeing[i+nx*j])*sqrt(1.-x*x-y*y); */
      } 
    } 
  }
  printf("Halfa en %f\n",halfa);
/*   printf("test %f %f\n",(1.<2.)*3,(1>2)*3); */
  cpgpanl(1,2);
/*   //cpgbeg(0,"?",2,2); */
  cpgswin(0,nx,0,ny);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  printf("Flujo del obejto  %e\n ",img[75+nx*70+nx*ny*40]);
  printf("Flujo del obejto  %e\n ",img[75+nx*71+nx*ny*40]);
  printf("Flujo del obejto  %e\n ",img[75+nx*72+nx*ny*40]);
  printf("Flujo del obejto  %e\n ",img[75+nx*73+nx*ny*40]); 
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  cpggray(img,nx,ny,1,nx,1,ny,1.e-21,0,tr);

  cpgpanl(1,2);
  cpgswin(O.ldomin, O.ldomax,0., fcont*G.ew/G.fwhm*2);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgmove(4000., 0);
/*   //fppo=fopen("espectro.dat","w"); */
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1);
    cpgdraw(ldo,img[(int)(nx/2)+nx*(int)(ny/2)+nx*ny*l]);
  }
  cpgswin(O.ldomin,O.ldomax,0.,1);
  cpgmove(4000., 0);
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1);
    cpgdraw(ldo,T(ldo)*QE(ldo));
  }

} 


void addstepgal(float *img,float *espec,int nx, int ny,float tamx,float tamy)
{ 

  /* De prueba para grabar la imagen directa */
/*   int status=0; */
/*   fitsfile *fits; */
/*   long naxes[2]; */
  /*         */

  int i,j,l; 
  float ldo,x,y; 
  float pi=4*atan(1.);
  float *imgdirect;
  float *imgseeing;
  float fcont;
  float halfa=6562.*(1+G.z);
  float r0;
  r0=G.reff/0.707;
/*   //r0=G.reff; */
  imgdirect=malloc(nx*ny*sizeof(float));
  imgseeing=malloc(nx*ny*sizeof(float));
/*   // Galaxia escalon	   ; */
/*   //Para normalizar hay que dividir por la integral de lo que va */
/*   // detras de imgdirect; pi*r2 */
  fcont=pow(10.,-0.4*(G.mag+zeropoint))/pi/r0/r0;
/*   printf("Continuum flux  %e W/m2 %e erg/cm2  log(erg/cm2) %f\n",pow(10.,-0.4*(G.mag+zeropoint)),pow(10.,-0.4*(G.mag+zeropoint))*1.e3,-0.4*(G.mag+zeropoint)+3); */
/*   printf(" r0 radius %f\n",r0); */
  printf("  Direct image size %d pix x %d pix\n",nx,ny);
/*   // fcont: W/(m·m·s·Angstrom·arcdeg·arcdeg) */
/*   // 39.39  Factor to convert to J/cm/cm/s... 1e4: to J/m/m/s... */

/*   //  tamx=60;tamy=60; */
/*   //  tamx=tamx/206265;  */
/*   //  tamy=tamy/206265;  */
/*   printf(" Perfil\n"); */
  for(j=0;j<ny;j++) { 
    y=tamy*(j+0.5)/ny-tamy/2; 
    for(i=0;i<nx;i++) { 
      x=tamx*(i+0.5)/nx-tamx/2;
/*       //imgdirect[i+nx*j]=((x*x+y*y)<( sin(r0/206265.)*sin(r0/206265.))); */
      imgdirect[i+nx*j]=((x*x+y*y)<(r0/206265.*r0/206265.));
    }
  }
  
/*   printf(" Seeing\n"); */
  
  Gausfilter2D(imgdirect,img,nx,ny,O.seeing/(tamx/nx*206265));
/*   //memcpy(img, imgdirect,nx*ny*sizeof(float)); */
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
/*     //espec[l]=(fcont*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-halfa)*(ldo-halfa)/2./G.fwhm/G.fwhm))); */
    espec[l]=(fcont*(1+G.ew/sqrt(2.*pi)/(G.fwhm/2.35)*exp(-(ldo-halfa)*(ldo-halfa)/2./(G.fwhm/2.35)/(G.fwhm/2.35))));
/*     //printf("voy por %d \n",l); */
  } 
/*   printf(" Finiquito\n"); */


  /* Para grabar la imagen directa   */
  /*

  ffinit(&fits,"direct.fits",&status);
  naxes[0]=nx;
  naxes[1]=ny;
  fits_create_img(fits,-32,2,naxes,&status);
  fits_write_img(fits,TFLOAT,1,nx*ny,img,&status);
  fits_close_file(fits,&status);
  if(status) fits_report_error(stderr,status);
  */
/*   //SaveOptions(I); */
  



}



void addgaussgal(float *img,float *espec,int nx, int ny,float tamx,float tamy)
{ 

  /* De prueba para grabar la imagen directa */
/*   int status=0; */
/*   fitsfile *fits; */
/*   long naxes[2]; */
  /*         */

/*   int i,j; */
  int l; 
  float ldo;
/*   int  x,y;  */
  float pi=4*atan(1.);
  float *imgdirect;
  float *imgseeing;
  float fcont;
  float halfa=6562.*(1+G.z);
  float r0;
/*   // Esto esta mal!! No seria G.reff/0.832 + O.seeing? */
  r0=(G.reff)/0.832+O.seeing/2.35;
  /* El 2.35 viene de que el seeing es el FWHM */ 
/*   Aqui sumo el seeing porque luego no lo hago pasar por la gaussiana, ya */
/*    que la convolucion de dos gaussianas es otra gaussiana */
  imgdirect=malloc(nx*ny*sizeof(float));
  imgseeing=malloc(nx*ny*sizeof(float));
  /*    Galaxia gaussiana	   ; */
  /*   Para normalizar hay que dividir por la integral de lo que va */
  /*    detras de imgdirect; pi*r2 para una exp(-r2) */
  fcont=pow(10.,-0.4*(G.mag+zeropoint)); /* El zeropoint calcula bien las cuentas */
  
/*   fcont=pow(10.,-0.4*(G.mag+38.52))*1.e4/pi/r0/r0; */ /* V de Jhonson */
/*   fcont=pow(10.,-0.4*(G.mag+38.66))*1.e4/pi/r0/r0; */ /* I de Jhonson */
/*   printf("Continuum flux  %e W/m2 %e erg/cm2  log(erg/cm2) %f\n",pow(10.,-0.4*(G.mag+zeropoint)),pow(10.,-0.4*(G.mag+zeropoint))*1.e3,-0.4*(G.mag+zeropoint)+3); */
/*   printf(" r0 radius %f\n",r0); */
  printf("  Direct image size %d x %d\n",nx,ny);
  /*    fcont: W/(m·m·s·Angstrom·arcdeg·arcdeg) */
  /*    39.39  Factor to convert to J/cm/cm/s... 1e4: to J/m/m/s... */
  
  /* Esto lo comento porque ya no lo uso */
/*   for(j=0;j<ny;j++) {  */
/*     y=tamy*(j+0.5)/ny-tamy/2;  */
/*     y=y/r0*206265; */
/*     for(i=0;i<nx;i++) {  */
/*       x=tamx*(i+0.5)/nx-tamx/2; */
/*       x=x/r0*206265; */
/*       imgdirect[i+nx*j]=exp(-(x*x+y*y))/pi/r0/r0; */
/*     } */
/*   } */
  
/*   printf(" Seeing\n"); */
  memcpy(img, imgdirect,nx*ny*sizeof(float));
/*   // Copio la imagen directa en la finall, porque ya tiene el seeing hecho */

  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    espec[l]=(fcont*(1+G.ew/sqrt(2.*pi)/(G.fwhm/2.35)*exp(-(ldo-halfa)*(ldo-halfa)/2./(G.fwhm/2.35)/(G.fwhm/2.35))));
  } 
/*   printf(" Finiquito\n"); */

  /* Para grabar la imagen directa   */
  /*
  ffinit(&fits,"direct.fits",&status);
  naxes[0]=nx;
  naxes[1]=ny;
  fits_create_img(fits,-32,2,naxes,&status);
  fits_write_img(fits,TFLOAT,1,nx*ny,img,&status);
  fits_close_file(fits,&status);
  if(status) fits_report_error(stderr,status);
  */
/*   //SaveOptions(I); */
  



}


void addexpgal(float *img,float *espec,int nx, int ny,float tamx,float tamy)
{ 
  /* De prueba para grabar la imagen directa */
/*   int status=0; */
/*   fitsfile *fits; */
/*   long naxes[2]; */
  /*         */






/*   //float tamx,tamy; */
  int i,j,l; 
  float ldo,x,y; 
  float pi=4*atan(1.);
  float *imgdirect;
  float *imgseeing;
  float fcont;
  float halfa=6562.*(1+G.z);
  float r0;
  r0=G.reff/0.31;
/*   //r0=G.reff; */
  imgdirect=malloc(nx*ny*sizeof(float));
  imgseeing=malloc(nx*ny*sizeof(float));
/*   // Galaxia exponencial  ; */
/*   //Para normalizar hay que dividir por la integral de lo que va */
/*   // detras de imgdirect; 2pi*r2 */
  fcont=pow(10.,-0.4*(G.mag+zeropoint))/2/pi/r0/r0;
  printf("Continuum flux  %e W/m2 %e erg/cm2  log(erg/cm2) %f\n",pow(10.,-0.4*(G.mag+zeropoint)),pow(10.,-0.4*(G.mag+zeropoint))*1.e3,-0.4*(G.mag+zeropoint)+3);
/*   // fcont: W/(m·m·s·Angstrom·arcsec·arcsec) */

/*   //tamx=60;tamy=60; */
/*   //tamx=tamx/206265;  */
/*   //tamy=tamy/206265;  */
  for(j=0;j<ny;j++) { 
    y=tamy*(j+0.5)/ny-tamy/2; 
    for(i=0;i<nx;i++) { 
      x=tamx*(i+0.5)/nx-tamx/2;
      imgdirect[i+nx*j]=(exp(-sqrt(x*x+y*y)/r0*206265.));    
    }
  }
  printf(" Performing seeing...");
/*   //printf("Anre \n"); */
  Gausfilter2D(imgdirect,img,nx,ny,O.seeing/(tamx/nx*206265));
/*   //memcpy(img, imgdirect,nx*ny*sizeof(float)); */
  printf("done\n");
/*   //printf("La caf\n"); */
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    
/*     //Aqui no es FWHM sino FWHM/2.35 */
    espec[l]=(fcont*(1+G.ew/sqrt(2.*pi)/(G.fwhm/2.35)*exp(-(ldo-halfa)*(ldo-halfa)/2./(G.fwhm/2.35)/(G.fwhm/2.35))));
    
  } 
/*   //printf("espexc %e\n",espec[20]); */
  /*
  ffinit(&fits,"direct.fits",&status);
  naxes[0]=nx;
  naxes[1]=ny;
  fits_create_img(fits,-32,2,naxes,&status);
  fits_write_img(fits,TFLOAT,1,nx*ny,img,&status);
  fits_close_file(fits,&status);
  if(status) fits_report_error(stderr,status);
  */

}


void addstar(float *img,float *espec,int nx, int ny)
{ 
  float tamx,tamy;
  int i,j,l; 
  float ldo,x,y; 
  float pi=4*atan(1.);
  float *imgdirect;
  float *imgseeing;
  float fcont;
  float lambdaeff=6538;
  float hck=1.434e8;
/*   //float teff=5500; */
  float dsigmad;
  float normaliza;
  tamx=60;tamy=60;
  tamx=tamx/206265; 
  tamy=tamy/206265; 
/*   //  dsigmad=2*O.seeing*O.seeing/(tamx/nx*206265)/(tamx/nx*206265); */
  dsigmad=2*O.seeing*O.seeing/206265/206265;
/*   //printf("dsginma %e see %e\n",dsigmad,O.seeing); */
  normaliza=1/(2*pi)/(dsigmad*206265*206265/2);
  imgdirect=malloc(nx*ny*sizeof(float));
  imgseeing=malloc(nx*ny*sizeof(float));
/*   // Estrella ; */
  fcont=pow(10.,-0.4*(E.mag+zeropoint))*(exp(hck/(lambdaeff*E.teff))-1);
  



  for(j=0;j<ny;j++) { 
    y=tamy*(j+0.5)/ny-tamy/2; 
    for(i=0;i<nx;i++) { 
      x=tamx*(i+0.5)/nx-tamx/2;
/*       // Perfil exponencial; */
      img[i+nx*j]=normaliza*(exp(-(x*x+y*y)/dsigmad));    
/*       //printf("normali %e img %e\n",normaliza,img[i+nx*j]); */
    }
  }
 
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    
    espec[l]=fcont*pow(lambdaeff/ldo,5)/(exp(hck/(ldo*E.teff))-1);
/*     //printf("con %e  %e   %e\n",pow(10.,-0.4*(E.mag+39.39))*1.e4,(exp(hck/(lambdaeff*E.teff))-1)*pow(lambdaeff/ldo,5)/(exp(hck/(ldo*E.teff))-1)); */
  } 
/*   //printf("espexc %e\n",espec[20]); */
  
} 







float addsky(void)
{
  float hc=1.979e-25;  /* // J·m Julio x metro */
  float pi=4*atan(1.);
  float skycont;
  int l;
/*   int nldo; */
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
  printf(" Contri %g \n",contri);


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

  printf(" Sky flux: %e  fot/s/pix\n",skycont);
  return(skycont);
}

void SaveOptions(struct instr I)
{
  FILE *fp;
/*   int i,j; */
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
  fprintf(fp,"MAGMIN  =%21f / Apparent   magnitude ( integrated)  Min       \n",FD.magmin  );
  fprintf(fp,"MAGMAX  =%21f / Apparent   magnitude ( integrated)  Max       \n",FD.magmax  );
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
  fprintf(fp,"COMMENT  Set INPUTCAT to 'NONE' if you want to use MAGMIN,..,PAMIN instead.    \n");

  
  fprintf(fp,"COMMENT  Algorithm parameters. Do not change if you are not sure what you do   \n");
  fprintf(fp,"LDOMIN  =%21f / Beginning wavelenght of the algorithm         \n",O.ldomin);
  fprintf(fp,"LDOMAX  =%21f / End wavelentgh of the algorithm               \n",O.ldomax);
  fprintf(fp,"NLDO    =%21d / Number of steps in wavelentgh calculation     \n",O.nldo);
  fprintf(fp,"NEW     =%21d / Number of steps in EW grid                    \n",FD.new);
  fprintf(fp,"NMAG    =%21d / Number of steps in apparente magnitude grid   \n",FD.nmag);
  fprintf(fp,"NZ      =%21d / Number of steps in redshift grid              \n",FD.nz);
  fprintf(fp,"NFWHM   =%21d / Number of steps in FWHM grid                  \n",FD.nfwhm);
  fprintf(fp,"NREFF   =%21d / Number of steps in effective radius grid      \n",FD.nsize);
  fprintf(fp,"NECCEN  =%21d / Number of steps in eccentrity grid            \n",FD.nexcen);
  fprintf(fp,"NPA     =%21d / Number of steps in position angle grid        \n",FD.nAP);
  fprintf(fp,"NGRID   =%21d / Number of galaxies in every dot of the grid   \n",FD.ngrid);
  fprintf(fp,"COMMENT  Output parameters                                                     \n");
  sprintf(ch21,"'%s'",imgfile);
  fprintf(fp,"FITSFILE= %-20s / Output FITS image  file                       \n",ch21); 
  sprintf(ch21,"'%s'",catfile);
  fprintf(fp,"CATFILE = %-20s / Log file with information about each galaxy   \n",ch21); 
  fprintf(fp,"END                                                           ");
  fclose(fp);



  
}
  
void LoadParam_file()
{
  int status=0;
  char comment[51];
  fitsfile *parfile;
  printf("Using parameter file <<%s>>\n",parfilename);
  if( ffopen2(&parfile,parfilename, READONLY, &status)) fits_report_error(stderr,status);
  if(status) {
    printf(" Error opening %s. Exiting\n",parfilename);
    exit(1);
  }

  ffgky(parfile,TFLOAT,"DIAMETER",&(I.diameter),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"FOCAL",&(I.focal),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"EXPTIME",&(I.texp),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"PRISMANG",&(I.alfa_prism),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"PIXSIZE",&(I.pixsize),comment,&status);
  fits_report_error(stderr,status);
  I.pixsize=I.pixsize/1e6;
  ffgky(parfile,TFLOAT,"GAIN",&(I.gain),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"NOISE_E",&(I.noise_e),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"BIAS",&(I.bias),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"DARK",&(I.dark),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"SATURA",&(I.satura),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"TOTEFF",&(I.eff),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NAXIS1",&(I.xnpix),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NAXIS2",&(I.ynpix),comment,&status);
  fits_report_error(stderr,status);

  ffgky(parfile,TSTRING,"FILTERFI",I.filter_file,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"QEFFFILE",I.qe_file,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"SKYFILE",O.sky_file,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"NFILE",I.prism_file,comment,&status);
  fits_report_error(stderr,status);

  ffgky(parfile,TFLOAT,"EWMIN",&(FD.ewmin),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"MAGMIN",&(FD.magmin),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"FWHMMIN",&(FD.fwhmmin),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"REFFMIN",&(FD.sizemin),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"ZMIN",&(FD.zmin),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"ECCMIN",&(FD.excenmin),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"PAMIN",&(FD.APmin),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"EWMAX",&(FD.ewmax),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"MAGMAX",&(FD.magmax),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"FWHMMAX",&(FD.fwhmmax),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"REFFMAX",&(FD.sizemax),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"ZMAX",&(FD.zmax),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"ECCMAX",&(FD.excenmax),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"PAMAX",&(FD.APmax),comment,&status);
  fits_report_error(stderr,status);

  ffgky(parfile,TFLOAT,"SPADIST",&(FD.spadist),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"PROFILE",&(FD.profile),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"PHOTBAND",O.photband,comment,&status);
  fits_report_error(stderr,status);

  ffgky(parfile,TFLOAT,"SEEING",&(O.seeing),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"LDOMAX",&(O.ldomax),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"LDOMIN",&(O.ldomin),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NLDO",&(O.nldo),comment,&status);
  fits_report_error(stderr,status);

  ffgky(parfile,TINT,"NEW",&(FD.new),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NMAG",&(FD.nmag),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NFWHM",&(FD.nfwhm),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NREFF",&(FD.nsize),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NZ",&(FD.nz),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NECCEN",&(FD.nexcen),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NPA",&(FD.nAP),comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"NGRID",&(FD.ngrid),comment,&status);
  fits_report_error(stderr,status);

  ffgky(parfile,TSTRING,"FITSFILE",imgfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"CATFILE",catfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"INPUTCAT",inputcat,comment,&status);
  fits_report_error(stderr,status);

  I.alfa_prism=I.alfa_prism/206265*3600;
  
  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file. Perhaps not enough keywords.\n");
    exit(1); 
  }
  fits_close_file(parfile,&status);
 
}

void LoadParam_kbd(void)
{
  printf("*****************\nParameters Input\n****************\n"); 
  printf("Telescope diameter (in meters): ");
  I.diameter=readf(1.);
  printf("Telescope focal (in meters): ");
  I.focal=readf(1.); 
  printf("Prism angle (in degrees): ");
  I.alfa_prism=readf(0.);
  I.alfa_prism=I.alfa_prism/206265*3600;
  printf("Exposure time (in seconds): ");
  I.texp=readf(1.);; 
  printf("pixel size (in microns): ");
  I.pixsize=readf(24.); 
  I.pixsize=I.pixsize/1e6;


  printf("CCD gain:");
  I.gain=readf(1.); 
  printf("CCD read-out noise (electrons / second): ");
  I.noise_e=readf(1.); 
  printf("CCD bias: ");
  I.bias=readf(0.); 
  printf("CCD dark: ");
  I.dark=readf(0.); 
  printf("CCD pixels in X direction: ");
  I.xnpix=readi(1000); 
  printf("CCD pixels in Y direction: ");
  I.ynpix=readi(1000); 
  printf("CCD saturation level (counts): ");
  I.satura=readi(65000); 

  printf("File with quantum efficiency: ");
  reads(I.qe_file,I.qe_file);
  

  printf("File with filter response: ");
  reads(I.filter_file,I.filter_file);
  printf("File with prism refraction index: ");
  reads(I.prism_file,I.prism_file);
  printf("Total efficiency of the system: ");
  I.eff=readf(1.); 
  printf("Seeing of the observations: ");
  O.seeing=readf(1.); 
  printf("File with sky brightness in W/(m·m·s·Angstrom·arc·arc): ");
  reads(O.sky_file,O.sky_file);


  printf("Minimum Apparent   magnitude of the galaxy integrated: ");
  FD.magmin=readf(10.); 
  printf("Maximum Apparent   magnitude of the galaxy integrated: ");
  FD.magmax=readf(20.); 
  printf("Minimum Equivalent width of the Halpha line ( Angstroms): ");
  FD.ewmin=readf(0.); 
  printf("Maxmum Equivalent width of the Halpha line ( Angstroms): ");
  FD.ewmin=readf(100.); 
  printf("Efective radius of the galaxy (arcseconds): ");
  G.reff=readf(O.seeing); 
  printf("Redshift of the galaxy: ");
  G.z=readf(0.); 


  printf("Beginning wavelenght of the algorithm: ");
  O.ldomin=readf(1000.); 
  printf("End  wavelenght of the algorithm: ");
  O.ldomax=readf(10000.); 
  printf("Number of steps in wavelentgh calculation: ");
  O.nldo=readf(1000.); 

  printf("Output image file: ");
  reads(imgfile,imgfile);
  printf("Output catalogue file: ");
  reads(catfile,catfile);
  
  
  /*
    printf("EW (Å): "); 
    scanf("%f",&(G.ew)); 
    printf("%f\n",G.ew);
    printf("G.Fwhmma (Å): "); 
    scanf("%f",&(G.fwhm)); 
    printf("Continuum flux (W/arcsec2/m2/Å): "); 
    scanf("%f",&(G.mag)); 
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
  *cosxs=cosx+alfa*(cocin-1);
  *cosys=cosy;
/*   printf("cosxs %f diff %f\n",*cosxs,alfa*(cocin-1)); */
/*   printf("ALFA %f\n",alfa); */
  return(0);
}




void Numbernet2param(int number,int *jew, int *jmag, int *jz, int *jfwhm, int *jsize,int  *jexcen,int *jAP,int *jgrid) {

  
  div_t divi;
  int nrest;
  int ntot;
  int nx,ny,i,j;
  ntot=(FD.ngrid*FD.new*FD.nmag*FD.nz*FD.nfwhm*FD.nsize*FD.nexcen*FD.nAP);
  
  divi=div(number,FD.ngrid*FD.new*FD.nmag*FD.nz*FD.nfwhm*FD.nsize*FD.nexcen);
  *jAP=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new*FD.nmag*FD.nz*FD.nfwhm*FD.nsize);
  *jexcen=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new*FD.nmag*FD.nz*FD.nsize);
  *jfwhm=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new*FD.nmag*FD.nsize);
  *jz=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new*FD.nmag);
  *jsize=divi.quot;
  nrest=divi.rem;
  divi=div(nrest,FD.ngrid*FD.new);
  *jmag=divi.quot;
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
  if(FD.nmag!=1) G.mag=FD.magmin+(FD.magmax-FD.magmin)**jmag/(FD.nmag-1);
  else G.mag=FD.magmin;
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

    printf(" Position %d %d    %d  %d %f %f \n",nx,ny,i,j,G.x,G.y);
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
/*   float fnul; */
  while(cnul=='#') {
    fgets(snul,2000,fileicat);
/*     //    printf("Leyendo dentro:  %s\n",snul); */
    cnul=snul[0];
/*     //printf(" cara %c\n",cnul); */
  }
/*   //printf("Leyendo:  %s\n",snul); */
  sscanf(snul," %d %f %f %f %f %f %f %f %f %f ",&inul,&G.x,&G.y,&G.mag,&G.ew,&G.fwhm,&G.z,&G.reff,&G.excen,&G.AP);

}
void Band2ZP(void) {

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
