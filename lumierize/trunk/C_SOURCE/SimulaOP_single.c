#include "modulos.h"
/* //#include "cpgplot.h" */
/* Esta version 8 lo que hace es calcular los pixelss del CCD que se ven contribuidos por el obejto. Solo integra en esa zona 
 */


struct instr {
  float    diameter;
  float    focal;
  char     filter_file[50];
  char     qe_file[50];
  char     prism_file[50];
  float    alfa_prism;
  float    texp;
  float    pixsize;
  float    gain;
  float    noise_e;
  float    bias;
  float    dark;
  float    eff;         /*  //Total efficiency of the system. Just a factor */
  int      xnpix;
  int      ynpix;
  int      satura;      /*   //Saturation level of CCD (o if no saturation) */
};

struct object {
  char esp_file[50];
  float *ldo;
  float *flux;
  int nldo;
  char spat_file[50];
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
struct object S;
struct star E;
struct other O;
struct prism P;


char imgfile[51];
/* //char catfile[51]; */
/* //char inputcat[51]; */
void addstepgal(float *direct,int nx, int ny,float tamangx,float tamangy);
void addexpgal(float *direct,int nx, int ny);
/* //void addgalaxy(float *img, struct gal G,int nx, int ny); */
void addgaussian(float *direct,int nx, int ny);
float addsky();
void SaveOptions();
void ReadOptions(FILE *filepar);
void InputOptions();
float QE(float ldo);
float T(float ldo);
float Skyesp(float ldo);
float n(float ldo);
float sed(float ldo);
void Readfilter( char file[100]);
void Readqe(char file[100]);
void Readsky(char file[100]);
void Readprism(char file[100]);
void Readspec(char file[100]);
int min(int x1,int x2);
int max(int x1,int x2);
float maxf(float x1,float x2);
int Disper2( float cosx, float cosy, float cocin, float alfa, float *cosxs, float *cosys );


int main(int argc, char **argv)
{
  float fnul,fnul2; 
/*   float alfa; */
  float xmin,xmax,xdesp,xs,ys,x,y; 

  float xpos,ypos;
/*   int *contador;  */
  char cnul;
/*   char snul[2000]; */
  float *imgccd; 
  float ldo,tamangx,tamangy;
  float psi,eta;
  int ipix,jpix;
  int ipix1,ipix2,jpix1,jpix2;
/*   int ic1,ic2,jc1,jc2; */
  int icmin,icmax,jcmin,jcmax;
  float hc=1.979e-25;  /* // J·m Julio x metro */
/*   //  float *img; */
  float *direct;
  float *espec;
/*   float *imgp; */
  int nx,ny;
/*   FILE *fimg; */
/*   //  char fileimg[50]="salida.fits"; */
  int i,j,l;
  int ii,jj;
  int i1,j1,i2,j2;
  float contri;
  float tr[6];
/*   float b=3.5e6; */
  float imin,imax;
  float pi=4*atan(1.);
  float skyfot;
  FILE *filepar;
/*   float *xesp; */
/*   float *yesp; */
  float *esp1d;
  /*   int number=0; */  /* Contador del numero de galaxia por el que vamos */ 
/*   int ntot;  */      /* // Numero total de galaxias simuladas. */
/*   int jew,jmr,jz,jfwhm,jsize,jexcen,jAP,jgrid; */

/*   int fromcat=0; */
/*   // Para saber si los parametros de las galaxias vendran del catalogo. */
/*   int pgid; */
  int status=0;
  fitsfile *fits;
  long naxes[2];
  /********************************************
  **********************************************
  LECTURA DE PARAMETROS
  *********************
  **********************/
  srandom((unsigned int)time(NULL)/2); 
  printf("\n Welcome to SimulaOP_field Version 4\n");
  if(argc < 2) InputOptions();
  else {
    printf("This program has been called with parameters file %s\n",argv[1]);
    if((filepar=fopen(argv[1],"r"))==NULL) {
      printf("ERROR: Can't open options file %s\n",argv[1]);
      exit(1);
    }
    ReadOptions(filepar);
  }


  Readfilter(I.filter_file);
  Readqe(I.qe_file);
  Readsky(O.sky_file);
  Readprism(I.prism_file);
  Readspec(S.esp_file);

  I.alfa_prism=-I.alfa_prism; /* // Tiene su razon de ser ; */
  skyfot=addsky();

  /**********************************************
    BUCLE
    *********************************************
    *********************************************
    ***********************************************/
  Disper(I.xnpix*I.pixsize/I.focal/2,0.,n(O.ldomax),I.alfa_prism,&xmax,&fnul);
/*   //printf("max %e  min %e\n",xmax,fnul); */
  Disper(-I.xnpix*I.pixsize/I.focal/2,0.,n(O.ldomin),I.alfa_prism,&xmin,&fnul);
/*   //printf("max %e  min %e\n",xmin,fnul); */
  Disper(I.xnpix*I.pixsize/I.focal/2,0.,n(O.ldomin),I.alfa_prism,&xmax,&fnul);  
/*   //printf("max %e  min %e\n",xmax,fnul); */
  Disper(-I.xnpix*I.pixsize/I.focal/2,0.,n(O.ldomax),I.alfa_prism,&xmin,&fnul);
/*   //printf("max %e  min %e\n",xmin,fnul); */
/*   //printf("entrad %e %e\n",I.xnpix*I.pixsize/I.focal/2,-I.xnpix*I.pixsize/I.focal/2); */
  xdesp=(xmax+xmin)/2.; 
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;


  
  printf("Field covered by the output image %f x %f arcsec\n",I.pixsize/I.focal*I.xnpix*206265,I.pixsize/I.focal*I.ynpix*206265); 
  printf("\nPlate scale in the output image %f arcsec/pix\n",I.pixsize/I.focal*206265); 

  setvbuf(stdin,"",_IOLBF,0);


  if((imgccd=malloc(I.xnpix*I.ynpix*sizeof(float)))==NULL) { printf("I cannot dimension imgccd of %d bytes",I.xnpix*I.ynpix*sizeof(float));exit(1);} 


  printf("ESTO ESW DIFERGE T\n");
  cpgopen("?");
  cpgswin(O.ldomin,O.ldomax,0.,1e-12);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgslw(2.);
  cpgmove(O.ldomin,0.);
  for(l=0;l<O.nldo;l++) { 
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    cpgdraw(ldo,sed(ldo));
    printf(" %f %e\n",ldo,sed(ldo));
  }
/*   //  exit(1); */
  cpgclos();
  cpgbeg(0,"?",2,2);

  
  printf("********************************************\n");
  printf("******** Beginning of computation **********\n");
  
  for(i=0;i<I.xnpix*I.ynpix;i++) imgccd[i]=0;
  
  
  
  /* Posicion del objeto */
  
  
  xpos=I.xnpix/2.;
  ypos=I.ynpix/2.;
  
  printf(" Galaxy at X=%f Y=%f \n",xpos,ypos);
  
  tamangx=maxf(I.pixsize/I.focal,3*(O.seeing)/206265);
  tamangy=maxf(I.pixsize/I.focal,3*(O.seeing)/206265);
  printf(" pixel %f seein %f\n",I.pixsize,O.seeing);
/*   // Estos tamaños van en radianes y son tamaños angulares!! */
  nx=(int)(tamangx*I.focal/I.pixsize*2);
  ny=(int)(tamangy*I.focal/I.pixsize*2);
/*   //  printf(" tamangx %f %f tamangy %f %f \n",I.pixsize/I.focal,3*(O.seeing)/206265,I.pixsize/I.focal,3*(O.seeing)/206265); */
/*   //  printf("Las dimesnciones %d %d\n",nx,ny); */
/*   //  printf("ER numero es %d\n",O.nldo); */
/*   //if((img=malloc(nx*ny*O.nldo*sizeof(float)))== NULL) { */
/*   //  printf("I can't dimension the matrix img of %d elements",nx*ny*O.nldo); */
/*   //  exit(1); */
/*   //}    */
  if((direct=malloc(nx*ny*sizeof(float)))== NULL) {
    printf("I can't dimension the matrix direct of %d elements",nx*ny);
    exit(1);
  }   
  if((espec=malloc(O.nldo*sizeof(float)))== NULL) {
    printf("I can't dimension the vector espec of %d elements",O.nldo);
    exit(1);
  }  
/*   //addgalaxy(img,G,nx,ny); */
/*   //addstepgal(direct,nx,ny,tamangx,tamangy); */
/*   //addexpgal(direct,nx,ny); */
  addgaussian(direct,nx,ny);  
  printf("End of generating direct image\n");
  
  
/*   //E.teff=4500; */
/*   //E.mr=G.mr; */
    
  if((tamangx/nx > I.pixsize /I.focal) || (tamangy/ny > I.pixsize/I.focal)) { 
    printf(" WARNING: Input image not enough discretized\n Do you want to continue anyway?  "); 
    fflush(NULL);
    scanf("%c",&cnul);
    if(cnul=='s' || cnul=='S');
    else exit(1); 
  }
  
  imin=1e38;
  imax=0;
  xs=-tamangx/2;
  ys=-tamangy/2;
  Disper(xs,ys,n(O.ldomin),I.alfa_prism,&x,&y);
  psi=y*I.focal;
  eta=(x-xdesp)*I.focal;
  icmin=(int)(eta/I.pixsize+xpos-.5)-1;
  jcmin=(int)(psi/I.pixsize+ypos-.5)-1;
/*   //printf(" desde %d %d \n",icmin,jcmin); */
  xs=tamangx/2;
  ys=tamangy/2;
  Disper(xs,ys,n(O.ldomax),I.alfa_prism,&x,&y);
  psi=y*I.focal;
  eta=(x-xdesp)*I.focal;
  icmax=(int)(eta/I.pixsize+xpos+.5)+1;
  jcmax=(int)(psi/I.pixsize+ypos+.5)+1;
/*   //printf(" hasta %d %d \n",icmax,jcmax);		   */
  icmin=max(0,icmin);
  icmax=min(I.xnpix,icmax);
  jcmin=max(0,jcmin);
  jcmax=min(I.ynpix,jcmax);
  
/*   //    printf(" Empieza lo bueno\n"); */
  
  printf(" This object extends in the range [%d:%d,%d:%d]\n",icmin,icmax,jcmin,jcmax);
  
  for(i=icmin;i<icmax;i++) { 
    printf(" x %e  x-des %e  xs %e   i %d\n",x,x-xdesp,xs,i);
    
    for(j=jcmin;j<jcmax;j++) {
      for(l=0;l<O.nldo;l++) { 
	ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
	eta=(i-.5-xpos)*I.pixsize;
	psi=(j-.5-ypos)*I.pixsize;
/* 	//x=eta/I.focal+xdesp; */
	x=eta/I.focal+xdesp;
	y=psi/I.focal;     
	Disper(x,y,n(ldo),-I.alfa_prism,&xs,&ys); /* //He cambiando -I.al */
/* 	//printf("x %e y %e xs %e ys %e,   ldo  %f desp %e\n",x,y,xs,ys,ldo,xdesp); */
/* 	//printf(" x %e  x-des %e",x,x-xdesp); */
	ipix1=(int)((xs+tamangx/2)*nx/tamangx);
	jpix1=(int)((ys+tamangy/2)*ny/tamangy);
	eta=(i+.5-xpos)*I.pixsize;
	psi=(j+.5-ypos)*I.pixsize;
	x=eta/I.focal+xdesp;
	y=psi/I.focal;
	Disper(x,y,n(ldo),-I.alfa_prism,&xs,&ys);
	ipix2=(int)((xs+tamangx/2)*nx/tamangx);
	jpix2=(int)((ys+tamangy/2)*ny/tamangy);
	
	contri=0;
	i1=min(ipix1,ipix2);
	i1=max(0,i1);
	i2=max(ipix1,ipix2);
	i2=min(nx,i2);
	j1=min(jpix1,jpix2);
	j1=max(0,j1);
	j2=max(jpix1,jpix2);
	j2=min(ny,j2);
/* 	//printf("Los iniced %d %d %d %d\n",ipix1,jpix1,ipix2,jpix2); */
/* 	//printf("Los iniced %d %d %d %d\n",i1,j1,i2,j2); */
	
	for(ii=i1;ii<i2+1;ii++) {
	  for(jj=j1;jj<j2+1;jj++) {	    
/* 	    //contri=contri+img[ii+jj*nx+l*nx*ny]; */
	    contri=contri+direct[ii+jj*nx];
/* 	    //printf("img %e imgold %e\n",direct[ii+jj*nx]*espec[l],img[ii+jj*nx+l*nx*ny]); */
/* 	    //printf("Estoy integrandao %d %d %d %d %d %d\n",ii,jj,i1,i2,j1,j2); */
/* 	    //printf("en i %d j %d\n",i,j); */
/* 	    //if(ii==20 && jj==20 && l==40) printf("En el 20, 20,50: %e %e\n",contri,img[ii+nx*jj+nx*ny*l]); */
	    
	  }
	}
	
/* 	//printf("el producot %d",(i2-i1+1)*(j2-j1+1)); */
	if((i2-i1+1)*(j2-j1+1)!=0) contri=contri/(i2-i1+1)/(j2-j1+1);
/* 	//printf("Los adas %e %e %d %d %d %d \n",contri,img[ipix1+jpix1*nx+l*nx*ny],ipix1,jpix1,ipix2,jpix2); */
	imgccd[i+j*I.xnpix]=imgccd[i+j*I.xnpix]+contri*sed(ldo)*ldo*T(ldo)*QE(ldo);
	
	
      }
      
      if(imgccd[i+j*I.xnpix]<imin) imin=imgccd[i+j*I.xnpix];
      if(imgccd[i+j*I.xnpix]>imax) imax=imgccd[i+j*I.xnpix];
      
    }

    
  }
  printf(" Acabao \n");
/*   //    free(img); */
  free(direct);
  printf(" Acabao 2 \n");
  free(espec);
  printf(" Acabao 3 \n");
  cpgpanl(1,1);
  cpgswin(0.,(float)I.xnpix,0.,(float)I.ynpix);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
/*   //		  printf("dibuho ccd\n"); */
  /* for(i=icmin;i<icmax;i++) { 
     for(j=jcmin;j<jcmax;j++) {
     imgccd[i+j*I.xnpix]=imgccd[i+j*I.xnpix]+imgp[i+j*I.xnpix];
     }
     }
  */
  
  cpggray(imgccd,I.xnpix,I.ynpix,1,I.xnpix,1,I.ynpix,imax ,0,tr);
    
  
  printf("*****************************************\n");
  printf("******* End of primary computation ******\n");
  
/*   //  printf("sal\n"); */
/*   //printf("imin %e imax %e\n",imin,imax); */
  cpgpanl(1,1);
  cpgswin(0.,(float)I.xnpix,0.,(float)I.ynpix);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
/*   //  printf("dibuho ccd\n"); */
  
  cpggray(imgccd,I.xnpix,I.ynpix,1,I.xnpix,1,I.ynpix,imax ,0,tr);
  cpgpanl(1,1);
/*   //  esp1d=malloc(I.xnpix*sizeof(float)); */
  printf("Pasa aqui papa\n");

/*   //for(i=0;i<I.xnpix;i++) { */
/*   //  esp1d[i]=0; */
/*   // for(j=50-3-1;j<50+3;j++) { */
/*   //    esp1d[i]=esp1d[i]+imgccd[i+I.xnpix*j]; */
/*   //  } */
/*   //} */
  printf("Pasa aqui\n");
/*   //  cpgswin(0.,(float)I.xnpix,esp1d[1]/1.2,esp1d[I.xnpix/2]*1.2); */
/*   //cpgbox("BCTNS",0,0,"BCTNS",0,0); */
/*   //cpgmove(0.,0.); */
/*   //for(i=0;i<I.xnpix;i++) { */
/*   //  cpgdraw(i,esp1d[i]); */
/*   //}  */
/*   //free(esp1d); */
  printf("********************\n");
  printf("********************\n");
  printf("********************\n");
  printf("Ahora vienen las contribuciones finales\n");
  printf("********************\n");
  printf("********************\n");
  printf("********************\n");
  printf("Contribucion del cielo %e fotones/s/pixel\n",skyfot);
  /* Aqui vienen ahora las contribuciones finales, Poisson, ADU,...,bias,..*/
  for(i=0;i<I.xnpix*I.ynpix;i++) {   
    imgccd[i]=imgccd[i]/O.nldo*(O.ldomax-O.ldomin)/hc/1.e10*(pi*I.diameter*I.diameter/4)*I.pixsize*I.pixsize/I.focal/I.focal*206265*206265;
    imgccd[i] +=skyfot;  /* //Sumamos el numero de fotones por s del cielo; */
    imgccd[i] *=I.texp*I.eff; /* //Esta ya son fotones por pixel; */
    if(imgccd[i]<1.e9)     imgccd[i]=Poidev((int)imgccd[i]);/* // Ruido poissoniano */
    else imgccd[i]=Gasdev()*sqrt(imgccd[i])+imgccd[i];  /* //Ruido gaussiano; */
    
    imgccd[i] +=I.noise_e*Gasdev();     /* //Anadimos ruido de lectura; */
    imgccd[i] /=I.gain;                  /* // Pasamos a ADU; */
    imgccd[i] +=I.bias;                  /* //Anadimos BIAS; */
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
    for(j=(int)(I.ynpix/2-3);j<(int)(I.ynpix/2+3);j++) {
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
  Disper(0,0.,n(O.ldomin),I.alfa_prism,&xmax,&fnul);
  Disper(0,0.,n(O.ldomax),I.alfa_prism,&xmin,&fnul);
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
    
  cpgswin(0.,(float)I.xnpix,0.,(float)I.ynpix);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpggray(imgccd,I.xnpix,I.ynpix,1,I.xnpix,1,I.ynpix,imax/fnul*I.ynpix/2 ,imin,tr);
  cpgend();
  /* Intentando salvar en FITS */
  printf(" Saving output fits image %s\n",imgfile);
  ffinit(&fits,imgfile,&status);
  naxes[0]=I.xnpix;
  naxes[1]=I.ynpix;
  fits_create_img(fits,-32,2,naxes,&status);
  fits_write_img(fits,TFLOAT,1,I.xnpix*I.ynpix,imgccd,&status);
  fits_close_file(fits,&status);
  if(status) fits_report_error(stderr,status);
  SaveOptions(I);
  return(0);
}


void addstepgal(float *img,int nx, int ny,float tamangx,float tamangy)
{ 
  int i,j; 
  float x,y; 
/*   float pi=4*atan(1.); */
  float *imgdirect;
  float *imgseeing;
/*   float fcont; */
  float r0;
/*   // ESto se queda asi de momento */
  r0=O.seeing/0.707;
/*   //r0=G.reff; */
  imgdirect=malloc(nx*ny*sizeof(float));
  imgseeing=malloc(nx*ny*sizeof(float));
/*   // Galaxia escalon	   ; */
/*   //  fcont=pow(10.,-0.4*(G.mr+39.39))*1.e4/pi/r0/r0; */
  
/*   //  tamangx=60;tamangy=60; */
/*   //  tamangx=tamangx/206265;  */
/*   //  tamangy=tamangy/206265;  */
  printf(" Perfil\n");
  for(j=0;j<ny;j++) { 
    y=tamangy*(j+0.5)/ny-tamangy/2; 
    for(i=0;i<nx;i++) { 
      x=tamangx*(i+0.5)/nx-tamangx/2;
/*       // Galaxia escalon; */
      imgdirect[i+nx*j]=((x*x+y*y)<( sin(r0/206265.)*sin(r0/206265.)));
/*       //printf("voy por %d %d\n",j,i); */
/*       //printf("sdsss %f   %f\n",((x*x+y*y)<( sin(r0/206265.)*sin(r0/206265.))),imgdirect[i+nx*j]); */
    }
  }
  
  printf(" Seeing\n");
  
  Gausfilter2D(imgdirect,img,nx,ny,O.seeing/(tamangx/nx*206265));
/*   //  for(l=0;l<O.nldo;l++) { */
/*   // ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1);  */
/*   //  espec[l]=(fcont*(1+G.ew/sqrt(2.*pi)/G.fwhm*exp(-(ldo-halfa)*(ldo-halfa)/2./G.fwhm/G.fwhm))); */
/*     //printf("voy por %d \n",l); */
/*   //}  */


}

void addexpgal(float *img,int nx, int ny)
{ 
  float tamangx,tamangy;
  int i,j; 
  float x,y; 
/*   float pi=4*atan(1.); */
  float *imgdirect;
  float *imgseeing;
/*   //  float fcont; */
/*   //float halfa=6562.*(1+G.z); */
  float r0;
  r0=O.seeing/0.31;
/*   //r0=G.reff; */
  imgdirect=malloc(nx*ny*sizeof(float));
  imgseeing=malloc(nx*ny*sizeof(float));
/*   // Galaxia exponencial  ; */
/*   //  fcont=pow(10.,-0.4*(G.mr+39.39))*1.e4/pi/r0/r0; */
  
  tamangx=tamangx/206265; 
  tamangy=tamangy/206265; 
  for(j=0;j<ny;j++) { 
    y=tamangy*(j+0.5)/ny-tamangy/2; 
    for(i=0;i<nx;i++) { 
      x=tamangx*(i+0.5)/nx-tamangx/2;
/*       // Galaxia escalon; */

      imgdirect[i+nx*j]=(exp(-sqrt(x*x+y*y)/r0*206265.));    
/*       //printf("Qe cas\n"); */

    }
  }
  
/*   //printf("Anre \n"); */
  Gausfilter2D(imgdirect,img,nx,ny,O.seeing/(tamangx/nx*206265));
/*   //printf("La caf\n"); */
/*   //  for(l=0;l<O.nldo;l++) { */
/*   //ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1);  */

/*     //Aqui no es FWHM sino FWHM/2.35 */
/*   //espec[l]=(fcont*(1+G.ew/sqrt(2.*pi)/(G.fwhm/2.35)*exp(-(ldo-halfa)*(ldo-halfa)/2./(G.fwhm/2.35)/(G.fwhm/2.35)))); */
    
/*   //}  */
/*   //printf("espexc %e\n",espec[20]); */

}


void addgaussian(float *img,int nx, int ny)
{ 
  float tamangx,tamangy;
  int i,j; 
  float x,y; 
  float pi=4*atan(1.);
  float *imgdirect;
  float *imgseeing;
  float fcont;
  float lambdaeff=6538;
  float hck=1.434e8;
/*   //float teff=5500; */
  float dsigmad;
  float normaliza;
  tamangx=60;tamangy=60;
  tamangx=tamangx/206265; 
  tamangy=tamangy/206265; 
/*   //  dsigmad=2*O.seeing*O.seeing/(tamangx/nx*206265)/(tamangx/nx*206265); */
  dsigmad=2*O.seeing*O.seeing/206265/206265;
/*   //printf("dsginma %e see %e\n",dsigmad,O.seeing); */
  normaliza=1/(2*pi)/(dsigmad*206265*206265/2);
  imgdirect=malloc(nx*ny*sizeof(float));
  imgseeing=malloc(nx*ny*sizeof(float));
/*   // Estrella ; */
  fcont=pow(10.,-0.4*(E.mr+39.39))*1.e4*(exp(hck/(lambdaeff*E.teff))-1);
  



  for(j=0;j<ny;j++) { 
    y=tamangy*(j+0.5)/ny-tamangy/2; 
    for(i=0;i<nx;i++) { 
      x=tamangx*(i+0.5)/nx-tamangx/2;
/*       // Perfil exponencial; */
      img[i+nx*j]=normaliza*(exp(-(x*x+y*y)/dsigmad));    
/*       //printf("normali %e img %e\n",normaliza,img[i+nx*j]); */
    }
  }
 
  
} 




float addsky()
{
  float hc=1.979e-25;  /* // J·m Julio x metro */
  float pi=4*atan(1.);
  float skycont;
  int l;
/*   int nldo; */
  double ldo;
  double contri;
  
  contri=0;
  for(l=0;l<O.nldo;l++) {
    ldo=O.ldomin+l*(O.ldomax-O.ldomin)/(O.nldo-1); 
    contri=contri+Skyesp(ldo)*ldo*T(ldo)*QE(ldo);
  }
  skycont=contri*(O.ldomax-O.ldomin)/O.nldo*1e-10;  /* //1e-10 From angstroms to m */
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
  printf("Name of output parameter file: ");
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
  fprintf(fp,"PIXSIZE =%21f / Pixel size of the CCD (in microns)            \n",I.pixsize*1e6);
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
  fprintf(fp,"COMMENT  Object parameters                                                     \n");
  sprintf(ch21,"'%s'",S.esp_file);
  fprintf(fp,"ESPFILE = %-20s / File with spectral energy distribution        \n",ch21); 
  fprintf(fp,"COMMENT  Algorithm parameters. Do not change if you are not sure what you do   \n");
  fprintf(fp,"LDOMIN  =%21f / Beginning wavelenght of the algorithm         \n",O.ldomin);
  fprintf(fp,"LDOMAX  =%21f / End wavelentgh of the algorithm               \n",O.ldomax);
  fprintf(fp,"NLDO    =%21d / Number of steps of the wavelength computation \n",O.nldo);
  fprintf(fp,"COMMENT  Output parameters                                                     \n");
  sprintf(ch21,"'%s'",imgfile);
  fprintf(fp,"FITSFILE= %-20s / Output FITS image  file                       \n",ch21); 
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
  I.pixsize=I.pixsize/1.e6;
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
  c+=f_kfls(filepar,"ESPFILE",S.esp_file);



/*   //printf("El z %f\n",G.z); */
  c+=f_kff(filepar,"SEEING",&(O.seeing));
  c+=f_kfi(filepar,"NAXIS1",&(I.xnpix));
  c+=f_kfi(filepar,"NAXIS2",&(I.ynpix));
  c+=f_kfi(filepar,"SATURA",&(I.satura));
  c+=f_kff(filepar,"LDOMAX",&(O.ldomax));
  c+=f_kff(filepar,"LDOMIN",&(O.ldomin));
  c+=f_kfi(filepar,"NLDO",&(O.nldo));



/*   //FD.ngrid=1; */

  c+=f_kfls(filepar,"FITSFILE",imgfile);
/*   //c+=f_kfls(filepar,"CATFILE",catfile); */
/*   //strcpy(inputcat,"NONE"); */
/*   //  c+=f_kfls(filepar,"INPUTCAT",inputcat); */

  
/*   //printf("El tm  %f \n",G.reff); */
/*   //printf("El tm  %f \n",G.mr); */
/*   //printf("El   %f \n",G.mr); */
  I.alfa_prism=I.alfa_prism/206265*3600;
  
  if(c!=45) {
    printf("Not enough parameters in parameters file\n EXITING\n");
    printf("La prxima vez........\n");
/*     //exit(1); */
  }

}

void InputOptions()
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
  getline(I.qe_file,50);
  

  printf("File with filter response: ");
  getline(I.filter_file,50);
  printf("File with prism refraction index: ");
  getline(I.prism_file,50);
  printf("Total efficiency of the system: ");
  I.eff=readf(1.); 
  printf("Seeing of the observations: ");
  O.seeing=readf(1.); 
  printf("File with sky brightness in W/(m·m·s·Angstrom·arc·arc): ");
  getline(O.sky_file,50);
  printf("File with spectral distribution: ");
  getline(S.esp_file,51); 



  printf("Beginning wavelenght of the algorithm: ");
  O.ldomin=readf(1000.); 
  printf("End  wavelenght of the algorithm: ");
  O.ldomax=readf(10000.); 
  printf("Number of steps in wavelentgh calculation: ");
  O.nldo=readf(1000.); 

  printf("Output image file: ");
  getline(imgfile,51);
  printf("Output catalogue file: ");
  getline(imgfile,51);
  
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

void Readspec( char file[100])
{
  FILE *fp;
  int i;       
  S.nldo=FileNLin(file);
  if((fp=fopen(file,"r"))==NULL) {
    printf("Cannot open file %s\n",file);
    exit(1);
  }
  if((S.ldo=malloc(S.nldo*sizeof(float)))== NULL) {
    printf("I can't dimension the vector P.ldo of %d elements",S.nldo);
    exit(1);
  }
  if((S.flux=malloc(S.nldo*sizeof(float)))== NULL) {
    printf("I can't dimension the vector P.y of %d elements",S.nldo);
    exit(1);
  }
  for (i=0;i<S.nldo;i++) {
    fscanf(fp," %f %f",S.ldo+i,S.flux+i);

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

float sed(float ldo)
{
  return(Lagr2(S.ldo,S.flux,S.nldo,ldo));
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
  return(0);
}




