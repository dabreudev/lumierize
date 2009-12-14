#include "modulos.h"



struct plot_main {
  float xmin,xmax,ymin,ymax;
  int limflag;
};



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


struct sel {
  int idx;
/*   //int numb; */
/*   //int numn; */
  double rab,ran;  /*      //Coordenadas ecuatoriales en  */
  double decb,decn; /*     //ambas imagenes */
  float magbroad;
  float errmagbroad;
  float magnarrow;
  float errmagnarrow;
  float ewwhole,errewwhole;
  float ewaper;
  float errewaper;
  float fluxline;
  float errfluxline;
  float fluxaperline,errfluxaperline;
  float z;
  float errz;
  float nsig;
  float xbroad,ybroad;
  float xnarrow,ynarrow;
  float pflag;
  float fluxbroad,errfluxbroad;
  float fluxnarrow,errfluxnarrow;
  float fluxaperbroad,errfluxaperbroad;
  float fluxapernarrow,errfluxapernarrow;
  float magaperbroad,errmagaperbroad;
  float magapernarrow,errmagapernarrow;


};

struct mcal {
  float Ab,Bb;  /*      //Constantes de la recta de Bouger */
  float An,Bn;  /*      //mag=A -2.5log(cts/s)+Bsec(z) */
  float texpbroad,texpnarrow;
  float seczb,seczn;
};

struct filter {
  float *ldo;
  float *y;
  int n;
};

struct param {
/*   // Parametros de entrada del programa */
  char catfile[51];              /*      //Fichero con objetos detectados */
  char narrowfile[51];         /*     //Imagen en banda estrecha */
  char broadfile[51];            /*   //Imagen en banda ancha */
  int colnumb,colnumn;
  int colxpb,colypb;
  int colxpn,colypn;
  int colfb,colfn;
  int colefb,colefn;
  int coleab,colebb,coletb;
  int colean,colebn,coletn;
  int colskyb,colskyn;
  float rfixap;
  char broadfilter[51],narrowfilter[51];
  int verbose;

};



/* //Variables para la astrometria */
struct WorldCoor *wcsbroad,*wcsnarrow;
char *headbroad,*headnarrow;
int lhead,nbfits;  
int wcsflag=1;


struct th thress;
struct plot_main pm;
struct param p;
struct sel *objs;
struct mcal magcal;
struct filter fbroad,fnarrow;

int nsel;
FILE *fc;
fitsfile *narrowfits;             
fitsfile *broadfits;              
float *imgbroad;
float *imgnarrow;
int nxb,nyb,nxn,nyn;
int thflag=0;
int magcalflag=0;
int aperflag=0;  /*   //0: Given fluxes. 1: Compute aperture fluxes */
float *numb,*numn,*xpb,*ypb,*xpn,*ypn,*fb,*fn,*efb,*efn;
float *fapb,*fapn,*efapb,*efapn,*skyb,*skyn; /* //Flujos de apertura */
float *elipab,*elipbb,*eliptb,*elipan,*elipbn,*eliptn;
float *nsig;
int *true;
float *x,*y;
int nobj;




int SeeObject(int i);
void SelectObject(int i);
/* //void GetOption_kbd(); */
/* //void GetOption_cbut(); */

void ReadMagCal();
void PlotThress();
void SaveSelected(char selfile[],char filepar[]);
void ReadAstrom();
void ReadImages();
void ReadBothCat();
void ReadOneCat();
void Readfilter( char file[100],struct filter F);
void ComputeAperFot();
void GetObjectOnClick(int *i);
void ThressholdPoly();
void ThressholdMean();
void ThressholdSpline();
void ThressholdTheor();
void SelectObjThres();
void PlotMain();
void ChangePlotMain();
void LoadParam_file(char file[100]);
void LoadParam_kbd();
void SaveParam();
void EW(float fb,float fn,float efb,float efn,float *ew,float *errew);
void MAG(float fb,float fn,float efb,float efn,float *mb,float *emb,float *mn,float *emn);
void FLUXLINE(float fb,float fn,float efb,float efn,float *fl,float *efl);

/* //void ReadParam_file(struct param p); */
/* //void ReadParam_kbd(struct param p); */
/* //void SaveParam(struct param p); */



int main(int argc, char **argv)
{
  int ii;
  char opt='A';
  char selfile[51];
/*   if((objs=malloc(66*sizeof(struct sel)))==NULL) printf("I cannot dimension objs   of %d elements \n",nobj); */

/*   aperflag=0; */
/*   SaveSelected("mierda","paramertros"); */
/*   exit(1); */

  
/*   //Inicializo algunas variables y flags */

  p.verbose=0;
  thflag=0;
  pm.limflag=1;
  thress.ngrad=1;
  thress.ngradsig=2;
  thress.nbin=10;

  thress.nsig=3;
  p.rfixap=0;
  if((objs=malloc(1*sizeof(struct sel)))==NULL) printf("I cannot dimension objs   of %d elements \n",1);
  if((thress.coef   =malloc(1*sizeof(float)))==NULL) printf("I cannot dimension coef     of %d elements \n",1);
  if((thress.coefsig=malloc(1*sizeof(float)))==NULL) printf("I cannot dimension coefsig  of %d elements \n",1);
  if((thress.x      =malloc(1*sizeof(float)))==NULL) printf("I cannot dimension coef     of %d elements \n",1   );
  if((thress.ym     =malloc(1*sizeof(float)))==NULL) printf("I cannot dimension coefsig  of %d elements \n",1      );
  if((thress.sig    =malloc(1*sizeof(float)))==NULL) printf("I cannot dimension coefsig  of %d elements \n",1 );

  if(argc<2) {
    LoadParam_kbd();
  }
  else LoadParam_file(argv[1]);
  ReadImages();
  ReadBothCat();
  ReadMagCal();
  Readfilter(p.broadfilter,fbroad);
  Readfilter(p.narrowfilter,fnarrow);
  rcpgbegok("?",0);
  cpgscf(2);
/*   //cpgopen("?"); */
  cpgask(0);


/*   //printf("Opt %c\n",opt); */
  while(opt!='E') {

    PlotMain();
    if(thflag) PlotThress();
    printf(" G See object nearest cursor\n");        
    printf(" N See object with a given number\n");       
    printf(" P Change main plot\n");
    printf(" O Select objects for given thresholds\n");
    printf(" T Compute thresholds statitically\n");
    printf(" C Compute thresholds with theorical curves\n");
    printf(" W Write file of selected objects\n");
    printf(" M Input magnitude calibration\n");
    printf(" S Save parameter file\n");
    printf(" E Exit\n");                            
    
    printf(" Choose your option: ");
    opt=readc('G');
/*     //printf(" No has ele\n"); */
    switch(opt) {
    case 'M':
    case 'm':
      ReadMagCal();
      break;
    case 'O':
    case 'o':
      SelectObjThres(thress);
      break;
    case 'G':
    case 'g':
      printf(" justo antes de getobject\n");
      GetObjectOnClick(&ii);
      if(SeeObject(ii)) {
	printf(" Object selected\n");
	SelectObject(ii);
      }
      break;
    case 'P':
    case 'p':
      ChangePlotMain();
      break;
    case 'N':
    case 'n':
      ii=readi(0);
      if(SeeObject(ii))  SelectObject(ii);
      break;
    case 'W':
    case 'w':
      printf(" Input output file with selected objects\n");
      reads(selfile,selfile);
      SaveSelected(selfile,argv[1]);
      break;
    case 'E':
    case 'e':
      cpgend();
      exit(1);
      break;
    case 'T':
    case 't':
/*       //ThressholdStat(&thress); */
      thflag=1;
      ThressholdMean(&thress);
/*       //ThressholdPoly(thress); */
      break;
    case 'C':
    case 'c':
/*       //ThressholdTheor(thress); */
      break;
    case 'S':
    case 's':
      SaveParam();
      break;
    }
  }
  exit(1);
}



int SeeObject(int i)
{

  char snul[200];
  int crossflag=1;
  float tr[6];
  float median,sigma;
  float bgb=0,fgb=0,bgn=0,fgn=0;
  float *buffb,*buffn;
  int j,k,lb,ln;
  int autoflag=1;
  int zoomflag=0;
/*   //Subrutina para ver toda la informacion disponible sobre un objeto */
/*   float fluxbroad,fluxnarrow; */
/*   float errfluxbroad,errfluxnarrow; */
  float magbroad,errmagbroad,magnarrow,errmagnarrow;
  float ewwhole,errewwhole;
  float ewaper,errewaper;
  float x1b,y1b,x2b,y2b;
  float x1n,y1n,x2n,y2n;
  float xbox,ybox;
  char opt='A',cnul;  
  float xc,yc;
  int nb;
  double xpm,ypm,xps,yps;
  double ra,dec;
  int off;
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  printf(" Numero %d\n",i);

  xbox=40;ybox=40;

  while(opt!='S' && opt!='E') {
    if(!zoomflag) {
      x1b=xpb[i]-xbox;
      x2b=xpb[i]+xbox;
      y1b=ypb[i]-ybox;
      y2b=ypb[i]+ybox;
    }

    if(aperflag) {
      xpm=(double)xpb[i];ypm=(double)ypb[i];
      off=0;
      if(wcsflag) {
	pix2wcs(wcsbroad,xpm-1, ypm-1 ,&ra, &dec);
	/*       //El -1 es por la cagada de las columnas */
	  wcs2pix(wcsnarrow,ra,dec,&xps,&yps,&off);
	  xps=xps+1;yps=yps+1; /* //Tambien por la cagada */
	  printf(" ASTROMETRY %f %f \n",ra,dec);
      }
      else {
	xps=xpb[i];
	yps=ypb[i];
      }
      if(off) {
	xps=0;
	yps=0;
      }

      printf(" Calculated coordinates in NARROW %f,%f\n",xps,yps);
      
    }
    else {
      xps=(double)xpn[i];
      yps=(double)ypn[i];
    }

    x1n=(float)xps-xbox;x2n=(float)xps+xbox;
    y1n=(float)yps-ybox;y2n=(float)yps+ybox;
    printf(" Broad  image region  xmin %f xmax %f ymin %f ymax %f\n",x1b,x2b,y1b,y2b);
    printf(" Narrow image region  xmin %f xmax %f ymin %f ymax %f\n",x1n,x2n,y1n,y2n);
    if((buffb=malloc(((int)x2b-(int)x1b)*((int)y2b-(int)y1b)*sizeof(float)))==NULL) printf("I cannot dimension buffb   of %d elements \n",((int)x2b-(int)x1b));
    if((buffn=malloc(((int)x2n-(int)x1n)*((int)y2n-(int)y1n)*sizeof(float)))==NULL) printf("I cannot dimension buffn   of %d elements \n",((int)x2b-(int)x1b));
    lb=0;
    for(j=(int)x1b;j<(int)x2b;j++) {
      for(k=(int)y1b;k<(int)y2b;k++) {
/* 	//buffb[lb]=imgbroad[j+k*((int)x2b-(int)x1b)]; */
	buffb[lb]=imgbroad[j+k*nxb];
/* 	//printf(" bif %d %d %f\n",j,k,buffb[lb]); */
/* 	//buffb[lb]=imgbroad[j*((int)y2b-(int)y1b)+k]; */
	lb++;
      }
    }
    ln=0;
    for(j=(int)x1n;j<(int)x2n;j++) {
      for(k=(int)y1n;k<(int)y2n;k++) {
/* 	//buffn[ln]=imgnarrow[j+k*((int)x2n-(int)x1n)]; */
	buffn[ln]=imgnarrow[j+k*nxn];
	ln++;
      }
    }
/*     //printf(" lb %d el orto %d\n",lb,((int)x2b-(int)x1b)*((int)y2b-(int)y1b)); */

/*     //printf(" Pasa\n"); */
    if(autoflag) {
      median=StMedia(lb-1,buffb,&sigma);
      bgb=median-sigma*4.5;
      fgb=median+sigma*2.5;
/*       //printf(" bro Median %f sigma %f\n",median,sigma); */
      median=StMedia(ln,buffn,&sigma);
      bgn=median-sigma*4.5;
      fgn=median+sigma*2.5;
/*       //printf(" nar Median %f sigma %f\n",median,sigma); */
    }
    cpgpage();
    cpgsvp(0.1,0.4,0.6,0.9);
/*     //cpgswin(x1b,x2b,y1b,y2b); */
    cpgwnad(x1b,x2b,y1b,y2b);
    cpggray(imgbroad,nxb,nyb,(int)x1b,(int)x2b,(int)y1b,(int)y2b,fgb,bgb,tr);
    strcpy(snul,"Broad band image \0");
    strcat(snul,p.broadfile);
    cpglab("X axis","Y axis",snul);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
/*     //cpglab("X pixel","Y pixel","Broad band image"); */
    cpgsci(2);
    printf(" fixap %f\n",p.rfixap);
    if(crossflag) {
      if(aperflag) cpgelip(xpb[i],ypb[i],p.rfixap,p.rfixap,0.);
      else         cpgelip(xpb[i],ypb[i],elipab[i],elipbb[i],eliptb[i]);
      cpgsch(3.); cpgsci(2);
      cpgpt1(xpb[i],ypb[i],2);
      cpgsci(1); cpgsch(1.);
    }
    cpgsci(1);
    cpgsvp(0.6,0.9,0.6,0.9);
/*     //cpgswin(x1n,x2n,y1n,y2n); */
    cpgwnad(x1n,x2n,y1n,y2n);
    cpggray(imgnarrow,nxn,nyn,(int)x1n,(int)x2n,(int)y1n,(int)y2n,fgn,bgn,tr);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
    strcpy(snul,"Narrow band image \0");
    strcat(snul,p.narrowfile);
    cpglab("X axis","Y axis",snul);
/*     //cpglab("X pixel","Y pixel","Narrow band image"); */
    cpgsci(2);
    if(crossflag) {
      if(aperflag) cpgelip((float)xps,(float)yps,p.rfixap,p.rfixap,0.);
      else         cpgelip((float)xps,(float)yps,elipan[i],elipbn[i],eliptn[i]);
      cpgsch(3.);
      if(aperflag) {
	cpgsci(2);
	cpgpt1(xps,yps,2);
      }
      cpgsci(4);
      cpgpt1(xpn[i],ypn[i],2);
      cpgsci(1);cpgsch(1.);
    }
    free(buffb);free(buffn);

    /*     cpgtext(); */
    /*     EW(fluxbroad,fluxnarrow,errfluxbroad,errfluxnarrow,&ew,&errew); */
    /*    MAG(fluxbroad,fluxnarrow,errfluxbroad,errfluxnarrow,&magbroad,&errmagbroad,&magnarrow,&errmagnarrow); */
    if(aperflag) EW(fapb[i],fapn[i],efapb[i],efapn[i],&ewaper,&errewaper);
    EW(fb[i]  ,fn[i]  ,efb[i]  ,efn[i]  ,&ewwhole,&errewwhole);
    MAG(fb[i],fn[i],efb[i],efn[i],&magbroad,&errmagbroad,&magnarrow,&errmagnarrow);
    printf(" Object number in broad image: %6.0f\n",numb[i]);
    printf(" Object number in narrow image: %6.0f\n",numn[i]);
    printf(" Position:  Broad band image : X= %f Y= %f\n",xpb[i],ypb[i]);
    printf("            Narrow band image: X= %f Y= %f\n",xpn[i],ypn[i]);

    if(aperflag) {
      printf(" Aperture flux in broad band image %f counts\n",fapb[i]);
      printf(" Aperture flux in narrow band image %f counts\n",fapn[i]);
    }
    printf(" Isophotal flux in broad band image %f counts\n",fb[i]);
    printf(" Isophotal flux in narrow band image %f counts\n",fn[i]);
    if(aperflag) printf(" Equivalent width inside aperture : %f +/- %f\n",ewaper,errewaper); 
    printf(" Equivalent width for the whole object : %f +/- %f\n",ewwhole,errewwhole); 
    printf(" Apparent broad magnitude: %f +/- %f\n",magbroad,errmagbroad); 
    printf(" Apparent narrow magnitude: %f +/- %f\n",magnarrow,errmagnarrow); 
    if(aperflag) {
      MAG(fapb[i],fapn[i],efapb[i],efapn[i],&magbroad,&errmagbroad,&magnarrow,&errmagnarrow);
      printf(" Apparent broad magnitude inside aperture phot: %f +/- %f\n",magbroad,errmagbroad); 
    }
    if(thflag) printf("\n Number of sigmas above median : %f\n",nsig[i]);


/*     //Botones */
    cbuttsbr(0.,1.,0.2,0.4);
    cbutton(1,"Automatic",0);
    cbutton(2,"Set BG & FG",0);
    cbutton(3,"Zoom in x2",0);
    cbutton(4,"Zoom out x2",0);
    cbutton(5,"Mouse zoom",0);
    cbutton(6,"Crosses ON/OFF",0);
    cbutton(7,"Delete",0);
    cbutton(8,"Select",0);
    cbutton(9,"EXIT",0);


/*     //Menu de seleccion */
    printf(" A Automatic grey scaling (bg & fg)\n");
    printf(" B Set background & foreground manually\n");
    printf(" I Zoom in x2\n");
    printf(" O Zoom out x2\n");
    printf(" Z Zoom with mouse\n");
    printf(" C Crosses ON/OFF\n");
    printf(" D Delete object for ever & Exit\n");
    printf(" S Select Object & Exit\n");
    printf(" E Exit without selecting Object\n");

    printf(" Broad band image cuts: %f -  %f\n",bgb,fgb);
    printf(" Narrow band image cuts: %f -  %f\n",bgn,fgn);
    
/*     //opt=readc('E'); */
    cpgband(0,0,0.,0.,&xc,&yc,&cnul);
    cifbutton(xc,yc,&nb);
    if(nb==1) opt='A';
    if(nb==2) opt='B';
    if(nb==3) opt='I';
    if(nb==4) opt='O';
    if(nb==5) opt='Z';
    if(nb==6) opt='C';
    if(nb==7) opt='D';
    if(nb==8) opt='S';
    if(nb==9) opt='E';
    switch(opt) {
    case 'C':
    case 'c':
      crossflag=1-1*crossflag;
      break;
    case 'D':
    case 'd':
      true[i]=0;
      return(0);
      break;
    case 'A':
    case 'a':
      autoflag=1;
      break;
    case 'B':
    case 'b':
      autoflag=0;
      printf(" Broad image\n");
      printf(" Background ");bgb=readf(bgb);
      printf(" Foreground ");fgb=readf(fgb);
      printf(" Narrow image\n");
      printf(" Background ");bgn=readf(bgn);
      printf(" Foreground ");fgn=readf(fgn);
      break;
    case 'I':
    case 'i':
      xbox/=2;ybox/=2;
      zoomflag=0;
/*       //y1b=(int)y-ybox;y2b=(int)y-ybox; */
/*       //x1b=(int)x-xbox;x2b=(int)x-xbox; */
/*       //y1n=(int)y-ybox;y2n=(int)y-ybox; */
/*       //x1n=(int)x-xbox;x2n=(int)x-xbox; */
      break;
    case 'O':
    case 'o':
      xbox*=2;ybox*=2;
      zoomflag=0;
/*       //y1b=(int)y-ybox;y2b=(int)y-ybox; */
/*       //x1b=(int)x-xbox;x2b=(int)x-xbox; */
/*       //y1n=(int)y-ybox;y2n=(int)y-ybox; */
/*       //x1n=(int)x-xbox;x2n=(int)x-xbox; */
      break;
    case 'Z':
    case 'z':
      zoomflag=1;
      printf(" Click on right image the bottom left square\n");
      cpgcurs(&x1b,&y1b,&cnul);
      printf(" Click on right image the upper right square\n");
      cpgsci(2);
      cpgband(2,1,x1b,y1b,&x2b,&y2b,&cnul);
      cpgsci(1);
      break;
    case 'S':
    case 's':
      return(1);
      break;
    case 'E':
    case 'e':
      return(0);
      break;
    }	
    

  
  }
  
  return(0);
  
}

void ReadImages() 
{
  
  /* Varaibales para leer el FITS */
  int status=0;
  int nfound, anynull;
/*   fitsfile *image; */
  long naxes[2], fpixel, npixels;
/*   long  ii,jj,nbuffer; */
  float nullval;
/*   float datamin, datamax; */
/*   float *buffer; */
  
  if(ffopen(&broadfits,p.broadfile , READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(broadfits, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  npixels=naxes[0]*naxes[1];
/*   //printf(" npixels %d\n",npixels); */
  if((imgbroad=malloc(npixels*sizeof(float)))==NULL) printf("I cannot dimension buffer of %ld elements \n",npixels);
  fpixel=1;
  nullval=0;
/*   //datamin=1.0e30;datamax=-1.0e30; */
  printf("\n...Reading image %s \n",p.broadfile);
  if(fits_read_img(broadfits, TFLOAT, fpixel, npixels, &nullval, imgbroad, &anynull, &status )) fits_report_error(stderr,status);
/*   //  printf("...Computing datamin and datamax \n"); */
/*   //for (ii=0 ;ii<npixels;ii++) { */
/*   //  if( buffer[ii]< datamin) datamin=buffer[ii]; */
/*   //  if(buffer[ii]> datamax) datamax=buffer[ii]; */
/*   //} */
/*   //printf(" Datamin %f Datamax %f \n",datamin,datamax); */
/*   //mean=StMedia(npixels,buffer,&sigma); */
/*   //datamin=mean-sigma*2; */
/*   //datamax=mean+sigma*2; */
  nxb=naxes[0];nyb=naxes[1];
 
  status=0;
  if(ffopen(&narrowfits,p.narrowfile , READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(narrowfits, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  npixels=naxes[0]*naxes[1];
/*   //printf(" npixels %d\n",npixels); */
  if((imgnarrow=malloc(npixels*sizeof(float)))==NULL) printf("I cannot dimension buffer of %ld elements \n",npixels);
  fpixel=1;
  nullval=0;
  printf("...Reading image %s \n",p.narrowfile);
  if(fits_read_img(narrowfits, TFLOAT, fpixel, npixels, &nullval, imgnarrow, &anynull, &status )) fits_report_error(stderr,status);
/*   //printf(" Finalizado\n"); */
  nxn=naxes[0];nyn=naxes[1];
}



void ReadBothCat()
{
/*   //int nobj; */
  int j;
  nobj=FileNLin(p.catfile);
  printf("...Reading catalogue file %s\n",p.catfile);
 
  if((numb=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension numb    of %d elements \n",nobj);
  if((numn=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension numn    of %d elements \n",nobj);

  if((xpb   =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xpb    of %d elements \n",nobj);
  if((ypb   =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension ypb    of %d elements \n",nobj);
  if((xpn   =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xpn    of %d elements \n",nobj);
  if((ypn   =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension ypn    of %d elements \n",nobj);
  if((fb    =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension fb     of %d elements \n",nobj);
  if((fn    =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension fn     of %d elements \n",nobj);
  if((efb   =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension efb    of %d elements \n",nobj);
  if((efn   =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension efn    of %d elements \n",nobj);
  if((fapb  =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension fapb   of %d elements \n",nobj);
  if((fapn  =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension fapn   of %d elements \n",nobj);
  if((efapb =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension efapb  of %d elements \n",nobj);
  if((efapn =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension efapn  of %d elements \n",nobj);
  if((elipab=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension elipab of %d elements \n",nobj);
  if((elipbb=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension elipbb of %d elements \n",nobj);
  if((eliptb=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension eliptb of %d elements \n",nobj);
  if((elipan=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension elipan of %d elements \n",nobj);
  if((elipbn=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension elipbn of %d elements \n",nobj);
  if((eliptn=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension eliptn of %d elements \n",nobj);
  if((skyb  =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension skyb   of %d elements \n",nobj);
  if((skyn  =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension skyn   of %d elements \n",nobj);
  if((true  =malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension true    of %d elements \n",nobj);
  if((x     =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension x    of %d elements \n",nobj);
  if((y     =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension y    of %d elements \n",nobj);
  if((nsig  =malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension nsig    of %d elements \n",nobj);

  for(j=0;j<nobj;j++) true[j]=0;


/*   //printf(" Antes delle¡ nln %d\n",nobj); */

  ReadNumcol(p.catfile, p.colnumb,numb,true,&nobj);
/*   //printf(" perofrms \n"); */
  ReadNumcol(p.catfile, p.colnumn,numn,true,&nobj);
  ReadNumcol(p.catfile, p.colxpb,xpb,true,&nobj);
  ReadNumcol(p.catfile, p.colypb,ypb,true,&nobj);
  ReadNumcol(p.catfile, p.colxpn,xpn,true,&nobj);
  ReadNumcol(p.catfile, p.colypn,ypn,true,&nobj);
  ReadNumcol(p.catfile, p.colfb ,fb ,true,&nobj);
  ReadNumcol(p.catfile, p.colfn ,fn ,true,&nobj);
  ReadNumcol(p.catfile, p.colefb,efb,true,&nobj);
  ReadNumcol(p.catfile, p.colefn,efn,true,&nobj);
  ReadNumcol(p.catfile, p.coleab,elipab,true,&nobj);
  ReadNumcol(p.catfile, p.colebb,elipbb,true,&nobj);
  ReadNumcol(p.catfile, p.coletb,eliptb,true,&nobj);
  ReadNumcol(p.catfile, p.colean,elipan,true,&nobj);
  ReadNumcol(p.catfile, p.colebn,elipbn,true,&nobj);
  ReadNumcol(p.catfile, p.coletn,eliptn,true,&nobj);
/*   //printf(" Todo igual sky %d %d\n",p.colskyb,p.colskyn); */
  ReadNumcol(p.catfile, p.colskyb,skyb,true,&nobj);
  ReadNumcol(p.catfile, p.colskyn,skyn,true,&nobj);
/*   //printf(" Ternima \n"); */
  ReadAstrom();
  if(p.rfixap!=0) {
    ComputeAperFot(); /* //Hacemos fotometria de apertura */
    aperflag=1;
  }

/*   //printf(" Aqui esta bien\n"); */
/*   //printf(" Operacion\n"); */
  for(j=0;j<nobj;j++) {
    if(true[j]) {
      if(aperflag) {
	x[j]=-2.5*log10(fapb[j]);
	y[j]=-2.5*log10(fapb[j]/fapn[j]);
      }
      else{
	x[j]=-2.5*log10(fb[j]);
	y[j]=-2.5*log10(fb[j]/fn[j]);
	exit(1);
      }
    }
  }


}


void GetObjectOnClick(int *i)
{
   float xcur,ycur;
   char cnul;
   float mindist=1.e15,dist;
   int j;
   printf(" Entra en esta getobecjtonclik\n");
   cpgcurs(&xcur,&ycur,&cnul);
   printf(" No pasa de qauiq\n");
   for(j=0;j<nobj;j++) {
     if(true[j]) {
       dist=((x[j]-xcur)*(x[j]-xcur)+(y[j]-ycur)*(y[j]-ycur));
       if(mindist*mindist>dist*dist) {
	 mindist=dist;
	 *i=j;
       }
     }
   }
   printf(" Cursor position %f %f\n",xcur,ycur);
   printf(" Nearest object positio %f %f\n",x[*i],y[*i]);
   cpgsci(2);
   cpgsch(2.);
   cpgpt1(x[*i],y[*i],3);
   cpgsci(1);
   cpgsch(1.);
   printf(" Press <CR> to continue\n");
   getline(&cnul,2);
}


void SelectObject(int i)
{
  double xp,yp;
/*   //  double ra,dec; */
  

  nsel++;
  if((objs=realloc(objs,nsel*sizeof(struct sel)))==NULL) printf("I cannot dimension objs   of %d elements \n",nsel);
  objs[nsel-1].idx=numb[i];
  if(aperflag) { 
    EW(fapb[i],fapn[i],efapb[i],efapn[i], &(objs[nsel-1].ewaper),&(objs[nsel-1].errewaper));
    MAG(fapb[i],fapn[i],efapb[i],efapn[i], &(objs[nsel-1].magaperbroad),&(objs[nsel-1].errmagaperbroad), &(objs[nsel-1].magnarrow),&(objs[nsel-1].errmagnarrow));
    FLUXLINE(fapb[i],fapn[i],efapb[i],efapn[i], &(objs[nsel-1].fluxaperline),&(objs[nsel-1].errfluxaperline));
  }
  EW(fb[i],fn[i],efb[i],efn[i], &(objs[nsel-1].ewwhole),&(objs[nsel-1].errewwhole));
  MAG(fb[i],fn[i],efb[i],efn[i], &(objs[nsel-1].magbroad),&(objs[nsel-1].errmagbroad), &(objs[nsel-1].magnarrow),&(objs[nsel-1].errmagnarrow));
  FLUXLINE(fb[i],fn[i],efb[i],efn[i], &(objs[nsel-1].fluxline),&(objs[nsel-1].errfluxline));
/*   //objs[nsel-1].magbroad=-2.5*log10(fb[i]);      //Esta mal */
/*   //objs[nsel-1].magnarrow=-2.5*log10(fn[i]);     //Esta mal */
/*   //objs[nsel-1].ew=(80.+fn[i])/(1000.);          //Esta mal */
  objs[nsel-1].xbroad=xpb[i];
  objs[nsel-1].ybroad=ypb[i];
  objs[nsel-1].xnarrow=xpn[i];
  objs[nsel-1].ynarrow=ypn[i];
  objs[nsel-1].z=0.24;                         /*  Esta mal */
  
  if(thflag) {
    objs[nsel-1].nsig=nsig[i];               /*  Esta mal */
  }
  else {
    objs[nsel-1].nsig=0.;
  }
  objs[nsel-1].fluxbroad=fb[i];          
  objs[nsel-1].fluxnarrow=fn[i];         
  
    
/*   //objs[nsel-1].numb=(int)numb[i]; */
/*   //objs[nsel-1].numn=(int)numn[i]; */
  if(wcsflag) {
    xp=(double)xpb[i];yp=(double)ypb[i];
    pix2wcs(wcsbroad ,xp-1, yp-1 ,&(objs[nsel-1].rab), &(objs[nsel-1].decb));
    xp=(double)xpn[i];yp=(double)ypn[i];
    pix2wcs(wcsnarrow,xp-1, yp-1 ,&(objs[nsel-1].ran), &(objs[nsel-1].decn));
  }
  else {
    objs[nsel-1].rab=0.;objs[nsel-1].decb=0.;objs[nsel-1].ran=0.;objs[nsel-1].decn=0.;
  }
/*   //El -1 es por la cagada */
/*   //Pongo WFC para acordame de arreglar eso. ES IMPORTANTE!! */

  
/*   //Esto es para los polinomios. */
/*   //objs[nsel-1].nsig=fabs(y[i]-poly(x[i],thress.ngrad,thress.coef)); */
/*   //objs[nsel-1].nsig/=(poly(x[i],thress.ngradsig,thress.coefsig)-poly(x[i],thress.ngrad,thress.coef)); */
  
  printf("Input oject priority (0-1): ");
  objs[nsel-1].pflag=readf(1.);
  printf(" Objeto seleccionado numero %d xpos %f\n",i,xpb[i]);
}



void PlotMain()
{
  char snul[200];
  int j;
  int ntrue=0;
  float ewmin,ewmax;
  float fnul;

  if(pm.limflag) pm.xmin=1.e30,pm.xmax=-1.e30,pm.ymin=1.e30,pm.ymax=-1.e30;
  for(j=0;j<nobj;j++) {
    if(true[j]) {

/*       //printf("Uno mas\n"); */
      if(aperflag) {
	x[j]=-2.5*log10(fapb[j]);
	y[j]=-2.5*log10(fapb[j]/fapn[j]);
      }
      else{
	x[j]=-2.5*log10(fb[j]);
	y[j]=-2.5*log10(fb[j]/fn[j]);
/* 	//exit(1); */
      }
      if(pm.limflag) {
	if(pm.xmin>x[j]) pm.xmin=x[j];
	if(pm.xmax<x[j]) pm.xmax=x[j];
	if(pm.ymin>y[j]) pm.ymin=y[j];
	if(pm.ymax<y[j]) pm.ymax=y[j];
      }
    }
  }
  cpgpage();
  cpgsvp(0.1,0.9,0.1,0.9);
  EW(powf(10.,-0.4*(pm.ymax)),1.  ,0.  ,0.  ,&ewmax,&fnul);
  EW(powf(10.,-0.4*(pm.ymin)),1.  ,0.  ,0.  ,&ewmin,&fnul);


/*   ewmax=(1-magcal.texpbroad/magcal.texpnarrow*(powf(10.,0.4*(pm.ymax+magcal.Ab-magcal.An+magcal.Bb*magcal.seczb-magcal.Bn*magcal.seczn))))*p.deltalambda; */
/*   ewmin=(1-magcal.texpbroad/magcal.texpnarrow*(powf(10.,0.4*(pm.ymin+magcal.Ab-magcal.An+magcal.Bb*magcal.seczb-magcal.Bn*magcal.seczn))))*p.deltalambda; */
/*   printf(" EWMIN EWMAX %f %f\n",ewmin,ewmax); */
/*   printf(" PMMIN PMMAX %f %f\n",pm.ymin,pm.ymax); */
  cpgswin(pm.xmin,pm.xmax,ewmax,ewmin);
  cpgbox("BCTNS",0,0,"CTMS",0,0);
/*   printf(" EWMIN EWMAX %f %f\n",ewmin,ewmax); */
/*   printf(" xmin %f max %f ymon %f max %f\n",pm.xmin,pm.xmax,pm.ymin,pm.ymax); */
  cpgswin(pm.xmin,pm.xmax,pm.ymin,pm.ymax);
  cpgbox("BCTNS",0,0,"BTNS",0,0);
  
  for(j=0;j<nobj;j++) {
    if(true[j]) {
      ntrue++;
/*       //printf(" j %d br %f na %f x %f y %f\n",j,fapb[j],fapn[j],x[j],y[j]); */
      cpgpt1(x[j],y[j],17);   /* //WFC esta era 4 */
    }
  }
  printf(" nsel %d\n",nsel);
  cpgsch(2.5);
  cpgsci(5);
  for(j=0;j<nsel;j++) {
/*     //printf(" Sleccion %d  xpos %f \n",objs[j].idx,xpb[objs[j].idx]); */
    cpgpt1(x[objs[j].idx],y[objs[j].idx],5);
  }
  cpgsch(1.);
  cpgsci(1);

/*   //cpglab("I magnitude","2.5log(Flux Narrow/Flux Broad)",""); //WFC */
  strcpy(snul,"Catalogue file \0");
  strcat(snul,p.catfile);
  cpglab("Apparent magnitude","2.5log(Flux Narrow/Flux Broad)",snul);

  printf(" Total number of objects= %d\n",nobj);  
  printf(" Number of objects in plot= %d\n",ntrue);  

}

void ThressholdPoly()
{
/*   int j; */
/*   thress.xmin=1.e30; */
/*   thress.xmax=-1.e30; */
/*   for(j=0;j<nobj;j++) { */
/*     if(true[j]) { */
/*       if(thress.xmin>x[j]) thress.xmin=x[j]; */
/*       if(thress.xmax<x[j]) thress.xmax=x[j]; */
/*     } */
/*   } */

/*   thress.nbin=10; */
  
/*   //De momento, por el metodo del polinomio */
  float *xpol, *ypol;
  int npol=0,j;
  char opt='K';
  float xdraw,ydraw;
  float nrej=1e30;
  float sigabo;
  if((xpol=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xpol   of %d elements \n",nobj);
  if((ypol=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension ypol   of %d elements \n",nobj);
  


   while(opt!='E') {

     if((thress.coef   =realloc(thress.coef   ,(thress.ngrad+1)*sizeof(float)))==NULL) printf("I cannot dimension coef     of %d elements \n",thress.ngrad+1);
     if((thress.coefsig=realloc(thress.coefsig,(thress.ngradsig+1)*sizeof(float)))==NULL) printf("I cannot dimension coefsig  of %d elements \n",thress.ngradsig+1);
     for(j=0;j<thress.ngrad;j++)        thress.coef[j]=0;
     for(j=0;j<thress.ngradsig;j++)     thress.coefsig[j]=1;

    
     printf(" Ya he hecho el realloc\n");
     npol=0;
     for(j=0;j<nobj;j++) {
       if(true[j]) {
	 xpol[npol]=x[j];
	 sigabo=0;
	 sigabo=fabs(y[j]-poly(x[j],thress.ngrad,thress.coef));
	 sigabo/=fabs((poly(x[j],thress.ngradsig,thress.coefsig)-poly(x[j],thress.ngrad,thress.coef)));
	 printf(" sigabo %f\n",sigabo);
	 if(sigabo<nrej) {
	   ypol[npol]=y[j];
	   npol++;
	 }
       }
     }
     printf(" He saignado\n");
     if(npol!=0) MCPN(npol,xpol,ypol,thress.ngrad,thress.coef);
     printf(" pooly %f\n",poly(-10.,thress.ngrad,thress.coef));
     printf(" He ajustado\n");
     for(j=0;j<thress.ngrad+1;j++) printf(" coef %d %f\n",j,thress.coef[j]);
     printf(" Ya estan los coef\n");

     printf(" UFF\n");
     npol=0;
     for(j=0;j<nobj;j++) {
       if(true[j]) {
	 xpol[npol]=x[j];
	 sigabo=0;
	 sigabo=fabs(y[j]-poly(x[j],thress.ngrad,thress.coef));
	 sigabo/=fabs(poly(x[j],thress.ngradsig,thress.coefsig)-poly(x[j],thress.ngrad,thress.coef));
	 if(sigabo<nrej) {
	   ypol[npol]=(y[j]-poly(x[j],thress.ngrad,thress.coef))*(y[j]-poly(x[j],thress.ngrad,thress.coef));
	   npol++;
	   printf(" x %f y %f yy %g pp %g\n",x[j],y[j],ypol[npol-1],poly(x[j],thress.ngrad,thress.coef));
	 }
       }
     }
     printf(" He vuelto a agisna,\n");
     if(npol!=0) MCPN(npol,xpol,ypol,thress.ngradsig,thress.coefsig);
     for(j=0;j<thress.ngradsig+1;j++) printf(" coef  otro %d %f\n",j,thress.coefsig[j]);

     printf(" Se ha acabdo tood\n");
/*      //Dibujo el spline */
     xdraw=pm.xmin;ydraw=poly(pm.xmin,thress.ngrad,thress.coef);cpgmove(xdraw,ydraw);
     for(j=0;j<100;j++) {
       xdraw=(pm.xmax-pm.xmin)/99.*j+pm.xmin;
       ydraw=poly(xdraw,thress.ngrad,thress.coef);
       cpgdraw(xdraw,ydraw);
     }
/*      //Dibujo +nsig */
     xdraw=pm.xmin;ydraw=poly(xdraw,thress.ngrad,thress.coef)+sqrt(fabs(poly(xdraw,thress.ngradsig,thress.coefsig)));
     cpgmove(xdraw,ydraw);
     for(j=0;j<100;j++) {
       xdraw=(pm.xmax-pm.xmin)/99.*j+pm.xmin;
       ydraw=poly(xdraw,thress.ngrad,thress.coef)+sqrt(fabs(poly(xdraw,thress.ngradsig,thress.coefsig)));
       cpgdraw(xdraw,ydraw);
     }
/*      //Dibujo -nsig */
     xdraw=pm.xmin;ydraw=poly(xdraw,thress.ngrad,thress.coef)-sqrt(fabs(poly(xdraw,thress.ngradsig,thress.coefsig)));
     cpgmove(xdraw,ydraw);
     for(j=0;j<100;j++) {
       xdraw=(pm.xmax-pm.xmin)/99.*j+pm.xmin;
       ydraw=poly(xdraw,thress.ngrad,thress.coef)-sqrt(fabs(poly(xdraw,thress.ngradsig,thress.coefsig)));
       cpgdraw(xdraw,ydraw);
     }


     printf(" P Change degree of polinomial\n");
     printf(" D Delete points above a given sigma\n");
     printf(" R Replot\n");
     printf(" E Exit\n");
     opt=readc('D');
     switch(opt) {
     case 'R':
     case 'r':
       PlotMain();
       break;
     case 'D':
     case 'd':
       printf(" New sigma threshold: ");
       nrej=readf(nrej);
       break;
     case 'P':
     case 'p':
       printf(" New polinomial degree: ");
       thress.ngrad=readi(thress.ngrad);
       printf(" New sigma polinomial degree: ");
       thress.ngradsig=readi(thress.ngradsig);
       break;
     case 'E':
     case 'e':
       break;
     }
   }
   
   free(xpol);free(ypol);
}


void ThressholdMean()
{
  float yobjth,sigobjth;
  int j,ib; 
  float *xbuf,*ybuf;
  int nbuf;
  float xdelt;
  char opt='D';
  float nrej=10;
  if((xbuf=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xbuf   of %d elements \n",nobj);
  if((ybuf=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension ybuf   of %d elements \n",nobj);
  for(j=0;j<nobj;j++) {
    if(true[j]) nsig[j]=0;
  }

  while(opt!='E') {
    if((thress.x    =realloc(thress.x   ,thress.nbin*sizeof(float)))==NULL) printf("I cannot dimension x        of %d elements \n",thress.nbin);
    if((thress.ym   =realloc(thress.ym  ,thress.nbin*sizeof(float)))==NULL) printf("I cannot dimension ym       of %d elements \n",thress.nbin);
    if((thress.sig  =realloc(thress.sig ,thress.nbin*sizeof(float)))==NULL) printf("I cannot dimension sig      of %d elements \n",thress.nbin);
    thress.xmin=pm.xmin; 
    thress.xmax=pm.xmax; 
/*     //printf(" cmin %f xmax %f\n",thress.xmin,thress.xmax); */
    for(ib=0;ib<thress.nbin;ib++) thress.x[ib]=(thress.xmax-thress.xmin)*(ib+.5)/thress.nbin+thress.xmin;
    xdelt=(thress.xmax-thress.xmin)/thress.nbin;
    
    for(ib=0;ib<thress.nbin;ib++) {
      nbuf=0;
      for(j=0;j<nobj;j++) {
	if(true[j]) {
	  if(x[j]>=thress.x[ib]-xdelt/2 && x[j]<thress.x[ib]+xdelt/2 && nsig[j]<nrej) {
	    xbuf[nbuf]=x[j];
	    ybuf[nbuf]=y[j];
	    nbuf++;
	  }
	}
      }
/*       //printf(" nbuf= %d\n",nbuf); */
      thress.ym[ib]=StMedia(nbuf,ybuf,&(thress.sig[ib]));
      printf(" x %f y %f sig %f\n",thress.x[ib],thress.ym[ib],thress.sig[ib]);
    }

/*     //Calculo nsig para todos lo objetos */
    for(j=0;j<nobj;j++) {
      if(true[j]) {
	yobjth=Lagr2(thress.x,thress.ym,thress.nbin,x[j]);
	sigobjth=Lagr2(thress.x,thress.sig,thress.nbin,x[j]);
	nsig[j]=(y[j]-yobjth)/sigobjth;
      }
    }
    PlotThress();
    
    printf(" B Change number of bins\n");
    printf(" S Change N sigma limit\n");
    printf(" D Delete points above a given sigma to compute mean and sigma\n");
    printf(" R Replot\n");
    printf(" E Exit\n");
    opt=readc('D');
    switch(opt) {
    case 'S':
    case 's':
      printf(" Number of N sigmas\n");
      thress.nsig=readf(thress.nsig);
      break;
    case 'B':
    case 'b':
      printf(" Number of bins: ");
      thress.nbin=readf(thress.nbin);
    case 'R':
    case 'r':
      PlotMain();
      break;
    case 'D':
    case 'd':
      printf(" New sigma threshold: ");
      nrej=readf(nrej);
      break;
    case 'E':
    case 'e':
      break;
    }
  }
  

  free(xbuf);free(ybuf);
}



void LoadParam_file(char file[100])
{


/*   int c,i,j; */
/*   int nlin; */
/*   char nul3[3],nul1[1]; */
/*   char keyf[9]="",key[9]=""; */
  int status=0;
/*   char string[51]; */
  char comment[51];
  fitsfile *parfile;
/*   int inte; */
/*   float crpix; */
  printf("Using parameter file <<%s>>\n",file);
  if( ffopen2(&parfile,file, READONLY, &status)) fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TSTRING,"BROADIMG",p.broadfile,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TSTRING,"NARROIMG",p.narrowfile,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TSTRING,"OBJCAT",p.catfile,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLIDX_B",&(p.colnumb),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLIDX_N",&(p.colnumn),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLXP_B",&(p.colxpb),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLXP_N",&(p.colxpn),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLYP_B",&(p.colypb),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLYP_N",&(p.colypn),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLFL_B",&(p.colfb),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLFL_N",&(p.colfn),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLEFL_B",&(p.colefb),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLEFL_N",&(p.colefn),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLELA_B",&(p.coleab),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLELB_B",&(p.colebb),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLELT_B",&(p.coletb),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLELA_N",&(p.colean),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLELB_N",&(p.colebn),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLELT_N",&(p.coletn),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLSKY_B",&(p.colskyb),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLSKY_N",&(p.colskyn),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"RFIXAP",&(p.rfixap),comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"RESBROAD",p.broadfilter,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"RESNARR",p.narrowfilter,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  fits_close_file(parfile,&status);

}


void SaveParam()
{
  FILE *fp;
/*   int i,j; */
  int nc=0,nt;
/*   char nul[2]; */
  char file[100];
  char ch51[51];
/*   char opt; */
  printf("Name of file: ");
  scanf("%s",file);
  if((fp=fopen(file,"w")) ==NULL) {
    printf("ERROR: Can't open file %s\n",file);
    return;
  }
  fprintf(fp,"COMMENT  Parameter file for NarrowSearch                                       \n");
  nc++;
  sprintf(ch51,"'%s'",p.broadfile);
  fprintf(fp,"BROADIMG= %-51.51s / Broad band img.\n",ch51);
  sprintf(ch51,"'%s'",p.narrowfile);
  fprintf(fp,"NARROIMG= %-51.51s / Narrow band img\n",ch51);
  sprintf(ch51,"'%s'",p.catfile);
  fprintf(fp,"OBJCAT  = %-51.51s / Object file    \n",ch51);
  fprintf(fp,"COLIDX_B=%21d / Column with identification in broad band image\n",p.colnumb);
  fprintf(fp,"COLIDX_N=%21d / Column with identification in narrow band img \n",p.colnumn);
  fprintf(fp,"COLXP_B =%21d / Column with X position in broad band image    \n",p.colxpb );
  fprintf(fp,"COLYP_B =%21d / Column with Y position in broad band image    \n",p.colypb );
  fprintf(fp,"COLXP_N =%21d / Column with X position in narrow band image   \n",p.colxpn );
  fprintf(fp,"COLYP_N =%21d / Column with Y position in narrow band image   \n",p.colypn );
  fprintf(fp,"COLFL_B =%21d / Column with flux in broad band image          \n",p.colfb  );
  fprintf(fp,"COLFL_N =%21d / Column with flux in narrow band image         \n",p.colfn  );
  fprintf(fp,"COLEFL_B=%21d / Column with error in flux broad band image    \n",p.colefb );
  fprintf(fp,"COLEFL_N=%21d / Column with error in flux narrow band image   \n",p.colefn );
  fprintf(fp,"COLELA_B=%21d / Column with A ellipse fit in broad band image \n",p.coleab );
  fprintf(fp,"COLELB_B=%21d / Column with B ellipse fit in broad band image \n",p.colebb );
  fprintf(fp,"COLELT_B=%21d / Column with angle ellipse in broad band image \n",p.coletb );
  fprintf(fp,"COLELA_N=%21d / Column with A ellipse fit in narrow band image\n",p.colean );
  fprintf(fp,"COLELB_N=%21d / Column with B ellipse fit in narrow band image\n",p.colebn );
  fprintf(fp,"COLELT_N=%21d / Column with angle ellipse in narrow band image\n",p.coletn );

  fprintf(fp,"COLSKY_B=%21d / Column with angle ellipse in narrow band image\n",p.colskyb);
  fprintf(fp,"COLSKY_N=%21d / Column with angle ellipse in narrow band image\n",p.colskyn);
  fprintf(fp,"RFIXAP  =%21f / Fixed aperture radius                         \n",p.rfixap );
  sprintf(ch51,"'%s'",p.broadfilter);
  fprintf(fp,"RESBROAD= %-51.51s / Broad response \n",ch51);
  sprintf(ch51,"'%s'",p.narrowfilter);
  fprintf(fp,"RESNARR = %-51.51s / Broad response \n",ch51);
  nc+=24;

  fprintf(fp,"COMMENT  Blank lines to assure the last block of 2880 caracters is complete    \n");
  fprintf(fp,"COMMENT                                                                        \n");
  nc+=2;
  
  for(nt=nc;nt<(1+(int)(nc/36))*36-1;nt++) {
    fprintf(fp,"COMMENT                                                                        \n");
  }
  fprintf(fp,"END                                                                            \n");
  fclose(fp);

}


void  EW(float fb,float fn,float efb,float efn,float *ew,float *errew)
{
  float fluxb,fluxn;
/*   float frt=1; */
  float transmitanciabroad=0;
  float transmitancianarrow=0;
  float deltabroad=0;
  float deltanarrow=0;
  float zapartirdelfiltroestrecho=0;
/*   float deltalambda=50; */

/*   //printf(" Flujo   br %f nar %f\n",fb,fn); */
  fluxb=fb/magcal.texpbroad *(powf(10.,0.4*(-magcal.Ab-magcal.Bb*magcal.seczb)));
  fluxn=fn/magcal.texpnarrow*(powf(10.,0.4*(-magcal.An-magcal.Bn*magcal.seczn)));
  printf(" fluxb %g fluxn %g\n",fluxb,fluxn);

  *ew=(fluxb*deltanarrow-fluxn*deltabroad)/(fluxn*transmitanciabroad-fluxb*transmitancianarrow);
  *ew=*ew/(1+zapartirdelfiltroestrecho);

  *errew=0. ;
}


void MAG(float fb,float fn,float efb,float efn,float *mb,float *emb,float *mn,float *emn)
{
/*   magcal.seczb=0; */
/*   magcal.seczn=0; */
  *mb=-2.5*log10(fb/magcal.texpbroad)+magcal.Ab+magcal.Bb*magcal.seczb;
  *emb=*mb;
  *mn=-2.5*log10(fn/magcal.texpnarrow)+magcal.An+magcal.Bn*magcal.seczn;
  *emn=*mn;
}



void FLUXLINE(float fb,float fn,float efb,float efn,float *fl,float *efl)
{

  float fluxb,fluxn;
  float transmitanciabroad=0;
  float transmitancianarrow=0;
  float deltabroad=0;
  float deltanarrow=0;

  fluxb=fb/magcal.texpbroad *(powf(10.,0.4*(-magcal.Ab-magcal.Bb*magcal.seczb)));
  fluxn=fn/magcal.texpnarrow*(powf(10.,0.4*(-magcal.An-magcal.Bn*magcal.seczn)));

  *fl=(fluxb*deltanarrow-fluxn*deltabroad)/(deltanarrow*transmitanciabroad-deltabroad*transmitancianarrow);

  *efl=0.;
}




void LoadParam_kbd()
{
  printf(" FITS broad band image: ");
  reads(p.broadfile,p.broadfile);
  printf(" FITS narrow image: ");
  reads(p.narrowfile,p.narrowfile);
  printf(" Input catalogue with detected objects in both images: ");
  reads(p.catfile,p.catfile);
  printf(" Column in %s with broad band image identification ",p.catfile);
  p.colnumb=readi(1);
  printf(" Column in %s with narrow band image identification ",p.catfile);
  p.colnumn=readi(23);

  printf(" Column in %s with broad band image X position ",p.catfile);
  p.colxpb=readi(14);
  printf(" Column in %s with broad band image Y position ",p.catfile);
  p.colypb=readi(15);
  printf(" Column in %s with narrow band image X position ",p.catfile);
  p.colxpn=readi(36);
  printf(" Column in %s with narrow band image Y position ",p.catfile);
  p.colypn=readi(37);
  printf(" Column in %s with broad band image object flux ",p.catfile);
  p.colfb=readi(4);
  printf(" Column in %s with broad band image object error flux ",p.catfile);
  p.colefb=readi(26);
  printf(" Column in %s with narrow band image object flux ",p.catfile);
  p.colfn=readi(4);
  printf(" Column in %s with narrow band image object error flux ",p.catfile);
  p.colefn=readi(26);
  printf(" Column in %s with A axis of the ellipse fit in broad band",p.catfile);
  p.coleab=readi(6);
  printf(" Column in %s with B axis of the ellipse fit in broad band",p.catfile);
  p.colebb=readi(7);
  printf(" Column in %s with angle of the ellipse fit in broad band",p.catfile);
  p.coletb=readi(8);
  printf(" Column in %s with A axis of the ellipse fit in narrow band",p.catfile);
  p.colean=readi(28);
  printf(" Column in %s with B axis of the ellipse fit in narrow band",p.catfile);
  p.colebn=readi(29);
  printf(" Column in %s with angle of the ellipse fit in narrow band",p.catfile);
  p.coletn=readi(30);
  printf(" Column in %s with background sky in broad band",p.catfile);
  p.colskyb=readi(12);
  printf(" Column in %s with background sky in narrow band",p.catfile);
  p.colskyn=readi(34);
  printf(" Fixed aperture photometry (0=use other flux");
  p.rfixap=readf(0);
}




void ChangePlotMain()
{
  char opt='A',cnul;

  while(opt!='E') {
    PlotMain();
    printf(" A Automatic limits\n");
    printf(" Z Zoom limits with cursor\n");
    printf(" L Set limits by keyboard\n");
    printf(" E Exit\n");                            
    
    printf(" Choose your option: ");
    opt=readc('Z');
/*     //printf(" No has ele\n"); */
    switch(opt) {
    case 'A': 
    case 'a':
      pm.limflag=1;
      break;
    case 'Z':
    case 'z':
      pm.limflag=0;
      printf(" Click on right image the bottom left square\n");
      cpgcurs(&pm.xmin,&pm.ymin,&cnul);
      printf(" Click on right image the upper right square\n");
      cpgsci(2);
      cpgband(2,1,pm.xmin,pm.ymin,&pm.xmax,&pm.ymax,&cnul);
      cpgsci(1);
      break;
    case 'S':
    case 's':
      printf(" Xmin: "); pm.xmin=readf(pm.xmin);
      printf(" Xmax: "); pm.xmax=readf(pm.xmax);
      printf(" Ymin: "); pm.ymin=readf(pm.ymin);
      printf(" Ymax: "); pm.ymax=readf(pm.ymax);
      break;
    }
  }

  
}


void ReadAstrom()
{
  
  if ((headbroad = fitsrhead (p.broadfile, &lhead, &nbfits)) == NULL) {
    fprintf (stderr, "Cannot read FITS header of file %s\n", p.broadfile);
    exit(1);
  }
  
  if((wcsbroad=wcsinit(headbroad))==NULL) {
    wcsflag=0;
    printf(" No WCS information found in header. Asuming pixel-to-pixel transformation.\n");  
    /*exit(1); */
    return;
  }
  else {
    printf(" WCS information from header:\n");
    PrintWCS(headbroad,1);
  }
  if ((headnarrow = fitsrhead (p.narrowfile, &lhead, &nbfits)) == NULL) {
    fprintf (stderr, "Cannot read FITS header of file %s\n", p.narrowfile);
    exit(1);
  }
  
  if((wcsnarrow=wcsinit(headnarrow))==NULL) {
    printf(" No WCS information found in header\n Exiting");  
    exit(1);
  }
  else {
    printf(" WCS information from header:\n");
    PrintWCS(headnarrow,1);
  }


}

void ComputeAperFot() 
{

/*   FILE *faperiraf; */
/*   char iraffile[51]; */
/*   //FILE *faperb,*fapern; */

/*   //Para el dibujo: */
  float tr[6];
/*   float meanb,sigmab; */
/*   float meann,sigman; */
  float bgb,fgb,bgn,fgn;
  float x1b,y1b,x2b,y2b;
  float x1n,y1n,x2n,y2n;
  double xpm,ypm;  /* //Posiciones en el fichero del catalogo  master */
  double xps,yps;  /* //Posiciones calculadas del catalogo  secundario (no master) */
  double ra,dec; /*  // Posicion astrometrica */
  int off;
  

  int npix;
/*   //float fb,fn; */
  int j;

  
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  
/*   //printf(" Input IRAF photometry file: "); */
/*   //scanf("%s",iraffile); */

/*   //faperiraf=fopen(iraffile,"r"); */

/*   //fapern=fopen("apernarrow.iraf","w"); */

  printf(" Phase 1: computing aperture photometry for the whole sample...\n");
  /* meanb=StMedia(nxb*nyb,imgbroad,&sigmab); */
/*   meann=StMedia(nxn*nyn,imgnarrow,&sigman); */
/*   bgb=meanb-sigmab*4.5; */
/*   fgb=meanb+sigmab*2.5; */
/*   bgn=meann-sigman*4.5; */
/*   fgn=meann+sigman*2.5; */

  bgb=4000;fgb=6500;
  bgn=90;fgn=300;
  
/*   //printf(" Using cuts: broad %f %f narrow %f %f\n",bgb,fgb,bgn,fgn);  */

  if(p.verbose)  cpgopen("?");
/*   //printf("fixap %f\n",p.rfixap); */
  printf(" Object 0000/4%d",nobj);
/*   //De momento cojo la broad como master file */
  for(j=0;j<nobj;j++) {
    if(true[j]) {
      printf("\b\b\b\b\b\b\b\b\b%4d/%4d",j,nobj);
/*       //printf(" Object %d in BROAD: %f,%f NARROW: %f,%f\n",j,xpb[j],ypb[j],xpn[j],ypn[j]); */
      fapb[j]=circ_aper(imgbroad,nxb,nyb,xpb[j]-1,ypb[j]-1,p.rfixap,&npix);
      fapb[j]-=npix*skyb[j];
      xpm=(double)xpb[j];
      ypm=(double)ypb[j];
      off=0;
/*       //printf(" x %f y %f\n",xpm,ypm); */
      pix2wcs(wcsbroad,xpm-1, ypm-1 ,&ra, &dec);
/*       //printf(" ASTROMETRY %f %f \n",ra,dec); */
      wcs2pix(wcsnarrow,ra,dec,&xps,&yps,&off);
/*       //printf(" Calculated coordinates in NARROW %f,%f\n",xps,yps); */
      xps=xps+1;yps=yps+1;
/*       //El -1 es por la cagada de las columnas. Hay que arreglarlo!! */

      if(!off) {
	fapn[j]=circ_aper(imgnarrow,nxn,nyn,(float)xps-1,(float)yps-1,p.rfixap,&npix);
	fapn[j]-=npix*skyn[j];
/* 	//printf(" Apertura %f  ",fapn[j]); */
/* 	//printf(" La otra %f\n",fn[j]); */
	  
      }
      else {
	true[j]=0;
/* 	//printf(" Object out of secondary image bounds\n"); */
      }

/*       //      fscanf(faperiraf," %f %f",&fapb[j],&fapn[j]);  */
/*       //printf("IRAF  %f %f\n",fapb[j],fapn[j]);  */

/*       //Para que no salgan tonterias: */
      if(fapb[j]<=0 || fapn[j]<=0) true[j]=0;
      if(xpn[j]<40 || xpn[j]>(nxn-40) || ypn[j]<40 || ypn[j]>(nyn-40)) true[j]=0;

/*       //De prueba para fotometria IRAF */
/*       //fprintf(faperb,"  %f %f\n",xpm,ypm); */
/*       //fprintf(fapern,"  %f %f\n",xps,yps); */
      

      if(p.verbose) {
	x1b=xpb[j]-p.rfixap*2.5;
	x2b=xpb[j]+p.rfixap*2.5;
	y1b=ypb[j]-p.rfixap*2.5;
	y2b=ypb[j]+p.rfixap*2.5;
	cpgsvp(0.05,0.45,0.05,0.95);
	
	cpgwnad(x1b,x2b,y1b,y2b);
	cpggray(imgbroad,nxb,nyb,(int)x1b,(int)x2b,(int)y1b,(int)y2b,fgb,bgb,tr);
	cpgbox("BCTNS",0,0,"BCTNS",0,0);
	cpglab("X pixel","Y pixel","Broad band image");
	cpgsci(2);
	cpgelip(xpb[j],ypb[j],p.rfixap,p.rfixap,0.); 
	cpgsci(1);
	
	printf(" Broad band aperture photometry  %f  ",fapb[j]);
	printf(" Narrow band aperture photometry  %f  ",fapn[j]);
	x1n=xps-p.rfixap*2.5;
	x2n=xps+p.rfixap*2.5;
	y1n=yps-p.rfixap*2.5;
	y2n=yps+p.rfixap*2.5;
	cpgsvp(0.55,0.95,0.05,0.95);
	cpgwnad(x1n,x2n,y1n,y2n);
	cpggray(imgnarrow,nxn,nyn,(int)x1n,(int)x2n,(int)y1n,(int)y2n,fgn,bgn,tr);
	cpgbox("BCTNS",0,0,"BCTNS",0,0);
	cpglab("X pixel","Y pixel","Narrow band image");
	cpgsci(2);
	cpgelip((float)xps,(float)yps,p.rfixap,p.rfixap,0.); 
	cpgsci(1);
	cpgsch(4.);
	cpgsci(2);
	cpgpt1(xps,yps,6);
	cpgsci(4);
	cpgpt1(xpn[j],ypn[j],6);
	cpgsci(1);
	cpgsch(1.);

	
	
      }



    }
/*     //cpgpage(); */
    if(p.verbose) cpgeras();
  }
  if (p.verbose) cpgclos();
/*   //exit(1); */
 
/*   //fclose(faperb);fclose(fapern); */
}

void SaveSelected(char selfile[],char filepar[])
{
  FILE *fs;
  int j;
  char sanul1[32],sanul2[32];
  char sdnul1[32],sdnul2[32];

  if((fs=fopen(selfile,"w")) == NULL) {
    fprintf(stderr,"NarrowSearch: ERROR. No such parameter file %s\n", selfile);
    exit(1);
  }

  printf(" Number of objects to save: %d\n",nsel);
  fprintf(fs,"#Selected objects for images %s %s\n# NarrowSearch called with parameter file %s\n",p.broadfile,p.narrowfile,filepar);
/*   Se prodria poner tanto como esto: */
/*   fprintf(fs,"#Num          RA_broad        DEC_narrow       Broad_mag    EW        */
/*              X_broad    Y_broad   X_narrow     Y_narrow   RA_narrow   DEC_narrow */
/*              Broad_counts Err_broad_counts  Narrow_counts Err_narrow_counts */
/*              Z  Err_Z   Nsig   Flux_line  Err_flux_line\n"); */
  if(aperflag) {
    fprintf(fs,"#Num_broad RA_broad     DEC_broad  Mag_iso_broad Flux_iso_broad");
    fprintf(fs," Mag_aper_broad Flux_aper_broad  Xpix_broad Ypix_broad");
    fprintf(fs,"    RA_narrow   DEC_narrow   Mag_iso_narrow Flux_iso_narrow     EW_iso     Flux_line_iso");
    fprintf(fs," Flux_aper_narrow  EW_aper Flux_line_aper  X_pix_narrow  Y_pix_narrow");
    fprintf(fs,"    Z      Nsig   Prio\n");
  }
  else {
    fprintf(fs,"#Num_broad RA_broad     DEC_broad      Mag_broad    Flux_broad  ");
    fprintf(fs,"  Xpix_broad  Ypix_broad");
    fprintf(fs,"      RA_narrow   DEC_narrow       Mag_narrow     Flux_narrow        EW         Flux_line  ");
    fprintf(fs," X_pix_narrow  Y_pix_narrow");
    fprintf(fs,"    Z      Nsig   Prio\n");
  }
  for(j=0;j<nsel;j++) {
    ra2str(sanul1,32,objs[j].rab,3);
    dec2str(sdnul1,32,objs[j].decb,2);
    ra2str(sanul2,32,objs[j].ran,3);
    dec2str(sdnul2,32,objs[j].decn,2);
    if(aperflag) { 
      fprintf(fs," %6d  %s %s       %7.2f    %9.2f  ",objs[j].idx,sanul1,sdnul1,objs[j].magbroad,objs[j].fluxbroad);
      fprintf(fs,"     %7.2f     %9.2f     %8.2f    %8.2f   ",objs[j].magaperbroad,objs[j].fluxaperbroad,objs[j].xbroad,objs[j].ybroad);
      fprintf(fs,"  %s  %s        %7.2f    %9.2f       %6.2f       %9.2f  ",sanul2,sdnul2,objs[j].magnarrow,objs[j].fluxnarrow,objs[j].ewwhole,objs[j].fluxline);
      fprintf(fs,"   %9.2f       %6.2f      %9.2f     %8.2f    %8.2f   ",objs[j].fluxapernarrow,objs[j].ewaper,objs[j].fluxaperline,objs[j].xnarrow,objs[j].ynarrow);
      fprintf(fs,"   %5.2f %8.2g  %4f  %5.3f\n",objs[j].z,objs[j].nsig,objs[j].fluxline,objs[j].pflag);
    }
    else {
      fprintf(fs," %6d  %s %s       %7.2f    %9.2f  ",objs[j].idx,sanul1,sdnul1,objs[j].magbroad,objs[j].fluxbroad);
      fprintf(fs,"     %8.2f    %8.2f   ",objs[j].xbroad,objs[j].ybroad);
      fprintf(fs,"  %s  %s        %7.2f    %9.2f       %6.2f       %9.2f  ",sanul2,sdnul2,objs[j].magnarrow,objs[j].fluxnarrow,objs[j].ewwhole,objs[j].fluxline);
      fprintf(fs,"   %8.2f    %8.2f   ",objs[j].xnarrow,objs[j].ynarrow);
      fprintf(fs,"   %5.2f %8.2g  %4f  %5.3f\n",objs[j].z,objs[j].nsig,objs[j].fluxline,objs[j].pflag);
    }
  }
  fclose(fs);
}

void SelectObjThres(struct th thrhol)
{
/*   int j; */
}
void PlotThress(){
  int ib;
  float *sigd,*sigu;
  
  if((sigd=malloc(100*sizeof(float)))==NULL) printf("I cannot dimension sigd   of %d elements \n",100);
  if((sigu=malloc(100*sizeof(float)))==NULL) printf("I cannot dimension sigu   of %d elements \n",100);
  
  for(ib=0;ib<thress.nbin;ib++) {
    sigu[ib]=thress.ym[ib]+thress.nsig*thress.sig[ib];
    sigd[ib]=thress.ym[ib]-thress.nsig*thress.sig[ib];  
/*     //printf(" x %f y %f sig %f\n",thress.x[ib],thress.ym[ib],thress.sig[ib]); */
  }
/*   //printf(" nbin %d\n",thress.nbin); */
  
  cpgsci(2);
  cpgline(thress.nbin,thress.x,thress.ym);
  cpgsci(4);
  cpgline(thress.nbin,thress.x,sigu);
  cpgline(thress.nbin,thress.x,sigd);
  cpgsci(1);
  free(sigu);free(sigd); 
}
void ReadMagCal() { 
  printf(" Magnitude Calibration--------------\n");
  printf(" Broad image: mag=-2.5log(cts/s)+A+Bsec(z)\n");
  printf(" Input A: ");
  magcal.Ab=readf(0);
  printf(" Input B: ");
  magcal.Bb=readf(0);
  printf(" Input sec z: ");
  magcal.seczb=readf(1);
  printf(" Input exposure time : ");
  magcal.texpbroad=readf(1);
  printf(" Narrow image: mag=-2.5log(cts/s)+A+Bsec(z)\n");
  printf(" Input A: ");
  magcal.An=readf(0);
  printf(" Input B: ");
  magcal.Bn=readf(0);
  printf(" Input sec z: ");
  magcal.seczn=readf(1);
  printf(" Input exposure time : ");
  magcal.texpnarrow=readf(1);
}


void Readfilter( char file[100],struct filter F)
{
  FILE *fp;
  int i;       
  if((fp=fopen(file,"r"))==NULL) {
    printf("Cannot open file %s\n",file);
    exit(1);
  }
  F.n=FileNLin(file);
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
  }  
}
