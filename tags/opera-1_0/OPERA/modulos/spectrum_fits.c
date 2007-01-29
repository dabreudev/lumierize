#include "modulos.h" 
#define DEBUG 0
#define NBUF 5


void ReadSpec(struct spectrum *sp) {
 
  int status=0; 
  int nfound, anynull;
  long fpixel;
  int nullval;
  long naxes;
  char comment[51];
  int i;
  int cal_flag;
  float cdelt1, crval1,crpix1;
  
  if(DEBUG) printf(" Entra a readspe\n");

  if(ffopen(&((*sp).filefits), (*sp).file, READONLY, &status)) fits_report_error(stderr,status);  
  if(fits_read_key((*sp).filefits,TLONG,"NAXIS", &naxes, comment, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng((*sp).filefits, "NAXIS", 1, 2, (*sp).naxes, &nfound, &status)) fits_report_error(stderr,status);
  if(status) {
    fits_report_error(stderr,status);
    exit(status);
  }
  ffgky((*sp).filefits,TFLOAT,"CDELT1",&cdelt1,comment,&status);
  if(status) {
    cal_flag=0;
    status=0;
  }
  else {
    ffgky((*sp).filefits,TFLOAT,"CRVAL1",&crval1,comment,&status);
    ffgky((*sp).filefits,TFLOAT,"CRPIX1",&crpix1,comment,&status);
    cal_flag=1;
    if(status) {
      fits_report_error(stderr,status);
      printf(" Something went wrong with wavelenght calibration\n");
      exit(1);
    }
  }
  
  if(naxes==1 || (naxes==2 && (*sp).naxes[1]==1))  {
    (*sp).npixels=(*sp).naxes[0];
    (*sp).naxes[1]=1;
  }
  else {
    if(DEBUG) printf(" naxes %ld na1 %ld na2 %ld\n",naxes,(*sp).naxes[0],(*sp).naxes[1]); 
    printf(" FITS file %s is not a 1-d image\n",(*sp).file);
    exit(1);
  }
  if(DEBUG) printf(" Aqui\n"); 
  if((*sp).aloc_flag) free((*sp).spec);
  if((*sp).alocldo_flag) free((*sp).ldo);
  if(DEBUG) printf(" Aqui mas\n"); 
  (*sp).spec=vector_f((*sp).npixels);
  (*sp).ldo=vector_f((*sp).npixels);
  (*sp).aloc_flag=1;
  (*sp).alocldo_flag=1;
  fpixel=1;
  nullval=0;
  (*sp).datamin=1.0e30;
  (*sp).datamax=-1.0e30;
  (*sp).nx=(*sp).naxes[0];
  if(cal_flag) {
    (*sp).ldomin=crval1-(crpix1-1.)*cdelt1;
    (*sp).deltaldo=cdelt1;
  }
  else {
    (*sp).ldomin=1;
    (*sp).deltaldo=1;
  }
  printf("...Reading spectrum %s \n",(*sp).file);
  if(fits_read_img((*sp).filefits, TFLOAT, fpixel, (*sp).npixels, &nullval, (*sp).spec, &anynull, &status )) fits_report_error(stderr,status);
  for(i=0;i<(*sp).nx;i++)   (*sp).ldo[i]=(*sp).ldomin+i*(*sp).deltaldo;
  MinMax((*sp).nx,(*sp).spec,&((*sp).datamin),&((*sp).datamax));  
  ffclos((*sp).filefits,&status);
  if(DEBUG) printf(" Finalizo leer\n");
  if(status) {
    fits_report_error(stderr,status);
    exit(status);
  }
}

void ReadSpec_buf(struct spectrum *sp) {

  static int nbuffer=0;
  static struct image imabuf[NBUF];
  char *ptr1;  
  char url[1000];
  char *imagefile;
  char *restname;
  char *number;
 
/*   int status=0;   */
/*   int nfound, anynull; */
/*   long fpixel; */
/*   int nullval; */
/*   long naxes; */
/*   char comment[51]; */
  int i;
  int cal_flag;
/*   float cdelt1, crval1,crpix1; */
  
  int j,ibuffer;
  static int ifree;
  int iline;

  strcpy(url,(*sp).file);
  ptr1=strchr(url, '[');
  if(!ptr1) {
    if(DEBUG) printf(" No hay [\n");
    ReadSpec(sp);
  }
  else {
    imagefile = strtok(url, "[");
    if(DEBUG) printf(" HHHHHHHAAAAYYYY hay iamagefile %s\n",imagefile);
    
    ibuffer=-1;
    for(j=0;j<NBUF;j++) {
      if(!strcmp(imagefile,imabuf[j].file))  ibuffer=j;
    }
    if(DEBUG) printf(" ibuffer %d\n",ibuffer);
    if(ibuffer==-1) {
      if(DEBUG) printf(" No image buffered\n");
      if(nbuffer==0) {
	imabuf[0].aloc_flag=0;
	ibuffer=0;
	nbuffer=1;
	ifree=1;
      }
      else if(nbuffer<NBUF) {
	if(DEBUG) printf(" Llevo %d buffers\n",nbuffer);
	nbuffer++;
	imabuf[nbuffer].aloc_flag=0;
	ibuffer=nbuffer-1;
	ifree=nbuffer;
	if(ifree==NBUF) ifree=0;
      }
      else if(nbuffer==NBUF) {
	if(DEBUG) printf(" He llegado all final del buffer\n");
	ibuffer=ifree;
	ifree++;
	if(ifree==NBUF) ifree=0;
      }
      if(DEBUG) printf(" Al final lo leo en el buffer %d\n",ibuffer); 
      strcpy(imabuf[ibuffer].file,imagefile);
      ReadImage(imabuf+ibuffer);
    }
    strcpy(url,(*sp).file);    
    restname=strstr(url,"[");
    if(DEBUG) printf(" restname %s u %s\n",restname,url);
    if(restname[1]!='*') {
      printf(" I cannot handle %s as a valid parse for a spectrum inside %s\n",restname,imagefile);
      exit(1);
    }
    else {
      number=strstr(restname,",")+1;
      if(DEBUG)  printf(" number %s\n",number);
      number=strtok(number,":");
      if(DEBUG)  printf(" number2 %s\n",number);
      iline=atoi(number);
    }
    if(DEBUG) printf(" iline %d\n",iline);

    /* Ahora que se que lo tengo en ibuffer lo copio en el espectro */
    (*sp).npixels=(imabuf[ibuffer]).naxes[0];
    (*sp).naxes[0]=(imabuf[ibuffer]).naxes[0];
    (*sp).naxes[1]=1;
    if((*sp).aloc_flag) free((*sp).spec);
    if((*sp).alocldo_flag) free((*sp).ldo);
    (*sp).spec=vector_f((*sp).npixels);
    (*sp).ldo=vector_f((*sp).npixels);
    (*sp).aloc_flag=1;
    (*sp).alocldo_flag=1;
    (*sp).datamin=1.0e30;
    (*sp).datamax=-1.0e30;
    (*sp).nx=(*sp).naxes[0];
    cal_flag=0;
    (*sp).ldomin=1;
    (*sp).deltaldo=1;
    for(i=0;i<(*sp).nx;i++)   (*sp).spec[i]=(imabuf[ibuffer]).array[i+(iline-1)*imabuf[ibuffer].nx];

    for(i=0;i<(*sp).nx;i++)   (*sp).ldo[i]=(*sp).ldomin+i*(*sp).deltaldo;
    MinMax((*sp).nx,(*sp).spec,&((*sp).datamin),&((*sp).datamax));     
    if(DEBUG) printf(" He terminado\n");
  }
}



void SaveSpec(struct spectrum sp) {
  
  int status=0;
  float crpix,crval,cdelt;
  

  sp.naxes[1]=1;
  crpix=1;
  crval=sp.ldomin;
  cdelt=sp.deltaldo;
  /* Salvando El espectro */
  printf(" Saving spectrum  %s \n",sp.file);
  if(DEBUG) printf(" NX %d %ld\n",sp.nx,sp.naxes[0]);
  if(ffinit(&(sp.filefits),sp.file,&status)) fits_report_error(stderr,status); 
  fits_create_img(sp.filefits, -32,1,sp.naxes,&status);
  if(status) fits_report_error(stderr,status);
  fits_write_img(sp.filefits,TFLOAT,1,sp.npixels,sp.spec,&status);
  if(status) fits_report_error(stderr,status); 
  fits_write_key(sp.filefits,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
  if(status) fits_report_error(stderr,status); 
  fits_write_key(sp.filefits,TFLOAT,"CRVAL1",&crval,"Coordinate at reference pixel",&status);
  if(status) fits_report_error(stderr,status); 
  fits_write_key(sp.filefits,TFLOAT,"CDELT1",&cdelt,"Coordinate increment per pixel",&status);
  if(status) fits_report_error(stderr,status); 
  fits_close_file(sp.filefits,&status);
  if(status) {
    fits_report_error(stderr,status);
    exit(status);
  }
}


void FillSpec(struct spectrum *sp, int npix, float *spec, float ldomin, float ldomax) {
  
  int i;
  int cal_flag=1;
  
  sprintf((*sp).file,"INDEF");
  (*sp).naxes[0]=npix;
  (*sp).naxes[1]=1;
  (*sp).npixels=(*sp).naxes[0];
  
  if(DEBUG) printf(" antes aki\n");
  
  if((*sp).aloc_flag) free((*sp).spec);
  if((*sp).alocldo_flag) free((*sp).ldo);
  (*sp).spec=vector_f((*sp).npixels);
  (*sp).ldo=vector_f((*sp).npixels);
  (*sp).aloc_flag=1;
  (*sp).alocldo_flag=1;
  (*sp).datamin=1.0e30;
  (*sp).datamax=-1.0e30;
  (*sp).nx=(*sp).naxes[0];
  if(ldomin==ldomax) cal_flag=0;
  if(cal_flag) {
    (*sp).ldomin=ldomin;
    (*sp).deltaldo=(ldomax-ldomin)/npix;
  }
  else {
    (*sp).ldomin=1;
    (*sp).deltaldo=1;
  }
  for(i=0;i<(*sp).nx;i++)  {
    (*sp).spec[i]=spec[i];
    (*sp).ldo[i]=(*sp).ldomin+i*(*sp).deltaldo;
  }
  
}

void CopySpec(struct spectrum *spdest, struct spectrum sporig) {
  
  int i;
/*   int cal_flag=1; */
  if(DEBUG) printf(" Entro en copy\n");
  
  sprintf((*spdest).file,"INDEF");
  (*spdest).naxes[0]=sporig.naxes[0];
  (*spdest).naxes[1]=sporig.naxes[1];
  (*spdest).npixels=sporig.npixels;
  
  if(DEBUG) printf(" antes aki\n");
  
  if((*spdest).aloc_flag) {
    if(DEBUG) printf(" LIBERO !!! %d\n",(*spdest).spec);
    free((*spdest).spec);
  }
  if((*spdest).alocldo_flag) {
    if(DEBUG) printf(" LIBERO !!! %d\n",(*spdest).ldo);
    free((*spdest).ldo);
  }
  (*spdest).spec=vector_f((*spdest).npixels);
  (*spdest).ldo=vector_f((*spdest).npixels);

  if(DEBUG) printf(" alocateos spec %d ldo %d\n",(*spdest).spec,(*spdest).ldo);

  (*spdest).aloc_flag=1;
  (*spdest).alocldo_flag=1;
  (*spdest).nx=(*spdest).naxes[0]; 
  (*spdest).ldomin=sporig.ldomin; 
  (*spdest).deltaldo=sporig.deltaldo; 
  (*spdest).datamin=sporig.datamin;  
  (*spdest).datamax=sporig.datamax;  
  if(DEBUG) printf(" Antes for nx %d nouxeks %ld\n",(*spdest).nx,(*spdest).npixels);
  if(DEBUG) printf(" nace[0] %ld orignax0 %ld\n",(*spdest).naxes[0],sporig.naxes[0]);
  for(i=0;i<(*spdest).nx;i++)  {
/*     printf(" %d  %f %f \n",i,sporig.spec[i],sporig.ldo[i]);  */
    (*spdest).spec[i]=sporig.spec[i];
    (*spdest).ldo[i]=sporig.ldo[i];
  }
  if(DEBUG) printf(" Al final\n");

  
}

void CloseSpec(struct spectrum *sp) {
  if((*sp).aloc_flag) {
    if(DEBUG) printf(" Si LIBERO %d\n",(*sp).spec);
    free((*sp).spec);
  }
  if((*sp).alocldo_flag) {
    if(DEBUG) printf(" Si LIBERO %d\n",(*sp).ldo);
    free((*sp).ldo);
  }
}


void PlotSpec(struct spectrum sp) {
  
  int i;
  
  if(sp.alocldo_flag!=1) {
    sp.ldo=vector_f(sp.nx);
    for(i=0;i<sp.nx;i++) {
      sp.ldo[i]=sp.ldomin+i*sp.deltaldo;
    }
  }
  pgLimits(sp.nx,sp.spec,&sp.datamin,&sp.datamax);
  if(DEBUG) printf(" min %f max %f\n",sp.datamin,sp.datamax);
  
  cpgask(0);
  cpgeras();
  cpgsch(1.5);
  cpgslw(2.);
  cpgvstd();
  cpgsch(1.);
  cpgswin(sp.ldomin,sp.ldomin+sp.nx*sp.deltaldo,sp.datamin-(sp.datamax-sp.datamin)*.2,sp.datamax+(sp.datamax-sp.datamin)*.45);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  cpgbin(sp.nx,sp.ldo,sp.spec,1);
  cpgscf(2);
  cpgsch(1.5); 
  cpglab("","",sp.file); 
  cpgscf(1); 

}


void PlotSpec_zoom(struct spectrum sp, float x1, float x2, float y1, float y2) {
  
  int i;
  
  if(sp.alocldo_flag!=1) {
    sp.ldo=vector_f(sp.nx);
    for(i=0;i<sp.nx;i++) {
      sp.ldo[i]=sp.ldomin+i*sp.deltaldo;
    }
  }
  
  cpgask(0);
  cpgeras();
  cpgsch(1.5);
  cpgslw(2.);
  cpgvstd();
  cpgsch(1.);
  cpgswin(x1,x2,y1,y2);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  cpgbin(sp.nx,sp.ldo,sp.spec,1);
  cpgscf(2);
  cpgsch(1.5); 
  cpglab("","",sp.file); 
  cpgscf(1); 

}

void PlotSpec_err(struct spectrum sp, struct spectrum errsp) {

  float *ldo; 
  float *upspec; 
  float *downspec; 
  int i; 
  float fnul; 
  float datamin,datamax; 
  char tittle[202]; 

  ldo     =vector_f(sp.nx); 
  upspec  =vector_f(sp.nx); 
  downspec=vector_f(sp.nx); 
  for(i=0;i<sp.nx;i++) { 
    ldo[i]=sp.ldomin+i*sp.deltaldo; 
    upspec[i]=sp.spec[i]+errsp.spec[i]; 
    downspec[i]=sp.spec[i]-errsp.spec[i]; 
  } 
  pgLimits(sp.nx,upspec  ,&fnul      ,&datamax); 
  pgLimits(sp.nx,downspec,&datamin,&fnul      );
  printf(" min %f max %f\n",datamin,datamax);
  
  cpgask(0);
  cpgeras();
  cpgvstd();
  cpgswin(sp.ldomin,sp.ldomin+sp.nx*sp.deltaldo,datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.45);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  cpgsci(1);
  cpgbin(sp.nx,ldo,sp.spec,1);
  cpgsci(2);
  cpgbin(sp.nx,ldo,upspec  ,1);
  cpgbin(sp.nx,ldo,downspec,1);
  sprintf(tittle," %s & %s",sp.file,errsp.file);
  cpgsci(1);
  cpglab("","",tittle);
  free(ldo);
  free(upspec);
  free(downspec);
}

void PlotSpec_err_zoom(struct spectrum sp, struct spectrum errsp, float x1, float x2, float y1, float y2) {

  float *ldo; 
  float *upspec; 
  float *downspec; 
  int i; 
/*   float fnul;  */
/*   float datamin,datamax;  */
  char tittle[202]; 

  ldo     =vector_f(sp.nx); 
  upspec  =vector_f(sp.nx); 
  downspec=vector_f(sp.nx); 
  for(i=0;i<sp.nx;i++) { 
    ldo[i]=sp.ldomin+i*sp.deltaldo; 
    upspec[i]=sp.spec[i]+errsp.spec[i]; 
    downspec[i]=sp.spec[i]-errsp.spec[i]; 
  } 
  
  cpgask(0);
  cpgeras();
  cpgvstd();
  cpgswin(x1,x2,y1,y2);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  cpgsci(1);
  cpgbin(sp.nx,ldo,sp.spec,1);
  cpgsci(2);
  cpgbin(sp.nx,ldo,upspec  ,1);
  cpgbin(sp.nx,ldo,downspec,1);
  sprintf(tittle," %s & %s",sp.file,errsp.file);
  cpgsci(1);
  cpglab("","",tittle);
  free(ldo);
  free(upspec);
  free(downspec);
}



void PlotSpec_pix(struct spectrum sp) {

  float *ldo;
  float datamin,datamax; 
  int i;

  if(DEBUG) printf(" El nx %d \n",sp.nx);

  ldo=vector_f(sp.nx);
  for(i=0;i<sp.nx;i++) {
    ldo[i]=1.+i;
  }
  
  if(DEBUG) printf(" Antes Minmax\n");

  pgLimits(sp.nx,sp.spec,&datamin,&datamax);
  printf(" min %f max %f\n",datamin,datamax);
  
  cpgask(0);
  cpgeras();
  cpgvstd();
  cpgswin(ldo[0],ldo[sp.nx-1],datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.45);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  cpgbin(sp.nx,ldo,sp.spec,1);
  cpglab("","",sp.file);
  free(ldo);
   

}



void PlotSpec_ov(struct spectrum sp) {

  float *ldo;
  int i;
  static int color=2;


  if(color>8) color=2;
  ldo=vector_f(sp.nx);
  for(i=0;i<sp.nx;i++) {
    ldo[i]=sp.ldomin+i*sp.deltaldo;
  }

  cpgsci(color);
  cpgbin(sp.nx,ldo,sp.spec,1);
  free(ldo);
  color++;
  cpgsci(1);
  cpglab("","",sp.file);
}

void PlotSpec_pix_ov(struct spectrum sp) {

  float *ldo;
  int i;
  static int color=2;


  if(color>8) color=2;
  ldo=vector_f(sp.nx);
  for(i=0;i<sp.nx;i++) {
    ldo[i]=1.+i;
  }

  cpgsci(color);
  cpgbin(sp.nx,ldo,sp.spec,1);
  free(ldo);
  color++;
  cpgsci(1);
  cpglab("","",sp.file);
}


void PlotSpec_pix_err(struct spectrum sp, struct spectrum errsp) {

  float *ldo;
  float *upspec;
  float *downspec;
  int i;
  float fnul;
  float datamin,datamax;
  char tittle[202];

  ldo     =vector_f(sp.nx);
  upspec  =vector_f(sp.nx);
  downspec=vector_f(sp.nx);
  for(i=0;i<sp.nx;i++) {
    ldo[i]=1.+i;
    upspec[i]=sp.spec[i]+errsp.spec[i];
    downspec[i]=sp.spec[i]-errsp.spec[i];
  }
  pgLimits(sp.nx,upspec  ,&fnul      ,&datamax);
  pgLimits(sp.nx,downspec,&datamin,&fnul      );
  if(DEBUG)   printf(" min %f max %f\n",datamin,datamax);

  cpgask(0);
  cpgeras(); 
  cpgvstd();
  cpgswin(ldo[0],ldo[sp.nx-1],datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.45);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  cpgsci(1);
  cpgbin(sp.nx,ldo,sp.spec,1);
  cpgsci(2);
  cpgbin(sp.nx,ldo,upspec  ,1);
  cpgbin(sp.nx,ldo,downspec,1);
  cpgsci(1);
  sprintf(tittle," %s & %s",sp.file,errsp.file);
  cpglab("","",tittle);
  free(ldo);
  free(upspec);
  free(downspec);
}


