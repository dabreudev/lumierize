#include "modulos.h"
 
struct coadd_info {
  
  

};

struct specitem {
  
  char rawspec[51];
  char errrawspec[51];
  float xp;
  float yp;
  double ra;
  double dec;
  float sky;
};

struct specfile {
  struct specitem *si;
  int nitems;
  char filename[51];
};


char dbfile[100];
char catfile[100];
int  colcatra;
int  colcatdec;
char pgdevice[100];
int  maxnbuffer;

struct SurveyDB sdb;


void LoadParam_kbd(void);

void ParseObjFile(void);
void FindFITS(double ra, double dec, struct SurveyDB surdb, char ***imalist, int *nima);
int  whithinimage(double ra, double dec, struct SurveyItem sitem);
void SpecExtractFromImages(char **imalist, int nlist, struct spectrum *spectra, int *nspec);
void SpecExtractFromList(double ra, double dec, char **specfilelist, int nlist, struct spectrum **spectra, struct spectrum **errspectra, int *nspec);
void ReadSpecFile(char *sfile, struct specfile  *sf);
void CoaddSpec(struct spectrum *spectra, int nspec, struct spectrum sumspec, struct coadd_info *ci);

double minequal(int n,double *val,int *nlog,int *eof);



int main()
{

 
  struct catalogue {
    FILE *fp;
    char catname[101];
    double *num;
    int    *log;
    int    nreg;
/*     //int    eof; */
  };
  struct catbody {
    double num;
    char   reg[1000];
  };
  int i,j;
  int ncat;
  int neof;
  int *irecord;
  int *eof;
  int *nlog;
  double *val;
  struct catalogue *cat;
  struct catbody   *bod;
  FILE *catout;
  char outfile[101];
/*   char  nul[1000]; */
/*   char word[1000]; */
/*   char c; */
  double min;


  int nobs1;
  int ccd1;
  int qnx1;
  int drift1;
  int nobs2;
/*   int ccd2; */
/*   int qnx2; */
/*   int drift2; */
  int sondos;

  char comando_iraf[10000];
  char specfile1[101],specfile2[101],specout[101];


  LoadParam_kbd();
  printf(" ANtes read\n");
  ReadSurDB(dbfile, &sdb);
  printf(" Despeus\n");

  ParseObjFile();
  
  
  

  ncat=2;
  cat= malloc(ncat*sizeof(struct catalogue));
  val= malloc(ncat*sizeof(float));
  nlog=malloc(ncat*sizeof(float));
  eof= malloc(ncat*sizeof(float));
  bod=malloc(ncat*sizeof(struct catbody));
  irecord=malloc(ncat*sizeof(int));

  printf(" FILE 1:\n");
  printf(" Input Number of observation: ");
  nobs1=readi(0);
  printf(" Input CCD: ");
  ccd1=readi(0);
  printf(" Input QNX: ");
  qnx1=readi(0);
  printf(" Input Drift: ");
  drift1=readi(0);
  printf(" FILE 2:\n");
  printf(" Input Number of observation: ");
  nobs2=readi(0);
/*   //exit(1); */
/*   printf(" Input CCD: "); */
/*   ccd2=readi(0); */
/*   printf(" Input QNX: "); */
/*   qnx2=readi(0); */
/*   printf(" Input Drift: "); */
/*   drift2=readi(0); */
  strcpy(cat[9].catname,"\0");
  sprintf(cat[0].catname,"%3d.%1d.%1d.%d.exusso",nobs1,qnx1,ccd1,drift1);
  printf(" ASN \n");
  sprintf(cat[1].catname,"%3d.%1d.%1d.%d.exusso",nobs2,qnx1,ccd1,drift1);
  sprintf(outfile,"%3d-%3d.%1d.%1d.%d.exusso",nobs1,nobs2,qnx1,ccd1,drift1);
  printf(" Aqui\n");
  
/*   //printf(" AUI\n"); */
/*   strcpy(cat[0].catname,"503.3.1.2.exusso"); */
/*   strcpy(cat[1].catname,"505.3.1.2.exusso"); */
  catout=fopen(outfile,"w");
/*   //printf(" AUSIA S\n"); */
  for(i=0;i<ncat;i++) {
    cat[i].fp=fopen(cat[i].catname,"r");
    cat[i].nreg=FileNLin(cat[i].catname);
    cat[i].num=malloc(cat[i].nreg*sizeof(double));
    cat[i].log=malloc(cat[i].nreg*sizeof(int));
    ReadDoublecol((cat[i]).catname,1,(cat[i]).num,(cat[i]).log,&((cat[i]).nreg));
/*     //cat[i].eof=0; */
    eof[i]=0;
    irecord[i]=0;
  }
  for(i=0;i<ncat;i++) {
    for(j=0;j<cat[i].nreg;j++) {
      printf(" %12.8f\n",cat[i].num[j]);
    }
    printf(" \n");
  }


/*   //printf(" SALLS\n"); */

/*   //printf(" KAJSD\n"); */
  neof=0;
  while(neof<ncat) {
    for(i=0;i<ncat;i++) {
/*       //printf(" KJH\n"); */
/*       fgets((bod[i]).reg,1000,cat[i].fp);  */
/*       //printf(" LEEE\n"); */
/*       c=nul[0];  */
/*       if(c=='#') { */
/* 	continue; */
/*       } */
/*       LeeWord(bod[i].reg,1,word); */
/*       val[i]=(double)atof(word); */
      val[i]=cat[i].num[irecord[i]];
/*       //printf("%d val %d %12.8f %d\n",i,irecord[i],cat[i].num[irecord[i]],cat[i].log[irecord[i]]); */
    }
/*     //printf(" ASDOSA\n"); */

    min=minequal(ncat,val,nlog,eof);
/*     //printf(" KLASJDL\n"); */
    fprintf(catout," %12.8f\n",min);
    sondos=0;
    for(i=0;i<ncat;i++) {
/*       //exit(1);       */

      if(!eof[i]) {
	if(nlog[i]) {
	  sondos++;


	  printf(" Numero %12.8f en fichero %d:  %s\n",val[i],i,bod[i].reg);
/* 	  //fprintf(catout," %s \n",bod[i].reg); */

/* 	  //printf(" irecord %d\n",irecord[i]); */
	  irecord[i]++;
	  while(!(cat[i].log[irecord[i]]))  {
/* 	    //printf(" CACUI\n"); */
	    irecord[i]++;
	  }
	  if(irecord[i]>=cat[i].nreg) {
	    neof++;
	    eof[i]=1;
	  }
/* 	  //printf(" asds\n"); */
	}
      }
      if(sondos==2) {
	sprintf(specfile1,"SP/%3d.%1d.%1d.%1d.wcsp%12.8f.fits",nobs1,qnx1,ccd1,drift1,min);
	sprintf(specfile2,"SP/%3d.%1d.%1d.%1d.wcsp%12.8f.fits",nobs2,qnx1,ccd1,drift1,min);
	sprintf(specout,"SC/%3d-%3d.%1d.%1d.%1d.wcsp%12.8f.fits",nobs1,nobs2,qnx1,ccd1,drift1,min);
	sprintf(comando_iraf,"set wd = `pwd`;cd /home/ceg/Iraf; cl << FIN\n\n\ncd /po2/ceg/UCM-CIDA/OBSERVACIONES/19_7_99/images3\nimarith (operand1=\"%s\",op=\"+\",operand2=\"%s\",result=\"TEMP.fits\",title=\"\",divzero=0.,hparams=\"\",pixtype=\"real\",calctyp=\"\",verbose+,noact-,mode=\"ql\")\nimarith (operand1=\"TEMP.fits\",op=\"/\",operand2=\"2.\",result=\"%s\",title=\"\",divzero=0.,hparams=\"\",pixtype=\"real\",calctyp=\"\",verbose+,noact-,mode=\"ql\")\n\nkeep\nlogout\nFIN\n\ncd /po2/ceg/UCM-CIDA/OBSERVACIONES/19_7_99/images3;rm TEMP.fits\n",specfile1,specfile2,specout);
/* 	//printf("COMANDO:\n%s\n",comando_iraf); */
	system(comando_iraf);
	system("cd /po2/ceg/UCM-CIDA/OBSERVACIONES/19_7_99/images3");
      }


/*       //fprintf(catout,"\n"); */
    }
  }
  fclose(catout);
  return(0);
}

double minequal(int n,double *val,int *nlog,int *eof)
{
  double min,max;
  int i;
  min=val[0];
  max=val[0];
  for(i=1;i<n;i++) {
    if(!eof[i]) {
      if(val[i] > max) max=val[i];
      if(val[i] < min) min=val[i];
    }
  }
  for(i=0;i<n;i++) {
    if(!eof[i] && val[i]==min) nlog[i]=1;
    else nlog[i]=0;
  }
  return(min);
}


void LoadParam_kbd(void) {

  printf(" Input file with objetcs to extract: ");
  reads("",catfile);

  printf(" Input column with RA: ");
  colcatra=readi(colcatra);
  printf(" Input column with DEC: ");
  colcatdec=readi(colcatdec);

  printf(" Input name of database file: ");
  reads("",dbfile);
  
  printf(" Input number of concurrent images in the buffer (as many as compute memory afford): ");
  maxnbuffer=readi(maxnbuffer);

  cpgopen("?");
  cpgask(1);

}

void FindFITS(double ra, double dec, struct SurveyDB surdb, char ***imalist, int *nima) {

  int i;
  char decstr[31];
  char rastr[31];

  *nima=0;

  ra2str (rastr ,31,ra *15,2);
  dec2str(decstr,31,dec,2);

  printf(" Looking images for RA = %s   DEC = %s\n",rastr,decstr);

  if(((*imalist)=(char **) malloc(1*sizeof(char *)))==NULL) {
    printf(" I cannot dimension imalist of 1 element\n");
    exit(1);
  }

  if(((*imalist)[0]=(char *) malloc(101*sizeof(char)))==NULL) {
    printf(" I cannot dimension imalist of 101 element\n");
    exit(1);
  }

  for(i=0;i<surdb.nitems;i++) {
    if(whithinimage(ra,dec,surdb.si[i])) { 
      ra2str (rastr ,31,surdb.si[i].alfac *15,2);
      dec2str(decstr,31,surdb.si[i].deltac,2);
       
      printf(" Match with image %s at RA = %s  DEC = %s  %f %f\n",surdb.si[i].image,rastr,decstr,surdb.si[i].alfac,surdb.si[i].deltac );
      
      (*nima)++;
      if(((*imalist)=(char **) realloc((*imalist), (*nima+1)*101*sizeof(char)))==NULL) {
	printf(" I cannot dimension imalist of %d element\n",*nima);
	exit(1);
      }
      if(((*imalist)[*nima]=(char *) malloc(101*sizeof(char)))==NULL) {
	printf(" I cannot dimension imalist of 101 element\n");
	exit(1);
      }
      strcpy((*imalist)[*nima-1],surdb.si[i].specfile);
    }
/*     printf(" NOOO in image %s at RA = %s  DEC = %s  %f %f\n",surdb.si[i].image,rastr,decstr,surdb.si[i].alfac,surdb.si[i].deltac );  */
  }
}



void SpecExtractFromImages(char **imalist, int nlist, struct spectrum *spectra, int *nspec) {
/*   static int nbuffer=0; */
/*   static struct image *imabuffer; */

  /* Este lo dejo aqui por si algun dia lo necesito.
     Seria extraer los espectros dadas las imagenes y  la ascension recta y declinacion */

}

void SpecExtractFromList(double ra, double dec, char **specfilelist, int nlist, struct spectrum **spectra, struct spectrum **errspectra, int *nspec) {
  static int nbuffer=0;
  static struct specfile *sf;

  int i,j,k;
  int ibuffer;
  int ifree=0;

  int ispec;

  *nspec=0;    


  for(i=0;i<nlist;i++) {
    /* Todo esto es aue viene es para averiguar en que buffer tengo el fichero de espectro.
       Si no lo tengo en ningun buffer, lo leo. Al final la variable que importa es ibuffer */

    ibuffer=-1;
    printf(" Busco si ya he leido %s\n",specfilelist[i]);
    for(j=0;j<nbuffer;j++) {
      if(!strcmp(sf[j].filename,specfilelist[i]))  ibuffer=j;
    }
    
    printf(" Esta en el buffer: %d\n",ibuffer);

    if(ibuffer==-1) {
      if(nbuffer==0) {
	printf(" alocateo por primera vez\n");
	if((sf=(struct specfile *)malloc(1*sizeof(struct specfile)))==NULL) {
	  printf(" Cannot allocate sf of 1 element\n");
	  exit(1);
	}
	ibuffer=0;
	nbuffer=1;
	ifree=1;
      }
      else if(nbuffer<maxnbuffer) {
	printf(" Llevo %d buffers\n",nbuffer);
	nbuffer++;
	if((sf=(struct specfile *)realloc(sf,nbuffer*sizeof(struct specfile)))==NULL) {
	  printf(" Cannot allocate sf of %d element\n",nbuffer);
	  exit(1);
	}
	ibuffer=nbuffer-1;
	ifree=nbuffer;
	if(ifree==maxnbuffer) ifree=0;
      }
      else if(nbuffer==maxnbuffer) {
	printf(" He llegado all final del buffer\n");
	free(sf[ifree].si);
	ibuffer=ifree;
	ifree++;
	if(ifree==maxnbuffer) ifree=0;
      }
      printf(" Al final lo leo en el buffer %d\n",ibuffer); 
      ReadSpecFile(specfilelist[i],sf+ibuffer);
/*       printf(" Ya fuera RA %f DEC %f xp %f yp %f\n",(sf[0]).si[0].ra,(sf[0]).si[0].dec,(sf[0]).si[0].xp,(sf[0]).si[0].yp); */
    }
    /* Hasta aqui  */
    /* Ahora recorro el fichero de espectros en busca del objeto */
    ispec=-1;
    for(k=0;k<sf[ibuffer].nitems;k++) {
/*       printf(" Me voy por %d de %d\n",k,sf[ibuffer].nitems); */
/*       printf(" Las coor %f %f en %s\n",sf[ibuffer].si[k].ra,sf[ibuffer].si[k].dec,sf[ibuffer].si[k].rawspec); */
/*       printf(" Comparo %f %f\n",ra,dec); */
      if(ra==sf[ibuffer].si[k].ra && dec==sf[ibuffer].si[k].dec) { 
/* 	printf(" Pues es %d\n",k); */
	ispec=k;
	break;
      }
    }
    printf(" Salio con ispec %d\n",ispec); 

    /* Ahora leo el espectro y lo guardo en spectra */
    if(ispec!= -1) {
      if((*nspec)==0) {
	printf(" Este es el primer espect\n");
	if(((*spectra)=(struct spectrum *)malloc(1*sizeof(struct spectrum)))==NULL) {
	  printf(" Cannot allocate spectra of 1 element\n");
	  exit(1);
	}
	if(((*errspectra)=(struct spectrum *)malloc(1*sizeof(struct spectrum)))==NULL) {
	  printf(" Cannot allocate errspectra of 1 element\n");
	  exit(1);
	}
	(*spectra)[0].aloc_flag=0;
	(*spectra)[0].alocldo_flag=0;
	(*errspectra)[0].aloc_flag=0;
	(*errspectra)[0].alocldo_flag=0;
	*nspec=1;
      }
      else {
	printf(" Este es el spec %d\n",*nspec);
	(*nspec)++;
	if(((*spectra)=(struct spectrum *)realloc((*spectra),(*nspec)*sizeof(struct spectrum)))==NULL) {
	  printf(" Cannot allocate spectra of %d element\n",*nspec);
	  exit(1);
	}
	if(((*errspectra)=(struct spectrum *)realloc((*errspectra),(*nspec)*sizeof(struct spectrum)))==NULL) {
	  printf(" Cannot allocate errspectra of %d element\n",*nspec);
	  exit(1);
	}
	(*spectra)[*nspec-1].aloc_flag=0;
	(*spectra)[*nspec-1].alocldo_flag=0;
	(*errspectra)[*nspec-1].aloc_flag=0;
	(*errspectra)[*nspec-1].alocldo_flag=0;

      }
      strcpy(((*spectra)[*nspec-1]).file,sf[ibuffer].si[ispec].rawspec);
      strcpy((*errspectra)[*nspec-1].file,sf[ibuffer].si[ispec].errrawspec);
      printf(" Antes de leer\n");
      ReadSpec(&((*spectra)[*nspec-1]));
      ReadSpec(&((*errspectra)[*nspec-1]));
    }
  }
  
}


void ReadSpecFile(char *sfile, struct specfile  *sf) {

  int colraw=1;
  int colerrraw=2;
  int colxp=3;
  int colyp=4;
  int colra=5;
  int coldec=6;
  int colsky=7;

  
  float *xp;
  float *yp;
  float *sky;

  double *ra;
  double *dec;
  char *rawspec;
  char *errrawspec;

  int *ilog;
  int nlin;
  int i;

  printf(" Reading specfile %s...",sfile);
  strcpy(sf->filename,sfile);

  nlin=FileNLin(sfile);
  xp=vector_f(nlin);
  yp=vector_f(nlin);
  sky=vector_f(nlin);
  ra=vector_d(nlin);
  dec=vector_d(nlin);

  if((rawspec   =malloc(51*nlin))==NULL) printf(" Cannot allocate rawspec of %d elements\n",nlin);
  if((errrawspec=malloc(51*nlin))==NULL) printf(" Cannot allocate errrawspec of %d elements\n",nlin);
  ilog=vector_i(nlin);


  ReadCharcol(sfile,colraw,rawspec,ilog,51,&nlin); 
  ReadCharcol(sfile,colerrraw,errrawspec,ilog,51,&nlin); 
  ReadWCScol(sfile,colra,ra,ilog,&nlin);
  ReadWCScol(sfile,coldec,dec,ilog,&nlin);
  ReadNumcol(sfile,colxp,xp,ilog,&nlin);
  ReadNumcol(sfile,colyp,yp,ilog,&nlin);
  ReadNumcol(sfile,colsky,sky,ilog,&nlin);

  (*sf).nitems=0;
  if(((*sf).si=(struct specitem *) malloc(1*sizeof(struct specitem)))==NULL) printf("I cannot dimension (*sf).si   of %d elements \n",1);
  for(i=0;i<nlin;i++) {
    if(ilog[i]) {
/*       printf(" AR %f DEC %f\n",alfac[i],deltac[i]); */
/*       printf(" %d RA %f DEC %f xp %f yp %f\n",(*sf).nitems,ra[i],dec[i],xp[i],yp[i]); */
      ((*sf).si[(*sf).nitems]).ra=ra[i];
      ((*sf).si[(*sf).nitems]).dec=dec[i];
      ((*sf).si[(*sf).nitems]).xp=xp[i];
      ((*sf).si[(*sf).nitems]).yp=yp[i];
      ((*sf).si[(*sf).nitems]).sky=sky[i];
      strcpy(((*sf).si[(*sf).nitems]).rawspec,rawspec+i*51);
      strcpy(((*sf).si[(*sf).nitems]).errrawspec,errrawspec+i*51);
/*       printf("    RA %f DEC %f xp %f yp %f\n",ra[i],dec[i],xp[i],yp[i]); */
      (*sf).nitems++;
      if(((*sf).si=(struct specitem *) realloc((*sf).si   ,((*sf).nitems+1)*sizeof(struct specitem)))==NULL) printf("I cannot dimension (*sf).si        of %d elements \n",(*sf).nitems);
    }
  }

/*   printf(" Primero RA %f DEC %f xp %f yp %f\n",(*sf).si[0].ra,(*sf).si[0].dec,(*sf).si[0].xp,(*sf).si[0].yp); */


  free(ra);free(dec);
  free(xp);free(yp);
  free(sky);
  free(rawspec);
  free(errrawspec);
  printf(" done\n");
}


void ParseObjFile(void) {
  
  int nobj;
  double *ra;
  double *dec;
  int *ilog;

  int i,j;

  char **speclist;
  int nlist;

  struct spectrum *spectra;
  struct spectrum *errspectra;
  int nspec;


  nobj=FileNLin(catfile);
  ra=vector_d(nobj);
  dec=vector_d(nobj);
  ilog=vector_i(nobj);

  ReadWCScol(catfile,colcatra ,ra, ilog,&nobj);
  ReadWCScol(catfile,colcatdec,dec,ilog,&nobj);

  for(i=0;i<nobj;i++) {
    if(ilog[i]){
/*       printf(" Las %d ra %f dec %f\n",i,ra[i],dec[i]); */
      FindFITS(ra[i],dec[i],sdb,&speclist,&nlist);
      for(j=0;j<nlist;j++) {
	printf(" Spec file %s \n",speclist[j]);
      }
      
      printf(" Nos vamos por i %d  \n",i);

      if(i==11) { /* Esto es para la leche en vinagre 666 */
	printf(" DENTRO!!!!\n");
	SpecExtractFromList(ra[i],dec[i],speclist,nlist,&spectra, &errspectra,&nspec);
	for(j=0;j<nspec;j++) {
	  printf(" Painting spectra  %s \n",spectra[j].file);
	  
	  cpgsls(1);
	  if(j==0) PlotSpec_pix(spectra[j]);
	  else PlotSpec_pix_ov(spectra[j]);
	  cpgsls(3);
	  j=readi(j);
	  PlotSpec_pix_ov(errspectra[j]);
	  j=readi(j);
	}
	
	
	if(nspec!=0) exit(1); 
      }

    }
  }
}


void CoaddSpec(struct spectrum *spectra, int nspec, struct spectrum sumspec, struct coadd_info *ci) {



}
