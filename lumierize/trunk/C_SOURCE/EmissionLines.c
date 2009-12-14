#include "modulos.h"

#define NBL 4   /* Number of Balmer lines registered */

#define DEBUG 0

/* This program is intented to compute redshifts, corrected extinction fluxes,
   etc... based on observed lambda, EW and fluxes for an object */


struct observedline {

  float lambda;
  float ew;
  float flux;
  float errlambda;
  float errew;
  float errflux;
  float weight;

  float ewrf;
  float errewrf;
  float ewrfcorabs;
  float errewrfcorabs;
  float fluxcorabs;
  float fluxcorext;
  float errfluxcorext;
  float restlambda;
  struct emissionline *assignedline;

};

struct emissionline {
  
  char label[20];
  float lambda;
  struct observedline *assignedline;

};

struct extlaw
{
  float *ldo;
  float *y;
  int n;
};

struct kinematics {

  float zobs,errzobs;
  float velobs,errvelobs;
  float ra,dec;
  
  float zgal,errzgal;
  float velgal,errvelgal;
  float d,errd;
  float dlum,errdlum;
  float mod,errmod;

};


struct reddening {

  float ebv,errebv;
  float c,errc;   /* logarithm 10 */
  float C,errC;   /* neperian log */



};

/* Global variables */

char parfilename[100];
char datafile[100];
int objcol,wlcol,errwlcol,ewcol,errewcol,fluxcol,errfluxcol;
float ebv;
char extinfile[100];
float finf=-1.09;
float ra,dec;
double JD;
float H0=60;
float q0=0.5;
char zfile[100];
char ewfile[100];
char fluxfile[100];
char reddfile[100];
char pgdevice[51];

/* Flags */
int extinflag;
int interact;
int errflag;


/* Line data variables */
int nlines;
char object[100];
struct observedline *OL;
int nel;
struct emissionline *EL;
int nassign;


/* Kinetics variables */

struct kinematics KG;



/*Variables for extinction */
struct extlaw EX;
struct reddening RE;


int iolcompare(const void *X1,const void *X2);
int ielcompare(const void *X1,const void *X2);
void RecordMainEmissionLines();
void LoadParam_kbd();
void LoadParam_file();
void SaveParam();

void MPA(float *data,float *weight,int ndata,float *reference,int nreference,float *a,float *erra,unsigned short int *match,unsigned short int *delete);
void ReadDataFile();
void Identify();
void TypeEL();
float Redshift(float *zerr);
void Velcorrect();
float Vannual(double ra,double dec);
float Vsolar(double l,double b);
float VLSR(double l,double b);
void EWCorrect();
void ReadExtlaw();
void Extinction();
void ExtinBalmer();
void ExtinExcess();
void WriteFiles();
void pol2(float x,float *p,float *y,float *dyda,int n);


float EXT(float ldo);
void Base2(unsigned short int binary[32],long int n,long int base);


int main(int argc, char **argv) 
{

/*   int i; */
  
  if(argc<2) {
    LoadParam_kbd();
    interact=1;
  }
  else {
    strcpy(parfilename,argv[1]);
    LoadParam_file(); 
  }
  
  interact=1;
  
  RecordMainEmissionLines();
  ReadDataFile();
  if(nlines<1) {
    WriteFiles();
    if(interact) SaveParam();
    exit(0);
  }

  Identify();

/*   printf(" PASA 0\n"); */

  Velcorrect();
  
/*   printf(" PASA 1\n"); */

  
  EWCorrect();

/*   printf(" PASA 2\n"); */
  
  Extinction();

  WriteFiles();
  
  if(interact) SaveParam();

  return 0;
}

void RecordMainEmissionLines() {


  nel=21;
  if((EL      =malloc(nel*sizeof(struct emissionline)))==NULL) printf("I cannot dimension EL    of %d elements \n",nel);

  strcpy(EL[0].label,"[MgII]2798");
  EL[0].lambda=2798.00;
  strcpy(EL[1].label,"[OII]3727");
  EL[1].lambda=3727.30;
  strcpy(EL[2].label,"[NeIII]3869");
  EL[2].lambda=3868.74;
  strcpy(EL[3].label,"[HeI]3889+H8");
  EL[3].lambda=3889.0;
  strcpy(EL[4].label,"[NeIII]3968+H7");
  EL[4].lambda=3970.07;
  strcpy(EL[5].label,"Hdelta");
  EL[5].lambda=4101.7;
  strcpy(EL[6].label,"Hgamma");
  EL[6].lambda=4340.5;
  strcpy(EL[7].label,"[OIII]4363");
  EL[7].lambda=4363.21;
  strcpy(EL[8].label,"[HeI]4471");
  EL[8].lambda=4471.47;
  strcpy(EL[9].label,"[HeII]4686");
  EL[9].lambda=4685.68;
  strcpy(EL[10].label,"Hbeta");
  EL[10].lambda=4861.33;
  strcpy(EL[11].label,"[OIII]4959");
  EL[11].lambda=4958.91;
  strcpy(EL[12].label,"[OIII]5007");
  EL[12].lambda=5006.84;
  strcpy(EL[13].label,"[HeI]5876");
  EL[13].lambda=5875.62;
  strcpy(EL[14].label,"[OI]6300");
  EL[14].lambda=6300.23;
  strcpy(EL[15].label,"[NII]6548");
  EL[15].lambda=6548.06;
  strcpy(EL[16].label,"Halfa");
  EL[16].lambda=6562.82;
  strcpy(EL[17].label,"[NII]6584");
  EL[17].lambda=6583.57;
  strcpy(EL[18].label,"[HeI]6678");
  EL[18].lambda=6678.0;
  strcpy(EL[19].label,"[SII]6716");
  EL[19].lambda=6717.00;
  strcpy(EL[20].label,"[SII]6731");
  EL[20].lambda=6731.30;

}



void LoadParam_file() {

  int status=0;
  char comment[51];
  fitsfile *parfile;
  printf("Using parameter file <<%s>>\n",parfilename);
  if( ffopen2(&parfile,parfilename, READONLY, &status)) fits_report_error(stderr,status);

  ffgky(parfile,TSTRING,"DATAFILE",datafile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"OBJCOL",&objcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"WLCOL",&wlcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"ERRWLCOL",&errwlcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EWCOL",&ewcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"ERREWCOL",&errewcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FLUXCOL",&fluxcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"ERRFLCOL",&errfluxcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"EBV",&ebv,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"EXTLAW",extinfile,comment,&status); 
  fits_report_error(stderr,status); 
  ffgky(parfile,TFLOAT,"FINF",&finf,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"RA",&ra,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"DEC",&dec,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TDOUBLE,"JD",&JD,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"H0",&H0,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"Q0",&q0,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"ZFILE",zfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"EWFILE",ewfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"FLUXFILE",fluxfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"REDDFILE",reddfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"DEVICE",pgdevice,comment,&status);
  fits_report_error(stderr,status);
  
  if(strcmp(extinfile,"NONE"))  extinflag=1;
  else     extinflag=0;
  
  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file. Not enough keywords.\n");
/*     exit(1); */
  }

  if(ewcol==0) errflag=0;

  fits_close_file(parfile,&status);
  
  cpgopen(pgdevice);
  cpgask(0);
}

void LoadParam_kbd() {

  objcol=1;
  wlcol=2;
  ewcol=3;
  fluxcol=4;
  errwlcol=0;
  errewcol=0;
  errfluxcol=0;
  ebv=0;
  finf=-1.09;

  printf(" Input file with data containing all information about lines: ");
  reads(datafile,datafile);
  printf(" Input column with object name in %s: ",datafile);
  objcol=readi(objcol);
  printf(" Input column with observed wavelenght of lines in %s: ",datafile);
  wlcol=readi(wlcol);
  printf(" Input column with error for observed wavelenght of lines in %s (0=no errors) : ",datafile);
  errwlcol=readi(errwlcol);
  printf(" Input column with observed EW of lines in %s: ",datafile);
  ewcol=readi(ewcol);
  printf(" Input column with error for EW of lines in %s (0=no errors) : ",datafile);
  errewcol=readi(errewcol);
  printf(" Input column with observed flux of lines in %s: ",datafile);
  fluxcol=readi(fluxcol);
  printf(" Input column with error for observed flux of lines in %s (0=no errors) : ",datafile);
  errfluxcol=readi(errfluxcol);
  printf(" Input E(B-V) extinction constant. If you say 0, then extinctions are computed using decrement in Balmer lines: ");
  ebv=readf(ebv);
  printf(" The extinction law is a table with this structure:\n");
  printf(" lambda     f(lambda)-f(H_beta)\n");
  printf(" Input file with extinction law (NONE=no extinction applied): ");
  reads(extinfile,extinfile);
  if(!strcmp(extinfile,"NONE")) {
    extinflag=0;
  }
  else {
    printf(" Input f(infinite)-f(H_beta): ");
    finf=readf(finf);
    extinflag=1;
  }

  printf(" Input RA for object (needed to compute redshift corrections) :\n");
  ra=readf(0);
  printf(" Input DEC for object:\n");
  dec=readf(0);
  printf(" Input Julian Day:\n");
  JD=readf(JD);
  printf(" Input Hubble constant H0\n");
  H0=readf(H0);
  printf(" Input deceleration parameter q0\n");
  q0=readf(q0);
  printf(" File with redshift information: ");
  reads(zfile,zfile);
  printf(" File with corrected EW: ");
  reads(ewfile,ewfile);
  printf(" File with extinction-corrected fluxes: ");
  reads(fluxfile,fluxfile);
  printf(" File with reddening information: ");
  reads(reddfile,reddfile);


  cpgopen("?");
  cpgask(1);


  if(ewcol==0) errflag=0;
  
  cpgopen("?");
  cpgask(0);

}



void SaveParam() {

  FILE *fp;
  int nc=0,nt;
  char ch51[51];
  char opt='n';
  printf(" Do you want to save parameters in a file?: ");
  opt=readc(opt);
  if(opt=='y') {
    printf("Name of parameter file: ");
    reads(parfilename,parfilename);
    if((fp=fopen(parfilename,"w")) ==NULL) {
      printf("ERROR: Can't open file\n");
      return;
    }
    fprintf(fp,"COMMENT  Parameter file for EmissionLines                                      \n");
    fprintf(fp,"COMMENT  If there are no errors, set those columns to 0                        \n");
    fprintf(fp,"COMMENT  If you provide constant c for extinction, set EBV                     \n");
    fprintf(fp,"COMMENT  setting  EBV=0, extinction is computed via Balmer lines if possible   \n");
    fprintf(fp,"COMMENT  set EXTLAW  to NONE if you don't want to apply extinction             \n");
    nc+=5;
    sprintf(ch51,"'%s'",datafile);
    fprintf(fp,"DATAFILE= %-51.51s / F. /w line data\n",ch51);
    fprintf(fp,"OBJCOL  =%21d / Column with object name in file DATAFILE      \n",objcol      );
    fprintf(fp,"WLCOL   =%21d / Column with observed wavelenght in DATAFILE   \n",wlcol       );
    fprintf(fp,"ERRWLCOL=%21d / Column with error of wavelenght in DATAFILE   \n",errwlcol    );
    fprintf(fp,"EWCOL   =%21d / Column with observed EW in file DATAFILE      \n",ewcol       );
    fprintf(fp,"ERREWCOL=%21d / Column with error of EW in file DATAFILE      \n",errewcol    );
    fprintf(fp,"FLUXCOL =%21d / Column with observed line flux in DATAFILE    \n",fluxcol     );
    fprintf(fp,"ERRFLCOL=%21d / Column with error of flux in file DATAFILE    \n",errfluxcol  );
    fprintf(fp,"EBV     =%21f / E(B-V) ext.  constant for reddening correction\n",ebv         );
    sprintf(ch51,"'%s'",extinfile);
    fprintf(fp,"EXTLAW  = %-51.51s / Extinction file\n",ch51);
    fprintf(fp,"FINF    =%21f / f(infinite)-f(H_beta) for extinction law      \n",finf        );
    fprintf(fp,"RA      =%21f / RA coordinate for galactic velocity correction\n",ra          );
    fprintf(fp,"DEC     =%21f / DEC coor.     for galactic velocity correction\n",dec         );
    fprintf(fp,"JD      =%21f / Julian Day of observations                    \n",JD          );
    fprintf(fp,"H0      =%21f / Hubble constant for distance calculation      \n",H0          );
    fprintf(fp,"Q0      =%21f / Deceleration parameter                        \n",q0          );
    sprintf(ch51,"'%s'",zfile);
    fprintf(fp,"ZFILE   = %-51.51s / Out file w/ Z  \n",ch51);
    sprintf(ch51,"'%s'",ewfile);
    fprintf(fp,"EWFILE  = %-51.51s / Out file w/ EWs\n",ch51);
    sprintf(ch51,"'%s'",fluxfile);
    fprintf(fp,"FLUXFILE= %-51.51s / Out file w/ fl.\n",ch51);
    sprintf(ch51,"'%s'",reddfile);
    fprintf(fp,"REDDFIEE= %-51.51s / Out file w/ red\n",ch51);
    sprintf(ch51,"'%s'","?" );
    fprintf(fp,"DEVICE  = %-51.51s / PGPLOT device  \n",ch51);
    nc+=19;
    fprintf(fp,"COMMENT  Blank lines to assure the last block of 2880 caracters is complete    \n");
    fprintf(fp,"COMMENT                                                                        \n");
    nc+=2;
    for(nt=nc;nt<(1+(int)(nc/36))*36-1;nt++) {
      fprintf(fp,"COMMENT                                                                        \n");
    }
    fprintf(fp,"END                                                                            \n");
    fclose(fp);
  }
}


void ReadDataFile() {

  float *wl,*ew,*flux,*errwl,*errew,*errflux;
  char *name;
  int *ilog;
  int nfile;
  int i;

  printf(" Reading file..."); 
  fflush(NULL);

  nfile=FileNLin(datafile);

/*   printf(" NOOO nfile %d\n",nfile);  */
  
  if((wl      =malloc(nfile*sizeof(float)))==NULL) printf("I cannot dimension wl      of %d elements \n",nfile);
  if((ew      =malloc(nfile*sizeof(float)))==NULL) printf("I cannot dimension ewo     of %d elements \n",nfile);
  if((flux    =malloc(nfile*sizeof(float)))==NULL) printf("I cannot dimension flux    of %d elements \n",nfile);
  if((errwl   =malloc(nfile*sizeof(float)))==NULL) printf("I cannot dimension errwl   of %d elements \n",nfile);
  if((errew   =malloc(nfile*sizeof(float)))==NULL) printf("I cannot dimension errew   of %d elements \n",nfile);
  if((errflux =malloc(nfile*sizeof(float)))==NULL) printf("I cannot dimension errflux of %d elements \n",nfile);
  if((ilog    =malloc(nfile*sizeof(int  )))==NULL) printf("I cannot dimension ilog    of %d elements \n",nfile);
  if((name    =malloc(nfile*51*sizeof(char)))==NULL) printf("I cannot dimension name    of %d elements \n",nfile);

  if((OL      =malloc(nfile*sizeof(struct observedline)))==NULL) printf("I cannot dimension OL    of %d elements \n",nfile);

/*    printf(" SSII\n");  */

  ReadCharcol(datafile,objcol   ,name   ,ilog,51,&nfile);
  ReadNumcol(datafile,wlcol     ,wl     ,ilog,&nfile);
  ReadNumcol(datafile,ewcol     ,ew     ,ilog,&nfile);
  ReadNumcol(datafile,fluxcol   ,flux   ,ilog,&nfile);
  if(errflag) {
    ReadNumcol(datafile,errwlcol  ,errwl  ,ilog,&nfile);
    ReadNumcol(datafile,errewcol  ,errew  ,ilog,&nfile);
    ReadNumcol(datafile,errfluxcol,errflux,ilog,&nfile);
  }
  nlines=0;
  strcpy(object ,name+51);

/*      printf(" MAS\n");  */

  for(i=0;i<nfile;i++) {
    if(ilog[i]) {
      OL[nlines].lambda=wl[i];
      OL[nlines].ew=ew[i];
      OL[nlines].flux=flux[i];
      if(errflag) {
	OL[nlines].errlambda=errwl[i];
	OL[nlines].errew=errew[i];
	OL[nlines].errflux=errflux[i];
      }
      nlines++;
    }
  }

  printf("OK\n");
  free(wl);free(ew);free(flux);free(errwl);free(errew);free(errflux);
  free(ilog);free(name);

}


void Identify() {

  float *data,*weight,*reference,*assigned;

  int ndata,nreference;

  int i,l;
  float a,erra;
  unsigned short int *match,*delete;

  static char agree='y';
  static int nsel1=1,nsel2=0;
  
  float x1,x2,y1,y2;

  float z,errz;
  float minsep=1e19;
  int iminsep=0;
  float *deltaldo;
  float meandeltaldo,sigmadeltaldo;
  int lminsep;
  
  printf(" Identifying lines...");
  fflush(NULL);

  ndata=nlines;
  nreference=nel;

/*   for(i=0;i<nreference;i++) { */
/*     printf("BEFO %i ldo %f\n",i+1,EL[i].lambda); */
/*   } */


  
  if((data      =malloc(ndata*sizeof(float)))==NULL) printf("I cannot dimension data of %d elements \n",ndata);
  if((weight    =malloc(ndata*sizeof(float)))==NULL) printf("I cannot dimension weight of %d elements \n",ndata);
  if((assigned  =malloc(ndata*sizeof(float)))==NULL) printf("I cannot dimension assigned of %d elements \n",ndata);
  if((delete    =malloc(ndata*sizeof(unsigned short)))==NULL) printf("I cannot dimension delete of %d elements \n",ndata);
  if((reference =malloc(nreference*sizeof(float)))==NULL) printf("I cannot dimension reference of %d elements \n",nreference);
  if((match     =malloc(nreference*sizeof(unsigned short)))==NULL) printf("I cannot dimension match of %d elements \n",nreference);
  if((deltaldo  =malloc(ndata*sizeof(float)))==NULL) printf("I cannot dimension data of %d elements \n",ndata);


  qsort(OL,nlines,sizeof(struct observedline),iolcompare);
  qsort(EL,nel,sizeof(struct emissionline),ielcompare);

  

  for(i=0;i<ndata;i++) {
    data[i]=OL[i].lambda;
    if(errflag) weight[i]=1/OL[i].errlambda/OL[i].errlambda;
    else        weight[i]=OL[i].ew*OL[i].ew;
    delete[i]=0;    
  }
  for(i=0;i<nreference;i++) reference[i]=EL[i].lambda;

/*    printf(" After sorting\n");  */
/*    for(i=0;i<nlines;i++) {  */
/*      printf(" %i ldo %f we %f ew %f fl %g \n",i+1,OL[i].lambda,OL[i].ew,OL[i].flux); */ 
/*   } */
/*   for(i=0;i<nreference;i++) { */
/*     printf(" %i ldo %f\n",i+1,EL[i].lambda); */
/*   } */
/*   for(i=0;i<nlines;i++) { */
/*     printf("DATA %i ld %g w %g\n",i+1,data[i],weight[i]); */
/*   } */
/*   for(i=0;i<nreference;i++) { */
/*     printf("REF %i ldo %f\n",i+1,reference[i]); */
/*   } */

  /*   exit(1); */
  
  if(ndata>1)   MPA(data,weight,ndata,reference,nreference,&a,&erra,match,delete);
  else {
    for(i=0;i<nreference;i++) {
      if(fabs(EL[i].lambda-OL[0].lambda)<minsep) {
	minsep=fabs(EL[i].lambda-OL[0].lambda);
	iminsep=i;
 	if(DEBUG) printf(" l1 %f l1 %f YES %d\n",EL[i].lambda,OL[0].lambda,i); 
      }
      if(DEBUG) printf(" minsep %f i %d imin %d \n",minsep,i,iminsep);
    }
    
    for(i=0;i<nreference;i++) {
      if(i==iminsep) match[i]=0;
      else           match[i]=1;
      if(DEBUG) printf(" %d match %d\n",i,match[i]);
    }
    delete[0]=0;
  }

  printf("OK\n");

/*   	for(i=0;i<nreference;i++) {   */
/* 	  if(match[i]) printf("T");   */
/* 	  else         printf("F");   */
/* 	  printf(" %d ",match[i]); */
/* 	}   */
/* 	printf("\n a %f erra %f \n",a,erra);  */

  l=0;
  for(i=0;i<nreference;i++) (EL[i]).assignedline=NULL;

  if(DEBUG) printf(" SEL: ");
  for(i=0;i<nreference;i++) {
    if(DEBUG) printf("%d",match[i]);
    if(!match[i]) {

      if(DEBUG) printf(" l1 %f l2 %f\n",OL[l].lambda,EL[i].lambda); 
      if(delete[l]) {
	(OL[l]).assignedline=NULL;
	(EL[i]).assignedline=NULL;
      }
      else          {
	(OL[l]).assignedline=EL+i;
	(EL[i]).assignedline=OL+l;
      }
      l++;
    }
  }

  z=Redshift(&errz);
  if(DEBUG) printf(" Redshift with these assignments:  z = %f +/- %f\n",z,errz);
  /* Try to reassign possible lines */
  
  if(errflag) {
    for(i=0;i<nlines;i++) {
      if((OL[i]).assignedline==NULL) {
	minsep=1e19;
	lminsep=-1;
	for(l=0;l<nreference;l++) {
	  if((EL[l]).assignedline==NULL) {
	    if(fabs(EL[l].lambda*(1+z)-OL[i].lambda)<OL[i].errlambda && fabs(EL[l].lambda*(1+z)-OL[i].lambda)<minsep) {
	      minsep=fabs(EL[l].lambda*(1+z)-OL[i].lambda);
	      lminsep=l;
	    }
	  }
	}
	if(lminsep!=-1) {
	  if(DEBUG) printf(" %g readmitido \n",OL[i].lambda);
	  (OL[i]).assignedline=EL+l;
	  (EL[l]).assignedline=OL+i;
	}
      }
    }
  }
  else {
    nassign=0;
    for(i=0;i<nlines;i++) {
      if((OL[i]).assignedline!=NULL) {
	deltaldo[nassign]=((OL[i].assignedline)->lambda)*(1+z)-OL[i].lambda;
	nassign++;
      }
    }
    meandeltaldo=StMedia(nassign,deltaldo,&sigmadeltaldo);
    if(DEBUG) printf(" La mean %f sigma %f\n",meandeltaldo,sigmadeltaldo);
    for(i=0;i<nlines;i++) {
      if((OL[i]).assignedline==NULL) {
	if(DEBUG) printf(" con %f\n",OL[i].lambda);
	minsep=1e19;
	lminsep=-1;
	for(l=0;l<nreference;l++) {
	  if((EL[l]).assignedline==NULL) {
	    if(fabs(EL[l].lambda*(1+z)-OL[i].lambda)<sigmadeltaldo*3. && fabs(EL[l].lambda*(1+z)-OL[i].lambda)<minsep) {
	      minsep=fabs(EL[l].lambda*(1+z)-OL[i].lambda);
	      lminsep=l;
	      if(DEBUG) printf(" min %f en %d \n",minsep,l);
	    }
	  }
	}
	if(lminsep!=-1) {
	  if(DEBUG) printf(" %g readmitido \n",OL[i].lambda);
	  (OL[i]).assignedline=EL+lminsep;
	  (EL[lminsep]).assignedline=OL+i;
	}
      }
    }
  }
  




  if(DEBUG) printf(" \n");

  do {
    nassign=0;
    for(i=0;i<nlines;i++) {
      if((OL[i]).assignedline!=NULL) {
	assigned[i]=((OL[i]).assignedline)->lambda;
	nassign++;
      }
      else {
	assigned[i]=assigned[0];
      }
    }
    z=Redshift(&errz);
    if(nassign==0) {
      z=0;errz=0;
    }
    if(DEBUG) {
      for(i=0;i<nel;i++) {
	printf(" %02d %9.3f %18s \n",i,EL[i].lambda,EL[i].label);
	if(EL[i].assignedline!=NULL) printf("    asignada %f \n",((EL[i]).assignedline)->lambda);
      }
    }

    printf(" Assignments:\n");
    printf("    Observed ldo            Assignment                  @ z = %f\n",z);
    for(i=0;i<nlines;i++) {
      if((OL[i]).assignedline!=NULL)  printf("%3i  %9.3f    ----->   %9.3f %18s   %9.3f\n",i+1,OL[i].lambda,((OL[i]).assignedline)->lambda,((OL[i]).assignedline)->label,(((OL[i]).assignedline)->lambda)*(z+1));
      else                            printf("%3i  %9.3f    ----->     No assigment for this line\n",i+1,OL[i].lambda);
    }
    if(nassign==0) printf(" Redshift with these assignments:  Not available\n");
    else           printf(" Redshift with these assignments:  z = %f +/- %f\n",z,errz);

    printf(" Number of lines used %d \n",nassign);
    cpgpage();
    cpgsci(1);
    
    pgLimits(nlines,data,&y1,&y2);
    pgLimits(nlines,assigned,&x1,&x2);
    
    cpgsch(2.);
    cpgenv(x1,x2,y1,y2,0,1);
    cpglab("Reference wavelenght","Input wavelenghts",object);
    cpgsch(1.);
    cpgsci(3);
    
    cpgmove(x1,(z+1)*x1);
    cpgdraw(x2,(z+1)*x2);
    cpgsci(2);
/*     printf(" KAKA FUTI\n"); */
    for(i=0;i<nlines;i++) {
/*       printf(" estoy en %d\n",i); */
/*       if((OL[i]).assignedline!=NULL) cpgpt1(((OL[i]).assignedline)->lambda,OL[i].lambda,9); */
      if((OL[i]).assignedline!=NULL) cpgpt1(((OL[i]).assignedline)->lambda,OL[i].lambda,9);
      else                           cpgpt1(x1,OL[i].lambda,29);
    }
/*     printf(" Sali\n"); */
    if(interact) {  
      printf(" Do you agree with these settings?: ");
      fflush(stdout);
      fflush(NULL);
      setvbuf(stdin,"",_IOLBF,0);
      setvbuf(stdout,"",_IOLBF,0);
      agree=readc(agree);
      if(agree!='y') {
	do {
	  printf(" Type number of line to change: ");
	  fflush(stdout);
	  nsel1=readi(nsel1);
	} while(nsel1<=0 || nsel1>nlines);
 	TypeEL();
	do {
	  printf(" Type number of line to assign to %9.3f (rest frame %9.3f with this redshift): ",(OL[nsel1-1]).lambda,(OL[nsel1-1]).lambda/(1+z));
	  fflush(stdout);
	  if(nsel2<0 || nsel2>nreference ) nsel2=0;
	  nsel2=readi(nsel2);
	} while(nsel2<0 || nsel2>nreference);
	
	if(nsel2==0) {
	  if((OL[nsel1-1]).assignedline==NULL); 
	  else {
	    ((OL[nsel1-1]).assignedline)->assignedline=NULL;
	    (OL[nsel1-1]).assignedline=NULL;
	  }
	}
	else         {
	  if((EL[nsel2-1]).assignedline!=NULL) {
	    printf(" WARNING: Line %s already assigned to %9.3f. Try again\n",(EL[nsel2-1]).label,((EL[nsel2-1]).assignedline)->lambda);
	    nsel2=-1;
	  }
	  else {
	    if((OL[nsel1-1]).assignedline==NULL) {
	      (OL[nsel1-1]).assignedline=EL+nsel2-1;
	      (EL[nsel2-1]).assignedline=OL+nsel1-1;
	    } 
	    else {
	      ((OL[nsel1-1]).assignedline)->assignedline=NULL;
	      (OL[nsel1-1]).assignedline=EL+nsel2-1;
	      (EL[nsel2-1]).assignedline=OL+nsel1-1;
	      if(DEBUG) printf(" Pasa por aqui\n");
	    }
	  }
	}
      }
    }
  }  while(agree!='y');  
  
/*   printf(" Redhift at last: %f\n",z); */
  free(data);free(reference);free(match);free(assigned);free(weight);
  
}






void MPA(float *data,float *weight,int ndata,float *reference,int nreference,float *a,float *erra,unsigned short int *match,unsigned short int *delete) {

  float *x,*y,*yy,*zz,*ww;
  int nx,ny;
  
  int nxreject;
  float sigma;

  float as,erras,ierras;

  unsigned short *bin,*del;

  int i,j,k,n1,l,f,m;
  
  float nsigma=2.5;

/*   float fnul,fnul2; */
  char cnul;

#undef  DEBUG
#define DEBUG 0

  /* Short abreviations */
  nx=ndata;ny=nreference;
  

  if((x     =malloc(nx*sizeof(float)))==NULL) printf("I cannot dimension x of %d elements \n",nx);
  if((ww    =malloc(nx*sizeof(float)))==NULL) printf("I cannot dimension ww of %d elements \n",nx);
  if((zz    =malloc(nx*sizeof(float)))==NULL) printf("I cannot dimension zz of %d elements \n",nx);
  if((y     =malloc(ny*sizeof(float)))==NULL) printf("I cannot dimension y of %d elements \n",ny);
  if((yy    =malloc(ny*sizeof(float)))==NULL) printf("I cannot dimension yy of %d elements \n",ny);
  if((bin   =malloc(ny*sizeof(float)))==NULL) printf("I cannot dimension bin of %d elements \n",ny);
  if((del   =malloc(nx*sizeof(float)))==NULL) printf("I cannot dimension del of %d elements \n",nx);

  /* Short abreviations */
  nx=ndata;ny=nreference;
  memcpy(x,data,nx*sizeof(float));
  memcpy(ww,weight,nx*sizeof(float));
  memcpy(y,reference,ny*sizeof(float));

  /* Sorting x and y. Needed in order for the algorithm to work*/
  /* data, weight and reference must be sorted */

  /*   qsort(x,nx,sizeof(float),icompare); */
  /*   qsort(y,ny,sizeof(float),icompare); */


  if(nx>ny) {
    printf(" Number of input data greater than reference data. Not implemented yet. Exiting\n");
    exit(1);
  }



  ierras=0;
  for(i=1;i<=pow(2,ny);i++) {
    Base2(bin,i,ny);
    n1=0;
    for(j=0;j<ny;j++) if(bin[j]) n1++;

    if(n1==(ny-nx)) {
      memcpy(x,data,nx*sizeof(float));
      memcpy(ww,weight,nx*sizeof(float));
      for(k=0;k<nx;k++) del[k]=0;
      l=0;
      for(k=0;k<ny;k++){
	if(!bin[k]) {
	  yy[l]=y[k];
	  zz[l]=x[l]/yy[l]-1;
	  ww[l]=yy[l]*yy[l]*yy[l]*yy[l]/x[l]/x[l]*ww[l];  /* Esto es como propagar los pesos */
	  /* sigma_z**2=x**2/y**4*sigma_y**2 */
	  /* ww_y=1/sigma_y**2 */
	  /* ww_z=(y**4/x**2)*ww_y */
	  l++;
	}
      }

      as=StWeightMedia(nx,zz,ww,&erras)+1; 
/*       MCLine(nx,yy,x,&as,&erras); */
      sigma=0;
      for(l=0;l<nx;l++)       sigma+=(x[l]-(as)*yy[l])*(x[l]-(as)*yy[l]);
      sigma=sqrt(sigma/(nx-1));

      m=0;
      nxreject=nx;
      if(nx>2) {
	for(l=0;l<nxreject;l++) {
	  if(DEBUG)  	printf(" %f %f %f %f %f\n",x[l],yy[l],(as)*x[l],sigma,((x[l]-(as)*yy[l])/sigma)); 
	  if(fabs((x[l]-(as)*yy[l])/sigma)>nsigma) {
 	    if(DEBUG) {
	      printf("%d-%f  SOL %f +/- %f\n",i,pow(2,ny),as-1,erras); 
	      cpgpage();
	      cpgswin(2000.,7000.,2000.,7000.);
	      cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
	      for(k=0;k<nx;k++) {
		cpgsci(3);
		cpgpt1(yy[k],x[k],9);
		cpgsci(1);
		cpgmove(2000.,2000*as);
		cpgdraw(7000.,7000*as);
		printf(" %f  --> %f %f %f %f\n",x[k],yy[k],(as)*x[k],sigma,((x[k]-(as)*yy[k])/sigma)); 
	      }
	    }
	    del[l]=1;
	    cpgsci(2);
	    cpgpt1(yy[l],x[l],8);
	    cpgsci(1);
	    
	    if(DEBUG) printf(" %d Borrando %d lambda %f\n",i,l,x[l]);
	    if(l<nxreject-1) {
	      memmove(yy+l,yy+l+1,(nxreject-l-1)*sizeof(float)); 
	      memmove(x +l,x +l+1,(nxreject-l-1)*sizeof(float)); 
	      memmove(ww+l,ww+l+1,(nxreject-l-1)*sizeof(float)); 
	      memmove(zz+l,zz+l+1,(nxreject-l-1)*sizeof(float)); 
	      nxreject--; 
	      l++;
	    }
	    if(DEBUG) printf(" Borrada. Quedan %d\n",nxreject);
	  }
	}
      }
      as=StWeightMedia(nxreject,zz,ww,&erras)+1; 

      sigma=0;
      for(l=0;l<nx;l++) 	sigma+=(x[l]-(as)*yy[l])*(x[l]-(as)*yy[l]);
      sigma=sqrt(sigma/(nx-1));
      m=0;
      if(nxreject>2) {
	for(l=0;l<nxreject;l++) {
	  if(del[l]) m++;
	  if(fabs((x[l]-(as)*yy[l])/sigma)>nsigma) {
	    del[l+m]=1;
	    if(l<nxreject-1) {
	      memmove(yy+l,yy+l+1,(nxreject-l-1)*sizeof(float)); 
	      memmove(x +l,x +l+1,(nxreject-l-1)*sizeof(float)); 
	      memmove(ww+l,ww+l+1,(nxreject-l-1)*sizeof(float)); 
	      memmove(zz+l,zz+l+1,(nxreject-l-1)*sizeof(float)); 
	      nxreject--; 
	      l++;
	    }
	  }
	}
      }
      as=StWeightMedia(nxreject,zz,ww,&erras)+1; 

      if((erras==0 || (ierras<(1/erras))) && as>=1) {
	if(DEBUG) {
	  printf(" %i del: ",i);
	  for(f=0;f<nx;f++) printf("%d",del[f]);
	  printf(" \n");
	  printf("%d-%f  SOL2 %f +/- %f\n",i,pow(2,ny),as-1,erras); 
 	  printf(" Anterior %f +/- %f\n",*a-1,*erra);
/* 	  while(cnul!='o') cpgcurs(&fnul,&fnul2,&cnul); */
	  cnul='l';
	  
	}
	ierras=1/erras;
	*a=as;
	*erra=erras;
	for(f=0;f<ny;f++) {
	  match[f]=bin[f];
	}
	for(f=0;f<nx;f++) {
	  delete[f]=del[f];
	}
      }

    }

  }
  free(x);free(y);free(yy);free(zz);free(ww);free(bin);

}

    

  



void Base2(unsigned short *binary,long int n,long int base) {
  
  /* binary must be allocated at least for base members */

  long int j,ii;
  float q;
  
  ii=n;

  if(n>pow(2,base)) exit(1);
  for(j=0;j<base;j++) {
    q=ii/2.;
    ii=ii/2;
    binary[j]=!(ii==q);
  }
}


void TypeEL() {

  
  int i,j;
  int nperline=2;

  printf("  0  NONE\n");

  for(i=0;i<(int)((nel/nperline)+1);i++) {
    for(j=0;j<nperline;j++) {
      if((i*nperline+j)<nel) printf(" %02d %9.3f %18s  ",i*nperline+j+1,EL[i*nperline+j].lambda,EL[i*nperline+j].label);
    }
    printf("\n");
  }
}

float Redshift(float *zerr) {

  float *z,*w;
  int i;
  int nused;

  float zmean;

  if((z     =malloc(nlines*sizeof(float)))==NULL) printf("I cannot dimension z of %d elements \n",nlines);
  if((w     =malloc(nlines*sizeof(float)))==NULL) printf("I cannot dimension w of %d elements \n",nlines);

  nused=0;
  for(i=0;i<nlines;i++) {
    if(((OL[i]).assignedline)!=NULL) {
      z[nused]=OL[i].lambda/(((OL[i]).assignedline)->lambda)-1;
      if(errflag) w[nused]=1/OL[i].errlambda/OL[i].errlambda;
      else        w[nused]=OL[i].ew*OL[i].ew;
      if(DEBUG) printf(" z %f w %f\n",z[nused],w[nused]);
      nused++;
    }
  }

  zmean=StWeightMedia(nused,z,w,zerr);
  if(nused==1) {
    zmean=z[0];
    *zerr=0;
  }

  for(i=0;i<nlines;i++) OL[i].restlambda=OL[i].lambda/(1+zmean);

  

  if(DEBUG) printf(" z %f zerr %f n %d\n",zmean,*zerr,nused);

  free(z);free(w);
  return(zmean);
}


void Velcorrect() {

  /* Aqui se hace solo la correccion por el paso del LSR al centro de la galaxia */
  /* Se han despreciado:
     VDIURNAL:  correccion diurna del orden de 0.46 km/s
     VLUNAR  :  correccion por paso al baricentro Tierra-Luna
     VANNUAL :  correcion por el mov. del baricentro alrededor del sol ~29km/s
     VLSR    :  Paso del sol al LSR.  */
  
  float zobs,velobs,zgal,velgal,d,mod;
  float errzobs,errvelobs,errzgal,errvelgal,errd,errmod; 
  
  double ra_d,dec_d;
  
  float dlum,errdlum;

  const float c=299792.46;        /* En km/s */

  double temp1,temp2,temp3;

  double l,b;

  char ch1[32],ch2[32];

  

  l=ra*15;
  b=dec;

  ra2str(ch1,32,ra*15,2);
  dec2str(ch2,32,dec,2);
  printf("\n Transforming RA = %s   DEC = %s  (J2000)\n",ch1,ch2);
  fk52gal(&l,&b);
  dec2str(ch1,32,l,2);
  dec2str(ch2,32,b,2);
  printf(" to galactic   l = %s    b = %s\n",ch1,ch2);

  zobs=Redshift(&errzobs);

  velobs=(((1+zobs)*(1+zobs))-1)/(((1+zobs)*(1+zobs))+1)*c;
  errvelobs=((4*zobs+4)/((zobs*zobs+2*zobs+2)*(zobs*zobs+2*zobs+2)))*c*errzobs;
  
  ra_d=(double)ra;
  dec_d=(double)dec;
  velgal=velobs+VLSR(l,b)+Vannual(ra_d,dec_d)+Vsolar(l,b);
  errvelgal=errvelobs;

/*   printf(" VLSR  %f  Vannual  %f  Vsolar  %f\n",VLSR(l,b),Vannual(ra_d,dec_d),Vsolar(l,b)); */
  
  zgal=sqrt((1+velgal/c)/(1-velgal/c))-1;
  errzgal=errvelgal*velgal/(c*c*zgal*((1-velgal/c)*(1-velgal/c)));

  d=c/H0*2*(2-2*q0*(1-zgal)-(2-2*q0)*sqrt(1+2*q0*zgal))/(4*q0*q0*(1+zgal));
  temp1=(2*q0-(1-q0)*2*q0/(sqrt(1+2*q0*zgal)));
  temp3=4*q0*q0*(1+zgal);
  temp1=temp1/temp3;
  temp2=((2-2*q0*(1-zgal)-(2-2*q0)*sqrt(1+2*q0*zgal))/((4*q0*q0*(1+zgal))*(1+zgal)));
  errd=fabs(c/H0*2*(temp1-temp2)*errzgal);

  dlum=(1+zgal)*d;
  errdlum=fabs(sqrt(d*d*errzgal*errzgal+(1+zgal)*(1+zgal)*errd*errd));
  mod=-5*log10(dlum/10e-6); /*  10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
  errmod=fabs(5/dlum/log(10))*errdlum;

  KG.zobs=zobs;KG.errzobs=errzobs;
  KG.velobs=velobs;KG.errvelobs=errvelobs;
  KG.zgal=zgal;KG.errzgal=errzgal;
  KG.velgal=velgal;KG.errvelgal=errvelgal;

  KG.ra=ra;KG.dec=dec;
  KG.d=d;KG.errd=errd;
  KG.dlum=dlum;KG.errdlum=errdlum;
  KG.mod=mod;KG.errmod=errmod;

  printf("\n");
  printf(">> Observed Z        : %f +/- %f\n",zobs,errzobs);
  printf(">> Observed Velocity : %f +/- %f  (km/s)\n",velobs,errvelobs);
  printf(">> Galactic Z        : %f +/- %f\n",zgal,errzgal);
  printf(">> Galactic Velocity : %f +/- %f  (km/s)\n",velgal,errvelgal);
  printf(">> Distance          : %f +/- %f  (Mpc)\n",d,errd);
  printf(">> Distance moduli   : %f +/- %f\n",mod,errmod);

  if(interact) printf(" Press Enter to continue\n");
  getchar();
  
}


void Extinction() {
  ReadExtlaw();
  ExtinBalmer();

}


void ExtinBalmer()
{
  struct balmerline {
    char label[20];
    float intens;
    float ewabsorption;
  };

  int i,j,k;
  struct balmerline BL[NBL];
  float x[NBL],y[NBL],erry[NBL];

  float xmin,xmax,ymin,ymax;


/*   float par[3],sigpar[3]; */
/*   int ipar[3]; */
/*   float covarpar[9]; */

  float c,errc;
  
  float a,erra,covac,chi2,q;

  strcpy(BL[0].label,"Halfa");
  BL[0].intens=2.86;
  BL[0].ewabsorption=3.;
  strcpy(BL[1].label,"Hbeta");
  BL[1].intens=1;
  BL[1].ewabsorption=3.;
  strcpy(BL[2].label,"Hgamma");
  BL[2].intens=0.466;
  BL[2].ewabsorption=3.;
  strcpy(BL[3].label,"Hdelta");
  BL[3].intens=0.256;
  BL[3].ewabsorption=3.;

  k=0;


  printf(" Antes ident\n");

  printf(" Identified Balmer lines: ");
  for(i=0;i<nlines;i++) {
    for(j=0;j<NBL;j++) {
/*       printf(" Comparing  %d %s %d %s\n",i,(OL[i].assignedline)->label,j,BL[j].label);  */
      if((OL[i].assignedline)!=NULL) {
/* 	printf(" Pero aui solo si...\n"); */
	if(!strcmp((OL[i].assignedline)->label,BL[j].label)) {
	  printf("  %s ",BL[j].label);
	  OL[i].ewrfcorabs=OL[i].ewrf+BL[j].ewabsorption;
	  OL[i].errewrfcorabs=OL[i].errewrf;
	  OL[i].fluxcorabs=OL[i].flux*(1+BL[j].ewabsorption/OL[i].ewrf);
	  /* 	printf(" %f %f\n",OL[i].flux,OL[i].fluxcorabs); */
	  
	  
	  /* 	printf(" o %d e %d\n",i,j); */
	  /* 	printf(" %s %s\n",(OL[i].assignedline)->label,BL[j].label); */
	  x[k]=EXT((OL[i].assignedline)->lambda)-EXT(EL[9].lambda);  /* EL[9] es Hbeta */
	  y[k]=log10(OL[i].fluxcorabs)-log10(BL[j].intens);
	  /* 	printf(" F %f BL %d %f y %f x %f\n",log10(OL[i].flux),j,log10(BL[j].intens),y[k],x[k]); */
	  if(errflag) erry[k]=OL[i].errflux/OL[i].fluxcorabs*log(10);
	  else        erry[k]=1/sqrt(OL[i].fluxcorabs)*log(10);    /* Poisson noise is assumed */
/* 	  printf(" x %f y %f erry %f\n",x[k],y[k],erry[k]); */
	  k++;
	}
	else {
	  /* 	OL[i].ewcorabs=OL[i].ew; */
	  /* 	OL[i].fluxcorabs=OL[i].flux; */
	}
      }
    }
  }
  
  if(k<2) {
    printf("\n Less than two Balmer lines, not possible to compute reddening\n");
/*     printf(" Computing extinction assuming a given E(B-V)\n"); */
    ExtinExcess();
    c=ebv/0.68;
    errc=0;
  }
  else {
    cpgsci(1);
    pgLimits(k,x,&xmin,&xmax);
    pgLimits(k,y,&ymin,&ymax);
    cpgenv(xmin,xmax,ymin,ymax,0,1);
    cpgsci(2);
    for(i=0;i<k;i++) cpgpt1(x[i],y[i],9);
    cpgsci(1);
    
    MCP1Weight_err(k,x,y,erry,&a,&c,&erra,&errc,&covac,&chi2,&q); 
/*     MCP1_err(k,x,y,&a,&c,&erra,&errc,&covac,&chi2,&q);  */
    if(!errflag) {
      if(k>2) {
	erra*=sqrt((chi2)/(k-2.0));  /* Esto es asi porque en realidad desconozco los errores */
	errc*=sqrt((chi2)/(k-2.0));
      }
      else {
	erra=0;errc=0;
      }
    }
    cpgmove(xmin,c*xmin+a);
    cpgdraw(xmax,c*xmax+a);
    c=-c;
  }
  /* C(logaritmos neperianos) = c(logaritmos decimales) / 0.4343 */
  /* c(logaritmos decimales)  = 1.47 E(B-V) Seaton MNRAS 187,73P */
  /* E(B-V) =  0.68 c(log dec.) */
  /* C(logaritmos neperianos) = 3.38 E(B-V)  */
  /* E(B-V) =  0.295 C(log nep.) */
  
  printf("\n Reddening constant: c= %f +/- %f  C= %f +/- %f E(B-V) = %f +/- %f \n",c,errc,c/0.4343,errc/0.4343,c/1.47,errc/1.47);
  RE.c=c;RE.errc=errc;
  RE.C=c/0.4343;RE.errC=errc/0.4343;
  RE.ebv=c/1.47;RE.errebv=errc/1.47;


  printf("  Obs_ldo        Assignment    Obs_Flux  Flux_cor_abs  Flux_cor_ext\n");

  for(i=0;i<nlines;i++) {
    OL[i].fluxcorext=OL[i].fluxcorabs*pow(10.,(c*(EXT(OL[i].restlambda)-finf)));
    OL[i].errfluxcorext=sqrt(OL[i].fluxcorext/OL[i].fluxcorabs*OL[i].errflux*OL[i].fluxcorext/OL[i].fluxcorabs*OL[i].errflux+OL[i].fluxcorext*errc*(EXT(OL[i].restlambda)-finf)*log(10)*OL[i].fluxcorext*errc*(EXT(OL[i].restlambda)-finf)*log(10));
    
    if((OL[i]).assignedline!=NULL) printf("%3i %9.3f  (%14s)  %10.3g  %10.3g  %10.3g\n",i+1,OL[i].lambda,((OL[i]).assignedline)->label,OL[i].flux,OL[i].fluxcorabs,OL[i].fluxcorext);
    else                           printf("%3i %9.3f  (      ---     )  %10.3g  %10.3g  %10.3g\n",i+1,OL[i].lambda,OL[i].flux,OL[i].fluxcorabs,OL[i].fluxcorext);
    
  }
  

/*   printf(" Computing extinction using Balmer lines\n");    */

}

void ExtinExcess() 
{
  printf(" Computing extinction using E(B-V)\n");
  printf(" A reference, mean extinctions for different object types:\n");
  printf(" SBN: 0.798 (Gallego 1995)    HII: 0.470 (Gallego 1995)  BCD: 0.096 (Gallego 1995)\n");
  if(interact) {
    printf(" Input E(B-V): ");
    ebv=readf(ebv);
  }
  
}


void ReadExtlaw() {
  FILE *fp;
  
  int i;
  EX.n=FileNLin(extinfile);
  if((fp=fopen(extinfile,"r"))==NULL) {
    printf("Cannot open file %s\n",extinfile);
    exit(1);
  }
  
  if((EX.ldo=malloc(EX.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector EX.ldo of %d elements",EX.n);
    exit(1);
  }
  if((EX.y=malloc(EX.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector EX.y of %d elements",EX.n);
    exit(1);
  }
  for (i=0;i<EX.n;i++) {
    fscanf(fp," %f %f",EX.ldo+i,EX.y+i);
  }
}

float EXT(float ldo)
{
  return(Lagr2(EX.ldo,EX.y,EX.n,ldo));
}


void pol2(float x,float *p,float *y,float *dyda,int n)
{
  /* Esta funcion es simplemente un polinomio de segundo grado, 
     para ajustarlo usando Mrq */
  *y=p[0]+p[1]*x;
  dyda[0]=1;
  dyda[1]=x;
}



float Vannual(double ra_d,double dec_d) {
  /* Sacada de /iraf/iraf/noao/astutil/asttools/astvorbit.x
     de IRAF. No funciona bien */

  double t;
  double lperi,eccen,manom,tanom,slong,oblq;
  double l,b;
  
  double x,y,z;

  double vearth;
  const float pi=4*atan(1.);

  float vannual;

  t=(JD-2415020)/36525.;
  manom = 358.47583+t*(35999.04975-t*(0.000150+t*0.000003)); 
  oblq = 23.452294-t*(0.0130125+t*(0.00000164-t*0.000000503)); 
  lperi = 101.22083+t*(1.7191733+t*(0.000453+t*0.000003));
  eccen = 0.01675104-t*(0.00004180+t*0.000000126);
  vearth = ((2*pi*149598500.)/(365.2564 * 86400.))/sqrt(1-eccen*eccen);

  
  while(lperi>360) lperi=lperi-360;
  while(manom>360) manom=manom-360;
  while(oblq<0)    oblq=oblq+360;



  tanom = manom/180*pi+(2*eccen-0.25*eccen*eccen*eccen)*sin(manom/180*pi)+1.25*eccen*eccen*sin(2*manom/180*pi)+13./12.*eccen*eccen*eccen*sin(3*manom/180*pi);
  while(tanom>2*pi) tanom=tanom-2*pi;
  slong = lperi + tanom*180/pi + 180;  
  while(slong>360) slong=slong-360;
/*   printf(" lperi %f manom %f oblq %f tanom %f t %f\n",lperi,manom,oblq,tanom,t); */

  x = cos (ra_d*15/180*pi) * cos (dec_d/180*pi);
  z = sin (ra_d*15/180*pi) * cos (dec_d/180*pi) * (- cos (pi/2-oblq/180*pi)) + sin(dec_d/180*pi) * sin (pi/2-oblq/180*pi);
  y = sin (ra_d*15/180*pi) * cos (dec_d/180*pi) * sin (pi/2-oblq/180*pi) -sin (dec_d/180*pi)  * (- cos (pi/2-oblq)/180*pi);
  l = atan2 (y, x)*180/pi;
  b = asin (z)*180/pi;
  while(l<0) l=l+360;


/*   printf(" l %f lperi %f b %f slong %f\n",l,lperi,b,slong); */
    
  vannual = vearth*cos(b/180*pi)*(sin((slong-l)/180*pi)-eccen*sin((lperi-l)/180*pi));
  
  return(vannual);
  
}

    
    
float Vsolar(double l,double b) {
  float vsolar;
  float bsun=18,lsun=30;   /* Solar apex */
  float U=20;           /* Solar velocity to the apex km/s */
  const float pi=4*atan(1.);
  vsolar=U*(sin(b/180*pi)*sin(bsun/180*pi)+cos(b/180*pi)*cos(bsun/180*pi)*cos((l-lsun)/180*pi));
  return(vsolar);
}

float VLSR(double l,double b) {
  float vlsr;
  const float pi=4*atan(1.);
  vlsr=300*sin(l/180.*pi)*cos(b/180.*pi);
  return(vlsr);
}


void WriteFiles() {

  FILE *few,*fflux,*fredd,*fz;
  int i;

  if((few=fopen(ewfile,"w"))==NULL) {
    printf(" ERROR: Can't open file %s for writing \n",ewfile);
    fclose(few);
    return; 
  }

  if((fflux=fopen(fluxfile,"w"))==NULL) {
    printf(" ERROR: Can't open file %s for writing \n",fluxfile);
    fclose(fflux);
    return; 
  }
  if((fredd=fopen(reddfile,"w"))==NULL) {
    printf(" ERROR: Can't open file %s for writing \n",fluxfile);
    fclose(fflux);
    return; 
  }

  if((fz=fopen(zfile,"w"))==NULL) {
    printf(" ERROR: Can't open file %s for writing \n",zfile);
    fclose(fz);
  }
  
  fprintf(fz,"#Object             Z_obs     err_Z_obs    V_obs      err_V_obs    Z_gal     err_Z_gal     V_gal      err_V_gal    Distance err_Distance  D_moduli    err_D_moduli\n");
  if(nlines<1 || nassign==0 ) fprintf(fz,"%-15s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n",object,"INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF");
  else         fprintf(fz,"%-15s  %10.7f  %10.7f  %10.3f  %10.3f  %10.7f  %10.7f  %10.3f  %10.3f  %10.3f  %10.3f  %10.5f  %10.5f\n",object,KG.zobs,KG.errzobs,KG.velobs,KG.errvelobs,KG.zgal,KG.errzgal,KG.velgal,KG.errvelgal,KG.d,KG.errd,KG.mod,KG.errmod);


  fprintf(fredd,"#Object               c      errc       C       errC      E(B-V)  errE(B-V)\n");
  if(nlines<1) fprintf(fredd,"%-14s %-9s %-9s %-9s %-9s %-9s %-9s\n",object,"INDEF","INDEF","INDEF","INDEF","INDEF","INDEF");
  else         fprintf(fredd,"%-14s %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",object,RE.c,RE.errc,RE.C,RE.errC,RE.ebv,RE.errebv);

  


  fprintf(few,"#Object       ");
  fprintf(fflux,"#Object       ");
  
  for(i=0;i<nel;i++) {
    fprintf(few," EW_%-14s ",EL[i].label);
    fprintf(few," errEW_%-14s ",EL[i].label);
    fprintf(fflux," Flux_%-14s ",EL[i].label);
    fprintf(fflux," errFlux_%-14s ",EL[i].label);
  }
  fprintf(few,"\n# Equivalent widths are corrected from stellar absorption in Balmer lines and from the factor (1+z)\n");
  if(errflag) fprintf(few,"# Errors include only propagated error in redshift\n");
  else        fprintf(few,"# Errors include both observational errors and  propagated error in redshift\n");
  fprintf(fflux,"\n# Fluxes are corrected from stellar absorption in Balmer lines");
  fprintf(fflux,"\n# Fluxes are also corrected from extinction using E(B-V) = %f\n",ebv);

  fprintf(few,"%-14s",object);
  fprintf(fflux,"%-14s",object);

  for(i=0;i<nel;i++) {
    if((EL[i].assignedline)==NULL) {
      fprintf(few," %-17s ","INDEF");
      fprintf(few," %-20s ","INDEF");
      fprintf(fflux," %-19s ","INDEF");
      fprintf(fflux," %-22s ","INDEF");
    }
    else {
      fprintf(few," %-17.5g ",(EL[i].assignedline)->ewrfcorabs);
      fprintf(few," %-20.7g ",(EL[i].assignedline)->errewrfcorabs);
      fprintf(fflux," %-19.5g ",(EL[i].assignedline)->fluxcorext);
      fprintf(fflux," %-22.7g ",(EL[i].assignedline)->errfluxcorext);
    }
  }

  fprintf(few,"\n");
  fprintf(fflux,"\n");
  
  fclose(few);
  fclose(fflux);
  fclose(fredd);
  fclose(fz);
}



void EWCorrect() {
  int i;
  
  for(i=0;i<nlines;i++) {  
    OL[i].ewrf=OL[i].ew/(1+KG.zobs);
    OL[i].errewrf=sqrt(OL[i].errew*OL[i].errew/(1+KG.zobs)/(1+KG.zobs)+OL[i].ew*OL[i].ew*KG.errzobs*KG.errzobs/(1+KG.zobs)/(1+KG.zobs)/(1+KG.zobs)/(1+KG.zobs));
    OL[i].ewrfcorabs=OL[i].ewrf;
    OL[i].errewrfcorabs=OL[i].errewrf;
    OL[i].fluxcorabs=OL[i].flux;
  }

}

int iolcompare(const void *X1,const void *X2)
{  
  struct observedline *x1,*x2; 
  x1=(struct observedline*) X1;
  x2=(struct observedline*) X2;
  if((*x1).lambda  < (*x2).lambda) return(-1);
  if((*x1).lambda == (*x2).lambda) return(0);
  if((*x1).lambda  > (*x2).lambda) return(+1);
  return(0);
}


int ielcompare(const void *X1,const void *X2)
{  
  struct emissionline *x1,*x2; 
  x1=(struct emissionline*) X1;
  x2=(struct emissionline*) X2;
  if((*x1).lambda  < (*x2).lambda) return(-1);
  if((*x1).lambda == (*x2).lambda) return(0);
  if((*x1).lambda  > (*x2).lambda) return(+1);
  return(0);
}


