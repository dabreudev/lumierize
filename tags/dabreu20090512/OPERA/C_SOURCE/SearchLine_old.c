#include "modulos.h"
#define NEW 10
#define NFWHM 5
#define NLDO 25
#define FTOL 0.00001



#define FACTORESP 1

/* //#define EWBUENA 26.8 */
#define EWBUENA 0
#define FWHMBUENA 14.3
#define LBUENA 6722.


/* Variables para la funcion respuesta */
struct spectrum response;


int whole_esp;
/*hasta aqui*/

/* Variables globales: los parametros de entrada */
char rootfile[51];
char inputfile[51];
char outfile[51];
char selecfile[51];
char respfile[51];
int namecol;
int errnamecol;
float minselew,maxselew,minselfwhm,maxselfwhm,minselldo,maxselldo;
/*Hasta aqui */

/* Variables para los espectros */
struct spectrum *spectra;
int nspectra;

void ReadResp(char respfile[]);
float resp(float ldo);
void LoadParam_file(char file[100]);
void LoadParam_kbd();
void SaveParam(void);
int Read_whole_esp(
float **esp, float *crval,float *crpix,float *cdelt,long *nldo,int nobj,double *iobj,int *log,int *ntot);
void PCA(float *esp,int nldo,int ntot,float *evec);
float Amoe_Funk(int n, float *x, float *y, float *p);
void FitSpec(struct spectrum spec);
void ParseSpectra(void);
 

int main(int argc, char **argv)
{
  char cnul='A';
  
  int pg1;
  float prob;
  float lineprob;
  float ewprob,fwhmprob;
  
  float probmax;

  int j,i,k;

  float especsig[7];
  int   issigma;
  float *espec;
  float *sigma;
  float *lambda;
  float *chi;
  long nldo;
  float ldomax,ldomin;
  float *sintec;
  float ew,fwhm;
  float *produc;
  float sume,sums,sumsigma,sumchi,sumguena,sumsprob;
  float *logL,*logL0,*logLmax;
  int nesp;
  
  float ldominsearch,ldomaxsearch;
  float halfa;
  /*   float halfaes; */
  int pixmax,pixmin;
  int iter;
  float xmin,xmax,ymin,ymax;
  
  fitsfile *fitsesp;
  FILE *finput;
  int status=0;
  long naxes[2], fpixel,  ii;
  float datamin, datamax, nullval;
  int nobj,ind;
  float *esp;
  int ntot;
  float *evec;
  float pi=4*atan(1);
  float par[6],sigpar[6];
  
  int nfound, anynull;
  float crpix,crval,cdelt;
  char comm[51];
  
  
  FILE *fout;
  FILE *fsel;
  int len=0;
  char record[1000],snul[1000];
  

  char especfile[100];
  
  
  /* Variables a leer del fichero de espectros */
  
  double *iobj;
  int *log;
   


  srandom((unsigned int)time(NULL)/2); 
  if(argc<2) {
    LoadParam_kbd();
  }
  else LoadParam_file(argv[1]);
  ReadResp(respfile);

  
  
  if((fout=fopen(outfile,"w"))==NULL){
    printf(" I could not open %s \n",outfile);
    exit(1);
  }
  if((fsel=fopen(selecfile,"w"))==NULL){
    printf(" I could not open %s \n",selecfile);
    exit(1);
  }
    
  nobj=FileNLin(inputfile);
  if((log=malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension log of %d elements \n",nobj);
  if((iobj=malloc(nobj*sizeof(double)))==NULL) printf("I cannot dimension iobj of %d elements \n",nobj);
  
  printf("Number of objects %d\n",nobj);
  
  ReadDoublecol(inputfile, 1,iobj,log,&nobj);
  /*   for(ind=0;ind<(nobj);ind++) { */
  /*     printf(" %12.8f\n",iobj[ind]); */
  /*   } */
  /*   exit(1); */
  
  fprintf(fout,"# Num     spectrum    ew      fwhm    line_pos     L      L0\n");
  setvbuf(stdin,"",_IOLBF,0);
  
  pg1=cpgopen("?");
  cpgask(0);
  
  if((finput=fopen(inputfile,"r"))==NULL){
    printf(" I could not open %s \n",inputfile);
    exit(1);
  }
  
  rewind(finput);
  fgets(snul,1000,finput); 
  len=strlen(snul);
  strncpy(record,snul,len-1);
  fprintf(fsel,"%s. The rest of the information is from SearchLine.\n",record);
  fgets(snul,1000,finput);
  len=strlen(snul);
  strncpy(record,snul,len-1);
  fprintf(fsel,"%s  EW_FIT FWHM_FIT LAMBDA_FIT  ",record);
  
  
  
  rewind(finput);
  for (i=0;i<nobj;i++) { 
    fgets(snul,1000,finput); 
    if(strlen(snul)>len) len=strlen(snul);
  }
  
  /*    Esto primero no sirve para nada en realidad porque luego calculo la longitud de registro para cada uno! */
  rewind(finput);
  
  whole_esp=Read_whole_esp(&esp,&crval,&crpix,&cdelt,&nldo,nobj,iobj,log,&ntot);
  PCA(esp,nldo,ntot,evec);

  /*   printf(" %d %f %f \n",esp,esp[0],esp[100]); */
  /*   exit(1); */
  /*   if(whole_esp) nobj=ntot; */
  
  /*   exit(1); */
  if((logL=malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension logL of %d elements \n",nobj);
  if((logL0=malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension logL0 of %d elements \n",nobj);
  if((logLmax=malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension logLmax of %d elements \n",nobj);
  
  nesp=0;
  for(ind=0;ind<(nobj);ind++) {
    strcpy(snul,"\0");
    fgets(snul,1000,finput); 
    /*     printf(" SNUL %s\n",snul); */
    len=strlen(snul);
    strncpy(record,snul,len-1);
    if(log[ind]) {
      
      if(!whole_esp) {
	
	status=0;
	strcpy(especfile,rootfile);
	sprintf(snul,"%06d.fits",(int)iobj[ind]);
	/* 	sprintf(snul,"%12.8f.fits",iobj[ind]); */
	strcat(especfile,snul);
	
	
	/* Leo la imagen FITS de entrada */
	printf("...Reading image %s \n",especfile);
	/*       printf(" Cos %g\n",iobj[ind]); */
	/*       exit(1); */
	
	if(ffopen(&fitsesp, especfile, READONLY, &status))     fits_report_error(stdout,status);
	if(ffgky(fitsesp,TFLOAT,"CRPIX1",&crpix,comm,&status)) fits_report_error(stdout,status);
	if(ffgky(fitsesp,TFLOAT,"CDELT1",&cdelt,comm,&status)) fits_report_error(stdout,status);
	if(ffgky(fitsesp,TFLOAT,"CRVAL1",&crval,comm,&status)) fits_report_error(stdout,status);
	if(fits_read_keys_lng(fitsesp, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stdout,status);
	nldo=naxes[0];
	if((espec=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension espec of %ld elements \n",nldo);
	if((sigma=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension sigma of %ld elements \n",nldo);
	if((chi=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension chi of %ld elements \n",nldo);
	if((sintec=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension sinte of %ld elements \n",nldo);
	if((lambda=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension lambda of %ld elements \n",nldo);
	if((produc=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension produc of %ld elements \n",nldo);
	fpixel=1;
	nullval=0;
	if(fits_read_img(fitsesp, TFLOAT, fpixel, nldo, &nullval, espec, &anynull, &status )) fits_report_error(stdout,status);
	if(fits_close_file(fitsesp,&status)) fits_report_error(stdout,status);
	
	printf("crpix %f cdelt %f crval %f \n",crpix,cdelt,crval);
	
      }
      else {
	if((espec=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension espec of %ld elements \n",nldo);
	if((sigma=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension sigma of %ld elements \n",nldo);
	if((chi=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension espec of %ld elements \n",nldo);
	if((sintec=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension sinte of %ld elements \n",nldo);
	if((lambda=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension lambda of %ld elements \n",nldo);
	if((produc=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension produc of %ld elements \n",nldo);
      }
      
      
      if((status==0 && log[ind]) || whole_esp) {
	
	
	/* 	  nldo=100; */
	/* 	printf(" COn el wjole  ind %d ntot %d ndlo %d\n",ind,ntot,nldo); */
	
	for(ii=0 ;ii<nldo;ii++) {
	  if(whole_esp) {
	    /* 	    printf(" nesp %d\n",nesp); */
	    espec[ii]=esp[ii+nldo*nesp];
	    /* 	    	    printf(" %f %f %f %d\n",espec[ii],(ii+1-crpix)*cdelt+crval,esp[ii+nldo*ind],ii+nldo*ind); */
	  }
	  
	  lambda[ii]=(ii+1-crpix)*cdelt+crval;
	  /* 	      espec[ii]=((1+30/sqrt(2.*pi)/20*exp(-(lambda[ii]-6500)*(lambda[ii]-6500)/2./20/20))); */
	  /* 	      espec[ii]=Gasdev()*.1+espec[ii]; */
	  if(ii<3) issigma=0;
	  if(nldo-ii<3) issigma=nldo-7;
	  issigma=ii-3;
	  especsig[0]=espec[issigma  ];
	  especsig[1]=espec[issigma+1];
	  especsig[2]=espec[issigma+2];
	  especsig[3]=espec[issigma+3];
	  especsig[4]=espec[issigma+4];
	  especsig[5]=espec[issigma+4];
	  especsig[6]=espec[issigma+4];
	  StMedia(7,especsig,sigma+ii);
	  sigma[ii]/=100.;
	  
	  /* 	  sigma[ii]=sqrt(fabs(espec[ii])); */
	  /* 	  printf(" ldo %f spec %f sig %f\n",lambda[ii],espec[ii],sigma[ii]); */
	  sigma[ii]=0.1;
	  /* 	  Este puede estar bien para un caso general: */
	  /* 	  sigma[ii]=resp(lambda[ii]); */
	  
	  
	}
	if(whole_esp) nesp++;
	
	/* 	exit(1); */
	
	ldomin=(1-crpix)*cdelt+crval;
	ldomax=(nldo-crpix)*cdelt+crval;
	printf(" Wavelength range %f - %f\n",ldomin,ldomax);
	
	
	datamin=1.0e30;
	datamax=-1.0e30;
	for (ii=0 ;ii<nldo;ii++) {
	  if(espec[ii]< datamin) datamin=espec[ii];
	  if(espec[ii]> datamax) datamax=espec[ii];
	  
	}
	/* 
	   cpgpage();
	   
	   cpgswin(ldomin,ldomax,datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.2);
	   cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
	   cpgline(nldo,lambda,espec);
	*/
	
	
	
	/* Mega bucle */
	
	/* Comento el espectro creado */
	
	datamin=1.0e30;
	datamax=-1.0e30;
	
	
	for (ii=0 ;ii<nldo;ii++) {
	  if( espec[ii]< datamin){
	    datamin=espec[ii];
	    pixmin=ii;
	  }
	  if(espec[ii]> datamax) {
	    datamax=espec[ii];
	    pixmax=ii;
	  }
	  lambda[ii]=(ii+1-crpix)*cdelt+crval;
	}
	/* 	printf(" Datamin %g Datamax %g \n",datamin,datamax); */
	
	/* Empieza el bucle para minimizar */
	probmax=0;
	ldominsearch=6300;
	ldomaxsearch=7200;
	for(k=0;k<NLDO;k++) {
	  /* 	  halfa=ldomin+15+(ldomax-30-ldomin)*k/(nldo/3.-1.); */
	  halfa=ldominsearch+(ldomaxsearch-ldominsearch)*k/(NLDO-1.);
	  for(i=0;i<NEW;i++) {
	    ew=-30+i*230/NEW;
	    ew=0+i*230/NEW;
	    for(j=0;j<NFWHM;j++) {
	      fwhm=5+j*200/NFWHM;
	      for(ii=0 ;ii<nldo;ii++) {
		sintec[ii]=((1+ew/sqrt(2.*pi)/fwhm*exp(-(lambda[ii]-halfa)*(lambda[ii]-halfa)/2./fwhm/fwhm)));
		/* Aqui meto lo de la funcion respuesta */
		if(strcmp(respfile,"NONE") || whole_esp) {
		  /* 		  printf(" sintec %f sintec2 %f\n",sintec[ii],sintec[ii]*resp(lambda[ii])); */
		  sintec[ii]=sintec[ii]*resp(lambda[ii]);
		}
	      }   
	      
	      sume=StSuma1(nldo,espec,1);
	      
	      /*  Esto es mejor, ya que suma aunque no este dentro del rango: */
	      sums=nldo+ew/cdelt;
	      sums=StSuma1(nldo,sintec,1);
	      sumsigma=StSuma1(nldo,sigma,1);
	      sumchi=0;
	      logLmax[ind]=0;
	      for (ii=0 ;ii<nldo;ii++) {
		/* 		 Aqui viene lo bueno del asunto.	 */
		
		chi[ii]=gaussian(espec[ii]/sume,sintec[ii]/sums,sigma[ii]/sumsigma);
		/* 		chi[ii]=gaussian(espec[ii]/sume,sintec[ii]/sums,0.1);	 */
		sumchi+=log10(chi[ii]);
		/* 			  printf(" espec  %g sintec %f sume %g sums %f\n",espec[ii],sintec[ii],sume,sums); */
		/* 		 printf("NORMA  espec  %g sintec %f\n",espec[ii]/sume,sintec[ii]/sums); */
		/* 		 printf(" chi %f sumchi %f\n",chi[ii],sumchi); */
		logLmax[ind]+=log10(gaussian(0.,0.,sigma[ii]/sumsigma));
	      }
	      
	      prob=sumchi;
	      
	      
	      
	      if(prob>probmax) {
		ewprob=ew;
		fwhmprob=fwhm;
		lineprob=halfa;
		probmax=prob;
		sumsprob=sums;
	      }
	      /*  printf("Probability %e with:\n EWCAL %f FWHMCAL %f LINECAL %f\n",prob,ew,fwhm,halfa); */
	      
	      
	    }
	  }
	  
	}
	printf("Maximum probability %e with:\n EWCAL %f FWHMCAL %f LINECAL %f\n",probmax,ewprob,fwhmprob,lineprob);
	
	par[0]=20.;
	par[1]=10.;
	par[2]=(ldomaxsearch+ldominsearch)/2;
	par[3]=crpix;
	par[4]=crval;
	par[5]=cdelt;
	sigpar[0]=20.;
	sigpar[1]=20.;
	sigpar[2]=ldomaxsearch/10.;
	sigpar[3]=0.;
	sigpar[4]=0.;
	sigpar[5]=0.;
	/*    De la otra solucion */
	par[0]=ewprob;
	par[1]=fwhmprob;
	par[2]=lineprob;
	
	iter=Amoeba(nldo,espec,sigma,6,par,sigpar,FTOL,1000,Amoe_Funk);
	printf(" %d \n",iter);
	if(iter==0) {
	  par[0]=ewprob;
	  par[1]=fwhmprob;
	  par[2]=lineprob;
	}
	
	logL[ind]=-Amoe_Funk(nldo,espec,sigma,par);
	
	/* 	Esto es para calcular sin linea */
	for(ii=0 ;ii<nldo;ii++) {
	  sintec[ii]=1.;
	  if(strcmp(respfile,"NONE") || whole_esp)   sintec[ii]=sintec[ii]*resp(lambda[ii]);
	}
	sums=StSuma1(nldo,sintec,1);
	sumguena=0;
	for (ii=0 ;ii<nldo;ii++) {
	  sumguena+=log10(gaussian(espec[ii]/sume,sintec[ii]/sums,sigma[ii]/sumsigma));
	}
	
	logL0[ind]=sumguena;
	/* 	printf("probabilidad con la buena %f\n",sumguena); */
	
	
	
	/* 	fprintf(fout,"%0.6d  %s   %e    %e   %e\n",(int)iobj[ind],especfile,ewprob,fwhmprob,lineprob); */
	fprintf(fout,"%06d  %s   %e    %e   %e   %e   %e   %e   %e   %e\n",(int)iobj[ind],especfile,par[0],par[1],par[2],logL[ind],logL0[ind],logLmax[ind],logL[ind]-logL0[ind],logL[ind]-logLmax[ind]);
	
	
	
	/* 	printf(" EW %f   FWHM %f     LINE %f\n",ewee,fwhmee,halfaes); */
	
	printf("AMOEBA: Maximum probability %e with:\n EWCAL %f FWHMCAL %f LINECAL %f\n\n",logL[ind],par[0],par[1],par[2]);
	if((logL[ind]-logL0[ind])<0) {
	  
	  printf("logL %e logL0 %e logLmax %e L-L0 %e L-Lmax  %e\n",logL[ind],logL0[ind],logLmax[ind],logL[ind]-logL0[ind],logL[ind]-logLmax[ind]);
	  par[0]=ewprob;
	  par[1]=fwhmprob;
	  par[2]=lineprob;
	  printf("With the other %e \n",-Amoe_Funk(nldo,espec,sigma,par));
	  
	  
	  /*     cpgcurs(&fnul,&fnul,&cnul); */
	}
	
	if(lineprob>minselldo && lineprob<maxselldo && ewprob>minselew && ewprob< maxselew && fwhmprob>minselfwhm && fwhmprob< maxselfwhm) {
	  printf(" OBJECT Selected\n");
	  /* 	  fprintf(fsel,"%s   %f  %f  %f\n",record,ewprob,fwhmprob,lineprob); */
	  /* 	  fprintf(fsel,"%12.8f   %f  %f  %f\n",iobj[ind],ewprob,fwhmprob,lineprob); */
	  fprintf(fsel,"%06d   %f  %f  %f\n",(int)iobj[ind],ewprob,fwhmprob,lineprob);
	  
	}
	
	for(ii=0 ;ii<nldo;ii++) {
	  sintec[ii]=((1+ewprob/sqrt(2.*pi)/fwhmprob*exp(-(lambda[ii]-lineprob)*(lambda[ii]-lineprob)/2./fwhmprob/fwhmprob)));
	  if(strcmp(respfile,"NONE") || whole_esp)   sintec[ii]=sintec[ii]*resp(lambda[ii]);
	  sintec[ii]=sintec[ii]/sumsprob;
	  espec[ii]=espec[ii]/sume;
	}
	
	datamin=1.0e30;
	datamax=-1.0e30;
	for (ii=0 ;ii<nldo;ii++) {
	  /* 	if(sigma[ii]< datamin) datamin=sigma[ii]; */
	  if(espec[ii]< datamin) datamin=espec[ii];
	  if(sintec[ii]< datamin) datamin=sintec[ii];
	  if(espec[ii]> datamax) datamax=espec[ii];
	  if(sintec[ii]> datamax)  datamax=sintec[ii];
	  /* 	if(sigma[ii]> datamax)  datamax=sigma[ii]; */
	  
	}
	cpgpage();
	cpgswin(ldomin,ldomax,datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.2);
	cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
	cpglab("Wavelength \\(2078)","Flux",especfile);
	cpgline(nldo,lambda,espec);
	cpgsci(3);
	cpgline(nldo,lambda,sintec);
	/* Ahora pinto el de Amoeba */
	sumsprob=0;
	for(ii=0 ;ii<nldo;ii++) {
	  sintec[ii]=((1+par[0]/sqrt(2.*pi)/par[1]*exp(-(lambda[ii]-par[2])*(lambda[ii]-par[2])/2./par[1]/par[1])));
	  if(strcmp(respfile,"NONE") || whole_esp)   sintec[ii]=sintec[ii]*resp(lambda[ii]);
	  sumsprob+=sintec[ii];
	}
	for(ii=0 ;ii<nldo;ii++) {
	  sintec[ii]=sintec[ii]/sumsprob;
	}
	cpgsci(4);
	cpgline(nldo,lambda,sintec);
	
	
	/*       cpgsci(4); */
	/*       cpgline(nldo,lambda,sigma); */
	cpgsci(1);
	if(lineprob>minselldo && lineprob<maxselldo && ewprob>minselew && ewprob< maxselew && fwhmprob>minselfwhm && fwhmprob< maxselfwhm) {
	  /* 	  cpgcurs(&fnul,&fnul,&cnul); */
	}
	/* 	cpgcurs(&fnul,&fnul,&cnul); */
	free(espec);free(sigma),free(chi);free(sintec);free(lambda);free(produc);
	/* 	exit(1); */
								      }
    }
  }
  xmin=45.;xmax=50.;ymin=-1.;ymax=1.;
  
  while(cnul=='m') {
    cpgpage();
    cpgswin(xmin,xmax,ymin,ymax);
    cpglab("S/N","logL-logL0","");
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    
    for(i=0;i<nobj;i++) {
      printf(" logL %f logL0 %f\n",logL[i],logL0[i]);
      cpgpt1(logL0[i],logL[i]-logL0[i],2);
      /*    cpgpt1((float)i,logL[i]-logL0[i],2); */
    }
    cpgcurs(&xmin,&ymin,&cnul);
    cpgband(2,1,xmin,ymin,&xmax,&ymax,&cnul);
  }
  cpgend();
  pg1=cpgopen("?");
  cpgpage();
  cpgswin(xmin,xmax,ymin,ymax);
  cpglab("S/N","logL-logL0","");
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  
  for(i=0;i<nobj;i++) {
    printf(" logL %f logL0 %f\n",logL[i],logL0[i]);
    cpgpt1(logL0[i],logL[i]-logL0[i],2);
    /*     cpgpt1((float)i,logL[i]-logL0[i],2); */
  }
  cpgpage();
  cpgswin(0.,(float)nobj,ymin,ymax);
  cpglab("S/N","logL-logL0","");
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  
  for(i=0;i<nobj;i++) {
    printf(" logL %f logL0 %f\n",logL[i],logL0[i]);
    cpgpt1((float)i,logL[i]-logL0[i],2);
    /*     cpgpt1((float)i,logL[i]-logL0[i],2); */
  }
  
  
  
  
  
  
  cpgend();
  
  if(argc<2) {
    SaveParam();
  }
  return(0);
}





void ReadResp(char respfile[])
{
  int i;
  if(strcmp(respfile,"NONE")) {
    response.aloc_flag=0;
    response.alocldo_flag=0;
    ReadSpec(&response);
  }
  for(i=0;i<response.nx;i++) {
    printf(" %d ldoresp %f resp %f resp int. %f %f\n",i,response.ldo[i],response.spec[i],resp(response.spec[i]),Lagr4(response.ldo,response.spec,response.nx,response.ldo[i])); 
  }
  
}


float resp(float ldo) {
  return(Lagr2(response.ldo,response.spec,response.nx,ldo));
}



void LoadParam_kbd() {
  printf("\n Input FITS root name with spectra: ");
  reads(rootfile,rootfile);
  printf("\n Input BUHO file: ");
  reads(inputfile,inputfile);
  printf(" Ouput fitting file: ");
  reads(outfile,outfile);

  
  printf(" Output buho file with candidates: ");
  reads(selecfile,selecfile);


  printf(" FITS file with response function of the system (NONE for no response function: ");
  reads(respfile,respfile);
  printf(" Minimum EW selected: ");
  minselew=readf(minselew);
  printf(" Maximum EW selected: ");
  maxselew=readf(maxselew);
  printf(" Minimum FWHM selected: ");
  minselfwhm=readf(minselfwhm);
  printf(" Maximum FWHM selected: ");
  minselfwhm=readf(minselfwhm);
  printf(" Minimum wavelenght selected: ");
  minselldo=readf(minselldo);
  printf(" Maximum wavelenght selected: ");
  maxselldo=readf(maxselldo);
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
/*    ffgky(parfile,TSTRING,"ROOTNAME",inputfile,comment,&status); */
/*    fits_report_error(stderr,status); */
/*    status=0; */
  ffgky(parfile,TINT,"COLNAME",&namecol,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"COLENAME",&errnamecol,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TSTRING,"OBJFILE",inputfile,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TSTRING,"OUTFILE",outfile,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TSTRING,"SELFILE",selecfile,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TSTRING,"RESPFILE",respfile,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"MINEW",&minselew,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"MAXEW",&maxselew,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"MINFWHM",&minselfwhm,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"MAXFWHM",&maxselfwhm,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"MINLDO",&minselldo,comment,&status);
  fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TFLOAT,"MAXLDO",&maxselldo,comment,&status);
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
  char parfile[100];
  char ch51[51];
/*   char opt; */
  printf("Name of parameter file (NONE if you do not want to save current parameters): ");
  reads(parfile,parfile);
  if(!strcmp(parfile,"NONE")) {
    if((fp=fopen(parfile,"w")) ==NULL) {
      printf("ERROR: Can't open file %s\n",parfile);
      return;
    }
    fprintf(fp,"COMMENT  Parameter file for NarrowSearch                                       \n");
    nc++;
    sprintf(ch51,"'%s'",rootfile);
    fprintf(fp,"ROOTNAME= %-51.51s / Root spec name \n",ch51);
    sprintf(ch51,"'%s'",inputfile);
    fprintf(fp,"OBJFILE = %-51.51s / File with objs.\n",ch51);
    sprintf(ch51,"'%s'",outfile  );
    fprintf(fp,"OUTFILE = %-51.51s / Fitting file   \n",ch51);
    sprintf(ch51,"'%s'",selecfile);
    fprintf(fp,"SELFILE = %-51.51s / Selection file \n",ch51);
    sprintf(ch51,"'%s'",respfile );
    fprintf(fp,"RESPFILE= %-51.51s / Response spec. \n",ch51);
    fprintf(fp,"MINEW   =%21f / Minimum selected equivalent width             \n",minselew    );
    fprintf(fp,"MAXEW   =%21f / Maximum selected equivalent width             \n",maxselew    );
    fprintf(fp,"MINFWHM =%21f / Minimum selected Full Width Half Maximum      \n",minselfwhm  );
    fprintf(fp,"MAXFWHM =%21f / Maximum selected Full Widt Half Maximum       \n",maxselfwhm  );
    fprintf(fp,"MINLDO  =%21f / Minimum selected wavelenght                   \n",minselldo   );
    fprintf(fp,"MAXLDO  =%21f / Minimum selected wavelenght                   \n",maxselldo   );
    nc+=11;
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


int Read_whole_esp(float **esp, float *crval,float *crpix,float *cdelt,long *nldo,int nobj,double *iobj,int *log,int *ntot) {
  fitsfile *fitsesp;
  char comm[51];
/*   FILE *finput; */
/*   int nobj; */
  int ind;
  char snul[1000];

/*   //float nldo=0; */
/*   //float crpix,cdelt,crval; */
  float acrpix,acdelt,acrval;
  float *espec;
  int status=0;
  char especfile[51];
  long naxes[2], fpixel,  ii;
  int nfound, anynull;
  float nullval;
  float norma;
  int nesp;
  
  nesp=0;
  *nldo=0;
/*   nobj=FileNLin(inputfile); */
/*   if((finput=fopen(inputfile,"r"))==NULL){ */
/*     printf(" I could not open %s \n",inputfile); */
/*     exit(1); */
/*   } */

  for(ind=0;ind<(nobj);ind++) {
/*     strcpy(snul,"\0"); */
/*     fgets(snul,1000,finput);  */
/*     //    printf(" SNUL %s\n",snul); */
/*     len=strlen(snul); */
/*     strncpy(record,snul,len-1); */
    
    if(log[ind]) {

      status=0;
      strcpy(especfile,rootfile);
      sprintf(snul,"%06d.fits",(int)iobj[ind]);
/*       //sprintf(snul,"%12.8f.fits",iobj[ind]); */
      strcat(especfile,snul);
      /* Leo la imagen FITS de entrada */
/*       //printf("INSIDE READ...Reading image %s \n",especfile); */
      /*       printf(" Cos %g\n",iobj[ind]); */
      /*       exit(1); */
      
      if(ffopen(&fitsesp, especfile, READONLY, &status)) fits_report_error(stdout,status);
      if(ffgky(fitsesp,TFLOAT,"CRPIX1",crpix,comm,&status)) fits_report_error(stdout,status);
      if(ffgky(fitsesp,TFLOAT,"CDELT1",cdelt,comm,&status)) fits_report_error(stdout,status);
      if(ffgky(fitsesp,TFLOAT,"CRVAL1",crval,comm,&status)) fits_report_error(stdout,status);
      if(fits_read_keys_lng(fitsesp, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stdout,status);
      if(*nldo==0) {
	*nldo=naxes[0];
	acrpix=*crpix;
	acdelt=*cdelt;
	acrval=*crval;
	if(((*esp)=malloc(*nldo*nobj*sizeof(float)))==NULL) printf("I cannot dimension esp of %ld elements \n",*nldo*nobj);;
/* 	//printf(" Dimensione %d\n",*nldo*nobj); */
/* 	//exit(1); */
      }
      else {
	if(!(*nldo==naxes[0] || acrpix==*crpix || acdelt==*cdelt || acrval==*crval)) {
	  printf(" Non-uniform sample of spectra.\n Giving up to compute PCA.\n");
	  return(0);
	}
      }
      if((espec=malloc(*nldo*sizeof(float)))==NULL) printf("I cannot dimension espec of %ld elements \n",*nldo);
      fpixel=1;
      nullval=0;
      if(fits_read_img(fitsesp, TFLOAT, fpixel, *nldo, &nullval, espec, &anynull, &status )) fits_report_error(stdout,status);
      if(fits_close_file(fitsesp,&status)) fits_report_error(stdout,status);
/*       //printf(" Status %d\n",status); */
/*       //  nldo=100; */
      if(!status) {
/* 	//(*ntot)++; */
	nesp++;
/* 	//printf(" ntot  %d\n",*ntot); */
	norma=0;
	for(ii=0 ;ii<*nldo;ii++) {
	  norma+=espec[ii];
	}
	for(ii=0 ;ii<*nldo;ii++) {
/* 	  //printf(" Espec %f\n",espec[ii]); */
/* 	  //printf(" COn el wjole ii %d ntot %d ndlo %d\n",ii,*ntot-1,*nldo); */
/* 	  //printf(" %f %f\n",espec[ii],(ii+1-*crpix)**cdelt+*crval); */
	  (*esp)[ii+*nldo*(nesp-1)]=espec[ii]/norma**nldo;
	  (*esp)[ii+*nldo*(nesp-1)]=espec[ii];
/* 	  //printf(" %d %f %f\n",ii+*nldo*(*ntot-1),esp[ii+*nldo*(*ntot-1)],(ii+1-*crpix)**cdelt+*crval); */
	}
/* 	//exit(1); */
      }
    }
  }
  *ntot=nesp;
/*   //printf(" %d %f %f \n",*esp,(*esp)[0],(*esp)[100]); */
  return(1);
}



void PCA(float *esp,int nldo,int ntot,float *evec) {
  float med[ntot];
  float *cov,*cor;
  float eval[nldo];
  int i,j;
  int nnllddoo;
  float *espec;
  float *lambda;
  float ldomin,ldomax,datamin,datamax;
/*   float fnul; */
  char cnul;
  float crpix=1.,crval=5870.4,cdelt=28.32;
  float proj1[ntot],proj2[ntot],proj3[ntot];
  float xmin,xmax,ymin,ymax;
  float mindist=1.e15,dist;
  float xcur,ycur;
  int jsel;

  printf(" Nldo vale %d\n",nldo);
  if((cov=malloc(nldo*nldo*sizeof(float)))==NULL) printf("I cannot dimension cov of %d elements \n",nldo*nldo);
  if((cor=malloc(nldo*nldo*sizeof(float)))==NULL) printf("I cannot dimension cor of %d elements \n",nldo*nldo);
  if((evec=malloc(nldo*nldo*sizeof(float)))==NULL) printf("I cannot dimension evec of %d elements \n",nldo*nldo);
  if((espec=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension espec of %d elements \n",nldo);
  if((lambda=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension lambda of %d elements \n",nldo);

  printf("Aqui\n");
  covcor(esp,ntot,nldo,med,cov,cor);
  printf("Aqui 2\n");
/*   for(i=0;i<nldo;i++) { */
/*     for(j=0;j<i+1;j++) { */
/*        printf(" %f ",cov[i+nldo*j]); */
/*     } */
/*     printf(" \n"); */
/*   } */
  eigen(cov,nldo,eval,evec);
/*   //eigen(cor,nldo,eval,evec); */
  printf(" Nldo vale %d\n",nldo);
  nnllddoo=nldo;
/*   //eigen(cov,50,eval,cov); */
  printf("Aqui 3\n");

  for(i=0;i<nldo;i++) {
    printf(" autoval %d: %g \n",i,eval[i]);
    printf(" autovec %d: ",i);
    for(j=0;j<nldo;j++) {
      espec[j]=evec[j+nldo*i];
      lambda[j]=(j+1-crpix)*cdelt+crval;
      printf(" %f ",evec[j+nldo*i]);
    }
    ldomin=(1-crpix)*cdelt+crval;
    ldomax=(nldo-crpix)*cdelt+crval;
    datamin=1.0e30;
    datamax=-1.0e30;
    for (j=0 ;j<nldo;j++) {
      if(espec[j]< datamin) datamin=espec[j];
      if(espec[j]> datamax) datamax=espec[j];
    }
    cpgpage();
    cpgswin(ldomin,ldomax,datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.2);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    cpglab("Wavelength \\(2078)","Flux","");
    cpgline(nldo,lambda,espec);
/*     //cpgcurs(&fnul,&fnul,&cnul); */

    
    printf("\n");
  }

  for(j=0;j<ntot;j++) {
    proj1[j]=0;
    proj2[j]=0;
    proj3[j]=0;
    for(i=0;i<nldo;i++) {
      proj1[j]+=evec[i+nldo*0]*esp[i+nldo*j];
      proj2[j]+=evec[i+nldo*1]*esp[i+nldo*j];
      proj3[j]+=evec[i+nldo*2]*esp[i+nldo*j];
    }
  }
  xmin=1.0e30;ymin=1.0e30;
  xmax=-1.0e30;ymax=-1.0e30;
  for (j=0 ;j<ntot;j++) {
    if(proj1[j]< xmin) xmin=proj1[j];
    if(proj1[j]> xmax) xmax=proj1[j];
    if(proj2[j]< ymin) ymin=proj2[j];
    if(proj2[j]> ymax) ymax=proj2[j];
  }

  while(!(cnul=='m')) {
    mindist=1.e15;
    cpgpage();
    cpgswin(xmin*1.1,xmax*1.1,ymin*1.1,ymax*1.1);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    cpglab("Projection on eigenvector 1","Projection on eigenvector 2","");
    cpgpt(ntot,proj1,proj2,1);
    cnul='m';
/*     //cpgcurs(&xcur,&ycur,&cnul); */
    for(j=0;j<ntot;j++) {
      dist=((proj1[j]-xcur)*(proj1[j]-xcur)+(proj2[j]-ycur)*(proj2[j]-ycur));
      if(mindist*mindist>dist*dist) {
	mindist=dist;
	jsel=j;
      }
    }
    
    for(j=0;j<nldo;j++) {
      espec[j]=esp[j+nldo*jsel];
    }
    
    datamin=1.0e30;
    datamax=-1.0e30;
    for (j=0 ;j<nldo;j++) {
      if(espec[j]< datamin) datamin=espec[j];
      if(espec[j]> datamax) datamax=espec[j];
    }
    cpgpage();
    cpgswin(ldomin,ldomax,datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.2);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    cpglab("Wavelength \\(2078)","Flux","");
    cpgline(nldo,lambda,espec);
    printf(" Proj1 %f Proj2 %f\n",proj1[jsel],proj2[jsel]);
    
    for(j=0;j<nldo;j++) {
      espec[j]=evec[j+nldo*0]*proj1[jsel];
    }
    cpgsci(2);
    cpgline(nldo,lambda,espec);
    for(j=0;j<nldo;j++) {
      espec[j]=evec[j+nldo*0]*proj1[jsel]+evec[j+nldo*1]*proj2[jsel];
    }
    cpgsci(3);
    cpgline(nldo,lambda,espec);
    for(j=0;j<nldo;j++) {
      espec[j]=evec[j+nldo*0]*proj1[jsel]+evec[j+nldo*1]*proj2[jsel]+evec[j+nldo*2]*proj3[jsel];
    }
    cpgsci(4);
    cpgline(nldo,lambda,espec);
    cpgsci(1);

/*     //cpgcurs(&fnul,&fnul,&cnul); */
  }
    
  if(!strcmp(respfile,"NONE")) {
    response.spec=vector_f(nldo);
    response.ldo=vector_f(nldo);
    response.nx=nldo;
    response.naxes[1]=nldo;
    response.naxes[0]=1;
    response.aloc_flag=1;
    response.alocldo_flag=1;
    for(i=0;i<response.nx;i++) {
      response.spec[i]=evec[i+nldo*0];
      response.ldo[i]=lambda[i];
    }
  }
  printf(" Ha salido de PCA\n");

}

float Amoe_Funk(int n, float *x, float *y, float *p)
{
  
  /* Esta funcion calcula el logaritmo de L */
  /* Los parametros que se pueden variar son la anchura equivalente, 
     la FWHM y la posicion de la linea */
  /* p[0]=EW      (angstroms)
     p[1]=FWHM    (angstroms)
     p[2]=LINEPOS (angstroms)
     p[3]=crpix     Truco para pasarle este parametro
     p[4]=crval     Truco para pasarle este parametro
     p[5]=cdelt     Truco para pasarle este parametro
  */
  
  int i;
  float ew,fwhm,halfa;
  float logL;
  float *sintec,*lambda,*chi;
  int nldo;
  float pi=4*atan(1);
  float sume,sums,sumsigma;
  float crpix,crval,cdelt;
  ew=p[0];
  fwhm=p[1];
  halfa=p[2];
  crpix=p[3];
  crval=p[4];
  cdelt=p[5];
  nldo=n;
  if(fwhm<0.) return(0.);

/*   //p[5] es cdelt */
/*   //printf(" cdelt= %f\n",cdelt); */
/*   //if(fwhm<cdelt && fwhm>0 ) fwhm=cdelt; */

  if((sintec=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension sinte of %d elements \n",nldo);
  if((lambda=malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension lambda of %d elements \n",nldo);
  if((chi   =malloc(nldo*sizeof(float)))==NULL) printf("I cannot dimension chi of %d elements \n",nldo);
  for(i=0 ;i<nldo;i++) {
    lambda[i]=(i+1-crpix)*cdelt+crval;
    sintec[i]=((1+ew/sqrt(2.*pi)/fwhm*exp(-(lambda[i]-halfa)*(lambda[i]-halfa)/2./fwhm/fwhm)));
/*     //printf(" lambda %f sintec %f\n",lambda[i],sintec[i]); */
    /* Aqui meto lo de la funcion respuesta */
    if(strcmp(respfile,"NONE") || whole_esp) {
/*       //printf(" sintec %f resp %f\n",sintec[i],resp(lambda[i])); */
/*       //printf(" sintec %f sintec2 %f\n",sintec[ii],sintec[ii]*resp(lambda[ii])); */
      sintec[i]=sintec[i]*resp(lambda[i]);
    }
  }
  if(halfa<lambda[0] || halfa>lambda[nldo-1]) {
    free(sintec);free(lambda);free(chi);
    return(0.);
  }

   
  sume=StSuma1(nldo,x,1);
  sums=StSuma1(nldo,sintec,1);
  sumsigma=StSuma1(nldo,y,1);

/*   if(halfa < lambda[0]   ) return(0.); */
/*   if(halfa > lambda[nldo]) return(0.); */

  logL=0.;
  for (i=0 ;i<nldo;i++) {
    chi[i]=gaussian(x[i]/sume,sintec[i]/sums,y[i]/sumsigma);	
/*     //printf(" i  %d  x[i] %f sume %f sintec[i] %f sums %f  y[i] %f suchi %f\n",i,x[i],sume,sintec[i],sums,y[i],logL); */
/*     //chi[ii]=gaussian(espec[ii]/sume,sintec[ii]/sums,0.1);	 */
    logL-=log10(chi[i]);
/*     //El - es para que al minimizar lo que en realidad haga es maximizar */
    
  }
  free(sintec);free(lambda);free(chi);
/*   //printf("   %5.1f        %5.1f        %5.1f      %10.7f\n",ew,fwhm,halfa,logL); */
    
  return(logL);
}
  
