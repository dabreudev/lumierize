#include "modulos.h"

#define DEBUG 0
 
struct param {

  /* Input parameters */
  char elgfile[51];
  char imafile[51];
  char errimafile[51];
  char respfile[51];
  char specfile[51];
  char regfile[51];
  struct disper_prism DP;
  float minselldo,maxselldo;
  float seeing, platescale;
  char reportprefix[101];
  struct cosmo_param cosmo;
  float Kcoc;
  float niifact;
  float meanebv;
  char fileprop[51];
  int Mabscol;
  int ewcol;
  int niicol;
  int ebvcol;
  double JD;

  /* Output parameters */
  char candidatefile[51];
  char outelgfile[51];

  int interactflag;
  int uacflag;
  
};

struct elgdata {
  char spectrum[101];
  char errspectrum[101];
  float ewobs;
  float errewobs;
  float ldo;;
  float errldo;
  float zobs,errzobs;
  float zgal,errzgal;
  float Klin;
  float x_ima,y_ima;
  float sky;
  double ra,dec;
  int   selected;
  double mag,errmag;
  double Mabs,errMabs;
  float ebv,errebv;
  double fluxline,errfluxline;
  double fluxline_corrext,errfluxline_corrext;
  double fluxline_corrext_corrnii,errfluxline_corrext_corrnii;
  double Lumline,errLumline;
  double Lumline_corrext,errLumline_corrext;
  double Lumline_corrext_corrnii,errLumline_corrext_corrnii;
  float mag_usnoB;
  float mag_usnoR;
  float Mabs_usnoB;
  float Mabs_usnoR;
  char  ucyname[100];
  char  image_sur[51];
  
  struct spectrum spec;
  struct spectrum errspec;
  struct spectrum fitspec;
};

float Vannual(double ra,double dec, double JD);
float Vsolar(double l,double b);
float VLSR(double l,double b);

void LoadParam_file(struct param *P, char file[100]);
void LoadParam_kbd(struct param *P);
void SaveParam(struct param P, char file[100]);

void ReadResp(struct spectrum *response, char respfile[]);
void ReadELGData(struct param P, struct elgdata **ELG, struct spectrum response, int *nelg);
void DoAllSpec(struct param P,struct elgdata *ELG, struct spectrum response);
void CreateSintec(struct spectrum *sintec, struct spectrum response, struct spectrum spec, float ew_pix, float lpix, float fwhm, float K, float minselldo, float maxselldo, struct disper_prism DP);
void ParseELG(struct param P,struct elgdata *ELG,int nelg,struct wcsimage *ima,struct image *errima,struct spectrum response);
void ShowReport(char pgdevice[100], struct elgdata E, struct wcsimage *ima,int dssflag);
void WriteCandidateFile(struct elgdata *ELG, int nelg, char candidatefile[51], char outelgfile[51] );
int main(int argc, char **argv) {

  struct param P;
  struct wcsimage ima;
  struct image errima;
  struct spectrum response;
  struct elgdata *ELG;
  int nelg;
  
  if(argc<2) {
    LoadParam_kbd(&P);
  }
  else LoadParam_file(&P,argv[1]);

  RNDInit();

  ima.image.aloc_flag=0;
  strcpy(ima.image.file,P.imafile);
  ReadWCSImage(&ima);
  
  errima.aloc_flag=0;
/*   strcpy(errima.file,P.errimafile); */
/*   ReadImage(&errima); */
  ReadResp(&response,P.respfile);
  printf("...Reading ELG data\n");
  if(DEBUG) printf(" Antes rld Ho %f \n",P.cosmo.H0);
  ReadELGData(P,&ELG,response, &nelg);

  if(DEBUG) printf(" Antes parse  el spect %s\n",ELG[0].spectrum);

  ParseELG(P,ELG,nelg,&ima,&errima,response);
  WriteCandidateFile(ELG,nelg,P.candidatefile, P.outelgfile);

  return 0;
}

void ParseELG(struct param P,struct elgdata *ELG,int nelg,struct wcsimage *ima,struct image *errima,struct spectrum response) {

  int i;
  char pginter[100]="?";
  char pgreport[100];
  char filereport[100];



  for(i=0;i<nelg;i++) {
    if(DEBUG) printf(" Voy por %d / %d\n",i,nelg);
    printf(" Making report for candidate %d\n",i);
    printf(" Interact %d \n",P.interactflag);
    if(P.interactflag) {
      if(DEBUG) printf(" Aqui el spect %s\n",ELG[i].spectrum);
      if(ELG[i].ewobs>3500) {
        ELG[i].selected=0;
      }
      else {
        ShowReport(pginter,ELG[i],ima,0);
        ELG[i].selected=readi(0);
        if(ELG[i].selected) {
	  if(P.uacflag)       ShowReport(pginter,ELG[i],ima,1);
	  else                ShowReport(pginter,ELG[i],ima,0);
          printf(" Are you sure?\n");
          ELG[i].selected=readi(ELG[i].selected);
        }
      }
    } else {
      ELG[i].selected=1;
    }
    if(ELG[i].selected) {
      sprintf(filereport,"%s%s.ps",P.reportprefix,ELG[i].ucyname);
      if(!access(filereport,F_OK)) {
        printf(" %s already exists. Trying other name\n",filereport);
        sprintf(filereport,"%s%sb.ps",P.reportprefix,ELG[i].ucyname);
        if(!access(filereport,F_OK)) {
          printf(" %s already exists. Trying other name\n",filereport);
          sprintf(filereport,"%s%sc.ps",P.reportprefix,ELG[i].ucyname);
        }
      }
      sprintf(pgreport,"%s/cps",filereport);
      if(P.uacflag)       ShowReport(pgreport ,ELG[i],ima,1);
      else                ShowReport(pgreport ,ELG[i],ima,0);
    }
  }

}

void ShowReport(char pgdevice[100], struct elgdata E,  struct wcsimage *ima, int dssflag) {

  char snul[300]; 
  int nsubx=100,nsuby=50;
  int nsubbx=55,nsubby=10;
  struct wcsimage dssimage;
  int  x1,x2,y1,y2;
  double ra1,ra2,dec1,dec2;
  char rastr[32],decstr[32];
  float first,median,third;
  float *buffer;
  int nimapix;
  int i,j,k;
  float fg,bg;

 
  cpgopen(pgdevice);

  cpgscf(2);
  cpgsch(1.7);
  cpgvstd();
  
/*   cpgmtxt("T",-3.0,0.5,0.5,E.ucyname);  */
  printf(" El nombre UCM %s\n",E.ucyname);
/*   cpglab("","",E.ucyname); */
  cpgsch(1.7);
  cpgsubp(2,3);
  cpgvstd();
  /* Draw Spectrum */
  cpgpanl(1,1);
  strcpy(E.spec.file,"");
  strcpy(E.errspec.file,"");
  strcpy(E.fitspec.file,"");
  PlotSpec_pix_err(E.spec,E.errspec);
  PlotSpec_pix_ov(E.fitspec);

  /* Display Info */
  cpgpage();
/*   cpgpanl(1,2); */
  cpgsch(2.5);
  cpgmtxt("T",0.0,0.5,0.5,E.ucyname);
  cpgsch(1.7);
  sprintf(snul," EW = %f +/- %f",E.ewobs,E.errewobs);
  cpgmtxt("T",-2.5,0.1,0.0,snul);
  sprintf(snul," z = %f +/- %f",E.zobs,E.errzobs);
  cpgmtxt("T",-4.0,0.1,0.0,snul);
  sprintf(snul,"mag\\dcal\\u = %f +/- %f",E.mag,E.errmag);
  cpgmtxt("T",-5.5,0.1,0.0,snul);
  sprintf(snul,"Mabs\\dcal\\u = %f +/- %f",E.Mabs,E.errMabs);
  cpgmtxt("T",-7.0,0.1,0.0,snul);
  sprintf(snul," B\\dUSNO\\u = %f",E.mag_usnoB);
  cpgmtxt("T",-8.5,0.1,0.0,snul);
  sprintf(snul," R\\dUSNO\\u = %f",E.mag_usnoR);
  cpgmtxt("T",-10.0,0.1,0.0,snul);
  ra2str(rastr,32,E.ra*15,3);
  sprintf(snul," RA\\dJ2000\\u = %s ",rastr);
  cpgmtxt("T",-11.5,0.1,0.0,snul);
  dec2str(decstr,32,E.dec,2);
  sprintf(snul," \\gd\\dJ2000\\u = %s ",decstr);
  cpgmtxt("T",-13.0,0.1,0.0,snul);

  /* Dibujar trozo imagen */
  
/*   cpgpanl(2,1); */
  cpgpage();
  


  x1=(int)(E.x_ima-nsubx/2);
  x2=(int)(E.x_ima+nsubx/2);
  y1=(int)(E.y_ima-nsuby/2);
  y2=(int)(E.y_ima+nsuby/2);

  if(DEBUG) printf(" Primero x1 %d x2 %d y1 %d y2 %d\n",x1,x2,y1,y2);

  if(x1<1) {
    x1=1;x2=nsubx;
  }
  if(y1<1) {
    y1=1;y2=nsuby;
  }
  if(x2>(*ima).image.naxes[0]-1) {
    x1=(*ima).image.naxes[0]-1-nsubx;x2=(*ima).image.naxes[0]-1;
  }
  if(y2>(*ima).image.naxes[1]-1) {
    y1=(*ima).image.naxes[1]-1-nsuby;y2=(*ima).image.naxes[1]-1; 
  } 
  nimapix=(x2-x1+1)*(y2-y1+1);
  buffer=vector_f(nimapix);

  k=0;
  for(i=x1;i<=x2;i++) {
    for(j=y1;j<=y2;j++) {
      buffer[k]=(ima->image).array[i+j*ima->image.nx];
      k++;
    }
  }
  Quartil(nimapix,buffer,&first,&median,&third);
  fg=median+(median-first)*8;
  bg=median-(median-first)*3;

  if((*ima).image.array[(int)E.x_ima+(int)E.y_ima*(*ima).image.nx] > median+(median-first)*8) fg = 1.2*(*ima).image.array[(int)E.x_ima+(int)E.y_ima*(*ima).image.nx];
 


  PlotImage_sec_cuts(&((*ima).image),x1,x2,y1,y2,bg,fg);
  cpglab("","",(*ima).image.file);
  cpgsci(3); 
  cpgsch(15.);
/*   cpgpt1(E.x_ima,E.y_ima,4); */
  cpgsci(1);
  cpgsch(1.);
  if(DEBUG) printf(" Image painted\n");
/*   pix2wcs((*ima).wcs,x1,y1,&ra1,&dnul); */
/*   pix2wcs((*ima).wcs,x2,y1,&ra2,&dnul); */
/*   pix2wcs((*ima).wcs,x1,y1,&dnul,&dec1); */
/*   pix2wcs((*ima).wcs,x1,y2,&dnul,&dec2); */
  if(DEBUG) printf(" Range %f %f  %f %f \n",ra1,ra2,dec1,dec2);
/*   cpgswin(ra1*3600,ra2*3600,dec1*3600,dec2*3600); */
/*   cpgtbox("ZXMH",0,0,"ZMD",0,0); */

   
  /* Dibujo barras 3-D */
  cpgpage();
  x1=(int)(E.x_ima-nsubbx/2);
  x2=(int)(E.x_ima+nsubbx/2);
  y1=(int)(E.y_ima-nsubby/2);
  y2=(int)(E.y_ima+nsubby/2);
  if(x1<0) {
    x1=1;x2=nsubbx;
  }
  if(y1<0) {
    y1=1;y2=nsubby;
  }
  if(x2>(*ima).image.naxes[0]-1) {
    x1=(*ima).image.naxes[0]-1-nsubbx;x2=(*ima).image.naxes[0]-2;
  }
  if(y2>(*ima).image.naxes[1]-1) {
    y1=(*ima).image.naxes[1]-1-nsubby;y2=(*ima).image.naxes[1]-2; 
  } 
  cpghist3d((*ima).image.array,(*ima).image.nx,(*ima).image.ny,x1,x2,y1,y2);

  /* Dibujar DSS */

/*   cpgpanl(2,1); */
  cpgpage();
  dssimage.image.aloc_flag=0;
  dssimage.wcsset=0;

  if(dssflag)  {
    retrieveDSS(&dssimage,E.ra,E.dec,6.,6.);
    PlotImage(&dssimage.image);
    cpgsci(2);
    cpgmove(dssimage.image.nx*0.3,dssimage.image.ny/2.);cpgdraw(dssimage.image.nx*0.4,dssimage.image.ny/2.);
    cpgmove(dssimage.image.nx*0.6,dssimage.image.ny/2.);cpgdraw(dssimage.image.nx*0.7,dssimage.image.ny/2.);
    cpgmove(dssimage.image.nx/2.,dssimage.image.ny*0.3);cpgdraw(dssimage.image.nx/2.,dssimage.image.ny*0.4);
    cpgmove(dssimage.image.nx/2.,dssimage.image.ny*0.6);cpgdraw(dssimage.image.nx/2.,dssimage.image.ny*0.7);
    cpgsci(1);
  }
  if(dssimage.image.aloc_flag) free(dssimage.image.array);
  if(dssimage.wcsset) wcsfree(dssimage.wcs);


  cpgclos();

  free(buffer); 

}
 
void WriteCandidateFile(struct elgdata *ELG, int nelg, char candidatefile[51], char outelgfile[51] ) {

  int i;
  int newelgflag;

  FILE *fcand,*felg;
  char rastr[32],decstr[32];

  if(strcmp(outelgfile,"NONE"))    newelgflag=1;
  if((fcand=fopen(candidatefile,"w"))==NULL) {
    printf(" I cannot open file %s for writing\n",candidatefile);
    exit(1);
  }
  if((felg=fopen(outelgfile,"w"))==NULL) {
    printf(" I cannot open file %s for writing\n",candidatefile);
    exit(1);
  }

  fprintf(fcand,"#%-20s  %-32s  %-32s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-15s  %-20s  %-20s  %-20s  %-20s  %-20s  %-20s  %-10s  %-51s\n","UCYname","RA_J2000","Dec_J2000","EW_obs","errEW_obs","ldo_line","errldo_line","z_obs","errz_obs","z_gal","errz_gal","mag_obs","errmag_obs","Mabs_obs","errMabs_obs","E(B-V)", "err_E(B-V)","flux_line_obs","err_fl_line_obs","fl_lin_corr_abs","err_fl_lin_corr_abs","fl_lin_corr_abs_nii","err_fl_lin_corr_abs_nii","color","Image_sur");

  fprintf(felg,"#Output candidate file from Select\n");
  fprintf(felg,"#%-101s  %-101s  %-15s %-15s %-15s %-15s %-15s\n","Spectrum","Err_spectrum","EW_fit","ERR_EW_fit","LDO_fit","ERR_LDO_fit","K_lin");
/*   fprintf(felg,"# Original spectra file: %s\n",specfile); */
/*   fprintf(fc,"# Selection cuts: EW %f - %f   LDO %f -%f\n",minselew,maxselew,minselldo,maxselldo); */


  for(i=0;i<nelg;i++) {
    if(ELG[i].selected) {
      ra2str(rastr,32,ELG[i].ra*15,3);
      dec2str(decstr,32,ELG[i].dec,2);
      fprintf(fcand,"%-20s  %-32s  %-32s  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-15.10g  %-20.15g  %-20.15g  %-20.15g  %-20.15g  %-20.15g  %-20.15g  %-10.5g  %-5s\n",ELG[i].ucyname,rastr,decstr,ELG[i].ewobs,ELG[i].errewobs,ELG[i].ldo,ELG[i].errldo,ELG[i].zobs,ELG[i].errzobs,ELG[i].zgal,ELG[i].errzgal,ELG[i].mag,ELG[i].errmag,ELG[i].Mabs,ELG[i].errMabs,ELG[i].ebv,ELG[i].errebv,ELG[i].fluxline,ELG[i].errfluxline,ELG[i].fluxline_corrext,ELG[i].errfluxline_corrext,ELG[i].fluxline_corrext_corrnii,ELG[i].errfluxline_corrext_corrnii,ELG[i].Klin,ELG[i].image_sur);
      if(newelgflag)       fprintf(felg,"%-101s  %-101s   %-15.5g %-15.5g %-15.5g %-15.5g %-15.10g\n",ELG[i].spectrum,ELG[i].errspectrum,ELG[i].ewobs,ELG[i].errewobs,ELG[i].ldo,ELG[i].errldo,ELG[i].Klin);
    }
  }

  fclose(fcand);
  fclose(felg);
  
}




void ReadResp(struct spectrum *response, char respfile[])
{
  if(strcmp(respfile,"NONE")) {
    (*response).aloc_flag=0;
    (*response).alocldo_flag=0;
    strcpy((*response).file,respfile);
    ReadSpec(response);
  }
}



void ReadELGData(struct param P, struct elgdata **ELG, struct spectrum response, int *nelg) {
  
  int i,j;

  char *spec;
  char *errspec;
  float *ew;
  float *errew;
  float *ldo;
  float *errldo;
  float *K;
  int *ilog;
  int speccol=1,errspeccol=2,ewcol=3,errewcol=4,ldocol=5,errldocol=6,Kcol=7;
  int nlines;

  float *x_ima,*y_ima;
  char *specspec;
  char *errspecspec;
  float *sky;
  double *ra,*dec;
  int *ilogspec;
  int nspec;
  int specspeccol=1,errspecspeccol=2,ximacol=3,yimacol=4,racol=5,deccol=6,skycol=7;

  float ldorest=6562.82;
  float gamma=0.6626;
  float delta=1315;
  float calzetti=0.44;
  float A_6563=2.53;

  double gnum[2000],gmag[2000],gmagb[2000];
  double gar[2000],gdec[2000];
  int platenum[2000];
  int nusnostar;

  int rah, ram;
  float ras;
  int decg, decm;
  float decs;
  char dsig;

  double l,b;
  double velobs,velgal;
  double errvelobs,errvelgal;
  const float c=299792.46;        /* En km/s */
  char rastr[32],decstr[32];




  nlines=FileNLin(P.elgfile);
  nspec=FileNLin(P.specfile);

  *ELG=malloc(nlines*sizeof(struct elgdata));
  spec=malloc(nlines*101*sizeof(char));
  errspec=malloc(nlines*101*sizeof(char));
  ew=vector_f(nlines);
  errew=vector_f(nlines);
  ldo=vector_f(nlines);
  errldo=vector_f(nlines);
  K=vector_f(nlines);
  ilog=vector_i(nlines);

  specspec=malloc(nspec*101*sizeof(char));
  errspecspec=malloc(nspec*101*sizeof(char));
  x_ima=vector_f(nspec);
  y_ima=vector_f(nspec);
  sky=vector_f(nspec);
  ra=vector_d(nspec);
  dec=vector_d(nspec);
  ilogspec=vector_i(nspec);

  for(i=0;i<nlines;i++) ilog[i]=0;

  ReadNumcol(P.elgfile,ewcol,ew,ilog,&nlines);
  ReadNumcol(P.elgfile,errewcol,errew,ilog,&nlines);
  ReadNumcol(P.elgfile,ldocol,ldo,ilog,&nlines);
  ReadNumcol(P.elgfile,errldocol,errldo,ilog,&nlines);
  ReadNumcol(P.elgfile,Kcol,K,ilog,&nlines);
  ReadCharcol(P.elgfile,speccol,spec,ilog,101,&nlines);
  ReadCharcol(P.elgfile,errspeccol,errspec,ilog,101,&nlines);

  ReadCharcol(P.specfile,specspeccol,specspec,ilogspec,101,&nspec);
  ReadCharcol(P.specfile,errspecspeccol,errspecspec,ilogspec,101,&nspec);
  ReadNumcol(P.specfile,ximacol,x_ima,ilogspec,&nspec);
  ReadNumcol(P.specfile,yimacol,y_ima,ilogspec,&nspec);
  ReadNumcol(P.specfile,skycol,sky,ilogspec,&nspec);
  ReadWCScol(P.specfile,racol,ra,ilogspec,&nspec);
  ReadWCScol(P.specfile,deccol,dec,ilogspec,&nspec);

  if(DEBUG) printf(" final read\n");

  *nelg=0;
  for(i=0;i<nlines;i++) {
    if(ilog[i]) {
      if(DEBUG) printf("n %d ew %f spec %s \n",i,ew[i],spec+i*101);
      (*ELG)[*nelg].ewobs=ew[i];
      (*ELG)[*nelg].errewobs=errew[i];
      (*ELG)[*nelg].ldo=ldo[i];
      (*ELG)[*nelg].errldo=errldo[i];
      (*ELG)[*nelg].Klin=K[i];
      strcpy((*ELG)[*nelg].spectrum,spec+i*101);
      strcpy((*ELG)[*nelg].errspectrum,errspec+i*101);
      if(DEBUG) printf(" Escrito spec %s\n",(*ELG)[*nelg].spectrum);
      for(j=0;j<nspec;j++) {
	if(ilogspec[j]) {
/* 	  printf(" Comparando %s %s\n",(*ELG)[*nelg].spectrum,specspec+j*101); */
	  if(!strcmp((*ELG)[*nelg].spectrum,specspec+j*101)) {
	    if(DEBUG) printf(" Macheado con %d\n",j);
	    (*ELG)[*nelg].x_ima=x_ima[j];
	    (*ELG)[*nelg].y_ima=y_ima[j];
	    (*ELG)[*nelg].sky=sky[j];
	    (*ELG)[*nelg].ra=ra[j];
	    (*ELG)[*nelg].dec=dec[j];
	  }
	}
      }
      if(P.uacflag) {
	nusnostar=uacread("UAC2",1,(*ELG)[*nelg].ra*15,(*ELG)[*nelg].dec,1/3600.,1/3600.,0.,WCS_J2000,2000.,0.0,0.,25.,200,gnum,gar,gdec,gmag,gmagb,platenum,1);
        printf(" He buscado en %f %f\n",(*ELG)[*nelg].ra,(*ELG)[*nelg].dec);
        ra2str(rastr,32,(*ELG)[*nelg].ra*15,3);
        dec2str(decstr,32,(*ELG)[*nelg].dec,2);
        printf("           RA %s DEC %s\n",rastr,decstr);
        printf(" La primera   %f %f \n",gar[0]/15.,gdec[0]);
        ra2str(rastr,32,gar[0],3);
        dec2str(decstr,32,gdec[0],2);
        printf("           RA %s DEC %s\n",rastr,decstr);
	
        if(fabs(gar[0]-(*ELG)[*nelg].ra*15) + fabs(gdec[0]-(*ELG)[*nelg].dec)  >0.0004 ) {
          if(nusnostar==2) {
            printf(" La segunda   %f %f \n",gar[1]/15.,gdec[1]);
            ra2str(rastr,32,gar[1],3);
            dec2str(decstr,32,gdec[1],2);
            printf("           RA %s DEC %s\n",rastr,decstr);
	    
            if(fabs(gar[1]-(*ELG)[*nelg].ra*15) + fabs(gdec[1]-(*ELG)[*nelg].dec)  >0.00004 ) {
	      printf(" ERROR!: Object not found in USNO catalog\n");
	      exit(1);
            }
          } else {
            printf(" ERROR!: Object not found in USNO catalog\n");
            exit(1);
          }
        }
      }
      else  {
        gnum[0]=0; gar[0]=0; gdec[0]=0; gmag[0]=0; gmagb[0]=0; platenum[0]=0;
      }
      
      r2hms((*ELG)[*nelg].ra*M_PI/12, &rah, &ram, &ras);
      r2gms((*ELG)[*nelg].dec*M_PI/180, &dsig,  &decg, &decm, &decs);
		  
		  
      sprintf((*ELG)[*nelg].ucyname,"UCY%02d%02d%c%02d%02d",rah,ram,dsig,decg,decm);

      

      (*ELG)[*nelg].zobs=(*ELG)[*nelg].ldo/ldorest-1;
      if((*ELG)[*nelg].zobs==0) (*ELG)[*nelg].zobs=0.00001;
      (*ELG)[*nelg].errzobs=(*ELG)[*nelg].errldo/ldorest;
      (*ELG)[*nelg].errmag=0.25;
      l=(*ELG)[*nelg].ra*15;
      b=(*ELG)[*nelg].dec;
      fk52gal(&l,&b);
      velobs=(((1+(*ELG)[*nelg].zobs)*(1+(*ELG)[*nelg].zobs))-1)/(((1+(*ELG)[*nelg].zobs)*(1+(*ELG)[*nelg].zobs))+1)*c;
      errvelobs=((4*(*ELG)[*nelg].zobs+4)/(((*ELG)[*nelg].zobs*(*ELG)[*nelg].zobs+2*(*ELG)[*nelg].zobs+2)*((*ELG)[*nelg].zobs*(*ELG)[*nelg].zobs+2*(*ELG)[*nelg].zobs+2)))*c*(*ELG)[*nelg].errzobs;
      velgal=velobs+VLSR(l,b)+Vannual((*ELG)[*nelg].ra,(*ELG)[*nelg].dec,P.JD)+Vsolar(l,b);
      errvelgal=errvelobs;
      (*ELG)[*nelg].zgal=sqrt((1+velgal/c)/(1-velgal/c))-1;
      (*ELG)[*nelg].errzgal=errvelgal*velgal/(c*c*(*ELG)[*nelg].zgal*((1-velgal/c)*(1-velgal/c)));
      (*ELG)[*nelg].mag_usnoR=gmag[0];
      (*ELG)[*nelg].mag_usnoB=gmagb[0];
      (*ELG)[*nelg].mag=gmag[0];
      Mag_err((double)(*ELG)[*nelg].zgal,(double)(*ELG)[*nelg].errzgal,(double)(*ELG)[*nelg].mag,(double)(*ELG)[*nelg].errmag,P.cosmo,&((*ELG)[*nelg].Mabs),&((*ELG)[*nelg].errMabs));
      printf(" Mag %g err %g\n",(*ELG)[*nelg].Mabs,(*ELG)[*nelg].errMabs);
      (*ELG)[*nelg].ebv=P.meanebv;
      (*ELG)[*nelg].errebv=0.;
      



      
      Flux_ew_mag_err((*ELG)[*nelg].ewobs, (*ELG)[*nelg].errewobs, (*ELG)[*nelg].mag,(*ELG)[*nelg].errmag, "R", gamma, delta, P.Kcoc, &((*ELG)[*nelg].fluxline), &((*ELG)[*nelg].errfluxline));
      (*ELG)[*nelg].fluxline_corrext=(*ELG)[*nelg].fluxline*pow(10.,0.4*A_6563*(*ELG)[*nelg].ebv); /* No hace falta Calzetti porque son flujos de líneas de emisión */
      (*ELG)[*nelg].errfluxline_corrext= sqrt( ((*ELG)[*nelg].errfluxline*pow(10.,0.4*calzetti*A_6563*(*ELG)[*nelg].ebv))*((*ELG)[*nelg].errfluxline*pow(10.,0.4*calzetti*A_6563*(*ELG)[*nelg].ebv)) + ((*ELG)[*nelg].fluxline*(*ELG)[*nelg].errebv*0.4*calzetti*A_6563*log(10)*pow(10.,0.4*calzetti*A_6563*(*ELG)[*nelg].ebv))*((*ELG)[*nelg].fluxline*(*ELG)[*nelg].errebv*0.4*calzetti*A_6563*log(10)*pow(10.,0.4*calzetti*A_6563*(*ELG)[*nelg].ebv)));
      (*ELG)[*nelg].fluxline_corrext_corrnii=(*ELG)[*nelg].fluxline_corrext/(1+P.niifact);
      (*ELG)[*nelg].errfluxline_corrext_corrnii=(*ELG)[*nelg].errfluxline_corrext/(1+P.niifact);
      
      (*ELG)[*nelg].Lumline=Lum((*ELG)[*nelg].zgal,(*ELG)[*nelg].fluxline,P.cosmo);
      (*ELG)[*nelg].Lumline_corrext=Lum((*ELG)[*nelg].zgal,(*ELG)[*nelg].fluxline_corrext,P.cosmo);
      (*ELG)[*nelg].Lumline_corrext_corrnii=Lum((*ELG)[*nelg].zgal,(*ELG)[*nelg].fluxline_corrext_corrnii,P.cosmo);

      strcpy((*ELG)[*nelg].image_sur,P.imafile);

      printf(" Antes doall\n");

      (*ELG)[*nelg].spec.aloc_flag=0;(*ELG)[*nelg].spec.alocldo_flag=0;
      (*ELG)[*nelg].errspec.aloc_flag=0;(*ELG)[*nelg].errspec.alocldo_flag=0;
      (*ELG)[*nelg].fitspec.aloc_flag=0;(*ELG)[*nelg].fitspec.alocldo_flag=0;


      DoAllSpec(P,&((*ELG)[*nelg]),response);
      
      printf(" Despues doall\n");

      (*ELG)[*nelg].selected=0;
      
      (*nelg)++;
    }
  }

  if(DEBUG) printf(" Al final read  el spect %s\n",(*ELG)[0].spectrum);
  

/*   ewOIIcorr[i]=ewOIIobs[i]*pow(10.,0.4*((1.-Calzet)*ebv[i]*A_3727)); */

  
  
  free(spec);
  free(errspec);
  free(ew);
  free(errew);
  free(ldo);
  free(errldo);
  free(K);
  free(ilog);
  free(x_ima);
  free(y_ima);
  free(sky);
  free(ra);
  free(dec);
  if(DEBUG) printf(" mas final read  el spect %s\n",(*ELG)[0].spectrum);
}

void CreateSintec(struct spectrum *sintec, struct spectrum response, struct spectrum spec, float ew_pix, float lpix, float fwhm, float K, float minselldo, float maxselldo, struct disper_prism DP) {


  float sigma;
  int i;
  float sumsintec,sumresp,sumspec;
  float lowpix, uppix;

  lowpix=pr_ldo2pix(minselldo,DP);
  uppix=pr_ldo2pix(maxselldo,DP);

  if(DEBUG) printf(" Los límites son lowpix %f uppix %f\n",lowpix,uppix);

  CopySpec(sintec,response);
  sigma=fwhm/2.35;
  for(i=0;i<response.nx;i++) {
    (*sintec).spec[i]=response.spec[i]*(1+K*i+ew_pix*gaussian((float)i,lpix,sigma));
  }
  sumresp  =StSuma1((int)uppix-(int)lowpix+1,response.spec+(int)lowpix,1);
  sumsintec=StSuma1((int)uppix-(int)lowpix+1,(*sintec).spec+(int)lowpix,1);
  sumspec  =StSuma1((int)uppix-(int)lowpix+1,spec.spec+(int)lowpix,1);
  for(i=0;i<response.nx;i++) (*sintec).spec[i]*=sumspec/sumsintec;

}



void LoadParam_kbd(struct param *P) {

  printf(" Input file with candidates \n");
  reads("",(*P).elgfile);



}

void LoadParam_file(struct param *P, char file[100]) {
  int status=0;
  char comment[51];
  fitsfile *parfile;
  
  printf("Using parameter file <<%s>>\n",file);
  if( ffopen2(&parfile,file, READONLY, &status)) fits_report_error(stderr,status);
  status=0;

  ffgky(parfile,TSTRING,"ELGFILE",(*P).elgfile,comment,&status);
  ffgky(parfile,TSTRING,"IMAGE",(*P).imafile,comment,&status);
  ffgky(parfile,TSTRING,"RESPFILE",(*P).respfile,comment,&status);
  ffgky(parfile,TSTRING,"SPECFILE",(*P).specfile,comment,&status);
/*   ffgky(parfile,TSTRING,"FILEREG",(*P).regfile,comment,&status); */
  ffgky(parfile,TSTRING,"CANDFILE",(*P).candidatefile,comment,&status);
  ffgky(parfile,TSTRING,"NELGFILE",(*P).outelgfile,comment,&status);
  ffgky(parfile,TSTRING,"REPPREF",(*P).reportprefix,comment,&status);
  ffgky(parfile,TFLOAT,"A_DISP",&((*P).DP.A),comment,&status);
  ffgky(parfile,TFLOAT,"B_DISP",&((*P).DP.B),comment,&status);
  ffgky(parfile,TFLOAT,"C_DISP",&((*P).DP.C),comment,&status);
  ffgky(parfile,TFLOAT,"PIXSIZE",&((*P).DP.tampix),comment,&status);
  ffgky(parfile,TFLOAT,"MINLDO",&((*P).minselldo),comment,&status);
  ffgky(parfile,TFLOAT,"MAXLDO",&((*P).maxselldo),comment,&status);
  ffgky(parfile,TFLOAT,"SEEING",&((*P).seeing),comment,&status);
  ffgky(parfile,TFLOAT,"SECPPIX",&((*P).platescale),comment,&status);
  ffgky(parfile,TDOUBLE,"H0",&(((*P).cosmo).H0),comment,&status);
  ffgky(parfile,TDOUBLE,"Q0",&(((*P).cosmo).q0),comment,&status);
  ffgky(parfile,TFLOAT,"COCCONT",&((*P).Kcoc),comment,&status);
  ffgky(parfile,TFLOAT,"NIIFACT",&((*P).niifact),comment,&status); 
  ffgky(parfile,TFLOAT,"MEANEBV",&((*P).meanebv),comment,&status); 
//  ffgky(parfile,TSTRING,"FILEPROP",(*P).fileprop,comment,&status); 
//  ffgky(parfile,TINT,"MABSCOL",&((*P).Mabscol),comment,&status); 
//  ffgky(parfile,TINT,"EWCOL",&((*P).ewcol),comment,&status); 
//  ffgky(parfile,TINT,"NIICOL",&((*P).niicol),comment,&status); 
//  ffgky(parfile,TINT,"EBVCOL",&((*P).ebvcol),comment,&status); 
  ffgky(parfile,TDOUBLE,"JD",&((*P).JD),comment,&status); 
  ffgky(parfile,TINT,"INTERACT",&((*P).interactflag),comment,&status); 
  ffgky(parfile,TINT,"UACFLAG",&((*P).uacflag),comment,&status); 

  printf(" ahora ho %f\n",(*P).cosmo.H0);
  printf(" nii %f\n",(*P).niifact);
  printf(" inter %d uac %d\n",((*P).interactflag),((*P).uacflag));
    
  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file.\nERROR: Not enough keywords or badformed parameter file\n");
    exit(1);
  }
  fits_close_file(parfile,&status);
}



void DoAllSpec(struct param P, struct elgdata *ELG,  struct spectrum response) {

  float ewpix,lpix;

  if(DEBUG) printf(" spectrum %s\n",(*ELG).spectrum);
  if(DEBUG) printf(" errspectrum %s\n",(*ELG).spectrum);
  
  strcpy((*ELG).spec.file,(*ELG).spectrum);
  strcpy((*ELG).errspec.file,(*ELG).errspectrum);

  printf(" Despues strc\n");


  ReadSpec(&((*ELG).spec));
  ReadSpec(&((*ELG).errspec));

  printf(" Despue read\n");

  ewpix=(*ELG).ewobs*pr_dpdl((*ELG).ldo,P.DP);
  lpix=pr_ldo2pix((*ELG).ldo,P.DP);

  CreateSintec(&((*ELG).fitspec),response,(*ELG).spec,ewpix,lpix,P.seeing/P.platescale,(*ELG).Klin,P.minselldo,P.maxselldo,P.DP); 

  printf(" Depsues create\n");

}


float Vannual(double ra_d,double dec_d, double JD) {
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


float MeanValue(float *Mabs, float *ew, float *value, int nsample, float Mabsgal, float errMabsgal, float ewgal, float errewgal) {

  int nminmean=10;
  int nused;
  float ewmin,ewmax;
  float Mabsmin,Mabsmax;
  float ewbox;
  float Mabsbox;
  int i;
  float sigma;
  float mean;
 
  float *values;

  values=vector_f(nminmean);

  MinMax(nsample,Mabs,&Mabsmin,&Mabsmax);
  MinMax(nsample,ew,&ewmin,&ewmax);

  Mabsbox=maxf(errMabsgal,fabs((Mabsmax-Mabsmin)/(float)nsample*sqrt(nminmean)));
  ewbox  =maxf(errewgal,  fabs((  ewmax-  ewmin)/(float)nsample*sqrt(nminmean)));


  while(nused<nminmean) {
    nused=0;
    for(i=0;i<nsample;i++) {
      if(fabs(Mabs[i]-Mabsgal)<Mabsbox && fabs(ew[i]-ewgal)<ewbox) {
        values[nused]=value[i];
        nused++;
      }
    }
    Mabsbox*=1.1;
    ewbox*=  1.1;
  }
  mean=StMedia(nused,values,&sigma);
  free(values);
  return(mean);
}
