#include "modulos.h"

#define DEBUG 1

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

  /* Output parameters */
  char candidatefile[51];
  char outelgfile[51];

  int interactflag;
  
};

struct elgdata {
  char spectrum[101];
  char errspectrum[101];
  float ew;
  float errew;
  float ldo;;
  float errldo;
  float z,errz;
  float Klin;
  float x_ima,y_ima;
  double ra,dec;
  int   selected;
  float mag,errmag;
  float Mabs,errMabs;
  double lineflux;
  double linelum;
  float mag_usnoB;
  float mag_usnoR;
  float Mabs_usnoB;
  float Mabs_usnoR;
  char  ucyname[100];
};


void LoadParam_file(struct param *P, char file[100]);
void LoadParam_kbd(struct param *P);
void SaveParam(struct param P, char file[100]);

void ReadResp(struct spectrum *response, char respfile[]);
void ReadELGData(struct elgdata *ELG, int *nelg, char elgfile[51]);
void DoAllSpec(struct param P,struct elgdata ELG, struct spectrum *spec, struct spectrum *errspec,struct spectrum *fitspec,struct spectrum response);
void CreateSintec(struct spectrum *sintec, struct spectrum response, float ew_pix, float lpix, float fwhm, float K, float minselldo, float maxselldo, struct disper_prism DP);
void ParseELG(struct param P,struct elgdata *ELG,int nelg,struct wcsimage *ima,struct image *errima,struct spectrum response);
void ShowReport(char pgdevice[100], struct elgdata E, struct spectrum spec, struct spectrum errspec, struct spectrum fitspec, struct wcsimage *ima);
int main(int argc, char **argv) {

  struct param P;
  struct wcsimage ima;
  struct image errima;
  struct spectrum response;
  struct elgdata *ELG=NULL;
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
  strcpy(errima.file,P.errimafile);
  ReadImage(&errima);
  ReadResp(&response,P.respfile);
  ReadELGData(ELG,&nelg,P.elgfile);

  ParseELG(P,ELG,nelg,&ima,&errima,response);

  return 0;
}

void ParseELG(struct param P,struct elgdata *ELG,int nelg,struct wcsimage *ima,struct image *errima,struct spectrum response) {

  int i;
  struct spectrum spec;
  struct spectrum errspec;
  struct spectrum fitspec;
  char pginter[100]="?";
  char pgreport[100];

  spec.aloc_flag=0;spec.alocldo_flag=0;
  errspec.aloc_flag=0;errspec.alocldo_flag=0;
  fitspec.aloc_flag=0;fitspec.alocldo_flag=0;

  for(i=0;i<nelg;i++) {
    if(P.interactflag) {
      DoAllSpec(P,ELG[i],&spec,&errspec,&fitspec,response);
      ShowReport(pginter,ELG[i],spec,errspec,fitspec,ima);
      ELG[i].selected=readi(1);
      if(ELG[i].selected) {
	sprintf(pgreport,"%s%s.ps/cps",P.reportprefix,ELG[i].ucyname);
	ShowReport(pgreport,ELG[i],spec,errspec,fitspec,ima);
      }
    }
  }

  CloseSpec(&spec);
  CloseSpec(&errspec);
  CloseSpec(&fitspec);


}

void ShowReport(char pgdevice[100], struct elgdata E, struct spectrum spec, struct spectrum errspec, struct spectrum fitspec, struct wcsimage *ima) {

  char snul[300];
  int nsubx=1,nsuby=1;
  struct wcsimage dssimage;
  float x1,x2,y1,y2;
  double ra1,ra2,dec1,dec2;
  char rastr[32],decstr[32];
  double dnul;

  cpgscf(2);
  cpgsch(1.7);
  cpgvstd();
  

  cpgmtxt("T",0.,0.5,0.5,E.ucyname);
  cpgsubp(2,2);
  /* Draw Spectrum */
  cpgpanl(1,1);
  PlotSpec_pix_err(spec,errspec);
  PlotSpec_pix_ov(fitspec);

  /* Display Info */
  cpgpanl(1,2);
  sprintf(snul," EW = %f +/- %f\n",E.ew,E.errew);
  cpgmtxt("T",.5,0.1,0.0,snul);
  sprintf(snul," z = %f +/- %f\n",E.z,E.errz);
  cpgmtxt("T",1.5,0.1,0.0,snul);
  sprintf(snul,"mag\\dcal\\u = %f +/- %f\n",E.mag,E.errmag);
  cpgmtxt("T",2.5,0.1,0.0,snul);
  sprintf(snul," B\\dUSNO\\u = %f\n",E.mag_usnoB);
  cpgmtxt("T",3.5,0.1,0.0,snul);
  sprintf(snul," R\\dUSNO\\u = %f\n",E.mag_usnoR);
  cpgmtxt("T",4.5,0.1,0.0,snul);
  ra2str(rastr,32,E.ra,3);
  sprintf(snul," RA\\dJ2000\\u = %s \n",rastr);
  dec2str(decstr,32,E.ra,2);
  sprintf(snul," \\gd\\dJ2000\\u = %s \n",decstr);
  cpgmtxt("T",5.5,0.1,0.0,snul);

  /* Dibujar trozo imagen */
  x1=(int)(E.x_ima-nsubx/2);
  x2=(int)(E.x_ima+nsubx/2);
  y1=(int)(E.y_ima-nsuby/2);
  y2=(int)(E.y_ima+nsuby/2);

  if(x1<0) {
    x1=0;x2=nsubx;
  }
  if(y1<0) {
    y1=0;y2=nsuby;
  }
  if(x1>(*ima).image.naxes[0]-1) {
    x1=(*ima).image.naxes[0]-1-nsubx;x2=(*ima).image.naxes[0]-1;
  }
  if(y1>(*ima).image.naxes[1]-1) {
    y1=(*ima).image.naxes[1]-1-nsuby;y2=(*ima).image.naxes[1]-1;
  }
  PlotImage_sec(&((*ima).image),x1,x2,y1,y2);
  wcspix((*ima).wcs,0.,(*ima).image.naxes[1]/2.,&ra1,&dnul);
  wcspix((*ima).wcs,(*ima).image.naxes[1],(*ima).image.naxes[1]/2.,&ra2,&dnul);
  wcspix((*ima).wcs,(*ima).image.naxes[1]/2.,0.,&dnul,&dec1);
  wcspix((*ima).wcs,(*ima).image.naxes[1]/2.,(*ima).image.naxes[1],&dnul,&dec2);
  cpgswin((float)ra1,(float)ra2,(float)dec1,(float)dec2);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);

  /* Dibujar DSS */

  retrieveDSS(&dssimage,E.ra,E.dec,10.,10.);
  PlotImage(&dssimage.image);
  
  if(dssimage.image.aloc_flag) free(dssimage.image.array);
  if(dssimage.wcsset) free(dssimage.wcs);

}

void WriteCandidateFile(struct elgdata *ELG, int nelg, char candidatefile[51], char outelgfile[51] ) {

  int i;
  int newelgflag=0;

  FILE *fcand,*felg;

  if(strcmp(outelgfile,"NONE"))    newelgflag=1;
  if((fcand=fopen(candidatefile,"w"))==NULL) {
    printf(" I cannot open file %s for writing\n",candidatefile);
    exit(1);
  }
  if((felg=fopen(outelgfile,"w"))==NULL) {
    printf(" I cannot open file %s for writing\n",candidatefile);
    exit(1);
  }

  fprintf(fcand,"#%-20s  %-32s  %-32s","UCYname","RA_J2000","Dec_J2000");

  fprintf(felg,"# Output candidate file from SearchLine\n");
/*   fprintf(felg,"# Original spectra file: %s\n",specfile); */
/*   fprintf(fc,"# Selection cuts: EW %f - %f   LDO %f -%f\n",minselew,maxselew,minselldo,maxselldo); */


  for(i=0;i<nelg;i++) {
    if(ELG[i].selected) {
      fprintf(fcand," ");
      if(newelgflag) fprintf(felg," ");
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



void ReadELGData(struct elgdata *ELG, int *nelg, char elgfile[51]) {
  
  int nlines;
  int i;

  char *spec;
  char *errspec;
  float *ew;
  float *errew;
  float *ldo;
  float *errldo;
  float *K;
  int *ilog;
  int speccol=1,errspeccol=2,ewcol=3,errewcol=4,ldocol=5,errldocol=6,Kcol=7;
  
  
  nlines=FileNLin(elgfile);

  ELG=malloc(nlines*sizeof(struct elgdata));
  spec=malloc(nlines*101*sizeof(char));
  errspec=malloc(nlines*101*sizeof(char));
  ew=vector_f(nlines);
  errew=vector_f(nlines);
  ldo=vector_f(nlines);
  errldo=vector_f(nlines);
  K=vector_f(nlines);
  ilog=vector_i(nlines);

  for(i=0;i<nlines;i++) ilog[i]=0;

  ReadNumcol(elgfile,ewcol,ew,ilog,&nlines);
  ReadNumcol(elgfile,errewcol,errew,ilog,&nlines);
  ReadNumcol(elgfile,ldocol,ldo,ilog,&nlines);
  ReadNumcol(elgfile,errldocol,errldo,ilog,&nlines);
  ReadNumcol(elgfile,Kcol,K,ilog,&nlines);
  ReadCharcol(elgfile,speccol,spec,ilog,101,&nlines);
  ReadCharcol(elgfile,errspeccol,errspec,ilog,101,&nlines);

  *nelg=0;
  for(i=0;i<nlines;i++) {
    if(ilog[i]) {
      ELG[*nelg].ew=ew[i];
      ELG[*nelg].errew=errew[i];
      ELG[*nelg].ldo=ldo[i];
      ELG[*nelg].errldo=errldo[i];
      ELG[*nelg].Klin=K[i];
      strcpy(ELG[*nelg].spectrum,spec+i*101);
      strcpy(ELG[*nelg].errspectrum,errspec+i*101);
      nelg++;
    }
  }
  
  
  free(spec);
  free(errspec);
  free(ew);
  free(errew);
  free(ldo);
  free(errldo);
  free(K);
  free(ilog);
}

void CreateSintec(struct spectrum *sintec, struct spectrum response, float ew_pix, float lpix, float fwhm, float K, float minselldo, float maxselldo, struct disper_prism DP) {


  float sigma;
  int i;
  float sumsintec,sumresp;
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
  for(i=0;i<response.nx;i++) (*sintec).spec[i]*=sumresp/sumsintec;

}



void GetDataFromSpec() {



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
  printf(" Por aqui %d \n",status);
  ffgky(parfile,TFLOAT,"MINLDO",&((*P).minselldo),comment,&status);
  ffgky(parfile,TFLOAT,"MAXLDO",&((*P).maxselldo),comment,&status);
  ffgky(parfile,TFLOAT,"SEEING",&((*P).seeing),comment,&status);
  ffgky(parfile,TFLOAT,"SECPPIX",&((*P).platescale),comment,&status);
  
  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file.\nERROR: Not enough keywords or badformed parameter file\n");
    exit(1);
  }
  fits_close_file(parfile,&status);
}



void DoAllSpec(struct param P, struct elgdata ELG, struct spectrum *spec, struct spectrum *errspec,struct spectrum *fitspec, struct spectrum response) {

  float ewpix,lpix;
  
  strcpy((*spec).file,ELG.spectrum);
  strcpy((*errspec).file,ELG.errspectrum);

  ReadSpec(spec);
  ReadSpec(errspec);

  ewpix=ELG.ew*pr_dpdl(ELG.ldo,P.DP);
  lpix=pr_ldo2pix(ELG.ldo,P.DP);

  CreateSintec(fitspec,response,ewpix,lpix,P.seeing/P.platescale,ELG.Klin,P.minselldo,P.maxselldo,P.DP);

}
