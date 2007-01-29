#include "modulos.h"
#define FTOL 1e-5
#define MAXITER 500
#define DEBUG 0
  
struct Polresp
{
  float p0;
  float p1;
  float p2;
  float p3;
  float p4;
  float ldoeff;
};


struct theo_resp
{
  float *ldo;
  float *y;
  int n;
};


/* Parametros del programa */
char respfile[101];
char errrespfile[101];
char calrespfile[101];
char errcalrespfile[101];
char disperfile[101];
char combrespfile[101];

char transfile[101];
int  ldocol;
int  transcol;

int flip;
float seeing;
float platescale;

/* Variables comunes */ 
struct spectrum resp;
struct spectrum errresp;
struct spectrum calresp;
struct spectrum errcalresp;
struct disper_prism DP;
struct theo_resp TR;
struct Polresp PR;
/* struct spectrum tr; */
struct spectrum tr_s;
/* struct spectrum tr_s_l; */

float pixip;
float fluxtr;
float fluxresp;

void LoadParam_kbd(void);
void ReadAll(void);
void ReadTheoResp( char file[100]);
void GuessInitDisper(float *A, float *B);
void ComputeDispersion(void);
void ComputeTR(void);
void SaveCalDisper(void);
void SaveCombResp(void);
float func_sp(int n, float *x, float *y, float *p);
float func_qe(int n, float *x, float *y, float *p);
float func_all(int n, float *x, float *y, float *p);
void funpoly(double x,double *p,int n);
float Chi2();
float TResp(float ldo);
float func_conv(float x);
float ldoeff();
float Polr(float ldo);

int main(void) {


  LoadParam_kbd();
  ReadAll();
  cpgpanl(1,1);
  PlotSpec(resp);
  printf(" Trying to compute dispersion using that data ");
  cpgpanl(2,1);
  ComputeDispersion();
  SaveCalDisper();
  SaveCombResp();

  exit(0);
}



void LoadParam_kbd(void) {

  printf(" FITS file with non-calibrated response ");
  reads("",respfile);
  printf(" FITS file with error in non-calibrated response ");
  reads("",errrespfile);
  printf(" File with theorical response ");
  reads("",transfile);
  printf(" Input column with ldo ");
  ldocol=readi(1);
  printf(" Input column with transmission ");
  transcol=readi(2);

  printf(" Input seeing ");
  seeing=readf(2);
  printf(" Input pixel size (microns) ");
  DP.tampix=readf(24.);
  printf(" Input plate scale (\"/pix) ");
  platescale=readf(1.);

  printf(" Is response flipped (0=no, 1=yes) ? ");
  flip=readi(0);
   
  printf(" Input parameter C for dispersion ");
  DP.C=readf(1606.);

  printf(" Output FITS file with calibrated response ");
  reads("",calrespfile);
  printf(" Output FITS file with error in calibrated response ");
  reads("",errcalrespfile);
  printf(" Output text file with dispersion solution ");
  reads("",disperfile);
  printf(" Output text file with combined response of filter + system ");
  reads("",combrespfile);


  cpgopen("?");
  cpgsubp(2,1);

}



void GuessInitDisper(float *A, float *B) {

  *B = DP.tampix*(1.-resp.nx)/(1./(TR.ldo[0]-DP.C)- 1./(TR.ldo[TR.n-1]-DP.C));

  *A = 1.*DP.tampix - *B / (TR.ldo[0]-DP.C);
  
  
}

void ComputeDispersion(void) {


  float par_sp[2];   /* lo dejo con dos, para minimizar A y B */
  float sigpar_sp[2];
  float par_qe[5];   
  float sigpar_qe[5];
  float par_all[7];   
  float sigpar_all[7];
  struct spectrum coccalresp;
  double *xfit,*yfit,*sigyfit;
  struct spectrum tmpldo1,tmpldo2;
  float min1,max1,min2,max2;
  int i;
  int pix_lowcut;
  int pix_upcut;
  double coefpol[8];
  double **covpar;
  double chi2;
   
  coccalresp.ldo=vector_f(resp.nx); 
  coccalresp.spec=vector_f(resp.nx);
  coccalresp.aloc_flag=1;
  coccalresp.alocldo_flag=1;
  coccalresp.nx=resp.nx;
  coccalresp.npixels=resp.nx;
  coccalresp.datamin=0; 
  coccalresp.datamax=0;

  GuessInitDisper(&(par_sp[0]),&(par_sp[1]));


  if(flip) {
    par_sp[0]=-1950;    /* Estos son los de la solucion de dispersion */
    par_sp[1]=+139e5;
  }
  else     {
    par_sp[0]=3089;    /* Estos son los de la solucion de dispersion */
    par_sp[1]=-139e5;
  }
  sigpar_sp[0]=par_sp[0]/10.;
  sigpar_sp[1]=par_sp[1]/10.;

  printf(" Computing dispersion...\n");
  Amoeba(resp.nx, resp.spec, errresp.spec, 2, par_sp, sigpar_sp, FTOL, MAXITER, func_sp); 
  DP.A=par_sp[0];
  DP.B=par_sp[1];

  tmpldo1.alocldo_flag=0;
  tmpldo1.aloc_flag=0;
  SP_pix2ldo(&tmpldo1,resp,DP);
  tmpldo2.alocldo_flag=0;
  tmpldo2.aloc_flag=0;
  SP_pix2ldo(&tmpldo2,tr_s,DP);
  coccalresp.ldomin=tmpldo1.ldomin;
  coccalresp.deltaldo=tmpldo1.deltaldo;
  MinMax(tmpldo1.nx,tmpldo1.spec,&min1,&max1);
  MinMax(tmpldo2.nx,tmpldo2.spec,&min2,&max2);
  for(i=0;i<coccalresp.nx;i++) {
    if(tmpldo1.spec[i]==0) coccalresp.spec[i]=0;
    else                   coccalresp.spec[i]=tmpldo1.spec[i]/tmpldo2.spec[i]*max2/max1;
    coccalresp.ldo[i]=tmpldo1.ldo[i];
  }
  pix_lowcut=(int)pr_ldo2pix(6600,DP);
  pix_upcut=(int)pr_ldo2pix(7200,DP);
  pix_lowcut=21;
  pix_upcut=42;
  printf(" low %d ip %d\n",pix_lowcut,pix_upcut);
  xfit=vector_d(pix_upcut-pix_lowcut+1);
  yfit=vector_d(pix_upcut-pix_lowcut+1);
  sigyfit=vector_d(pix_upcut-pix_lowcut+1);
  for(i=pix_lowcut;i< pix_upcut+1;i++) {
    xfit[i-pix_lowcut]=(float)i-30.;
/*     xfit[i-pix_lowcut]=tmpldo1.ldo[i]; */
    yfit[i-pix_lowcut]=coccalresp.spec[i];
    sigyfit[i-pix_lowcut]=coccalresp.spec[i]/10.;
  }

      

  covpar=matrix_d(8,8);
  svdfit_d(xfit,yfit,sigyfit,pix_upcut-pix_lowcut+1,coefpol,5,covpar,&chi2,&funpoly);  
  PR.ldoeff=tmpldo1.ldo[30];

  printf(" Ajuesta SVD p0 %f p1 %f p2 %f p3 %f p4 %f\n",coefpol[0],coefpol[1],coefpol[2],coefpol[3],coefpol[4]);   
  cpgpanl(2,1);
  cpgswin((float)pix_lowcut-30,(float) pix_upcut-30, 0.,2.);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgline_d(pix_upcut-pix_lowcut+1,xfit,yfit);
  cpgsci(2);
  for(i=0;i< pix_upcut-pix_lowcut+1;i++) {
    cpgdraw(xfit[i],poly_d(xfit[i],4,coefpol));
  }
  cpgsci(1);
  
  if(DEBUG) {
    
    printf(" Ajuesta SVD p0 %f p1 %f p2 %f p3 %f\n",coefpol[0],coefpol[1],coefpol[2],coefpol[3]); 
    printf(" Chi2 %f\n",chi2);
    /*   MCPN(pix_upcut-pix_lowcut+1,coccalresp.ldo+pix_lowcut,coccalresp.spec+pix_lowcut,3,coefpol); */
    cpgpanl(2,1);
    /*   PlotSpec(coccalresp); */
    cpgswin((float)pix_lowcut-30,(float) pix_upcut-30, 0.,2.);
    /*   cpgswin(6500.,7300., 0.,2.); */
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
    cpgline_d(pix_upcut-pix_lowcut+1,xfit,yfit);
    
    cpgsci(2);
    for(i=0;i< pix_upcut-pix_lowcut+1;i++) {
      cpgdraw(xfit[i],poly_d(xfit[i],4,coefpol));
      /*     printf("x %f y %f\n",xfit[i],poly_d(xfit[i],2,coefpol)); */
    }
    for(i=0;i< 1000000000;i++) ;
    cpgpage();
    cpgsci(1);
    cpgswin(6500.,7300., 0.,2.); 
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
    for(i=0;i< pix_upcut-pix_lowcut+1;i++) xfit[i]=tmpldo1.ldo[i+pix_lowcut];
    cpgline_d(pix_upcut-pix_lowcut+1,xfit,yfit);
    coefpol[1]=coefpol[1]/tmpldo1.deltaldo;
    coefpol[2]=coefpol[2]/tmpldo1.deltaldo/tmpldo1.deltaldo;
    coefpol[3]=coefpol[3]/tmpldo1.deltaldo/tmpldo1.deltaldo/tmpldo1.deltaldo;
    coefpol[4]=coefpol[4]/tmpldo1.deltaldo/tmpldo1.deltaldo/tmpldo1.deltaldo/tmpldo1.deltaldo;
    printf(" delta %f\n",tmpldo1.deltaldo);
    printf(" ldoeff %f\n",PR.ldoeff);  
    cpgsci(2);
    for(i=0;i< pix_upcut-pix_lowcut+1;i++) {
      cpgdraw(xfit[i],poly_d(xfit[i]-PR.ldoeff,4,coefpol));
      printf("x %f y %f\n",xfit[i],poly_d(xfit[i]-PR.ldoeff,4,coefpol));
    }
    cpgpanl(2,1);
    for(i=0;i<tmpldo2.nx;i++) tmpldo2.spec[i]=tmpldo2.spec[i]/max2*max1;
    cpgsci(1);
    PlotSpec(tmpldo1);
    PlotSpec_ov(tmpldo2);
    exit(1);
    i=readi(i);
    cpgpanl(1,1);
    cpgsci(1);
    PlotSpec(resp);
    cpgpanl(2,1);
  }

  printf(" Computing global response...\n");
  par_qe[0]=coefpol[0];  /* Estos son los coeficientes de un polinomio que ajuste la respues total */
  par_qe[1]=coefpol[1]/tmpldo1.deltaldo;
  par_qe[2]=coefpol[2]/tmpldo1.deltaldo/tmpldo1.deltaldo;
  par_qe[3]=coefpol[3]/tmpldo1.deltaldo/tmpldo1.deltaldo/tmpldo1.deltaldo;
  par_qe[4]=coefpol[4]/tmpldo1.deltaldo/tmpldo1.deltaldo/tmpldo1.deltaldo/tmpldo1.deltaldo;
  sigpar_qe[0]=par_qe[0]/30.;
  sigpar_qe[1]=par_qe[1]/20.;
  sigpar_qe[2]=par_qe[2]/20.;
  sigpar_qe[3]=par_qe[3]/20.;
  sigpar_qe[4]=par_qe[4]/20.;
  Amoeba(resp.nx, resp.spec, errresp.spec, 5, par_qe, sigpar_qe, FTOL, MAXITER, func_qe);  

  printf(" Recomputing everything...\n");
  par_all[0]=par_sp[0];  /* Estos son los coeficientes de un polinomio que ajuste la respues total */
  par_all[1]=par_sp[1];
  par_all[2]=par_qe[0];
  par_all[3]=par_qe[1];
  par_all[4]=par_qe[2];
  par_all[5]=par_qe[3];
  par_all[6]=par_qe[4];
  sigpar_all[0]=par_sp[0]/40.;
  sigpar_all[1]=par_sp[1]/40.;
  sigpar_all[2]=par_qe[0]/100.;
  sigpar_all[3]=par_qe[1]/10.;
  sigpar_all[4]=par_qe[2]/10.;
  sigpar_all[5]=par_qe[3]/10.;
  sigpar_all[6]=par_qe[4]/10.;
  Amoeba(resp.nx, resp.spec, errresp.spec, 7, par_all, sigpar_all, FTOL, MAXITER, func_all);  


  printf(" Salgo anoemba\n");


}


void SaveCalDisper(void) {
  
  int i;
  float std;
/*   char inter='y'; */
  float fnul;
  char cnul;

  FILE *fp;
  
  printf(" Final solution: A %f B %f C %f\n",DP.A,DP.B,DP.C);
  
  ComputeTR();
  fluxresp=resp.nx*StMedia(resp.nx,resp.spec,&std);
  fluxtr  =tr_s.nx*StMedia(tr_s.nx,tr_s.spec,&std);
  
  for(i=0;i<tr_s.nx;i++) tr_s.spec[i]=tr_s.spec[i]*fluxresp/fluxtr;
  cpgsubp(1,1);
  PlotSpec(resp);
  PlotSpec_ov(tr_s);
  while(!(cnul=='m') && 0) {
    ComputeTR();
    fluxresp=resp.nx*StMedia(resp.nx,resp.spec,&std);
    fluxtr  =tr_s.nx*StMedia(tr_s.nx,tr_s.spec,&std);
    for(i=0;i<tr_s.nx;i++) tr_s.spec[i]=tr_s.spec[i]*fluxresp/fluxtr;
    PlotSpec(resp);
    PlotSpec_ov(tr_s);
    cpgcurs(&fnul,&fnul,&cnul); 
    if(cnul=='A') PR.p2*=1.03;
    if(cnul=='X') PR.p2*=0.97;
    printf("  p0 %f p1 %f p2 %f p3 %f p4 %f\n",PR.p0,PR.p1,PR.p2,PR.p3,PR.p4);
  }

/*   printf(" Are you satisfied? "); */
/*   inter=readc(inter); */
/*   if(inter!='n') exit(1); */
               
  printf(" Despue spintar\n");

  calresp.alocldo_flag=0;
  calresp.aloc_flag=0;
  errcalresp.alocldo_flag=0;
  errcalresp.aloc_flag=0;
  SP_pix2ldo(&calresp,resp,DP);
  printf(" Transformada la primera\n");
  SP_pix2ldo(&errcalresp,errresp,DP);  
  printf(" Aqui chungo\n");
  
  strcpy(calresp.file,calrespfile);
  strcpy(errcalresp.file,errcalrespfile);

  printf(" ANtes salvar\n");

  SaveSpec(calresp);
  SaveSpec(errcalresp);

  if((fp=fopen(disperfile,"w"))==NULL) {
    printf("ERROR: Can't open output file %s\n",disperfile);
    exit(1);
  }
  
  fprintf(fp,"# Dispersion solution for response %s (seeing = %f)\n# using theorical file %s\n",respfile,seeing,transfile);
  fprintf(fp,"#  A      B     C     pixsize\n");
  fprintf(fp," %f  %f  %f  %f\n",DP.A,DP.B,DP.C,DP.tampix);
  fclose(fp);

}


float func_sp(int n, float *x, float *y, float *p) {
  
  int npix;
  float chi2;
  static int iter=0;

  iter++;
  if(DEBUG) printf(" Estoy dentro de func_sp\n");
  
  npix=n;
  
  DP.A=p[0];
  DP.B=p[1];
  PR.p0=1;
  PR.p1=0;
  PR.p2=0;
  PR.p3=0;
  PR.p4=0;
  ComputeTR();
  if(DEBUG) printf(" Despues CPTRM\n");
  cpgpanl(2,1);
  if(DEBUG) PlotSpec(tr_s);
  if(DEBUG) printf(" Antes chi2\n");
  chi2 = Chi2();
  if(DEBUG)  printf(" A %f B %f  p1 %f p2 %f chi2 %f\n",DP.A,DP.B,PR.p1,PR.p2,chi2); 

/*   chi2=readf(chi2); */
  return(chi2);
}

float func_qe(int n, float *x, float *y, float *p) {
  
  int npix;
  float chi2;
  static int iter=0;

  iter++;
  if(DEBUG) printf(" Estoy dentro de func_sp\n");
  
  npix=n;
  PR.p0=p[0];
  PR.p1=p[1];
  PR.p2=p[2];
  PR.p3=p[3];
  PR.p4=p[4];

  ComputeTR();
  if(DEBUG) printf(" Despues CPTRM\n");
  if(DEBUG) cpgpanl(2,1);
  if(DEBUG) PlotSpec(tr_s);
  if(DEBUG) printf(" Antes chi2\n");
  chi2 = Chi2();
  if(DEBUG)  printf(" A %f B %f  p1 %f p2 %f chi2 %f\n",DP.A,DP.B,PR.p1,PR.p2,chi2); 

/*   chi2=readf(chi2);   */
  return(chi2);
}

float func_all(int n, float *x, float *y, float *p) {
  
  int npix;
  float chi2;
  static int iter=0;

  iter++;
  if(DEBUG) printf(" Estoy dentro de func_sp\n");
  
  npix=n;
  DP.A=p[0];
  DP.B=p[1];
  PR.p0=p[2];
  PR.p1=p[3];
  PR.p2=p[4];
  PR.p3=p[5];
  PR.p4=p[6];
  ComputeTR();
  if(DEBUG) printf(" Despues CPTRM\n");
  cpgpanl(2,1);
  if(DEBUG) PlotSpec(tr_s);
  if(DEBUG) printf(" Antes chi2\n");
  chi2 = Chi2();
  if(DEBUG)  printf(" A %f B %f  p1 %f p2 %f chi2 %f\n",DP.A,DP.B,PR.p1,PR.p2,chi2); 

/*   chi2=readf(chi2);  */
  return(chi2);
}



float Chi2() {
  
  int i;
  double chi2=0;
  float std;
  float w;
  fluxresp=resp.nx*StMedia(resp.nx,resp.spec,&std);
  fluxtr  =tr_s.nx*StMedia(tr_s.nx,tr_s.spec,&std);
  for(i=0;i<resp.nx;i++) {
    w=resp.spec[i]/fluxresp;
    w=1;   
    if(errresp.spec[i]==0);
    else    chi2+=w*w*((tr_s.spec[i]/fluxtr-resp.spec[i]/fluxresp)/errresp.spec[i]*fluxresp)*((tr_s.spec[i]/fluxtr-resp.spec[i]/fluxresp)/errresp.spec[i]*fluxresp);
/*     printf(" i %d sp %f sr %f esr %f  chi %f\n",i,tr_s.spec[i]/fluxtr,resp.spec[i]/fluxresp,errresp.spec[i],chi2);   */
  }
  return((float)chi2);
}


void ReadTheoResp( char file[100]) {
  
  int i;       
  
  float *ldo,*y;
  int *ilog;
  int n;
  
  n=FileNLin(file);

  TR.ldo=vector_f(n);
  TR.y=vector_f(n);
  ldo=vector_f(n);
  y=vector_f(n);
  ilog=vector_i(n);
  
  ReadNumcol(file,ldocol,ldo,ilog,&n);
  ReadNumcol(file,transcol,y,ilog,&n);

  TR.n=0;
  for (i=0;i<n;i++) {
    if(ilog[i]) {
      TR.ldo[TR.n]=ldo[i];
      TR.y[TR.n]=y[i];
      TR.n++;
    }
  }
  
  free(ldo);
  free(y);
  free(ilog);
}


void ReadAll(void) {
  
  strcpy(resp.file,respfile);
  ReadSpec(&resp);
  strcpy(errresp.file,errrespfile);
  ReadSpec(&errresp);

  ReadTheoResp(transfile);

  PR.ldoeff=ldoeff();
}


void ComputeTR(void) {
  int i,j;
  int nip=2;
  struct spectrum *s;
  float pix1,pix2;
  float conv;
  float intpix;

  if(DEBUG) printf(" Dentro de Comptetr resp.nx %d\n",resp.nx);

  s=&tr_s;

  if(s->aloc_flag==1) free(s->spec);
  s->spec=vector_f(resp.nx);
  if(s->alocldo_flag==1) free(s->ldo);
  s->ldo=vector_f(resp.nx);
  s->aloc_flag=1;
  s->alocldo_flag=1;
  s->nx=resp.nx;
  s->npixels=resp.nx;
  s->ldomin=1;
  s->deltaldo=1;

  if(DEBUG) printf(" Estoy por aui\n");

  
  for(i=0;i<resp.nx;i++) {
    pix1=(float)(i+1-0.5);
    pix2=(float)(i+1+0.5);
    intpix=0;
    for(j=0;j<nip;j++) {
      pixip=pix1+j*(pix2-pix1)/(nip-1);
      conv=gaussinther(func_conv,pixip,seeing/2.35/sqrt(2.),13);
      intpix+=conv;
    }
    intpix=intpix/(float)nip;  
    /* Estro de aqui abajo esta mal */
/*     intpix=TResp(pr_pix2ldo(pixip,DP))*pr_dldp(pixip,DP); */
/*     printf(" pix %f ldo %f TR %f dldp %f\n",pixip,pr_pix2ldo(pixip,DP),TResp(pr_pix2ldo(pixip,DP)),intpix); */
 
    (*s).spec[i]=intpix;
    (*s).ldo[i]=i+1.;
  }

  if(DEBUG) printf(" Salgo \n");

}



float func_conv(float x) {
  
  float pixconv;
  float ldoconv;
  float especreal;

  pixconv=x;
  ldoconv=pr_pix2ldo(pixconv,DP);
  if(flip)  especreal=-TResp(ldoconv)*Polr(ldoconv-PR.ldoeff);  
  else      especreal=+TResp(ldoconv)*Polr(ldoconv-PR.ldoeff);  
  return(gaussian(x,pixip,seeing/2.35)*especreal*pr_dldp(pixconv,DP));
}

float TResp(float ldo)
{
  return(Lagr2(TR.ldo,TR.y,TR.n,ldo));
}


float Polr(float ldo)
{
  double coef[5];
  coef[0]=PR.p0;
  coef[1]=PR.p1;
  coef[2]=PR.p2;
  coef[3]=PR.p3;
  coef[4]=PR.p4;
  return((float)poly_d(ldo,4,coef));
}

float ldoeff() 
{
  int   l; 
  float ldo; 
  double filt=0.,lfilt=0.;
  
  for(l=0;l<TR.n;l++) {
    ldo=TR.ldo[0]+l*(TR.ldo[TR.n-1]-TR.ldo[0])/(TR.n-1); 
    filt+=TResp(ldo);
    lfilt+=ldo*TResp(ldo);
  } 

  return(lfilt/filt);
}


void SaveCombResp() 
{

  int l;
  float ldo;
  int npt=1000;
  FILE *of;
  
  if((of=fopen(combrespfile,"w")) ==NULL) {
    printf(" I could not open %s \n",combrespfile);
    exit(1);
  }

  fprintf(of,"#ldo     QE\n");

  for(l=0;l<npt;l++) {
    ldo=6000+l*(8000-6000)/(npt-1); 
    fprintf(of," %10.3f  %10.5f\n",ldo,TResp(ldo)*Polr(ldo-PR.ldoeff));
  }
  fclose(of);
}

void funpoly(double x,double *p,int n) {
  int j;
  p[0]=1.0;
  for(j=1;j<n;j++) p[j]=p[j-1]*x;
}
