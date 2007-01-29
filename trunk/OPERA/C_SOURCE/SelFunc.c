#include "modulos.h" 
 
#define DEBUG 1
#define DEBUG2 1

struct objsim {
  float x,y;
  double alfa, delta;
  float ew;
  float z;
  float mag;

  float seeing;
  float sky;
  float transparency;
  char  instsetup[200];

  int selected;
  struct objdet *od; 
};

struct objdet {
  double alfa, delta;
  double x,y;
  float ew;
  float z;
  float mag;

  float seeing;
  float sky;
  float transparency;

  struct objsim *os;
};

struct simcat {
  char file[101];
  float seeing;
  float sky;
  float transparency;
  char setup[200];
  
  /* Ficticius Astrometry */
  struct WorldCoor *wcsim;
  double alfac,deltac;
  float platescale,epoch,rotang,flip;
  int nx,ny;

  int x1,x2,y1,y2;

  char filedet[101];
};



void DoPOCompute(void);
void DoDICompute(void);
void DoPOPlot(void);
void DoDIPlot(void);
void ReadCatAndFillStructFunSel_PO(struct poselfunc *SF,struct simcat *SC, int nfiles, struct objsim **OS, int *nsim, struct objdet **OD, int *ndet);
void MatchSimDet(struct objsim *OS, int nsim, struct objdet *OD, int ndet);
void GetIntervals(double *values, int nvalues, double *samplevalues, int nsample);
void GetStrings(char **values, int nvalues, char **samplevalues, int *nsample);
void ComputePOFunSel( struct poselfunc *SF,struct objsim *OS,int nsim,struct objdet *OD,int ndet);
int main() {
  
  char opt='E';
  do{
    printf("\n P Compute selection function from simulations for PO images\n"); 
    printf(" D Compute selection function from simulations for direct images\n"); 
    printf(" C Cleaning and smoothing PO selection function\n"); 
    printf(" S Plot PO selection function\n"); 
    printf(" V Plot direct selection function\n"); 
    printf(" E Exit\n");
    opt=readc(opt); 
    switch (opt) { 
    case 'P':
    case 'p':
      DoPOCompute();
      break;
    case 'D':
    case 'd':
      DoDICompute();
      break;
    case 'S':
    case 's':
      DoPOPlot();
      break;
    case 'V':
    case 'v':
      DoDIPlot();
      break;
    }
  }while(opt!='E' && opt!='e');
  
  return(0);
}



void DoPOCompute() {
  
  struct poselfunc SF;
  struct objsim *OS;
  int nsim;
  struct objdet *OD;
  int ndet;
  struct simcat *SC;
  int ncats;
  int i;
  
  char selfile[101];
  char rastr[32],decstr[32];

  printf(" Number of catalags to compute selection function from ");
  ncats=readi(1);

  if((SC=(struct simcat*) malloc( ncats*sizeof(struct simcat)))==NULL) {
    printf("I cannot dimension OS of %d bytes \n",(ncats+1)*sizeof(struct simcat));
  }


  for(i=0;i<ncats;i++) {
    printf(" Name of catalog: ");
    reads("",SC[i].file);
    printf(" Identification of instrumental setup: ");
    reads("",SC[i].setup);
    printf(" Seeing for this catalog: ");
    SC[i].seeing=readf(1);
    printf(" Sky background for this catalog: ");
    SC[i].sky=readf(100);
    printf(" Transparency for this catalog: ");
    SC[i].transparency=readf(1);
    printf(" Central RA for this catalog: ");
    SC[i].alfac=readd(1);
    SC[i].alfac*=15.;
    printf(" Central Dec for this catalog: ");
    SC[i].deltac=readd(1);
    printf(" Input plate scale for this catalog: ");
    SC[i].platescale=readf(1);
    printf(" Input epoch for this catalog: ");
    SC[i].epoch=readf(2000);
    printf(" Input rotation for this catalog: ");
    SC[i].rotang=readf(0);
    printf(" Input flipping for this catalog: ");
    SC[i].flip=readf(0);
    printf(" Input X dimension for this catalog: ");
    SC[i].nx=readf(2048);
    printf(" Input Y dimension for this catalog: ");
    SC[i].ny=readf(2048);
    printf(" Input X lower limit of range of useful image: ");
    SC[i].x1=readf(1);
    printf(" Input X upper limit of range of useful image: ");
    SC[i].x2=readf(SC[i].nx);
    printf(" Input Y lower limit of range of useful image: ");
    SC[i].y1=readf(1);
    printf(" Input Y upper limit of range of useful image: ");
    SC[i].y2=readf(SC[i].ny);
    printf(" Name of file with detected ELGs for this catalog: ");
    reads("",SC[i].filedet);
    SC[i].wcsim=wcsxinit((SC[i].alfac),(SC[i].deltac),(double)(SC[i].platescale),SC[i].nx/2.,SC[i].ny/2.,SC[i].nx,SC[i].ny,0.,2000,SC[i].epoch,"TAN");
    if(!SC[i].flip) wcsdeltset(SC[i].wcsim, SC[i].wcsim->cdelt[0],SC[i].wcsim->cdelt[1],(double)SC[i].rotang);
    else      wcsdeltset(SC[i].wcsim,-SC[i].wcsim->cdelt[0],SC[i].wcsim->cdelt[1],(double)SC[i].rotang);
  }

  printf(" Number of intervals in seeing ");
  SF.nseeing=readi(3);
  printf(" Number of intervals in sky ");
  SF.nsky=readi(3);
  printf(" Number of intervals in transparency ");
  SF.ntransparency=readi(3);

  printf(" Number of intervals in equivalent width ");
  SF.nEW=readi(3);
  printf(" Number of intervals in z ");
  SF.nz=readi(3);
  printf(" Number of intervals in mag");
  SF.nmag=readi(40);

  
  ReadCatAndFillStructFunSel_PO(&SF,SC,ncats,&OS,&nsim,&OD,&ndet);
  
  if(DEBUG2) {
    for(i=0;i<nsim;i++) {
      ra2str(rastr,32,OS[i].alfa,3);  dec2str(decstr,32,OS[i].delta,2);
      printf(" objsim %d ew %f mag %f z %f ra %s dec %s\n",i,OS[i].ew,OS[i].mag,OS[i].z,rastr,decstr);
    }
    for(i=0;i<ndet;i++) {
      ra2str(rastr,32,OD[i].alfa,3);  dec2str(decstr,32,OD[i].delta,2);
      printf(" objdet %d ew %f mag %f z %f ra %s dec %s\n",i,OD[i].ew,OD[i].mag,OD[i].z,rastr,decstr); 
    }
  }

  printf(" Number of useful objects in catalogue %d. Detected ones: %d\n",nsim,ndet);

  MatchSimDet(OS,nsim,OD,ndet);

  ComputePOFunSel(&SF,OS,nsim,OD,ndet);
  
  

  printf(" Name of output selection file\n");
  reads("",selfile);

  writeposelfunc(SF, selfile);  
}

void DoPOPlot() {

  struct poselfunc SF;

  int i;

  char selfile[101];
  char opt='E';

  double *p;
  double *x,*y;
  float   min,max;
  int nbin,iproj=7;
  int ifix1=5;
/*   int ifix2=6; */
  double valuefix1=0;
  /*   double valuefix2=0; */
  double **p2;
  int iprojx=7,iprojy=5;
  int nx,ny;
  int ls;
  char xlabel[100],ylabel[100];
  float xpl=0, ypl=0;
  double xmin,xmax;
  char leyend[100];
  /*   int n_col=256; */
  /*   int rr[256],gg[256],bb[256]; */
  double az=30,al=60;

  /*   int nmovie=50,imovie; */
  double alt[50],azi[50];
  double az1=10,az2=80,al1=10,al2=80;
  /*   char rootname[200],moviename[200]; */
  double ew=0,z=0,magn=0;
  struct SurveyItem si;
  double pd;
     
  printf(" Name of selection file: ");
  reads("",selfile);
  
  readposelfunc(&SF, selfile);

  do{
    printf(" S Projection of SF in one parameter\n");
    printf(" P Projection of SF in one parameter. Other parameter  fixed\n"); 
    printf(" O Projection of SF in one parameter. Diferent figures in other parameter\n"); 
    printf(" L Projection of SF in one parameter. Other two parameters  fixed\n"); 
    printf(" D Projection of SF in two parameters\n");
    printf(" M Movie of projection of SF in two parameters\n");
    printf(" V Value of selection function for a given set of values\n");
    printf(" E Exit\n");
    opt=readc(opt); 
    switch (opt) { 
    case 'S':
    case 's':
      printf(" 1 Project in seeing\n");
      printf(" 2 Project in sky brightness\n");
      printf(" 3 Project in system efficiency\n");
      printf(" 4 Project in setups\n");
      printf(" 5 Project in EW\n");
      printf(" 6 Project in z\n");
      printf(" 7 Project in mag\n");
      printf(" 8 Project in line flux\n");
      printf(" Type variable to which to project: ");
      iproj=readi(iproj);
      projectposf(SF,&p,&x,&nbin,iproj);
      
      for(i=0;i<nbin;i++) 
	printf(" x %f p %f\n",x[i],p[i]);

      pgLimits_d(nbin,p,&min,&max);
      cpgopen("?");
      cpgscf(2);
      cpgsch(1.7);
      cpgvstd();
/*       cpghist3d((*ima).image.array,(*ima).image.nx,(*ima).image.ny,x1,x2,y1,y2); */
      cpgswin(x[0]-0.1*(x[nbin-1]-x[0]),x[nbin-1]+0.1*(x[nbin-1]-x[0]),min,max);
      cpgbox("BCTNS",0,0,"BCTNS",0,0);
      cpgline_d(nbin,x,p);
      if(iproj==1 ) cpglab("Seeing","P\\dDet","");
      if(iproj==2 ) cpglab("Sky level","P\\dDet","");
      if(iproj==3 ) cpglab("Efficiency of the system","P\\dDet","");
      if(iproj==4 ) cpglab("Instrumental Setup","P\\dDet","");
      if(iproj==5 ) cpglab("EW","P\\dDet","");
      if(iproj==6 ) cpglab("z","P\\dDet","");
      if(iproj==7 ) cpglab("m","P\\dDet","");
      if(iproj==8 ) cpglab("m\\dL","P\\dDet","");
      
      cpgclos();     
      
      free(x);free(p); 
      break;
    case 'D':
    case 'd':
      printf(" 1 Project in seeing\n");
      printf(" 2 Project in sky brightness\n");
      printf(" 3 Project in system efficiency\n");
      printf(" 4 Project in setups\n");
      printf(" 5 Project in EW\n");
      printf(" 6 Project in z\n");
      printf(" 7 Project in mag\n");
      printf(" 8 Project in line flux \n");
      printf(" Type variable to which to project in X: ");
      iprojx=readi(iprojx);
      printf(" Type variable to which to project in Y: ");
      iprojy=readi(iprojy);
      printf(" Type azimut viewing: ");
      az=readd(az);
      printf(" Type altitud viewing: ");
      al=readd(al);
      project2posf(SF,&p2,&x,&y,&nx,&ny,iprojx,iprojy);
      
      if(iprojx==1 ) strcpy(xlabel,"Seeing");
      if(iprojx==2 ) strcpy(xlabel,"Sky level");
      if(iprojx==3 ) strcpy(xlabel,"Efficiency of the system");
      if(iprojx==4 ) strcpy(xlabel,"Instrumental Setup");
      if(iprojx==5 ) strcpy(xlabel,"EW");
      if(iprojx==6 ) strcpy(xlabel,"z");
      if(iprojx==7 ) strcpy(xlabel,"m");                         
      if(iprojx==8 ) strcpy(xlabel,"m\\dL");                         
      if(iprojy==1 ) strcpy(ylabel,"Seeing");
      if(iprojy==2 ) strcpy(ylabel,"Sky level");
      if(iprojy==3 ) strcpy(ylabel,"Efficiency of the system");
      if(iprojy==4 ) strcpy(ylabel,"Instrumental Setup");
      if(iprojy==5 ) strcpy(ylabel,"EW");
      if(iprojy==6 ) strcpy(ylabel,"z");
      if(iprojy==7 ) strcpy(ylabel,"m");                         
      if(iprojy==8 ) strcpy(ylabel,"m\\dL");                         

/*       pgLimits(nbin,p,&min,&max); */

#ifdef OPERA_HAVE_PLPLOT
      printf(" Using PlPlot libraries\n");
      plinit();
      pllightsource(1.,1.,1.);
      for (i=0;i<n_col;i++)
        rr[i] = gg[i] = bb[i] = i*256/n_col;
      plscmap1(rr,gg,bb,n_col);
      pladv(0);
      plvpor(0.0, 0.9, 0.0, 0.9); 
      plwind(-1.0, 1.0, -0.9, 1.4);
      plcol0(1);
      plw3d(1.0,1.0, 1.0, x[0], x[nx-1], y[0], y[ny-1], 0, 1.0, al, az);
      plbox3("bnstu", xlabel, 0.0, 0,"bnstu", ylabel, 0.0, 0,"bcdmnstuv", "p#dDet", 0.0, 0);
      plcol0(2);
      plot3d(x,y,p2,nx,ny,3,1);
      plend();
       
#else 
      cpgopen("?");
      cpgscf(2);
      cpgsch(1.7);
      cpgvstd();
      cpgswin(x[0]-0.1*(x[nbin-1]-x[0]),x[nbin-1]+0.1*(x[nbin-1]-x[0]),min,max);
      cpgbox("BCTNS",0,0,"BCTNS",0,0);
      cpgline_d(nbin,x,p);
      cpgclos();     
#endif
      
      free(x);free_matrix_d(p2,nx,ny); 
      break;
    case 'M':
    case 'm':
      printf(" 1 Project in seeing\n");
      printf(" 2 Project in sky brightness\n");
      printf(" 3 Project in system efficiency\n");
      printf(" 4 Project in setups\n");
      printf(" 5 Project in EW\n");
      printf(" 6 Project in z\n");
      printf(" 7 Project in mag\n");
      printf(" 8 Project in line flux\n");
      printf(" Type variable to which to project in X: ");
      iprojx=readi(iprojx);
      printf(" Type variable to which to project in Y: ");
      iprojy=readi(iprojy);
      printf(" Type azimut viewing range: ");
      az1=readd(az1);az2=readd(az2);
      printf(" Type altitud viewing range: ");
      al1=readd(al1);al2=readd(al2);
      project2posf(SF,&p2,&x,&y,&nx,&ny,iprojx,iprojy);
      
      if(iprojx==1 ) strcpy(xlabel,"Seeing");
      if(iprojx==2 ) strcpy(xlabel,"Sky level");
      if(iprojx==3 ) strcpy(xlabel,"Efficiency of the system");
      if(iprojx==4 ) strcpy(xlabel,"Instrumental Setup");
      if(iprojx==5 ) strcpy(xlabel,"EW");
      if(iprojx==6 ) strcpy(xlabel,"z");
      if(iprojx==7 ) strcpy(xlabel,"m");                         
      if(iprojx==8 ) strcpy(xlabel,"m\\dL");                         
      if(iprojy==1 ) strcpy(ylabel,"Seeing");
      if(iprojy==2 ) strcpy(ylabel,"Sky level");
      if(iprojy==3 ) strcpy(ylabel,"Efficiency of the system");
      if(iprojy==4 ) strcpy(ylabel,"Instrumental Setup");
      if(iprojy==5 ) strcpy(ylabel,"EW");
      if(iprojy==6 ) strcpy(ylabel,"z");
      if(iprojy==7 ) strcpy(ylabel,"m");                         
      if(iprojy==8 ) strcpy(ylabel,"m\\dL");                         

/*       pgLimits(nbin,p,&min,&max); */

#ifdef OPERA_HAVE_PLPLOT
      printf(" Using PlPlot libraries\n");
      pllightsource(1.,1.,1.);
      printf(" root name of JPEG images\n");
      reads("",rootname);
      for(imovie=0;imovie<nmovie;imovie++) {
	sprintf(moviename,"%s%03d.png",rootname,imovie);
	plsdev("png");
	plsfnam(moviename);
	plinit();
	alt[imovie]=al1+((float)imovie/(float)nmovie)*((float)imovie/nmovie)*(al2-al1);
	azi[imovie]=az1+sqrt((float)imovie/nmovie)*(az2-az1);
	for (i=0;i<n_col;i++)
	  rr[i] = gg[i] = bb[i] = i*256/n_col;
	plscmap1(rr,gg,bb,n_col);
	pladv(0);
	plvpor(0.0, 0.9, 0.0, 0.9); 
	plwind(-1.0, 1.0, -0.9, 1.4);
	plcol0(1);
	plw3d(1.0, 1.0, 1.0, x[0], x[nx-1], y[0], y[ny-1], 0, 1.0, alt[imovie], azi[imovie]);
	plbox3("bnstu", xlabel, 0.0, 0,"bnstu", ylabel, 0.0, 0,"bcdmnstuv", "p#dDet", 0.0, 0);
	plcol0(2);
	plot3d(x,y,p2,nx,ny,3,1);
	printf(" Con alt %f azi %f\n",alt[imovie],azi[imovie]);
	plend();
      }
      
#else 
      cpgopen("?");
      cpgscf(2);
      cpgsch(1.7);
      cpgvstd();
      cpgswin(x[0]-0.1*(x[nbin-1]-x[0]),x[nbin-1]+0.1*(x[nbin-1]-x[0]),min,max);
      cpgbox("BCTNS",0,0,"BCTNS",0,0);
      cpgline_d(nbin,x,p);
      cpgclos();     
#endif
      
      free(x);free_matrix_d(p2,nx,ny); 
      break;
    case 'P':
    case 'p':
      printf(" 1 Project in seeing\n");
      printf(" 2 Project in sky brightness\n");
      printf(" 3 Project in system efficiency\n");
      printf(" 4 Project in setups\n");
      printf(" 5 Project in EW\n");
      printf(" 6 Project in z\n");
      printf(" 7 Project in mag\n");
      printf(" 8 Project in line flux\n");
      printf(" Type variable to which to project: ");
      iproj=readi(iproj);
      printf(" Type variable to be fixed: ");
      ifix1=readi(ifix1);
      printf(" Input fixed value: ");
      valuefix1=readf(valuefix1);
      projectposf_fixone(SF,&p,&x,&nbin,iproj,ifix1,valuefix1);
      
      for(i=0;i<nbin;i++) 
	printf(" x %f p %f\n",x[i],p[i]);

      pgLimits_d(nbin,p,&min,&max);
      cpgopen("?");
      cpgscf(2);
      cpgsch(1.7);
      cpgvstd();
      cpgswin(x[0]-0.1*(x[nbin-1]-x[0]),x[nbin-1]+0.1*(x[nbin-1]-x[0]),min,max);
      cpgbox("BCTNS",0,0,"BCTNS",0,0);
      cpgline_d(nbin,x,p); 
      if(iproj==1 ) cpglab("Seeing","P\\dDet","");
      if(iproj==2 ) cpglab("Sky level","P\\dDet","");
      if(iproj==3 ) cpglab("Efficiency of the system","P\\dDet","");
      if(iproj==4 ) cpglab("Instrumental Setup","P\\dDet","");
      if(iproj==5 ) cpglab("EW","P\\dDet","");
      if(iproj==6 ) cpglab("z","P\\dDet","");
      if(iproj==7 ) cpglab("m","P\\dDet","");
      if(iproj==8 ) cpglab("m\\dL","P\\dDet","");
        
      cpgclos();     
      
      free(x);free(p); 
      break;
    case 'O':
    case 'o':
      printf(" 1 Project in seeing\n");
      printf(" 2 Project in sky brightness\n");
      printf(" 3 Project in system efficiency\n");
      printf(" 4 Project in setups\n");
      printf(" 5 Project in EW\n");
      printf(" 6 Project in z\n");
      printf(" 7 Project in mag\n");
      printf(" 8 Project in line flux\n");
      printf(" Type variable to which to project: ");
      iproj=readi(iproj);
      printf(" Type variable to be fixed: ");
      ifix1=readi(ifix1);
      cpgopen("?");
      cpgscf(2);
      cpgsch(1.7);
      cpgvstd();
      min=0;max=1;
      if(iproj==1 ) cpglab("Seeing","P\\dDet","");
      if(iproj==2 ) cpglab("Sky level","P\\dDet","");
      if(iproj==3 ) cpglab("Efficiency of the system","P\\dDet","");
      if(iproj==4 ) cpglab("Instrumental Setup","P\\dDet","");
      if(iproj==5 ) cpglab("EW","P\\dDet","");
      if(iproj==6 ) cpglab("z","P\\dDet","");
      if(iproj==7 ) cpglab("m","P\\dDet","");
      if(iproj==8 ) cpglab("m\\dL","P\\dDet","");
      ls=1;
      while(projectposf_fixonerot(SF,&p,&x,&nbin,iproj,ifix1)) {
	cpgsls(ls);
	xmin=x[0]-0.1*(x[nbin-1]-x[0]);
	xmax=x[nbin-1]+0.1*(x[nbin-1]-x[0]);
	cpgswin(xmin,xmax,min,max);
	cpgbox("BCTNS",0,0,"BCTNS",0,0);
	cpgline_d(nbin,x,p); 
 	printf(" Leyend (empty = leyend): ");
	reads("",leyend);
	if(strcmp("",leyend)) {
	  printf(" X position: ");
	  xpl=readf(xpl);
	  printf(" Y position: ");
	  ypl=readf(ypl);
	  cpgptxt(xpl+(xmax-xmin)/35,ypl,0.,0.,leyend);
	  cpgmove(xpl-(xmax-xmin)/80,ypl+0.01);
	  cpgdraw(xpl-(xmax-xmin)/30,ypl+0.01);       
	} 
	ls++;
      }
      cpgclos();     
      
      free(x);free(p); 
      break;
    case 'V':
    case 'v':
      printf(" Seeing: ");
      si.seeing=readd(si.seeing);
      printf(" Sky brightness: ");
      si.sky=readd(si.sky);
      printf(" System efficiency: ");
      si.transparency=readd(si.transparency);
      printf(" Setup: ");
      reads(si.instsetup,si.instsetup);
      printf(" EW: ");
      ew=readd(ew);
      printf(" z: ");
      z=readd(z);
      printf(" mag: ");
      magn=readd(magn);
      pd=prob_poselfunc_scale(SF,si,ew,z,magn);
      printf(" Probabiliy: %f\n",pd);


      break;
    }
  }while(opt!='E' && opt!='e');


}

void DoDICompute() {
  printf(" Not implemented yet\n");
}

void DoDIPlot() {
  printf(" Not implemented yet\n");
}

void ReadCatAndFillStructFunSel_PO(struct poselfunc *SF,struct simcat *SC,int ncats, struct objsim **OS, int *nsim, struct objdet **OD, int *ndet) {
  int i;

  double *seeing;
  double *sky;
  double *transparency;
  char  **setup;
        
  double *readx;
  double *ready;
  double *readew;
  double *readz;
  double *readmag;
  int *readlog;
  int nlines;
  int j;
 
  double *candew,*candmag,*candz;
  double *candalfa,*canddelta;
  int *candlog;
  int ncand;

  int off;

  SF->seeing      =vector_d(SF->nseeing      );
  SF->sky         =vector_d(SF->nsky         );
  SF->transparency=vector_d(SF->ntransparency);
  SF->instsetup    =vector_s(ncats,200);

  seeing=vector_d(ncats);sky=vector_d(ncats);transparency=vector_d(ncats);
  setup=vector_s(ncats,200);
  for(i=0;i<ncats;i++) {
    seeing[i]=SC[i].seeing;sky[i]=SC[i].sky;transparency[i]=SC[i].transparency;
    strcpy(setup[i],SC[i].setup);
  }

  GetIntervals(seeing      ,ncats,SF->seeing      ,SF->nseeing      );  
  GetIntervals(sky         ,ncats,SF->sky         ,SF->nsky         );  
  GetIntervals(transparency,ncats,SF->transparency,SF->ntransparency);  
  GetStrings(  setup       ,ncats,SF->instsetup    ,&(SF->ninstsetup) );



  *nsim=0;
  *ndet=0;
  for(i=0;i<ncats;i++) {
    nlines=FileNLin(SC[i].file);
    readx=vector_d(nlines);ready=vector_d(nlines);
    readew=vector_d(nlines);readz=vector_d(nlines);readmag=vector_d(nlines);
    readlog=vector_i(nlines);
    for(j=0;j<nlines;j++) readlog[j]=0;

    printf(" Number of lines in catalogue %s : %d\n",SC[i].file,nlines);
    ReadDoublecol(SC[i].file,2,readx,readlog,&nlines);
    ReadDoublecol(SC[i].file,3,ready,readlog,&nlines);
    ReadDoublecol(SC[i].file,5,readew,readlog,&nlines);
    ReadDoublecol(SC[i].file,7,readz,readlog,&nlines);
    ReadDoublecol(SC[i].file,4,readmag,readlog,&nlines);
    for(j=0;j<nlines;j++) {
      if(readlog[j] && readx[j]>SC[i].x1 && readx[j]<SC[i].x2 && ready[j]>SC[i].y1 && readx[j]<SC[i].y2) {
	if(*nsim==0) {
	  if(((*OS)=(struct objsim*) malloc( (1)*sizeof(struct objsim)))==NULL) {
	    printf("I cannot dimension OS of %d bytes \n",(*nsim+1)*sizeof(struct objsim));
	  }
	}
	else {
	  if(((*OS)=(struct objsim*) realloc(*OS, (*nsim+1)*sizeof(struct objsim)))==NULL) {
	    printf("I cannot dimension OS of %d bytes \n",(*nsim+1)*sizeof(struct objsim));
	  }
	}
	(*OS)[*nsim].x=readx[j];(*OS)[*nsim].y=ready[j];
	(*OS)[*nsim].ew=readew[j];(*OS)[*nsim].z=readz[j];(*OS)[*nsim].mag=readmag[j];

	pix2wcs(SC[i].wcsim,(*OS)[*nsim].x,(*OS)[*nsim].y,&((*OS)[*nsim].alfa),&((*OS)[*nsim].delta));
	/* (*OS)[*nsim].alfa=666.;	(*OS)[*nsim].delta=666.; */

	if(DEBUG2) printf(" %d ew %f m %f %d ew %f m %f al %f delta %f\n",j,readew[j],readmag[j],*nsim,(*OS)[*nsim].ew,(*OS)[*nsim].mag,(*OS)[*nsim].alfa,(*OS)[*nsim].delta);
	
	(*OS)[*nsim].seeing=SC[i].seeing;(*OS)[*nsim].sky=SC[i].sky;(*OS)[*nsim].transparency=SC[i].transparency;
	strcpy((*OS)[*nsim].instsetup,SC[i].setup);
	(*OS)[*nsim].selected=0;
	(*OS)[*nsim].od=NULL;
	(*nsim)++;
      }
    }
    ncand=FileNLin(SC[i].filedet);
    candew=vector_d(ncand);candmag=vector_d(ncand);candz=vector_d(ncand);
    candalfa=vector_d(ncand);canddelta=vector_d(ncand);
    candlog=vector_i(ncand);
    for(j=0;j<ncand;j++) candlog[j]=0;
    
    printf(" Number of lines in candidate file %s: %d\n",SC[i].filedet,ncand);
    ReadWCScol(SC[i].filedet,2, candalfa ,candlog,&ncand);
    ReadWCScol(SC[i].filedet,3, canddelta,candlog,&ncand);
    ReadDoublecol(SC[i].filedet,4, candew   ,candlog,&ncand);
    ReadDoublecol(SC[i].filedet,8, candz    ,candlog,&ncand);
    ReadDoublecol(SC[i].filedet,12, candmag ,candlog,&ncand);
    if(DEBUG) printf(" Despues de leer todo\n");
    for(j=0;j<ncand;j++) {
      if(candlog[j]) {
	if(*ndet==0) {
	  if(((*OD)=(struct objdet*) malloc( (1)*sizeof(struct objdet)))==NULL) {
	    printf("I cannot dimension OD of %d bytes \n",(*ndet+1)*sizeof(struct objdet));
	  }
	}
	else {
	  if(((*OD)=(struct objdet*) realloc(*OD, (*ndet+1)*sizeof(struct objdet)))==NULL) {
	    printf("I cannot dimension OD of %d bytes \n",(*ndet+1)*sizeof(struct objdet));
	  }
	}
	(*OD)[*ndet].alfa=candalfa[j];
	(*OD)[*ndet].alfa*=15.;;
	(*OD)[*ndet].delta=canddelta[j];
	(*OD)[*ndet].ew=candew[j];
	(*OD)[*ndet].z=candz[j];
    	(*OD)[*ndet].mag=candmag[j];
	off=0;
	wcs2pix(SC[i].wcsim,(*OD)[*ndet].alfa,(*OD)[*ndet].delta,&((*OD)[*ndet].x),&((*OD)[*ndet].y),&off);
	if(DEBUG2) printf(" alfa %f delta %f x %f y %f\n",(*OD)[*ndet].alfa,(*OD)[*ndet].delta,(*OD)[*ndet].x,(*OD)[*ndet].y);
	if(off) {printf(" Somewhat went wrong with astrometry. Exiting\n");exit(1); }
	
	(*OD)[*ndet].seeing=SC[i].seeing;(*OD)[*ndet].sky=SC[i].sky;(*OD)[*ndet].transparency=SC[i].transparency;
        if( (*OD)[*ndet].x>SC[i].x1 && (*OD)[*ndet].x<SC[i].x2 && (*OD)[*ndet].y>SC[i].y1 && (*OD)[*ndet].y<SC[i].y2 ) (*ndet)++;
      }
    }
    free(readx);free(ready);free(readew);free(readz);free(readmag);free(readlog);
    free(candalfa);free(canddelta);free(candew);free(candz);free(candmag);
  }
}



void GetIntervals(double *values, int nvalues, double *samplevalues, int nsample) {

  double min,max;
  double delta;
  int i;
  if(DEBUG) printf(" Antes min max n %d\n",nvalues);
  MinMax_d(nvalues,values, &min, &max);
  if(DEBUG) printf(" min %f max %f\n",min,max);
  delta=(max-min)/(nsample-1.);
  samplevalues[0]=min;
  samplevalues[nsample-1]=max;
  if(DEBUG) printf(" ANtes for \n");
  for(i=1;i<nsample-1;i++)   samplevalues[i]=min+i*delta;
  if(DEBUG) printf("Salio\n");
}


void GetStrings(char **values, int nvalues, char **samplevalues, int *nsample) {
  
  int i,j;
  int foundflag;

  *nsample=0;
  for(i=0;i<nvalues;i++) {
    foundflag=0;
    if(DEBUG) printf(" Voy por valor %d  %s\n",i,values[i]);
    for(j=0;j<*nsample;j++) {
      if(DEBUG) printf(" Comparing <%s> <%s>\n",values[i],samplevalues[j]);
      if(!strcmp(values[i],samplevalues[j])) { 
	if(DEBUG) printf(" Son iguales!!!\n");
	foundflag=1;
	break;
      }
    }
    if(!foundflag) {
      strcpy(samplevalues[*nsample],values[i]);
      (*nsample)++;
      if(DEBUG) printf(" Lo añado. Llevo %d\n",*nsample);
    }
  }
}


void MatchSimDet(struct objsim *OS, int nsim, struct objdet *OD, int ndet) {


  int i,j;
  float tol=1.5; 
  float tolwcs=1.5; 
  double mindist,mindistwcs,distwcs;
  char rastr[32],decstr[32];
  int jmin=0;

  for(i=0;i<ndet;i++) {
    printf("%05d/%05d\b\b\b\b\b\b\b\b\b\b\b",i,ndet);
    fflush(NULL);
    if(DEBUG2) printf(" Matching %d x %f y %f\n",i,OD[i].x,OD[i].y);
    mindist=1e39;
    mindistwcs=1e39;
    for(j=0;j<nsim;j++) {
/*       dist=sqrt((OD[i].x-OS[j].x)*(OD[i].x-OS[j].x)+(OD[i].y-OS[j].y)*(OD[i].y-OS[j].y)); */
      distwcs=wcsdist(OD[i].alfa,OD[i].delta,OS[j].alfa,OS[j].delta);
      if(distwcs<mindistwcs) {
/*       if(dist<mindist) { */
/* 	mindist=dist; */
	mindist=sqrt((OD[i].x-OS[j].x)*(OD[i].x-OS[j].x)+(OD[i].y-OS[j].y)*(OD[i].y-OS[j].y));
	mindistwcs=distwcs;
	jmin=j;
      }
    }
    if(mindistwcs*3600.<tolwcs) {
      OS[jmin].selected=1;
      OS[jmin].od=&(OD[i]);
      OD[i].os=&(OS[jmin]);
      if(DEBUG2) printf(" Match OS %d x %f y %f OD %d x %f y %f  dist %f distwcs %g\n",jmin,OS[jmin].x,OS[jmin].y,i,OD[i].x,OD[i].y,mindist,mindistwcs*3600.); 
      ra2str(rastr,32,OD[i].alfa,3);  dec2str(decstr,32,OD[i].delta,2);
      if(DEBUG2) printf("       OD alfa %f %s  delta %f %s\n",OD[i].alfa,rastr,OD[i].delta,decstr);
      ra2str(rastr,32,OS[jmin].alfa,3);  dec2str(decstr,32,OS[jmin].delta,2);
      if(DEBUG2) printf("       OS alfa %f %s  delta %f %s\n",OS[jmin].alfa,rastr,OS[jmin].delta,decstr);
      
    }
    else {
      printf("ERROR: No detected matched object \n");
      if(DEBUG) printf(" Match OS %d x %f y %f OD %d x %f y %f\n",jmin,OS[jmin].x,OS[jmin].y,i,OD[i].x,OD[i].y);
      ra2str(rastr,32,OD[i].alfa,3);  dec2str(decstr,32,OD[i].delta,2);
      if(DEBUG) printf("       OD alfa %f %s  delta %f %s\n",OD[i].alfa,rastr,OD[i].delta,decstr);
      ra2str(rastr,32,OS[jmin].alfa,3);  dec2str(decstr,32,OS[jmin].delta,2);
      if(DEBUG) printf("       OS alfa %f %s  delta %f %s\n",OS[jmin].alfa,rastr,OS[jmin].delta,decstr);
      exit(1);
    }
  }
  if(DEBUG ) printf(" Salio de match\n");
}


void ComputePOFunSel( struct poselfunc *SF,struct objsim *OS,int nsim,struct objdet *OD,int ndet) {


  double *objew;
  double *objz;
  double *objmag;
  int i;

  int ntmp;
  double ptmp;

  int iew,iz,imag,iseeing,isky,itrans,isetup;

  int   *******nbin;
  int   *******nbindet;

  if(DEBUG) printf(" Entro en el calculo\n");
  
  SF->ewbin=vector_d(SF->nEW);
  SF->zbin=vector_d(SF->nz);
  SF->magbin=vector_d(SF->nmag);

  if(DEBUG) printf(" Despues alloc\n");

  objew=vector_d(nsim);objz=vector_d(nsim);objmag=vector_d(nsim);
  for(i=0;i<nsim;i++) {
    objew[i]=OS[i].ew;objz[i]=OS[i].z;objmag[i]=OS[i].mag;
  }

  if(DEBUG) printf(" Antes interv\n");

  GetIntervals(objew,nsim,SF->ewbin,SF->nEW);
  GetIntervals(objz,nsim,SF->zbin,SF->nz);
  GetIntervals(objmag,nsim,SF->magbin,SF->nmag);

  if(DEBUG) printf(" Despu inter\n");

  nbin   =tensor7_i(SF->nseeing,SF->nsky,SF->ntransparency,SF->ninstsetup,SF->nEW,SF->nz,SF->nmag);
  nbindet=tensor7_i(SF->nseeing,SF->nsky,SF->ntransparency,SF->ninstsetup,SF->nEW,SF->nz,SF->nmag);
  (*SF).p   =tensor7_d(SF->nseeing,SF->nsky,SF->ntransparency,SF->ninstsetup,SF->nEW,SF->nz,SF->nmag);
  (*SF).errp=tensor7_d(SF->nseeing,SF->nsky,SF->ntransparency,SF->ninstsetup,SF->nEW,SF->nz,SF->nmag);

  if(DEBUG) printf("antes buvl\n");

  if(DEBUG) printf(" KAKA\n");

  for(iseeing=0;iseeing<SF->nseeing;iseeing++) {  
    for(isky=0;isky<SF->nsky;isky++) {  
      for(itrans=0;itrans<SF->ntransparency;itrans++) {  
	for(isetup=0;isetup<SF->ninstsetup;isetup++) {  
	  for(iew=0;iew<SF->nEW;iew++) {
	    for(iz=0;iz<SF->nz;iz++) {
	      for(imag=0;imag<SF->nmag;imag++) {
		nbin[iseeing][isky][itrans][isetup][iew][iz][imag]=0;
		nbindet[iseeing][isky][itrans][isetup][iew][iz][imag]=0;
		//		p[iseeing][isky][itrans][isetup][iew][iz][imag]=45.;
	      }
	    }
	  }
	}
      }
    }
  }

  	  
  if(DEBUG) {
    printf(" EW intervals \n");
    for(i=0;i<SF->nEW;i++) printf(" %d %f\n",i,SF->ewbin[i]);
    printf(" z intervals \n");
    for(i=0;i<SF->nz;i++) printf(" %d %f\n",i,SF->zbin[i]);
    printf(" Mag intervals \n");
    for(i=0;i<SF->nmag;i++) printf(" %d %f\n",i,SF->magbin[i]);
  }


  for(i=0;i<nsim;i++) {
    iew     =sf_selectbin(OS[i].ew,SF->ewbin,SF->nEW);
    iz      =sf_selectbin(OS[i].z,SF->zbin,SF->nz);
    imag    =sf_selectbin(OS[i].mag,SF->magbin,SF->nmag);
    iseeing =sf_selectbin(OS[i].seeing,SF->seeing,SF->nseeing);
    isky    =sf_selectbin(OS[i].sky,SF->sky,SF->nsky);
    itrans  =sf_selectbin(OS[i].transparency,SF->transparency,SF->ntransparency);
    isetup  =sf_selectbin_char(OS[i].instsetup,SF->instsetup,SF->ninstsetup);


    if(DEBUG2) printf(" Parametros galaxia  EW %f z %f mag %f\n",OS[i].ew,OS[i].z,OS[i].mag);
    if(DEBUG2) printf(" Asignados  iew %d iz %d imag %d \n",iew,iz,imag);
    if(DEBUG2) printf("       see %f sky %f tra %f  setup %s\n",OS[i].seeing,OS[i].sky,OS[i].transparency,OS[i].instsetup);
    if(DEBUG2) printf(" Asignados  iseeing %d isky %d itrans %d iset %d\n",iseeing,isky,itrans,isetup);
    if(DEBUG2) printf("      setup %s Asignado %d\n",OS[i].instsetup,isetup);
    
    //    if(DEBUG) i=readi(i); 
    (nbin[iseeing][isky][itrans][isetup][iew][iz][imag])++;
    if(OS[i].selected)     (nbindet[iseeing][isky][itrans][isetup][iew][iz][imag])++;

    if(DEBUG2) {if(OS[i].selected) printf(" ESTE SELECCCCC\n");}
    if(DEBUG2)printf(" Este bin %d sel %d\n",nbin[iseeing][isky][itrans][isetup][iew][iz][imag],nbindet[iseeing][isky][itrans][isetup][iew][iz][imag]);

  }



  printf(" Objects assigned to bins.\n");
/*   i=readi(i); */
        
  printf(" New %d\n",SF->nEW);
  for(iseeing=0;iseeing<SF->nseeing;iseeing++) {  
    for(isky=0;isky<SF->nsky;isky++) {   
      for(itrans=0;itrans<SF->ntransparency;itrans++) {  
	for(isetup=0;isetup<SF->ninstsetup;isetup++) {  
	  for(iew=0;iew<SF->nEW;iew++) {
	    for(iz=0;iz<SF->nz;iz++) {
	      for(imag=0;imag<SF->nmag;imag++) {
    		ntmp=(nbin[iseeing][isky][itrans][isetup][iew][iz][imag]);
		if(ntmp==0) {
		  printf(" ERROR: bin with 0 galaxies. Exiting\n");
		  printf(" ise %d isk %d it %d ise %d iew %d iz %d im %d\n",iseeing,isky,itrans,isetup,iew,iz,imag); 
		  exit(1);
		}
		if(DEBUG) printf(" ise %d isk %d it %d ise %d iew %d iz %d im %d\n",iseeing,isky,itrans,isetup,iew,iz,imag); 
 		if(DEBUG2) printf(" ntmp %d \n",ntmp); 
 		ptmp=(float)(nbindet[iseeing][isky][itrans][isetup][iew][iz][imag]); 
		ptmp/=ntmp;
		if(DEBUG) printf(" ntmp %d ptmp %f\n",ntmp,ptmp);
 		SF->p[iseeing][isky][itrans][isetup][iew][iz][imag]=ptmp; 
  		SF->errp[iseeing][isky][itrans][isetup][iew][iz][imag]=sqrt(ptmp*(1-ptmp)/ntmp); 
	      }
	    }
	  }
	}
      }
    } 
  }
  printf(" Probability of detection computed\n");

  free(objew);
  free(objz);
  free(objmag);
}
