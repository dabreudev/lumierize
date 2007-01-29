#include "modulos.h"
#define DEBUG 0


void ReadSurDB(char *dbfile, struct SurveyDB  *sdb) {

  int colimage=1;
  int colRA=2;
  int colDEC=3;
  int colRAsize=4;
  int colDECsize=5;
  int colrot=6;
  int colepoch=7;
  int colstdx=8;
  int colstdy=9;
  int colcalfile=10;
  int colnobj=11;
  int colcalA=12;
  int colerrcalA=13;
  int colcalB=14;
  int colerrcalB=15;
  int colcovcalAB=16;
  int colcalrms=17;
  int colmaxmag=18;
  int colmodemag=19;
  int colmfermi=20;
  int colerrmfermi=21;
  int coldelfermi=22;
  int colerrdelfermi=23;
  int colalpha=24;
  int colerralpha=25;
  int colseeing=26;
  int colstdseeing=27;
  int colsky=28;
  int colsigsky=29;
  int colspecfile=30;
  int colrespfile=31;
  int colinstsetup=32;
  
  double *alfac,*deltac;
  float *xdim,*ydim;
  float *rot;
  float *epoch;
  float *stdx,*stdy;

  float *a,*b;
  float *erra,*errb,*covab;
  float *rms;

  float *mmaximum;
  float *mmode;
  float *mfermicut;
  float *deltafermi;
  float *gammapowerlaw;
  float *errmfermicut;
  float *errdeltafermi;
  float *errgammapowerlaw;
  
  float *seeing;
  float *stdseeing;
  float *sky;
  float *sigsky;

  char *imagename;
  char *specfile;
  char *respfile;
  char *calfile;
  char *instsetup;
  float * nobj;

  int *ilog;

  int nlin;

  int i;

  printf(" Reading database...");
  if(DEBUG) printf(" Al principio\n");

  nlin=FileNLin(dbfile);

  alfac=vector_d(nlin);
  deltac=vector_d(nlin);
  xdim=vector_f(nlin);
  ydim=vector_f(nlin);
  rot=vector_f(nlin);
  epoch=vector_f(nlin);
  stdx=vector_f(nlin);
  stdy=vector_f(nlin);
  
  a=vector_f(nlin);
  b=vector_f(nlin);
  erra=vector_f(nlin);
  errb=vector_f(nlin);
  covab=vector_f(nlin);
  rms=vector_f(nlin);
  
  mmaximum=vector_f(nlin);
  mmode=vector_f(nlin);
  mfermicut=vector_f(nlin);
  deltafermi=vector_f(nlin);
  gammapowerlaw=vector_f(nlin);
  errmfermicut=vector_f(nlin);
  errdeltafermi=vector_f(nlin);
  errgammapowerlaw=vector_f(nlin);
  
  seeing=vector_f(nlin);
  stdseeing=vector_f(nlin);

  sky=vector_f(nlin);
  sigsky=vector_f(nlin);

  nobj=vector_f(nlin);

  if((imagename=malloc(101*nlin))==NULL) printf(" Cannot allocate imagename of %d elements\n",nlin);
  if((specfile =malloc(101*nlin))==NULL) printf(" Cannot allocate specfile of %d elements\n",nlin);
  if((respfile =malloc(101*nlin))==NULL) printf(" Cannot allocate respfile of %d elements\n",nlin);
  if((calfile  =malloc(101*nlin))==NULL) printf(" Cannot allocate catfile of %d elements\n",nlin);
  if((instsetup=malloc(101*nlin))==NULL) printf(" Cannot allocate instsetup of %d elements\n",nlin);

  ilog=vector_i(nlin);
 
  if(DEBUG) printf(" Termino alloc \n");

  ReadCharcol(dbfile,colimage,imagename,ilog,101,&nlin); 

  ReadWCScol(dbfile,colRA,alfac,ilog,&nlin);
  ReadWCScol(dbfile,colDEC,deltac,ilog,&nlin);
  ReadNumcol(dbfile,colRAsize,xdim,ilog,&nlin);
  ReadNumcol(dbfile,colDECsize,ydim,ilog,&nlin);
  ReadNumcol(dbfile,colrot,rot,ilog,&nlin);
  ReadNumcol(dbfile,colepoch,epoch,ilog,&nlin);
  ReadNumcol(dbfile,colstdx,stdx,ilog,&nlin);
  ReadNumcol(dbfile,colstdy,stdy,ilog,&nlin);

  ReadNumcol(dbfile,colcalA,a,ilog,&nlin);
  ReadNumcol(dbfile,colcalB,b,ilog,&nlin);
  ReadNumcol(dbfile,colerrcalA,erra,ilog,&nlin); 
  ReadNumcol(dbfile,colerrcalB,errb,ilog,&nlin);
  ReadNumcol(dbfile,colcovcalAB,covab,ilog,&nlin);
  ReadNumcol(dbfile,colcalrms,rms,ilog,&nlin);
  ReadNumcol(dbfile,colmaxmag,mmaximum,ilog,&nlin);
  ReadNumcol(dbfile,colmodemag,mmode,ilog,&nlin);
  ReadNumcol(dbfile,colmfermi,mfermicut,ilog,&nlin);
  ReadNumcol(dbfile,colerrmfermi,errmfermicut,ilog,&nlin); 
  ReadNumcol(dbfile,coldelfermi,deltafermi,ilog,&nlin);
  ReadNumcol(dbfile,colerrdelfermi,errdeltafermi,ilog,&nlin);
  ReadNumcol(dbfile,colalpha,gammapowerlaw,ilog,&nlin);
  ReadNumcol(dbfile,colerralpha,errgammapowerlaw,ilog,&nlin);
  ReadNumcol(dbfile,colseeing,seeing,ilog,&nlin);
  ReadNumcol(dbfile,colstdseeing,stdseeing,ilog,&nlin);
  ReadNumcol(dbfile,colsky,sky,ilog,&nlin);
  ReadNumcol(dbfile,colsigsky,sigsky,ilog,&nlin);


  if(DEBUG) printf(" Antes spec\n");

  ReadCharcol(dbfile,colspecfile,specfile,ilog,101,&nlin); 
  ReadCharcol(dbfile,colrespfile,respfile,ilog,101,&nlin); 
  ReadCharcol(dbfile,colcalfile ,calfile ,ilog,101,&nlin); 
  ReadCharcol(dbfile,colinstsetup ,instsetup ,ilog,101,&nlin); 
  ReadNumcol(dbfile,colnobj,nobj,ilog,&nlin);

  if(DEBUG) printf(" Termino Read\n");

  (*sdb).nitems=0;
  if(((*sdb).si=(struct SurveyItem *) malloc(1*sizeof(struct SurveyItem)))==NULL) printf("I cannot dimension (*sdb).si   of %d elements \n",1);
  for(i=0;i<nlin;i++) {
    if(ilog[i]) {
      /* printf(" AR %f DEC %f\n",alfac[i],deltac[i]); */
      if(DEBUG) printf(" %d yd %f xd %f\n",i,xdim[i],ydim[i]); 
      ((*sdb).si[(*sdb).nitems]).alfac=alfac[i];
      ((*sdb).si[(*sdb).nitems]).deltac=deltac[i];
      ((*sdb).si[(*sdb).nitems]).xdim=xdim[i];
      ((*sdb).si[(*sdb).nitems]).ydim=ydim[i];
      ((*sdb).si[(*sdb).nitems]).rot=rot[i];
      ((*sdb).si[(*sdb).nitems]).epoch=epoch[i];
      ((*sdb).si[(*sdb).nitems]).stdx=stdx[i];
      ((*sdb).si[(*sdb).nitems]).stdy=stdy[i];
      ((*sdb).si[(*sdb).nitems]).a=a[i];
      ((*sdb).si[(*sdb).nitems]).b=b[i];
      ((*sdb).si[(*sdb).nitems]).erra=erra[i];
      ((*sdb).si[(*sdb).nitems]).errb=errb[i];
      ((*sdb).si[(*sdb).nitems]).covab=covab[i];
      ((*sdb).si[(*sdb).nitems]).rms=rms[i];
      ((*sdb).si[(*sdb).nitems]).mmaximum=mmaximum[i];
      ((*sdb).si[(*sdb).nitems]).mmode=mmode[i];
      ((*sdb).si[(*sdb).nitems]).mfermicut=mfermicut[i];
      ((*sdb).si[(*sdb).nitems]).deltafermi=deltafermi[i];
      ((*sdb).si[(*sdb).nitems]).gammapowerlaw=gammapowerlaw[i];
      ((*sdb).si[(*sdb).nitems]).errmfermicut=errmfermicut[i];
      ((*sdb).si[(*sdb).nitems]).errdeltafermi=errdeltafermi[i];
      ((*sdb).si[(*sdb).nitems]).errgammapowerlaw=errgammapowerlaw[i];
      ((*sdb).si[(*sdb).nitems]).seeing=seeing[i];
      ((*sdb).si[(*sdb).nitems]).stdseeing=stdseeing[i]; 
      ((*sdb).si[(*sdb).nitems]).sky=sky[i]; 
      ((*sdb).si[(*sdb).nitems]).sigsky=sigsky[i]; 
      ((*sdb).si[(*sdb).nitems]).nobj=(int)(nobj[i]); 
      strcpy(((*sdb).si[(*sdb).nitems]).image,imagename+i*101);
      strcpy(((*sdb).si[(*sdb).nitems]).specfile,specfile+i*101);
      strcpy(((*sdb).si[(*sdb).nitems]).respfile,respfile+i*101);
      strcpy(((*sdb).si[(*sdb).nitems]).calfile,calfile+i*101);
      strcpy(((*sdb).si[(*sdb).nitems]).instsetup,instsetup+i*101);
      (*sdb).nitems++;
      if(((*sdb).si=(struct SurveyItem *) realloc((*sdb).si   ,((*sdb).nitems+1)*sizeof(struct SurveyItem)))==NULL) printf("I cannot dimension (*sdb).si        of %d elements \n",(*sdb).nitems);
    }
  }

  free(alfac);free(deltac);
  free(xdim);free(ydim);
  free(rot);
  free(epoch);
  
  free(a);free(b);
  free(erra);free(errb);free(covab);
  free(rms);
  
  free(mmaximum);
  free(mmode);
  free(mfermicut);
  free(deltafermi);
  free(gammapowerlaw);
  free(errmfermicut);
  free(errdeltafermi);
  free(errgammapowerlaw);
  
  free(seeing);
  free(stdseeing);
  free(sky);
  free(sigsky);

  free(imagename);
  free(specfile);
  free(respfile);
  free(instsetup);
 
  free(ilog);
  
  printf(" done\n");
}

int  whithinimage(double ra, double dec, struct SurveyItem sitem) {

  /* RA en horas 
     DEC en grados */

  double px[4],py[4];
  double pxrot[4],pyrot[4];
  
  double vec0x, vec0y;
  double vec1x, vec1y;
  double vec2x, vec2y;

  double alpha, beta;
  double vec1vec1;
  double vec2vec2;
  double vec1vec2;
  double vec0vec1;  
  double vec0vec2;

  double ptvec1;
  double ptvec2;

  double xplac,yplac;
  double offsetx;
  double offsety;


  px[0]=-sitem.xdim/3600./2./180*M_PI;py[0]=-sitem.ydim/3600./2./180*M_PI;
  px[1]=+sitem.xdim/3600./2./180*M_PI;py[1]=-sitem.ydim/3600./2./180*M_PI;
  px[2]=+sitem.xdim/3600./2./180*M_PI;py[2]=+sitem.ydim/3600./2./180*M_PI;
  px[3]=-sitem.xdim/3600./2./180*M_PI;py[3]=+sitem.ydim/3600./2./180*M_PI;

  offsetx=sitem.alfac /180*M_PI*15.;
  offsety=sitem.deltac /180*M_PI;
  pxrot[0]= (px[0]*cos(sitem.rot/180*M_PI)-py[0]*sin(sitem.rot/180*M_PI));
  pyrot[0]= px[0]*sin(sitem.rot/180*M_PI)+py[0]*cos(sitem.rot/180*M_PI);
  	   
  pxrot[1]= (px[1]*cos(sitem.rot/180*M_PI)-py[1]*sin(sitem.rot/180*M_PI));
  pyrot[1]= px[1]*sin(sitem.rot/180*M_PI)+py[1]*cos(sitem.rot/180*M_PI);
  	   
  pxrot[2]= (px[2]*cos(sitem.rot/180*M_PI)-py[2]*sin(sitem.rot/180*M_PI));
  pyrot[2]= px[2]*sin(sitem.rot/180*M_PI)+py[2]*cos(sitem.rot/180*M_PI);
  	   
  pxrot[3]= (px[3]*cos(sitem.rot/180*M_PI)-py[3]*sin(sitem.rot/180*M_PI));
  pyrot[3]= px[3]*sin(sitem.rot/180*M_PI)+py[3]*cos(sitem.rot/180*M_PI);

  
  vec1x=pxrot[0]-pxrot[1];
  vec1y=pyrot[0]-pyrot[1];
  vec2x=pxrot[2]-pxrot[1];
  vec2y=pyrot[2]-pyrot[1];

  vec0x=pxrot[1];
  vec0y=pyrot[1];
  
  vec1vec1=vec1x*vec1x+vec1y*vec1y;
  vec2vec2=vec2x*vec2x+vec2y*vec2y;
  vec1vec2=vec1x*vec2x+vec1y*vec2y;
  vec0vec1=vec0x*vec1x+vec0y*vec1y;
  vec0vec2=vec0x*vec2x+vec0y*vec2y;

/*   printf(" ra %f dec %f\n",ra,dec); */

  Ecu2Plac(ra/180*M_PI*15.,dec/180*M_PI,offsetx,offsety,&xplac,&yplac); 


  ptvec1=xplac*vec1x+yplac*vec1y;
  ptvec2=xplac*vec2x+yplac*vec2y;

/*   printf(" xplac %f ra-ra %f yplac %f dec-dec %f\n",xplac,ra/180*M_PI-offsetx,yplac,dec/180*M_PI-offsety); */


  beta= ((ptvec1-vec0vec1)*vec1vec2-(ptvec2-vec0vec2)*vec1vec1)/(vec1vec2*vec1vec2-vec2vec2*vec1vec1);
  alpha=((ptvec1-vec0vec1)*vec2vec2-(ptvec2-vec0vec2)*vec1vec2)/(vec1vec1*vec2vec2-vec1vec2*vec1vec2);

  if(alpha>0 && alpha<1 && beta>0 && beta<1) {

/*     printf(" ra %f dec %f p0 %f %f p1 %f %f p2 %f %f\n",ra,dec,pxrot[0],pyrot[0],pxrot[1],pyrot[1],pxrot[2],pyrot[2]); */
/*     printf(" vec1 %f %f vec2 %f %f\n",vec1x,vec1y,vec2x,vec2y); */
/*     printf(" alpha %f beta %f\n",alpha,beta); */
    
    return 1;
  } 
  else return 0;

}


float Surveyrad(struct SurveyDB sdb) {

  /* Devuelve el resultado en radianes cuadrados */

  int i,j,k;
  int nsub=60;
  int icont=0;

  int nra=100;
  int ndec=100;
  double area;
  int covered;
  float ra_,dec_;
  

  float ramin,ramax,decmin,decmax,minxdim,maxxdim,minydim,maxydim;
  float *ra,*dec,*xdim,*ydim;

  printf(" Computing area... ");

  ra =vector_f(sdb.nitems);
  dec=vector_f(sdb.nitems);
  xdim =vector_f(sdb.nitems);
  ydim=vector_f(sdb.nitems);
  for(i=0;i<sdb.nitems;i++) {
    ra[i] =(sdb.si[i]).alfac;
    dec[i]=(sdb.si[i]).deltac;
    xdim[i] =(sdb.si[i]).xdim;
    ydim[i]=(sdb.si[i]).ydim;
    if((sdb.si[i]).xdim==0 || (sdb.si[i]).ydim==0) {
      xdim[i]=xdim[i-1];
      ydim[i]=ydim[i-1];
      if((sdb.si[i]).alfac==0 && (sdb.si[i]).deltac==0) {
	ra[i] =ra[i-1];
	dec[i]=dec[i-1];
      }
    }
  }
  MinMax(sdb.nitems,ra,&ramin,&ramax);
  MinMax(sdb.nitems,dec,&decmin,&decmax);
  MinMax(sdb.nitems,xdim,&minxdim,&maxxdim);
  MinMax(sdb.nitems,ydim,&minydim,&maxydim);
  ramin=ramin-maxxdim/1./3600./15.;
  ramax=ramax+maxxdim/1./3600./15.;
  decmin=decmin-maxydim/1./3600.;
  decmax=decmax+maxydim/1./3600.;    

  if(DEBUG) {
    printf(" minxdim %f max %f\n",minxdim,maxxdim);
    printf(" minydim %f max %f\n",minydim,maxydim);
  }
  
  nra =(int)((ramax-ramin)/(minxdim/3600./15./nsub));
  ndec=(int)((decmax-decmin)/(minydim/3600./nsub));
  area=0;
  if(DEBUG) {
    printf(" DESDE ra %f %f \n",ramin,ramax);
    printf(" DESDE dec %f %f \n",decmin,decmax);
    printf(" CON ASO %d %d \n",nra,ndec);
  }
  printf("   %%");
  printf("\b\b\b\b");
  for(i=0;i<nra;i++) {
    if((i*100./(nra-1))>icont) {
      printf("%3d\b\b\b",icont);
      icont++;
      fflush(NULL);
    }
    ra_=ramin+i*(ramax-ramin)/(nra-1);
    for(j=0;j<ndec;j++) {
      dec_=decmin+j*(decmax-decmin)/(ndec-1);
      covered=0;
      for(k=0;k<sdb.nitems;k++) {
	if(whithinimage(ra_,dec_,sdb.si[k])) {
	  covered=1;
	  break;
	}
      }
      area+=covered*cos(dec_/180*M_PI);
    }
  }
  area=area*(ramax-ramin)*M_PI/12*(decmax-decmin)*M_PI/180/nra/ndec;

  if(DEBUG) printf(" area1 %f 2 %f\n",area*180/M_PI*180/M_PI,area*180/M_PI*180/M_PI);

  printf("                   ");
  printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
  printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
  free(ra);
  free(dec);
  free(xdim);
  free(ydim);
  return((float)area);
}

void SaveSurDB(char *dbfile, struct SurveyDB sdb) {

  int i;
  int asstatus=1,mcstatus=1,mlstatus=1,sestatus=1,specstatus=1;

  FILE *fr;
  char wcsar[16],wcsdec[16];
  if((fr=fopen(dbfile,"w")) == NULL) {
    printf(" Couldn't create file %s. Exiting\n",dbfile);
    exit(1);
  }
  fprintf(fr,"#%-51s %-16s  %-16s  %-9s  %-9s  %-8s  %10s  %10s  %10s  %-51s  %-10s  %-9s  %-9s  %-9s  %-9s  %-9s  %-6s  %-7s  %-7s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-7s  %-7s  %-51s  %-51s  %-51s\n","Image","Center_RA","Center_DEC","RA_size","DEC_size","Rot(ang)","Epoch","Std_X","Std_Y","Calibration_file","Num_obj","Mag_cal_A","Err_Mag_A","Mag_cal_B","Err_Mag_B","Cov_AB","RMS_MC","max_mag","mode_ma","m_fermi","errm_fermi","del_fer","errdel_fer","Alpha_PL","errAlpha_PL","Seeing","Std_See","Spectra_file","Response_file","Instrumental_setup");
  
  for(i=0;i<sdb.nitems;i++) {
    
    ra2str(wcsar,16,15*sdb.si[i].alfac,3);
    dec2str(wcsdec,16,sdb.si[i].deltac,3);
  if(asstatus)
    fprintf(fr,"%-51s  %-16s  %-16s  %9.2g  %9.2g  %8.3f  %10.3f  %10.5f  %10.5f  ",sdb.si[i].image,wcsar,wcsdec,sdb.si[i].xdim,sdb.si[i].ydim,sdb.si[i].rot,sdb.si[i].epoch,sdb.si[i].stdx,sdb.si[i].stdy);
  else 
    fprintf(fr,"%-51s  %-16s  %-16s  %-9s  %-9s  %-8s  %-10s  %10s  %10s  ",sdb.si[i].image,"INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF");
  if(mcstatus)
    fprintf(fr,"%-51s  %10d  %9.4f  %9.4g  %9.4f  %9.4g  %9.4g  %6.3f  ",sdb.si[i].calfile,sdb.si[i].nobj,sdb.si[i].a,sdb.si[i].erra,sdb.si[i].b,sdb.si[i].errb,sdb.si[i].covab,sdb.si[i].rms);
  else 
    fprintf(fr,"%-51s  %10d  %-9s  %-9s  %-9s  %-9s  %-9s  %-6s  ",sdb.si[i].calfile,sdb.si[i].nobj,"INDEF","INDEF","INDEF","INDEF","INDEF","INDEF");
  if(mlstatus)
    fprintf(fr,"%7.3f  %7.3f  %9.4f  %9.4f  %9.4g  %9.4f  %9.4f  %9.4g  ",sdb.si[i].mmaximum,sdb.si[i].mmode,sdb.si[i].mfermicut,sdb.si[i].errmfermicut,sdb.si[i].deltafermi,sdb.si[i].errdeltafermi,sdb.si[i].gammapowerlaw,sdb.si[i].errgammapowerlaw);
  else 
    fprintf(fr,"%-7s  %-7s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  ","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF");
  if(sestatus)
    fprintf(fr,"%7.3f  %7.3f  ",sdb.si[i].seeing,sdb.si[i].stdseeing);
  else 
    fprintf(fr,"%-7s  %-7s  ","INDEF","INDEF");
  if(specstatus)
    fprintf(fr,"%-51s  %-51s  ",sdb.si[i].specfile,sdb.si[i].respfile);
  else 
    fprintf(fr,"%-51s  %-51s  ","INDEF","INDEF");
  fprintf(fr,"%-51s\n",sdb.si[i].instsetup);

  }
  fclose(fr);
}


void PlotSurDB(struct SurveyDB sdb, int labflag) {

  float *ra,*dec;
  int i;

  double px[4],py[4];
  double pxrot[4],pyrot[4];
  double raproj[4],decproj[4];
  float ramin,ramax,decmin,decmax;

  
  ra =vector_f(sdb.nitems);
  dec=vector_f(sdb.nitems);

  for(i=0;i<sdb.nitems;i++) {
    ra[i] =(sdb.si[i]).alfac;
    dec[i]=(sdb.si[i]).deltac;
  }
  MinMax(sdb.nitems,ra,&ramin,&ramax);
  MinMax(sdb.nitems,dec,&decmin,&decmax);
  
  if(DEBUG)   printf(" decmin %f decmax %f\n",decmin,decmax);
  if(DEBUG)   printf(" ramin %f ramax %f\n",ramin,ramax);
  
  ramin=ramin-(sdb.si[0]).xdim/1./3600./15.;
  ramax=ramax+(sdb.si[0]).xdim/1./3600./15.;
  decmin=decmin-(sdb.si[0]).ydim/1./3600.;
  decmax=decmax+(sdb.si[0]).ydim/1./3600.;    

  if(DEBUG) printf(" decmin %f decmax %f\n",decmin,decmax);

  cpgpage();
  cpgvstd();
  cpgwnad(ramax*3600*15.,ramin*3600*15.,decmin*3600,decmax*3600);
  cpgswin(ramax*3600,ramin*3600,decmin*3600,decmax*3600);
  cpgswin(ramax*3600,ramin*3600,decmin*3600,decmax*3600);
  cpgtbox("ZXBCTNSHG",0.0,0,"ZBCTNSDGV",0.0,0);
  cpgswin(ramax,ramin,decmin,decmax);
  
  for(i=0;i<sdb.nitems;i++) {
    
    px[0]=-(sdb.si[i]).xdim/3600./2./180*M_PI;py[0]=-(sdb.si[i]).ydim/3600./2./180*M_PI;
    px[1]=+(sdb.si[i]).xdim/3600./2./180*M_PI;py[1]=-(sdb.si[i]).ydim/3600./2./180*M_PI;
    px[2]=+(sdb.si[i]).xdim/3600./2./180*M_PI;py[2]=+(sdb.si[i]).ydim/3600./2./180*M_PI;
    px[3]=-(sdb.si[i]).xdim/3600./2./180*M_PI;py[3]=+(sdb.si[i]).ydim/3600./2./180*M_PI;
    pxrot[0]= (px[0]*cos((sdb.si[i]).rot/180*M_PI)-py[0]*sin((sdb.si[i]).rot/180*M_PI));
    pyrot[0]= px[0]*sin((sdb.si[i]).rot/180*M_PI)+py[0]*cos((sdb.si[i]).rot/180*M_PI);
    pxrot[1]= (px[1]*cos((sdb.si[i]).rot/180*M_PI)-py[1]*sin((sdb.si[i]).rot/180*M_PI));
    pyrot[1]= px[1]*sin((sdb.si[i]).rot/180*M_PI)+py[1]*cos((sdb.si[i]).rot/180*M_PI);
    pxrot[2]= (px[2]*cos((sdb.si[i]).rot/180*M_PI)-py[2]*sin((sdb.si[i]).rot/180*M_PI));
    pyrot[2]= px[2]*sin((sdb.si[i]).rot/180*M_PI)+py[2]*cos((sdb.si[i]).rot/180*M_PI);
    pxrot[3]= (px[3]*cos((sdb.si[i]).rot/180*M_PI)-py[3]*sin((sdb.si[i]).rot/180*M_PI));
    pyrot[3]= px[3]*sin((sdb.si[i]).rot/180*M_PI)+py[3]*cos((sdb.si[i]).rot/180*M_PI);
    Plac2Ecu(pxrot[0],pyrot[0],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+0,decproj+0);
    Plac2Ecu(pxrot[1],pyrot[1],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+1,decproj+1);
    Plac2Ecu(pxrot[2],pyrot[2],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+2,decproj+2);
    Plac2Ecu(pxrot[3],pyrot[3],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+3,decproj+3);
    raproj[0]*=180/M_PI/15; decproj[0]*=180/M_PI;
    raproj[1]*=180/M_PI/15; decproj[1]*=180/M_PI;
    raproj[2]*=180/M_PI/15; decproj[2]*=180/M_PI;
    raproj[3]*=180/M_PI/15; decproj[3]*=180/M_PI;
    
    
    cpgsfs(1);
    cpgsci(0);
    cpgpoly_d(4,raproj,decproj);
    cpgsci(2);
    cpgsfs(2);
    cpgpoly_d(4,raproj,decproj);
    if(DEBUG) printf(" Image %d ra %f %f %f %f dec %f %f %f %f\n",i,raproj[0],raproj[1],raproj[2],raproj[3],decproj[0],decproj[1],decproj[2],decproj[3]);
    
    cpgsci(1);
    if(labflag) cpgtext((pxrot[0]+pxrot[1]+pxrot[2]+pxrot[3])/4,(pyrot[0]+pyrot[1]+pyrot[2]+pyrot[3])/4,sdb.si[i].image);
    
  }
  
  cpglab("RA","","");
  
  cpgptxt(ramin-(ramax-ramin)*0.08,decmin+(decmax-decmin)/2.,90,0.5,"Dec");
  
  
  free(ra);
  free(dec);



}

void PlotSurDB_zoom(struct SurveyDB sdb, float ramin,float ramax, float decmin, float decmax, int labflag) {

  float *ra,*dec;
  int i;

  double px[4],py[4];
  double pxrot[4],pyrot[4];

  double raproj[4],decproj[4];

  
  ra =vector_f(sdb.nitems);
  dec=vector_f(sdb.nitems);

  for(i=0;i<sdb.nitems;i++) {
    ra[i] =(sdb.si[i]).alfac;
    dec[i]=(sdb.si[i]).deltac;
  }
  if(ramin==ramax && decmin==decmax) {
    MinMax(sdb.nitems,ra,&ramin,&ramax);
    MinMax(sdb.nitems,dec,&decmin,&decmax);
    ramin=ramin-(sdb.si[0]).xdim/1./3600./15.;
    ramax=ramax+(sdb.si[0]).xdim/1./3600./15.;
    decmin=decmin-(sdb.si[0]).ydim/1./3600.;
    decmax=decmax+(sdb.si[0]).ydim/1./3600.;    
  }
  
  cpgpage();
  cpgvstd();
  cpgwnad(ramax*3600*15.,ramin*3600*15.,decmin*3600,decmax*3600);
  cpgswin(ramax*3600,ramin*3600,decmin*3600,decmax*3600);
  cpgswin(ramax*3600,ramin*3600,decmin*3600,decmax*3600);
  cpgtbox("ZXBCTNSHG",0.0,0,"ZBCTNSDGV",0.0,0);
  cpgswin(ramax,ramin,decmin,decmax);
  
  for(i=0;i<sdb.nitems;i++) {
    
    px[0]=-(sdb.si[i]).xdim/3600./2./180*M_PI;py[0]=-(sdb.si[i]).ydim/3600./2./180*M_PI;
    px[1]=+(sdb.si[i]).xdim/3600./2./180*M_PI;py[1]=-(sdb.si[i]).ydim/3600./2./180*M_PI;
    px[2]=+(sdb.si[i]).xdim/3600./2./180*M_PI;py[2]=+(sdb.si[i]).ydim/3600./2./180*M_PI;
    px[3]=-(sdb.si[i]).xdim/3600./2./180*M_PI;py[3]=+(sdb.si[i]).ydim/3600./2./180*M_PI;
    pxrot[0]= (px[0]*cos((sdb.si[i]).rot/180*M_PI)-py[0]*sin((sdb.si[i]).rot/180*M_PI));
    pyrot[0]= px[0]*sin((sdb.si[i]).rot/180*M_PI)+py[0]*cos((sdb.si[i]).rot/180*M_PI);
    pxrot[1]= (px[1]*cos((sdb.si[i]).rot/180*M_PI)-py[1]*sin((sdb.si[i]).rot/180*M_PI));
    pyrot[1]= px[1]*sin((sdb.si[i]).rot/180*M_PI)+py[1]*cos((sdb.si[i]).rot/180*M_PI);
    pxrot[2]= (px[2]*cos((sdb.si[i]).rot/180*M_PI)-py[2]*sin((sdb.si[i]).rot/180*M_PI));
    pyrot[2]= px[2]*sin((sdb.si[i]).rot/180*M_PI)+py[2]*cos((sdb.si[i]).rot/180*M_PI);
    pxrot[3]= (px[3]*cos((sdb.si[i]).rot/180*M_PI)-py[3]*sin((sdb.si[i]).rot/180*M_PI));
    pyrot[3]= px[3]*sin((sdb.si[i]).rot/180*M_PI)+py[3]*cos((sdb.si[i]).rot/180*M_PI);
    Plac2Ecu(pxrot[0],pyrot[0],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+0,decproj+0);
    Plac2Ecu(pxrot[1],pyrot[1],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+1,decproj+1);
    Plac2Ecu(pxrot[2],pyrot[2],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+2,decproj+2);
    Plac2Ecu(pxrot[3],pyrot[3],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+3,decproj+3);
    raproj[0]*=180/M_PI/15; decproj[0]*=180/M_PI;
    raproj[1]*=180/M_PI/15; decproj[1]*=180/M_PI;
    raproj[2]*=180/M_PI/15; decproj[2]*=180/M_PI;
    raproj[3]*=180/M_PI/15; decproj[3]*=180/M_PI;
    
    
    cpgsfs(1);
    cpgsci(0);
    cpgpoly_d(4,raproj,decproj);
    cpgsci(2);
    cpgsfs(2);
    cpgpoly_d(4,raproj,decproj);
    cpgsci(1);
    if(labflag) cpgtext((pxrot[0]+pxrot[1]+pxrot[2]+pxrot[3])/4,(pyrot[0]+pyrot[1]+pyrot[2]+pyrot[3])/4,sdb.si[i].image);
    
  }
  
  cpglab("RA","","");
  
  cpgptxt(ramin-(ramax-ramin)*0.08,decmin+(decmax-decmin)/2.,90,0.5,"Dec");
  
  
  free(ra);
  free(dec);
}
