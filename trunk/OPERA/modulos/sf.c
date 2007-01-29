#include "modulos.h"



void writeposelfunc(struct poselfunc SF, char selfile[101]) {
  
  int i;
  FILE *fp;
  int iew,iz,imag,iseeing,isky,itrans,isetup;

  if((fp=fopen(selfile,"w"))==NULL) {
    printf("ERROR: Can't open output file %s\n",selfile);
    exit(1);
  }


  fprintf(fp,"# Selection function file\n");
  fprintf(fp,"#nseeing  nsky    ntrans  nsetup  new    nz    nmag \n");
  fprintf(fp," %d %d %d %d %d %d %d\n",SF.nseeing,SF.nsky,SF.ntransparency,SF.ninstsetup,SF.nEW,SF.nz,SF.nmag);
  fprintf(fp,"#Seeing bins\n");
  for(i=0;i<SF.nseeing;i++)    fprintf(fp," %f\n",SF.seeing[i]);
  fprintf(fp,"#Sky bins\n");
  for(i=0;i<SF.nsky;i++)    fprintf(fp," %f\n",SF.sky[i]);
  fprintf(fp,"#Transparency bins\n");
  for(i=0;i<SF.ntransparency;i++)    fprintf(fp," %f\n",SF.transparency[i]);
  fprintf(fp,"#Setups\n");
  for(i=0;i<SF.ninstsetup;i++)    fprintf(fp," %s\n",SF.instsetup[i]);
  fprintf(fp,"#EW bins\n");
  for(i=0;i<SF.nEW;i++)    fprintf(fp," %f\n",SF.ewbin[i]);
  fprintf(fp,"#Redshift bins\n"); 
  for(i=0;i<SF.nz;i++)    fprintf(fp," %f\n",SF.zbin[i]);
  fprintf(fp,"#apparent magnitude bins\n");
  for(i=0;i<SF.nmag;i++)    fprintf(fp," %f\n",SF.magbin[i]);

  fprintf(fp,"#Selection Function\n");
  for(isetup=0;isetup<SF.ninstsetup;isetup++) {
    for(iseeing=0;iseeing<SF.nseeing;iseeing++) {
      for(isky=0;isky<SF.nsky;isky++) {
	for(itrans=0;itrans<SF.ntransparency;itrans++) {
	  for(iew=0;iew<SF.nEW;iew++) {
	    for(iz=0;iz<SF.nz;iz++) {
	      for(imag=0;imag<SF.nmag;imag++) {
		fprintf(fp," %f  %f\n",SF.p[iseeing][isky][itrans][isetup][iew][iz][imag],SF.errp[iseeing][isky][itrans][isetup][iew][iz][imag]);
	      }
	    }
	  }
	}
      }
    }
  }
	  
  fclose(fp);

}


void readposelfunc(struct poselfunc *SF, char selfile[101]) {
  
  int i;
  FILE *fp;
  int iew,iz,imag,iseeing,isky,itrans,isetup;

  char snul[1000];

  if((fp=fopen(selfile,"r"))==NULL) {
    printf("ERROR: Can't open input file %s\n",selfile);
    exit(1);
  }


  fgetline(fp,snul,1000);
  fgetline(fp,snul,1000);
  
  fscanf(fp," %d %d %d %d %d %d %d\n",&(*SF).nseeing,&(*SF).nsky,&(*SF).ntransparency,&(*SF).ninstsetup,&(*SF).nEW,&(*SF).nz,&(*SF).nmag);
/*   printf(" Ns %d %d %d %d %d %d %d\n", (*SF).nseeing, (*SF).nsky, (*SF).ntransparency, (*SF).ninstsetup, (*SF).nEW, (*SF).nz, (*SF).nmag); */
  fgetline(fp,snul,1000); 

  (*SF).seeing=vector_d((*SF).nseeing);
  (*SF).sky   =vector_d((*SF).nsky   );
  (*SF).transparency=vector_d((*SF).ntransparency);
  (*SF).instsetup=vector_s((*SF).ninstsetup,200);
  (*SF).ewbin=vector_d((*SF).nEW);
  (*SF).zbin=vector_d((*SF).nz);
  (*SF).magbin=vector_d((*SF).nmag);

  (*SF).p   =tensor7_d(SF->nseeing,SF->nsky,SF->ntransparency,SF->ninstsetup,SF->nEW,SF->nz,SF->nmag);
  (*SF).errp=tensor7_d(SF->nseeing,SF->nsky,SF->ntransparency,SF->ninstsetup,SF->nEW,SF->nz,SF->nmag);


  for(i=0;i<(*SF).nseeing;i++)    fscanf(fp," %lf\n",&(*SF).seeing[i]);
  fgetline(fp,snul,1000); 
  for(i=0;i<(*SF).nsky;i++)    fscanf(fp," %lf\n",&(*SF).sky[i]);
  fgetline(fp,snul,1000);
  for(i=0;i<(*SF).ntransparency;i++)    fscanf(fp," %lf\n",&(*SF).transparency[i]);
  fgetline(fp,snul,1000); 
  for(i=0;i<(*SF).ninstsetup;i++)    fscanf(fp," %s\n",(*SF).instsetup[i]);
  fgetline(fp,snul,1000); 
  for(i=0;i<(*SF).nEW;i++)    fscanf(fp," %lf\n",&(*SF).ewbin[i]);
  fgetline(fp,snul,1000); 
  for(i=0;i<(*SF).nz;i++)    fscanf(fp," %lf\n",&(*SF).zbin[i]);
  fgetline(fp,snul,1000); 
  for(i=0;i<(*SF).nmag;i++)    fscanf(fp," %lf\n",&(*SF).magbin[i]);

/*   for(i=0;i<(*SF).nmag;i++)    printf("MAg %f\n",(*SF).magbin[i]); */

  fgetline(fp,snul,1000); 
  for(isetup=0;isetup<(*SF).ninstsetup;isetup++) {  
    for(iseeing=0;iseeing<(*SF).nseeing;iseeing++) {  
      for(isky=0;isky<(*SF).nsky;isky++) {   
	for(itrans=0;itrans<(*SF).ntransparency;itrans++) {  
	  for(iew=0;iew<(*SF).nEW;iew++) {
	    for(iz=0;iz<(*SF).nz;iz++) {
	      for(imag=0;imag<(*SF).nmag;imag++) {
		fscanf(fp," %lf  %lf\n",&(*SF).p[iseeing][isky][itrans][isetup][iew][iz][imag],&(*SF).errp[iseeing][isky][itrans][isetup][iew][iz][imag]);
/* 		printf(" p %f\n",(*SF).p[iseeing][isky][itrans][isetup][iew][iz][imag]); */
	      }
	    }
	  }
	}
      }
    }
  }
	  
  fclose(fp);

}



void projectposf(struct poselfunc SF, double **p, double **x, int *nbin, int iproj) {
  
  int *nn;
  int i;
  int ibin;
  int iew,iz,imag,iseeing,isky,itrans,isetup;
  double delta,min,max;
  double cval;

  if(iproj==1 ) { *nbin=SF.nseeing;        min=SF.seeing[0];max=SF.seeing[SF.nseeing-1];}
  if(iproj==2 ) { *nbin=SF.nsky   ;        min=SF.sky[0];max=SF.sky[SF.nsky-1];}
  if(iproj==3 ) { *nbin=SF.ntransparency;  min=SF.transparency[0];max=SF.transparency[SF.ntransparency-1];}
  if(iproj==4 ) { *nbin=SF.ninstsetup ;    min=1;max=SF.ninstsetup;}
  if(iproj==5 ) { *nbin=SF.nEW;            min=SF.ewbin[0];max=SF.ewbin[SF.nEW-1];}
  if(iproj==6 ) { *nbin=SF.nz;             min=SF.zbin[0];max=SF.zbin[SF.nz-1];}
  if(iproj==7 ) { *nbin=SF.nmag;           min=SF.magbin[0];max=SF.magbin[SF.nmag-1];}
  if(iproj==8 ) { *nbin=SF.nmag;           min=SF.magbin[0]-2.5*log10(SF.ewbin[SF.nEW-1]);max=SF.magbin[SF.nmag-1];}

  printf(" min %f max %f\n",min,max);
/*   exit(1); */

  nn=vector_i(*nbin);
  *p=vector_d(*nbin);
  *x=vector_d(*nbin);


  for(i=0;i<*nbin;i++) {
    (*p)[i]=0;
    nn[i]=0;
  }

    
  for(iseeing=0;iseeing<SF.nseeing;iseeing++) {  
    for(isky=0;isky<SF.nsky;isky++) {   
      for(itrans=0;itrans<SF.ntransparency;itrans++) {  
	for(isetup=0;isetup<SF.ninstsetup;isetup++) {  
	  for(iew=0;iew<SF.nEW;iew++) {
	    for(iz=0;iz<SF.nz;iz++) {
	      for(imag=0;imag<SF.nmag;imag++) {
		if(iproj==1 ) ibin=iseeing;
		if(iproj==2 ) ibin=isky;
		if(iproj==3 ) ibin=itrans;
		if(iproj==4 ) ibin=isetup;
		if(iproj==5 ) ibin=iew;
		if(iproj==6 ) ibin=iz;
		if(iproj==7 ) ibin=imag;
		if(iproj==8 ) {
		  cval=(SF.magbin[0]+imag*(SF.magbin[SF.nmag-1]-SF.magbin[0])/(SF.nmag-1))-2.5*log10((SF.ewbin[0]+iew*(SF.ewbin[SF.nEW-1]-SF.ewbin[0])/(SF.nEW-1)));
		  ibin=(int)((*nbin-1)*(cval-min)/(max-min)+0.5);
		}
/* 		printf(" ibin %d\n",ibin); */
/* 		kak=SF.p[iseeing][isky][itrans][isetup][iew][iz][imag]; */
		((*p)[ibin])+=SF.p[iseeing][isky][itrans][isetup][iew][iz][imag];
/* 		printf(" kak %f\n",kak); */
		(nn[ibin])++;
	      }
	    }
	  }
	}
      }
    }
  }  

  
  for(i=0;i<*nbin;i++) {
    (*p)[i]=(*p)[i]/nn[i]; 
    if(nn[i]==0) (*p)[i]=0;
  }


  
  delta=(max-min)/(*nbin-1);
  
  for(i=0;i<*nbin;i++)  (*x)[i]=min+i*delta;
  
}


void projectposf_fixone(struct poselfunc SF, double **p, double **x, int *nbin, int iproj, int ifix, double valuefix) {
  
  int *nn;
  int i;
  int ibin;
  int iew,iz,imag,iseeing,isky,itrans,isetup;
  int iew_f=0,iz_f=0,imag_f=0,iseeing_f=0,isky_f=0,itrans_f=0,isetup_f=0;
  int ivalue,  nbin_f;
  double delta,min,max;
  double min_f, max_f;
  double cval;

  int trueflag=0;

  if(iproj==1 ) { *nbin=SF.nseeing;       min=SF.seeing[0];max=SF.seeing[SF.nseeing-1];}
  if(iproj==2 ) { *nbin=SF.nsky   ;       min=SF.sky[0];max=SF.sky[SF.nsky-1];}
  if(iproj==3 ) { *nbin=SF.ntransparency; min=SF.transparency[0];max=SF.transparency[SF.ntransparency-1];}
  if(iproj==4 ) { *nbin=SF.ninstsetup ;   min=1;max=SF.ninstsetup;}
  if(iproj==5 ) { *nbin=SF.nEW;           min=SF.ewbin[0];max=SF.ewbin[SF.nEW-1];}
  if(iproj==6 ) { *nbin=SF.nz;            min=SF.zbin[0];max=SF.zbin[SF.nz-1];}
  if(iproj==7 ) { *nbin=SF.nmag;          min=SF.magbin[0];max=SF.magbin[SF.nmag-1];}
  if(iproj==8 ) { *nbin=SF.nmag;           min=SF.magbin[0]-2.5*log10(SF.ewbin[SF.nEW-1]);max=SF.magbin[SF.nmag-1];}

  if(ifix==1 ) { nbin_f=SF.nseeing;       iseeing_f=1;  min_f=SF.seeing[0];max_f=SF.seeing[SF.nseeing-1];}
  if(ifix==2 ) { nbin_f=SF.nsky   ;       isky_f=1;     min_f=SF.sky[0];max_f=SF.sky[SF.nsky-1];}
  if(ifix==3 ) { nbin_f=SF.ntransparency; itrans_f=1;   min_f=SF.transparency[0];max_f=SF.transparency[SF.ntransparency-1];}
  if(ifix==4 ) { nbin_f=SF.ninstsetup ;   isetup_f=1;   min_f=1;max_f=SF.ninstsetup;}
  if(ifix==5 ) { nbin_f=SF.nEW;           iew_f=1;      min_f=SF.ewbin[0];max_f=SF.ewbin[SF.nEW-1];}
  if(ifix==6 ) { nbin_f=SF.nz;            iz_f=1;       min_f=SF.zbin[0];max_f=SF.zbin[SF.nz-1];}
  if(ifix==7 ) { nbin_f=SF.nmag;          imag_f=1;     min_f=SF.magbin[0];max_f=SF.magbin[SF.nmag-1];}


  ivalue=(int)((nbin_f-1)*(valuefix-min_f)/(max_f-min_f)+0.5);
  if(ivalue<0) ivalue=0;
  if(ivalue>nbin_f-1) ivalue=nbin_f-1;
  

  printf(" Me quedo solo con el %d en el rango %f-%f\n",ivalue,min_f,max_f);

  nn=vector_i(*nbin);
  *p=vector_d(*nbin);
  *x=vector_d(*nbin);


  for(i=0;i<*nbin;i++) {
    (*p)[i]=0;
    nn[i]=0;
  }

    
  for(iseeing=0;iseeing<SF.nseeing;iseeing++) {  
    for(isky=0;isky<SF.nsky;isky++) {   
      for(itrans=0;itrans<SF.ntransparency;itrans++) {  
	for(isetup=0;isetup<SF.ninstsetup;isetup++) {  
	  for(iew=0;iew<SF.nEW;iew++) {
	    for(iz=0;iz<SF.nz;iz++) {
	      for(imag=0;imag<SF.nmag;imag++) {
		if(iproj==1 ) ibin=iseeing;
		if(iproj==2 ) ibin=isky;
		if(iproj==3 ) ibin=itrans;
		if(iproj==4 ) ibin=isetup;
		if(iproj==5 ) ibin=iew;
		if(iproj==6 ) ibin=iz;
		if(iproj==7 ) ibin=imag;
		if(iproj==8 ) {
		  cval=(SF.magbin[0]+imag*(SF.magbin[SF.nmag-1]-SF.magbin[0])/(SF.nmag-1))-2.5*log10((SF.ewbin[0]+iew*(SF.ewbin[SF.nEW-1]-SF.ewbin[0])/(SF.nEW-1)));
		  ibin=(int)((*nbin-1)*(cval-min)/(max-min)+0.5);
		}
		trueflag=1;
		if(iseeing_f) {if(!(iseeing==ivalue)) trueflag=0;}
		if(isky_f) {if(!(isky==ivalue)) trueflag=0;}
		if(itrans_f) {if(!(itrans==ivalue)) trueflag=0;}
		if(isetup_f) {if(!(isetup==ivalue)) trueflag=0;}
		if(iew_f) { if(!(iew==ivalue)) trueflag=0;}
		if(iz_f) {if(!(iz==ivalue)) trueflag=0;}
		if(imag_f) {if(!(imag==ivalue)) trueflag=0;}

		if(trueflag) {
		  ((*p)[ibin])+=SF.p[iseeing][isky][itrans][isetup][iew][iz][imag];
		  (nn[ibin])++;
		} else {
/* 		  printf(" He excludio \n"); */
/* 		  printf(" ise %d isk %d it %d ise %d iew %d iz %d im %d\n",iseeing,isky,itrans,isetup,iew,iz,imag); */
		}
	      }
	    }
	  }
	}
      }
    }
  }  

  
  for(i=0;i<*nbin;i++) {
    (*p)[i]=(*p)[i]/nn[i]; 
    if(nn[i]==0) (*p)[i]=0;
  }


  
  delta=(max-min)/(*nbin-1);
  
  for(i=0;i<*nbin;i++)  (*x)[i]=min+i*delta;
  
}

int projectposf_fixonerot(struct poselfunc SF, double**p, double**x, int *nbin, int iproj, int ifix) {
  
  int *nn;
  int i;
  int ibin;
  int iew,iz,imag,iseeing,isky,itrans,isetup;
  int iew_f=0,iz_f=0,imag_f=0,iseeing_f=0,isky_f=0,itrans_f=0,isetup_f=0;
  int  nbin_f;
  double delta,min,max;
  double min_f, max_f;
  int static ivalue=0;
  double cval;

  int trueflag=0;

  if(iproj==1 ) { *nbin=SF.nseeing;       min=SF.seeing[0];max=SF.seeing[SF.nseeing-1];}
  if(iproj==2 ) { *nbin=SF.nsky   ;       min=SF.sky[0];max=SF.sky[SF.nsky-1];}
  if(iproj==3 ) { *nbin=SF.ntransparency; min=SF.transparency[0];max=SF.transparency[SF.ntransparency-1];}
  if(iproj==4 ) { *nbin=SF.ninstsetup ;   min=1;max=SF.ninstsetup;}
  if(iproj==5 ) { *nbin=SF.nEW;           min=SF.ewbin[0];max=SF.ewbin[SF.nEW-1];}
  if(iproj==6 ) { *nbin=SF.nz;            min=SF.zbin[0];max=SF.zbin[SF.nz-1];}
  if(iproj==7 ) { *nbin=SF.nmag;          min=SF.magbin[0];max=SF.magbin[SF.nmag-1];}
  if(iproj==8 ) { *nbin=SF.nmag;          min=SF.magbin[0]-2.5*log10(SF.ewbin[SF.nEW-1]);max=SF.magbin[SF.nmag-1];}
  if(ifix==1 )  { nbin_f=SF.nseeing;      iseeing_f=1;  min_f=SF.seeing[0];max_f=SF.seeing[SF.nseeing-1];}
  if(ifix==2 )  { nbin_f=SF.nsky   ;      isky_f=1;     min_f=SF.sky[0];max_f=SF.sky[SF.nsky-1];}
  if(ifix==3 )  { nbin_f=SF.ntransparency;itrans_f=1;   min_f=SF.transparency[0];max_f=SF.transparency[SF.ntransparency-1];}
  if(ifix==4 )  { nbin_f=SF.ninstsetup ;  isetup_f=1;   min_f=1;max_f=SF.ninstsetup;}
  if(ifix==5 )  { nbin_f=SF.nEW;          iew_f=1;      min_f=SF.ewbin[0];max_f=SF.ewbin[SF.nEW-1];}
  if(ifix==6 )  { nbin_f=SF.nz;           iz_f=1;       min_f=SF.zbin[0];max_f=SF.zbin[SF.nz-1];}
  if(ifix==7 )  { nbin_f=SF.nmag;         imag_f=1;     min_f=SF.magbin[0];max_f=SF.magbin[SF.nmag-1];}

  if(ivalue<0) ivalue=0;
  if(ivalue>nbin_f-1) {
    ivalue=0;
    return(0);
  }
  

  printf(" Me quedo solo con el %d en el rango %f-%f\n",ivalue,min_f,max_f);

  nn=vector_i(*nbin);
  *p=vector_d(*nbin);
  *x=vector_d(*nbin);


  for(i=0;i<*nbin;i++) {
    (*p)[i]=0;
    nn[i]=0;
  }

    
  for(iseeing=0;iseeing<SF.nseeing;iseeing++) {  
    for(isky=0;isky<SF.nsky;isky++) {   
      for(itrans=0;itrans<SF.ntransparency;itrans++) {  
	for(isetup=0;isetup<SF.ninstsetup;isetup++) {  
	  for(iew=0;iew<SF.nEW;iew++) {
	    for(iz=0;iz<SF.nz;iz++) {
	      for(imag=0;imag<SF.nmag;imag++) {
		if(iproj==1 ) ibin=iseeing;
		if(iproj==2 ) ibin=isky;
		if(iproj==3 ) ibin=itrans;
		if(iproj==4 ) ibin=isetup;
		if(iproj==5 ) ibin=iew;
		if(iproj==6 ) ibin=iz;
		if(iproj==7 ) ibin=imag;
		if(iproj==8 ) {
		  cval=(SF.magbin[0]+imag*(SF.magbin[SF.nmag-1]-SF.magbin[0])/(SF.nmag-1))-2.5*log10((SF.ewbin[0]+iew*(SF.ewbin[SF.nEW-1]-SF.ewbin[0])/(SF.nEW-1)));
		  ibin=(int)((*nbin-1)*(cval-min)/(max-min)+0.5);
		}
		trueflag=1;
		if(iseeing_f) {if(!(iseeing==ivalue)) trueflag=0;}
		if(isky_f) {if(!(isky==ivalue)) trueflag=0;}
		if(itrans_f) {if(!(itrans==ivalue)) trueflag=0;}
		if(isetup_f) {if(!(isetup==ivalue)) trueflag=0;}
		if(iew_f) { if(!(iew==ivalue)) trueflag=0;}
		if(iz_f) {if(!(iz==ivalue)) trueflag=0;}
		if(imag_f) {if(!(imag==ivalue)) trueflag=0;}

		if(trueflag) {
		  ((*p)[ibin])+=SF.p[iseeing][isky][itrans][isetup][iew][iz][imag];
		  (nn[ibin])++;
		} else {
/* 		  printf(" He excludio \n"); */
/* 		  printf(" ise %d isk %d it %d ise %d iew %d iz %d im %d\n",iseeing,isky,itrans,isetup,iew,iz,imag); */
		}
	      }
	    }
	  }
	}
      }
    }
  }  

  
  for(i=0;i<*nbin;i++) {
    (*p)[i]=(*p)[i]/nn[i]; 
    if(nn[i]==0) (*p)[i]=0;
  }


  
  delta=(max-min)/(*nbin-1);
  
  for(i=0;i<*nbin;i++)  (*x)[i]=min+i*delta;

  ivalue++;
  return(1);
}


void project2posf(struct poselfunc SF, double***p, double**x, double**y, int *nx, int *ny, int iprojx, int iprojy) {
  
  int **nn;
  int i,j;
  int ibinx,ibiny;
  int iew,iz,imag,iseeing,isky,itrans,isetup;
  double deltax,minx,maxx;
  double deltay,miny,maxy;
  double cval;


  if(iprojx==1 ) { *nx=SF.nseeing;       minx=SF.seeing[0];maxx=SF.seeing[SF.nseeing-1];}
  if(iprojx==2 ) { *nx=SF.nsky   ;       minx=SF.sky[0];maxx=SF.sky[SF.nsky-1];}
  if(iprojx==3 ) { *nx=SF.ntransparency; minx=SF.transparency[0];maxx=SF.transparency[SF.ntransparency-1];}
  if(iprojx==4 ) { *nx=SF.ninstsetup ;   minx=1;maxx=SF.ninstsetup;}
  if(iprojx==5 ) { *nx=SF.nEW;           minx=SF.ewbin[0];maxx=SF.ewbin[SF.nEW-1];}
  if(iprojx==6 ) { *nx=SF.nz;            minx=SF.zbin[0];maxx=SF.zbin[SF.nz-1];}
  if(iprojx==7 ) { *nx=SF.nmag;          minx=SF.magbin[0];maxx=SF.magbin[SF.nmag-1];}
  if(iprojx==8 ) { *nx=SF.nmag;          minx=SF.magbin[0]-2.5*log10(SF.ewbin[SF.nEW-1]);maxx=SF.magbin[SF.nmag-1];}
  if(iprojx==8 ) { *nx=SF.nmag;          minx=SF.magbin[0]-2.5*log10(SF.ewbin[SF.nEW-1]);maxx=11;}
  if(iprojy==1 ) { *ny=SF.nseeing;       miny=SF.seeing[0];maxy=SF.seeing[SF.nseeing-1];}
  if(iprojy==2 ) { *ny=SF.nsky   ;       miny=SF.sky[0];maxy=SF.sky[SF.nsky-1];}
  if(iprojy==3 ) { *ny=SF.ntransparency; miny=SF.transparency[0];maxx=SF.transparency[SF.ntransparency-1];}
  if(iprojy==4 ) { *ny=SF.ninstsetup ;   miny=1;maxy=SF.ninstsetup;}
  if(iprojy==5 ) { *ny=SF.nEW;           miny=SF.ewbin[0];maxy=SF.ewbin[SF.nEW-1];}
  if(iprojy==6 ) { *ny=SF.nz;            miny=SF.zbin[0];maxy=SF.zbin[SF.nz-1];}
  if(iprojy==7 ) { *ny=SF.nmag;          miny=SF.magbin[0];maxy=SF.magbin[SF.nmag-1];}
  if(iprojy==8 ) { *ny=SF.nmag;          miny=SF.magbin[0]-2.5*log10(SF.ewbin[SF.nEW-1]);maxy=SF.magbin[SF.nmag-1];}
  if(iprojy==8 ) { *ny=SF.nmag;          miny=SF.magbin[0]-2.5*log10(SF.ewbin[SF.nEW-1]);maxy=11;}


  nn=matrix_i(*nx,*ny);
  *p=matrix_d(*nx,*ny);
  *x=vector_d(*nx);
  *y=vector_d(*ny);


  for(i=0;i<*nx;i++) {
    for(j=0;j<*ny;j++) {
      (*p)[i][j]=0;
      nn[i][j]=0;
    }
  }

    
  for(iseeing=0;iseeing<SF.nseeing;iseeing++) {  
    for(isky=0;isky<SF.nsky;isky++) {   
      for(itrans=0;itrans<SF.ntransparency;itrans++) {  
	for(isetup=0;isetup<SF.ninstsetup;isetup++) {  
	  for(iew=0;iew<SF.nEW;iew++) {
	    for(iz=0;iz<SF.nz;iz++) {
	      for(imag=0;imag<SF.nmag;imag++) {
		if(iprojx==1 ) ibinx=iseeing;
		if(iprojx==2 ) ibinx=isky;
		if(iprojx==3 ) ibinx=itrans;
		if(iprojx==4 ) ibinx=isetup;
		if(iprojx==5 ) ibinx=iew;
		if(iprojx==6 ) ibinx=iz;
		if(iprojx==7 ) ibinx=imag;
		if(iprojx==8 ) {
		  cval=(SF.magbin[0]+imag*(SF.magbin[SF.nmag-1]-SF.magbin[0])/(SF.nmag-1))-2.5*log10((SF.ewbin[0]+iew*(SF.ewbin[SF.nEW-1]-SF.ewbin[0])/(SF.nEW-1)));
		  ibinx=(int)((*nx-1)*(cval-minx)/(maxx-minx)+0.5);
		}
		if(iprojy==1 ) ibiny=iseeing;
		if(iprojy==2 ) ibiny=isky;
		if(iprojy==3 ) ibiny=itrans;
		if(iprojy==4 ) ibiny=isetup;
		if(iprojy==5 ) ibiny=iew;
		if(iprojy==6 ) ibiny=iz;
		if(iprojy==7 ) ibiny=imag;
		if(iprojy==8 ) {
		  cval=(SF.magbin[0]+imag*(SF.magbin[SF.nmag-1]-SF.magbin[0])/(SF.nmag-1))-2.5*log10((SF.ewbin[0]+iew*(SF.ewbin[SF.nEW-1]-SF.ewbin[0])/(SF.nEW-1)));
		  ibiny=(int)((*ny-1)*(cval-minx)/(maxx-minx)+0.5);
		}
/* 		printf(" ibin %d\n",ibin); */
/* 		kak=SF.p[iseeing][isky][itrans][isetup][iew][iz][imag]; */
		((*p)[ibinx][ibiny])+=SF.p[iseeing][isky][itrans][isetup][iew][iz][imag];
/* 		printf(" kak %f\n",kak); */
		(nn[ibinx][ibiny])++;
	      }
	    }
	  }
	}
      }
    }
  }  

  
  for(i=0;i<*nx;i++) {
    for(j=0;j<*ny;j++) {
      (*p)[i][j]=(*p)[i][j]/nn[i][j]; 
      if(nn[i][j]==0) (*p)[i][j]=0;
    }
  }


  deltax=(maxx-minx)/(*nx-1);
  for(i=0;i<*nx;i++)  (*x)[i]=minx+i*deltax;
  deltay=(maxy-miny)/(*ny-1);
  for(i=0;i<*ny;i++)  (*y)[i]=miny+i*deltay;
  
}


double prob_poselfunc_scale(struct poselfunc SF, struct SurveyItem si, double ew, double z, double magn) {

  int iew,iz,imag,iseeing,isky,itrans,isetup;
  double trans_sf,sky_sf;
  double facttrans,factsky;
  double magoffset;

  /* Parametros observacionales */
  iseeing =sf_selectbin(si.seeing      ,SF.seeing      ,SF.nseeing      );
  isky    =sf_selectbin(si.sky         ,SF.sky         ,SF.nsky         );
  itrans  =sf_selectbin(si.transparency,SF.transparency,SF.ntransparency);

  sky_sf  =SF.sky[isky];
  trans_sf=SF.transparency[itrans];

  facttrans=trans_sf/si.transparency;
  factsky=sqrt(sky_sf/(si.sky*facttrans));



  magoffset=-2.5*log10(facttrans*factsky);
  


  

  /* Parametros de las galaxias */
  iew =sf_selectbin(ew ,SF.ewbin ,SF.nEW );
  iz  =sf_selectbin(z  ,SF.zbin  ,SF.nz  );
  imag=sf_selectbin(magn - magoffset,SF.magbin,SF.nmag);

  if(iseeing==-2 || isky==-2 || itrans==-2 || iew==-2 || iz==-2 || imag==-2) {
    printf(" ft %f tr %f fs %f sk %f mo %f\n",facttrans,si.transparency,factsky,si.sky,magoffset);
    exit(1);
    printf(" iew %d iz %d imag %d\n",iew,iz,imag);
    printf(" iseeing %d isy %d itrans %d\n",iseeing,isky,itrans);
  }
  
  isetup  =sf_selectbin_char(si.instsetup  ,SF.instsetup   ,SF.ninstsetup   );

  return(SF.p[iseeing][isky][itrans][isetup][iew][iz][imag]);

}




double prob_poselfunc(struct poselfunc SF, struct SurveyItem si, double ew, double z, double magn) {

  int iew,iz,imag,iseeing,isky,itrans,isetup;

  /* Parametros de las galaxias */
  iew =sf_selectbin(ew  ,SF.ewbin ,SF.nEW );
  iz  =sf_selectbin(z   ,SF.zbin  ,SF.nz  );
  imag=sf_selectbin(magn,SF.magbin,SF.nmag);

  /* Parametros observacionales */
  iseeing =sf_selectbin(si.seeing      ,SF.seeing      ,SF.nseeing      );
  isky    =sf_selectbin(si.sky         ,SF.sky         ,SF.nsky         );
  itrans  =sf_selectbin(si.transparency,SF.transparency,SF.ntransparency);
  
  isetup  =sf_selectbin_char(si.instsetup  ,SF.instsetup   ,SF.ninstsetup   );

  return(SF.p[iseeing][isky][itrans][isetup][iew][iz][imag]);

}


double errprob_poselfunc_scale(struct poselfunc SF, struct SurveyItem si, double ew, double z, double magn) {

  int iew,iz,imag,iseeing,isky,itrans,isetup;
  double trans_sf,sky_sf;
  double facttrans,factsky;
  double magoffset;

  /* Parametros observacionales */
  iseeing =sf_selectbin(si.seeing      ,SF.seeing      ,SF.nseeing      );
  isky    =sf_selectbin(si.sky         ,SF.sky         ,SF.nsky         );
  itrans  =sf_selectbin(si.transparency,SF.transparency,SF.ntransparency);

  sky_sf  =SF.sky[isky];
  trans_sf=SF.transparency[itrans];

  facttrans=trans_sf/si.transparency;
  factsky=sqrt(sky_sf/(si.sky*facttrans));

  magoffset=-2.5*log10(facttrans*factsky);
  

  /* Parametros de las galaxias */
  iew =sf_selectbin(ew ,SF.ewbin ,SF.nEW );
  iz  =sf_selectbin(z  ,SF.zbin  ,SF.nz  );
  imag=sf_selectbin(magn - magoffset,SF.magbin,SF.nmag);

  
  isetup  =sf_selectbin_char(si.instsetup  ,SF.instsetup   ,SF.ninstsetup   );

  return(SF.errp[iseeing][isky][itrans][isetup][iew][iz][imag]);

}




double errprob_poselfunc(struct poselfunc SF, struct SurveyItem si, double ew, double z, double magn) {

  int iew,iz,imag,iseeing,isky,itrans,isetup;

  /* Parametros de las galaxias */
  iew =sf_selectbin(ew  ,SF.ewbin ,SF.nEW );
  iz  =sf_selectbin(z   ,SF.zbin  ,SF.nz  );
  imag=sf_selectbin(magn,SF.magbin,SF.nmag);

  /* Parametros observacionales */
  iseeing =sf_selectbin(si.seeing      ,SF.seeing      ,SF.nseeing      );
  isky    =sf_selectbin(si.sky         ,SF.sky         ,SF.nsky         );
  itrans  =sf_selectbin(si.transparency,SF.transparency,SF.ntransparency);
  
  isetup  =sf_selectbin_char(si.instsetup  ,SF.instsetup   ,SF.ninstsetup   );

  return(SF.errp[iseeing][isky][itrans][isetup][iew][iz][imag]);

}



int sf_selectbin(double value, double *bins, int nbin) {
  int i;
  
  if(nbin==1) return(0);
  if(value>=(bins[nbin-1]+bins[nbin-2])/2.) return(nbin-1);

  for(i=0;i<nbin;i++) {
    if(value<(bins[i]+bins[i+1])/2.) return(i);
  }
  
  printf(" sf_selectbin: Never should reach this point.\n Value %g  bin range %g-%g\n Exiting\n",value,bins[0],bins[nbin-1]);
  return(-2);
  exit(1);
}


int sf_selectbin_char(char value[], char **bins, int nbin) {
  int i;
  
  for(i=0;i<nbin;i++) {
    if(!strcmp(value,bins[i])) return(i);
  }
  printf(" sf_selectbin_char: No correspondency found for %s\n",value);
  exit(1);
}
  
