#include "modulos.h"


void Plotpts();

struct cat {
  char catfile[51];
  int colra,coldec;
  double *ar,*dec;
  int  nobj;
  int *logi;
  int initflag;
};

int pg1,pg2;
long s_naxes[1], s_pixels;
float especmin,especmax;
float closest;
char rootspec[51];
fitsfile *spec;
char specfile[51];

float *spectrum,*x;
float fnul;
float xcu,ycu;

int num=0;
double armin,armax,decmin,decmax;
char snul[50];
FILE *filecat;
float xbmin,xbmax,ybmin,ybmax;
float xfmin,xfmax,yfmin,yfmax;
char opt[10];
char option='C';
float flux;
char cnul;
float imin,imax,jmin,jmax;
/* Variables del dibujo */
float factorarr=1;
float trans[6];
float mean, sigma;

float fnul1,fnul2;

struct cat *pc;
int ncats=0;

void InitCat(int i);

int main()
{
  int i;
  trans[0]=0;trans[1]=1;trans[2]=0;trans[3]=0;trans[4]=0;trans[5]=1;
  
  ncats=1;
  if((pc=(struct cat*) malloc(1*sizeof(struct cat)))==NULL) {
    printf("I cannot dimension plotcats       of %d bytes \n",11*sizeof(struct cat));
  }
  pc[0].initflag=0;
  pc[0].colra=1;
  pc[0].coldec=2;
  pc[0].nobj=0;

  InitCat(0);

  printf(" ...searching limits\n");

  pgLimits_d(pc[0].nobj,pc[0].ar,&fnul1,&fnul2);
  armin=fnul1;armax=fnul2;
  pgLimits_d(pc[0].nobj,pc[0].dec,&fnul1,&fnul2);
  decmin=fnul1;decmax=fnul2;
  

  pg1=cpgopen("?");
  cpgscf(2);
  cpgask(0);
  

  printf(" Drawing box [%f:%f,%f:%f]\n",armin,armax,decmin,decmax);
  
  cpgswin(armax*3600,armin*3600,decmin*3600,decmax*3600);
  cpgtbox("ZXBCTNSH",0.0,0,"ZBCTNSD",0.0,0);
  cpgswin((float)armax,(float)armin,(float)decmin,(float)decmax);
  

  /* Termino de leerlo */
  cpgsci(1);
  num=0;
  Plotpts(0);
  cpgsci(1);


  
  
/*   //exit(1); */
 /* Menu principal */

  while(opt[0]!='E') {
    
    printf("\n\n  Z Zoom image with cursor\n");
    printf("  X Zoom image with keyboard\n");
    printf("  O Zoom out\n");
    printf("  N Write numbers [%d]\n",num);
    printf("  A Add other catalog\n");
    printf("  E Exit\n");
    
    fflush(NULL);
    fflush(stdin);
    fflush(stdout);
    scanf("%s",opt);
    switch (opt[0]) {
    case 'A':
    case 'a':
      printf(" Current Buffers: \n");
      for(i=0;i<ncats;i++) {
	if(pc[i].initflag) printf(" %3d  File: %s  Ra col: %d Dec col: %d\n",i+1,pc[i].catfile,pc[i].colra,pc[i].coldec);
	else printf(" %3d  Not Loaded \n",i+1);
      }
      printf(" Input buffer ");
      i=readi(1);
      InitCat(i-1);
      break;
    case 'Z' :
    case 'z':
      cpgsci(2);
      printf(" Press bottom left square with mouse...\n");
      cpgcurs(&xbmin,&ybmin,&cnul);
      printf(" Press top right square with mouse...\n");
      cpgband(2,1,xbmin,ybmin,&xbmax,&ybmax,&cnul);
      armin=xbmax;armax=xbmin;
      decmin=ybmin;decmax=ybmax;
      cpgsci(1);
      break;
    case 'O' :
    case 'o':
      pgLimits_d(pc[0].nobj,pc[0].ar,&fnul1,&fnul2);
      armin=fnul1;armax=fnul2;
      pgLimits_d(pc[0].nobj,pc[0].dec,&fnul1,&fnul2);
      decmin=fnul1;decmax=fnul2;
      break;
    case 'X' :
    case 'x':
      printf("RA min: ");
      armin=readd(armin);
      printf("RA max: ");
      armax=readd(armax);
      printf("DEC min: ");
      decmin=readd(decmin);
      printf("DEC max: ");
      decmax=readd(decmax);
      break;
    case 'E' :
    case 'e':
      cpgend();
      exit(0);
      break;
    case 'N' :
    case 'n':
      if(num) num=0;	
      else num=1;
      break;
    }
    cpgpage();
    cpgsvp(0.1,0.9,0.1,0.9);
/*     //    cpgswin(0.,1.,0.,1.); */
/*     //    cpgwnad(xbmin,xbmax,ybmin,ybmax); */
    cpgsch(1.);
/*     //cpgswin(armax*3600,armin*3600,decmin*3600,decmax*3600); */
    cpgsvp(0.1,0.9,0.1,0.9);
/*     //cpgswin(armax,armin,decmin,decmax); */
    cpgwnad(armax*3600.*15.,armin*3600.*15.,decmin*3600.,decmax*3600);
    cpgswin(armax*3600.,armin*3600.,decmin*3600.,decmax*3600);
    cpgtbox("ZXBCTNSH",0.0,0,"ZBCTNSD",0.0,0);
    cpgswin((float)armax,(float)armin,(float)decmin,(float)decmax);



    /*    printf(" Drawing subimage [%f:%f,%f:%f]\n",xbmin,xbmax,ybmin,ybmax);


    cpggray(buffer,naxes[0],naxes[1],xbmin,xbmax,ybmin,ybmax,datamax,datamin,trans)
;
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);    
    cpgsci(2);
    */
    printf(" Drawing box [%f:%f,%f:%f]\n",armin,armax,decmin,decmax);

    for(i=0;i<ncats;i++) {
      if(pc[i].initflag)    Plotpts(i);
      else printf(" %3d  Not Loaded \n",i+1);
    }

    
  }
  

  return(0);


  

}


void Plotpts(int icat) {
  int j;
  
  if(pc[icat].initflag) {
    cpgsci(1);
    
    for(j=0;j<pc[icat].nobj;j++) {
      if(pc[icat].logi[j]) {
	cpgsch(0.8);
	cpgsci(icat+1);
	cpgpt1((float)(pc[icat].ar[j]),(float)(pc[icat].dec[j]),4+icat);
	/*       sprintf(snul,"%d",(int)iobj[i]); */
	sprintf(snul,"%d",j);
	if(num) cpgtext((float)(pc[icat].ar[j]),(float)(pc[icat].dec[j]),snul);
	
      }
    }
  }
  cpgsci(1);
}

void InitCat(int i) {
  int j;

  if(i>ncats-1) {
    printf(" Aqui debiste i %d\n",i);
    if((pc=(struct cat*) realloc(pc, (i+1)*sizeof(struct cat)))==NULL) {
      printf("I cannot dimension pc of %d bytes \n",(i+1)*sizeof(struct cat));
    }
    for(j=ncats;j<i+1;j++) {
      pc[j].initflag=0;
      pc[j].colra=1;
      pc[j].coldec=2;
      pc[j].nobj=0;
    }
    ncats=i+1;
  }
    
  printf(" Input catalogue file: ");
  reads(pc[i].catfile,pc[i].catfile);
  printf(" Column with RA: ");
  pc[i].colra=readi(pc[i].colra);
  printf(" Column with DEC: ");
  pc[i].coldec=readi(pc[i].coldec);
  pc[i].nobj=FileNLin(pc[i].catfile);
  if(pc[i].initflag) {
    free(pc[i].dec);
    free(pc[i].ar);
    free(pc[i].logi);
  }
  if((pc[i].dec=malloc(pc[0].nobj*sizeof(double)))==NULL) printf("I cannot dimension dec of %d elements \n",pc[i].nobj);
  if((pc[i].ar=malloc(pc[0].nobj*sizeof(double)))==NULL) printf("I cannot dimension ar of %d elements \n",pc[i].nobj);
  if((pc[i].logi=malloc(pc[0].nobj*sizeof(int)))==NULL) printf("I cannot dimension logi of %d elements \n",pc[i].nobj);
  pc[i].initflag=1;
  
  printf("Number of objects %d\n",pc[i].nobj);
  
  ReadWCScol(pc[i].catfile, pc[i].colra,pc[i].ar,pc[i].logi,&(pc[i].nobj));
  
  ReadWCScol(pc[i].catfile, pc[i].coldec,pc[i].dec,pc[i].logi,&(pc[i].nobj));
  
}
