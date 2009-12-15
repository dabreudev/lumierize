#include "modulos.h"
void V(float v,float *z,float *d, float *m);
void Z(float *v,float z,float *d, float *m);
 
float H0;
float q0;

float c=299792.46;        /* En km/s */
  
int main()
{
  int col;
  int colc;
  char opt1[10];
  char opt2[10];
  char opt3='Y';
  char infile[100],outfile[100];
  float *zz=NULL,*vv=NULL,*dd=NULL,*in,*modd=NULL;
  float z,v;
  int *tesin;
  int n=0,i;
  int nlin;
  FILE *outfp;
  char *object=NULL;
  char *inc;
  printf(" Input H0 : ");
  H0=readf(50.);
  printf(" Input q0 : ");
  q0=readf(0.5);
/*   //scanf("%f",&H0); */
  printf(" output file: ");
  reads("/dev/null",outfile);
/*   //scanf("%s",outfile); */
  outfp=fopen(outfile,"w");
  printf(" 1  Calculate v,d from z\n");
  printf(" 2  Calculate z,d from v\n");
  printf(" 3  Calculate z,v from d\n\n");
  printf(" Input option :");
/*   //  scanf("%s",opt1);  */
  opt1[0]=readc('1');

  printf(" Input data by [F]ile or by [K]eyboard? ");
  opt2[0]=readc('K');
/*   //scanf("%s",opt2); */
  fflush(stdin);
  switch (opt2[0]){
  case 'F' :
    printf(" Input file: ");
    reads(infile,infile);
/*     //scanf("%s",infile); */
    printf(" Input column with object name: ");
    colc=readi(1);
/*     //scanf("%d",&colc); */
    printf(" Input column with data: ");
    col=readi(2);
/*     //scanf("%d",&col); */
/*     //printf("%d",FileNLin(infile)); */
    nlin=FileNLin(infile);
    in=malloc(nlin*sizeof(float));
    tesin=malloc(nlin*sizeof(int));
    for(i=0;i<nlin;i++) tesin[i]=0;
    printf("Reading data...\n");
    ReadNumcol(infile,col,in,tesin,&nlin);
    inc=malloc(nlin*51);
    printf(" Reading object names..\n");
    ReadCharcol(infile,colc,inc,tesin,51,&nlin);

    vv=malloc(nlin*sizeof(float));
    zz=malloc(nlin*sizeof(float));
    dd=malloc(nlin*sizeof(float));
    modd=malloc(nlin*sizeof(float));
    object=malloc(nlin*51);
    n=-1;
    switch (opt1[0]) {
    case '1' :
      for(i=0;i<nlin;i++) {
	if(tesin[i]) {
	  n++;
	  strcpy(object+n*51,inc+i*51);
	  zz[n]=in[i];
	  Z(vv+n,in[i],dd+n,modd+n);
	}
      }
      break;
    case '2' :
      for(i=0;i<nlin;i++) {
	if(tesin[i]) {
	  n++;
	  strcpy(object+n*51,inc+i*51);
	  vv[n]=in[i];
	  V(in[i],zz+n,dd+n,modd+n);
	}
      }
 
      break;
    case '3' :
      printf("Input d: ");
      z=readf(0);
/*       //scanf("%f",&z); */
      exit(1);
      break;
    }
    break;
  case 'K' :
   
    vv=malloc(1000*sizeof(float));
    zz=malloc(1000*sizeof(float));
    dd=malloc(1000*sizeof(float));
    modd=malloc(1000*sizeof(float));
    object=malloc(1000*51);
    n=-1;
    while(!(opt3=='N') && !(opt3=='n')) {
      n++;
      if(n!=0) strcpy(object+n*51,object+(n-1)*51);
      printf(" Input Object name: ");
      reads(object+n*51,object+n*51);
/*       //scanf("%s",object+n*51); */
      switch (opt1[0]) {
      case '1' :
	printf(" Input z: ");
	z=readf(0);
/* 	//scanf("%f",&z); */
	zz[n]=z;
/* 	//printf("De aquipasa\n"); */
	Z(vv+n,zz[n],dd+n,modd+n);
	break;
      case '2' :
	printf(" Input v: ");
	v=readf(0);
/* 	//scanf("%f",&v); */
	vv[n]=v;
	V(vv[n],zz+n,dd+n,modd+n);
	
	break;
      case '3' :
	printf(" Input d: ");
	z=readf(0);
/* 	//scanf("%f",&z); */
	exit(1);
	break;
      }
      printf("\n>> Z                : %f\n",zz[n]);
      printf(">> Velocity(km/s)   : %f\n",vv[n]);
      printf(">> Distance(Mpc)    : %f\n",dd[n]);
      printf(">> Distance moduli  : %f\n",modd[n]);


      printf("\n Again (Y/N) [Y]? ");
      opt3=readc('Y');
/*       //scanf("%s",&opt3); */
    }
    break;
  }
  printf("#Object                                                  z    velocity(km/s) distance(Mpc)  Dis_moduli\n");
  fprintf(outfp,"#Object                                                  z    velocity(km/s) distance(Mpc)  Dist_moduli\n");
  for(i=0;i<n+1;i++) {
    printf("%-51s %8.5f   %10.4f   %10.3f   %10.4f\n",object+i*51,zz[i],vv[i],dd[i],modd[i]);
    fprintf(outfp,"%-51s %8.5f  %10.4f  %10.3f  %10.4f\n",object+i*51,zz[i],vv[i],dd[i],modd[i]);
  }

  return(0);
}
void V(float v,float *z,float *d, float *m) 
{
  float dlum;
  *z=sqrt((1+v/c)/(1-v/c))-1;
  /*   *d=v/H0; */ /*Esto es la aproximación clásica */
  
  *d=c/H0*2*(2-2*q0*(1-*z)-(2-2*q0)*sqrt(1+2*q0**z))/(4*q0*q0*(1+*z));
  dlum=(1+*z)**d;
  *m=-5*log10(dlum/10e-6); /*  //10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
/*   printf(" Hogg %f\n",*m); */
/*   *m=-5*log10((1+*z)-sqrt(1+*z))+5*log10(H0)-53.89; */
/*   printf(" Salzer %f\n",*m); */
/*    Esta segunda es la formula de Salzer, pero dan lo mismo para q0=0.5 */
}
void Z(float *v,float z,float *d,float *m)
{
  float dlum;
  *v=(((1+z)*(1+z))-1)/(((1+z)*(1+z))+1)*c;
  *d=*v/H0; /* //Esto es la aproximación clásica */
/*   //printf(" Dista clasica: %f\n",*d); */
  *d=c/H0*2*(2-2*q0*(1-z)-(2-2*q0)*sqrt(1+2*q0*z))/(4*q0*q0*(1+z));
/*   //printf(" Dista relat: %f\n",*d); */
/*   //printf("a puecto %f\n",-5*log10((1+z)-sqrt(1+z))+5*log10(H0)-53.89); */
  dlum=(1+z)**d;
  *m=-5*log10(dlum/10e-6); /*  //10e6 son 10pc en megaparsecs.Segun Hogg astroph9905116 */
  printf(" Hogg %f\n",*m);
/*   *m=-5*log10((1+z)-sqrt(1+z))+5*log10(H0)-53.89; */
/*   printf(" Salzer %f\n",*m); */
/*   // Esta segunda es la formula de Salzer, y dan lo mismo para q0=0.5 */
}
