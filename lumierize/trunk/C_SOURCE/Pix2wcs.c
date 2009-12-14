#include "modulos.h"



//Variables para la astrometria



//
void ReadAstrom();
void DoAstro();
void ReadCat();
void WriteFile();

struct WorldCoor *wcsim;
int lhead,nbfits;  
char *header;
char image[51]="";
char outfile[51]="";
char infile[51]="";
float *xp,*yp;
double *ra,*dec;
int nobj;
int *logi;
int colx=14,coly=15;

int main(int argc, char **argv)
{
  
  printf(" Input file with X, Y coordinates: ");
  reads(infile,infile);
  printf(" Input column with X coordinate: ");
  colx=readi(colx);
  printf(" Input column with Y coordinate: ");
  coly=readi(coly);
  printf(" Input image with astrometric solution: ");
  reads(image,image);
  printf(" Output file with RA, DEC coordinates: ");
  reads(outfile,outfile);

  ReadAstrom();
  ReadCat();
  DoAstro();
  WriteFile();

  return 0;
}


void ReadAstrom()
{
  
  if ((header = fitsrhead (image, &lhead, &nbfits)) == NULL) {
    fprintf (stderr, "Cannot read FITS header of file %s\n", image);
    exit(1);
  }
  
  if((wcsim=wcsinit(header))==NULL) {
    printf(" No WCS information found in header\n Exiting");  
    exit(1);
  }
  else {
    printf(" WCS information from header:\n");
    PrintWCS(header,1);
  }

}

void DoAstro() {
  int j;
  double xpm,ypm;
  int off;
  

  off=0;

  
  for(j=0;j<nobj;j++) {
    if(logi[j]) {
      xpm=(double)xp[j];ypm=(double)yp[j];
      pix2wcs(wcsim,xpm-1, ypm-1 ,(ra+j), (dec+j));
    }
  }
}

void WriteFile() {
  
  int j;
  char sanul[32],sdnul[32];
  FILE *fout;

  if((fout=fopen(outfile,"w")) == NULL) {
    fprintf(stderr,"Pix2wcs: ERROR. No such parameter file %s\n", outfile);
    exit(1);
  }
  for(j=0;j<nobj;j++) {
    if(logi[j]) {
      ra2str(sanul,32,ra[j],3);
      dec2str(sdnul,32,dec[j],2);
      fprintf(fout," %8.1f %8.1f %s %s\n",xp[j],yp[j],sanul,sdnul);
      
    }
  }
  fclose(fout);
}








void ReadCat() {

  int j;  
  nobj=FileNLin(infile);
  if((xp=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((yp=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension yp of %d elements \n",nobj);
  if((logi=malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension logii of %d elements \n",nobj);
  if((ra=malloc(nobj*sizeof(double)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((dec=malloc(nobj*sizeof(double)))==NULL) printf("I cannot dimension yp of %d elements \n",nobj);

  for(j=0;j<nobj;j++)     logi[j]=0;
    
  ReadNumcol(infile,colx,xp,logi,&nobj);
  ReadNumcol(infile,coly,yp,logi,&nobj);

}
