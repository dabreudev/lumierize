#include "modulos.h"
#include <sys/stat.h>


#define DEBUG 0
#define INTERACT 0

/* Formato de los ficheros .cat de USNO */
/* Cada estrella ocupa 12 bytes  */
/*     int sysref=WCS_J2000;  */     /* Catalog coordinate system */
/*     double eqref=2000.0;   */     /* Catalog equinox */

#define ABS(a) ((a) < 0 ? (-(a)) : (a))

typedef struct {
    int rasec, decsec, magetc;
} UACstar;

int nstars=0;

void TestRead(char file[1000]);
static int my_uacstar (int istar, UACstar  *star, FILE *fcat);
static void my_uacswap (char *string);
static double my_uacra (int rasec);
static double my_uacdec (int decsec);
static double my_uacmagb (int magetc);
static double my_uacmagr (int magetc);
static int my_uacplate (int magetc);
static int my_uacsra (double rax0, char usnofile[1000]);

static int writeuacstar ( int iuacstar, UACstar  star, char usnofile[1000]);
void ReadCat(char catfile[1000],int colx,int coly,int colmag,float **xp,float **yp,float **mag,int *nobj);
void DoAstro(  struct WorldCoor *wcsim, float *xp,float *yp,double **ra,double **dec,int nobj);
void WriteUSNOFile(char usnofile[1000],double *ra,double *dec,float *mag,int nobj);

int main() {


  char usnofile[1000];
  double alfac=0,deltac=0;
  double platescale;
  double rotang;
  int flip=0;
  double epoch=2000.0;
  int colx=2,coly=3,colmag=4;
  char catfile[1000];
  int naxes0,naxes1;
  float *xp,*yp,*mag;
  double *ra,*dec;
  int nobj;

  struct WorldCoor *wcsim=0;      /* World coordinate system structure */


  printf(" Output USNO catalog ");
  reads("",usnofile);

  printf(" Input RA central coordinate in hours.\n");
  alfac=readf(alfac/15.);
  alfac=alfac*15;
  printf(" Input DEC central coordinate in degrees.\n");
  deltac=readf(deltac);
  printf(" Input plate scale in arcsec/pix: ");
  platescale=readf(1);
  printf(" Input observation epoch (ej: 1998.5): ");
  epoch=(double)readf(epoch);
  printf(" Rotation ");
  rotang=readf(0.);
  printf(" Flipping?: ");
  flip=readi(flip);
  printf(" Input file with X, Y coordinates: ");
  reads("",catfile);
  printf(" Input column with X coordinate: ");
  colx=readi(colx);
  printf(" Input column with Y coordinate: ");
  coly=readi(coly);
  printf(" Input column with magnitude: ");
  colmag=readi(colmag);
  printf(" naxes X ");
  naxes0=readi(2048);
  printf(" naxes Y ");
  naxes1=readi(2048);

  wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(platescale),naxes0/2.,naxes1/2.,naxes0,naxes1,0.,2000,epoch,"TAN");
  if(!flip) wcsdeltset(wcsim, wcsim->cdelt[0],wcsim->cdelt[1],(double)rotang);
  else      wcsdeltset(wcsim,-wcsim->cdelt[0],wcsim->cdelt[1],(double)rotang);

  //TestRead(usnofile);  

  ReadCat(catfile,colx,coly,colmag,&xp,&yp,&mag,&nobj);
  DoAstro(wcsim,xp,yp,&ra,&dec,nobj);
  WriteUSNOFile(usnofile,ra,dec,mag,nobj);
  
  return 0;
}


static int my_uacstar (int istar, UACstar  *star, FILE *fcat)
{
  int nbs, nbr, nbskip;
  
  if (istar < 1 || istar > nstars) {
    printf ("UACstar %d is not in catalog\n",istar);
    return (-1);
  }
  nbskip = 12 * (istar - 1);
  if (fseek (fcat,nbskip,SEEK_SET))
    return (-1);
  nbs = sizeof (UACstar);
  nbr = fread (star, nbs, 1, fcat) * nbs;
  if (nbr < nbs) {
    printf ("UACstar %d / %d bytes read\n",nbr, nbs);
    return (-2);
  }
  my_uacswap ((char *)star);
  return (0);
}

static int writeuacstar ( int iuacstar, UACstar  star, char usnofile[1000])
{
  int nbs, nbr,nbskip;
  FILE *fu;
  struct stat statinfo;
  int sizebuf,sizew;
  char *buf;


  nbskip = 12 * (iuacstar - 1);
  stat(usnofile,&statinfo);
  fu=fopen(usnofile,"r+");

  fseek(fu,nbskip,SEEK_SET);
  sizebuf=statinfo.st_size-nbskip;
  buf=malloc(sizebuf*sizeof(char));
  fread(buf,1,sizebuf,fu);
  fseek(fu,nbskip,SEEK_SET);

  my_uacswap ((char *)(&star));
  nbs = sizeof (UACstar);
  nbr = fwrite(&star, nbs, 1, fu) *nbs;
  if (nbr < nbs) {
    fprintf (stderr, "UACstar %d / %d bytes read\n",nbr, nbs);
    return (-2);
  }
  
  sizew = fwrite(buf,1,sizebuf,fu);
  if (sizew < sizebuf) {
    fprintf (stderr, "Not all file has been written\n");
    return (-2);
  }
  
  fclose(fu);
  free(buf);
  return (0);
}


static void my_uacswap (char *string)
{
  char *sbyte, *slast;
  char temp0, temp1, temp2, temp3;
  int nbytes = 12; /* Number of bytes to reverse */
  
  slast = string + nbytes;
  sbyte = string;
  while (sbyte < slast) {
    temp3 = sbyte[0];
    temp2 = sbyte[1];
    temp1 = sbyte[2];
    temp0 = sbyte[3];
    sbyte[0] = temp0;
    sbyte[1] = temp1;
    sbyte[2] = temp2;
    sbyte[3] = temp3;
    sbyte = sbyte + 4;
  }
  return;
}


void TestRead(char file[100]) {

  FILE *fcat;
  struct stat statbuff;
  int i;
  UACstar star;
  int plate;
  double mag, magb;
  double ra,dec;
  char rastr[50],decstr[50];

/*   FILE *fout; */
 
  if (stat (file, &statbuff)) {
    fprintf (stderr,"UA zone catalog %s has no entries\n",file);
    return ;
  }
  else
    nstars = (int) statbuff.st_size / 12;
  
  /* Open zone catalog */
  if (!(fcat = fopen (file, "r"))) {
    fprintf (stderr,"UA zone catalog %s cannot be read\n",file);
    return ;
  }
/*   if (!(fout = fopen ("salida", "w"))) { */
/*     fprintf (stderr,"Wrting error in %s\n","salida"); */
/*     return ; */
/*   } */
  
  for(i=1;i<=nstars;i++) {
    my_uacstar(i,&star,fcat);
    mag = my_uacmagr (star.magetc); /* Red magnitude */
    plate = my_uacplate (star.magetc);
    ra = my_uacra (star.rasec);
    dec = my_uacdec (star.decsec);
    magb = my_uacmagb (star.magetc);
    ra2str(rastr,32,ra,3);
    dec2str(decstr,32,dec,2);
    printf (" %04d.%08d %s %s %s %5.2f %5.2f %d\n",825,i,rastr,decstr,"J2000",magb,mag,plate);
    printf(" ma %d\n",star.magetc); 
  }
  
  fclose(fcat);

}


static int my_uacplate (int magetc)
{
    if (magetc < 0)     return ((-magetc / 1000000) % 1000);
    else                return ((magetc / 1000000) % 1000);
}
static double my_uacmagb (int magetc)
{
    if (magetc < 0)     return ((double) ((-magetc / 1000) % 1000) * 0.1);
    else                return ((double) ((magetc / 1000) % 1000) * 0.1);
}
static double my_uacmagr (int magetc)
{
    if (magetc < 0)     return ((double) (-magetc % 1000) * 0.1);
    else                return ((double) (magetc % 1000) * 0.1);
}
static double my_uacdec (int decsec)
{
    return ((double) (decsec - 32400000) / 360000.0);
}
static double my_uacra (int rasec)
{
    return ((double) (rasec) / 360000.0);
}


void ReadCat(char catfile[1000],int colx,int coly,int colmag,float **xp,float **yp,float **mag,int *nobj) {

  int *logi;
  int i;

  *nobj=FileNLin(catfile);
  *xp=vector_f(*nobj);
  *yp=vector_f(*nobj);
  *mag=vector_f(*nobj);
  logi=vector_i(*nobj);

  for(i=0;i<*nobj;i++)     logi[i]=0;
  
  ReadNumcol(catfile,colx  ,*xp ,logi,nobj);
  ReadNumcol(catfile,coly  ,*yp ,logi,nobj);
  ReadNumcol(catfile,colmag,*mag,logi,nobj);


  for(i=0;i<*nobj;i++) {
    if(!logi[i]) {
      memmove(*xp  +i,*xp  +i+1,(*nobj-i-1)*sizeof(float));
      memmove(*yp  +i,*yp  +i+1,(*nobj-i-1)*sizeof(float));
      memmove(*mag +i,*mag +i+1,(*nobj-i-1)*sizeof(float));
      memmove( logi+i, logi+i+1,(*nobj-i-1)*sizeof(int));
      (*nobj)--;
      i--;
    }
  }

}


void DoAstro(  struct WorldCoor *wcsim, float *xp,float *yp,double **ra,double **dec,int nobj) {

  int j;
  double xpm,ypm;
  int off;
  off=0;
  
  *ra=vector_d(nobj);
  *dec=vector_d(nobj);
  printf(" Aqui nobjs %d\n",nobj);
  for(j=0;j<nobj;j++) {
    xpm=(double)xp[j];ypm=(double)yp[j];
    pix2wcs(wcsim,xpm, ypm ,(*ra+j), (*dec+j));
  }
  
  printf(" Salgo do\n");
}


void WriteUSNOFile(char usnofile[1000],double *ra,double *dec,float *mag,int nobj) {

  int i;
  int plate=666;
  int iuacstar;

  UACstar star;


  for(i=0;i<nobj;i++) {
    iuacstar=my_uacsra(ra[i],usnofile);
    if(DEBUG) printf(" SRA devuelve %d\n",iuacstar);
    star.rasec=(int)(ra[i]*360000.0);
    star.decsec=(int)(dec[i]*360000.0+32400000.);
    star.magetc=plate*1000000+((int) (mag[i]*10) % 1000)*1000 + ((int)(mag[i]*10) % 1000);
    //printf(" ra %d dec %d maget %d\n",star.rasec,star.decsec,star.magetc);
    if(DEBUG) printf(" %d  / %d \n",i,nobj);
    writeuacstar(iuacstar,star,usnofile);
    if(INTERACT) i=readi(i);
    nstars++;
  }

}


static int my_uacsra (double rax0, char usnofile[1000])
{
  int istar, istar1, istar2, nrep;
  double rax, ra1, ra, rau, ram, rdiff, rdiff1, rdiff2, sdiff;
  UACstar star;       /* UA catalog entry for one star */
  char rastrx[16];
  int debug = DEBUG;
  FILE *fu;
  struct stat statbuff;
  int istarmin,istarmax;
  static long idum=-1;

  if (stat (usnofile, &statbuff)) {
    fprintf (stderr,"UA zone catalog %s has no entries\n",usnofile);
    return 0;
  }
  else
    nstars = (int) statbuff.st_size / 12;


  if (!(fu = fopen (usnofile, "r"))) {
    fprintf (stderr,"UA zone catalog %s cannot be read\n",usnofile);
    return 0;
  }
  
  rax = rax0;
  if (debug)
    ra2str (rastrx, 16, rax, 3);
  istar1 = 1;
  if (my_uacstar (istar1, &star,fu)) {
    if(debug) printf("Salgo Con el primer error\n");
    return (0);
  }
  ra1 = my_uacra (star.rasec);
  if(rax<ra1) {
    fclose(fu);
    if(debug) printf("Salgo primero ra1 %f rax %f\n",ra1,rax);
    return(1);
  }
  if (my_uacstar (nstars, &star,fu)) {
    if(debug) printf("Salgo Con el segundo error\n");
    return (0);
  }
  rau = my_uacra (star.rasec);
  if(rax>rau) {
    fclose(fu);
    if(debug) printf("Salgo ultimo  rau %f rax %f\n",rau,rax);
    return(nstars+1);
  }
  istar = nstars;
  nrep = 0;
  istarmin=1;
  istarmax=nstars;
  while (  nrep < 50 && istarmax-istarmin > 1) {
    if (my_uacstar (istar, &star, fu)) {
      printf(" Salgo por este fallo nrep %d\n",nrep);
      break;
    }
    else {
      ra = my_uacra (star.rasec);
      if(ra<rax && istar > istarmin) istarmin=istar;
      if(ra>rax && istar < istarmax) istarmax=istar;
      if (debug) {
	char rastr[16];
	ra2str (rastr, 16, ra, 3);
	printf ("UACSRA %d %d: %s (%s)\n",
		 nrep,istar,rastr,rastrx);
      }
      rdiff = ra1 - ra;
      rdiff1 = ra1 - rax;
      rdiff2 = ra - rax;
      sdiff = (double)(istar - istar1) * rdiff1 / rdiff;
      if(debug)printf(" istar %d istar1 %d sfidd %f (%d)\n",istar,istar1,sdiff,(int) (sdiff + 0.5));
      istar2 = istar1 + (int) (sdiff + 0.5);
      if(rdiff == 0 && nrep == 1 ) istar2 = istar1 - 1;
      if(rdiff == 0 && nrep == 2 ) istar2 = istar1 + 1;
      if(rdiff == 0 && nrep > 15 ) istar2 = istarmin + (istarmax-istarmin)*ran2(&idum);
      ra1 = ra;
      istar1 = istar;
      istar = istar2;
      if (debug) {
	printf (" ra1=    %.5f ra=     %.5f rax=    %.5f\n",
		 ra1,ra,rax);
	printf (" rdiff=  %.5f rdiff1= %.5f rdiff2= %.5f\n",
		 rdiff,rdiff1,rdiff2);
	printf (" istar1= %d istar= %d istar2= %d\n",
		 istar1,istar,istar2);
      }
      if (istar < istarmin)
	istar = istarmin + (istarmax-istarmin)*ran2(&idum);
      if (istar > istarmax)
	istar = istarmin + (istarmax-istarmin)*ran2(&idum);
      nrep++;
    }
  }
  my_uacstar (istar, &star, fu);
  ra = my_uacra (star.rasec);
  if(my_uacstar (istar-1, &star, fu)) ram=-1;
  else                                ram = my_uacra (star.rasec);
  fclose(fu);
  if(debug) printf(" ra %f rax %f ra-1 %f\n",ra,rax,ram);
  if(debug) printf(" istarmin %d istarmax %d\n",istarmin,istarmax);
  if(ra>rax) {
    if(ram<rax) return(istar);
    else        return(istar-1);
  }
  else       return (istar+1);
}
