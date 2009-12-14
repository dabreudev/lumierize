#include "modulos.h"
#include <sys/stat.h>


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
/*   double alfac=0,deltac=0; */
/*   double platescale; */
/*   double rotang; */
/*   int nobj; */

/*   struct WorldCoor *wcsim=0;  */     /* World coordinate system structure */


  printf(" Output USNO catalog ");
  reads("",usnofile);


  TestRead(usnofile);  

  
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
  
  printf("\n");
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
    //printf(" ma %d\n",star.magetc); 
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
    star.rasec=(int)(ra[i]*360000.0);
    star.decsec=(int)(dec[i]*360000.0+32400000.);
    star.magetc=plate*1000000+((int) (mag[i]*10) % 1000)*1000 + ((int)(mag[i]*10) % 1000);
    //printf(" ra %d dec %d maget %d\n",star.rasec,star.decsec,star.magetc);
    //printf(" %d  / %d \n",i,nobj);
    writeuacstar(iuacstar,star,usnofile);
    nstars++;
  }

}


static int my_uacsra (double rax0, char usnofile[1000])
{
  int istar, istar1, istar2, nrep;
  double rax, ra1, ra, ram, rdiff, rdiff1, rdiff2, sdiff;
  UACstar star;       /* UA catalog entry for one star */
  char rastrx[16];
  int debug = 0;
  FILE *fu;

  if (!(fu = fopen (usnofile, "r"))) {
    fprintf (stderr,"UA zone catalog %s cannot be read\n",usnofile);
    return 0;
  }
  
  rax = rax0;
  if (debug)
    ra2str (rastrx, 16, rax, 3);
  istar1 = 1;
  if (my_uacstar (istar1, &star,fu))
    return (0);
  ra1 = my_uacra (star.rasec);
  if(rax<ra1) {
    fclose(fu);
    return(1);
  }
  istar = nstars;
  nrep = 0;
  while (istar != istar1 && nrep < 50) {
    if (my_uacstar (istar, &star, fu))
      break;
    else {
      ra = my_uacra (star.rasec);
      if (ra == ra1)
	break;
      if (debug) {
	char rastr[16];
	ra2str (rastr, 16, ra, 3);
	printf ("UACSRA %d %d: %s (%s)\n",
		 nrep,istar,rastr,rastrx);
      }
      rdiff = ra1 - ra;
      rdiff1 = ra1 - rax;
      rdiff2 = ra - rax;
      if (nrep > 40 && ABS(rdiff2) > ABS(rdiff1)) {
	istar = istar1;
	break;
      }
      nrep++;
      sdiff = (double)(istar - istar1) * rdiff1 / rdiff;
      if(debug)printf(" istar %d istar1 %d sfidd %f (%d)\n",istar,istar1,sdiff,(int) (sdiff + 0.5));
      istar2 = istar1 + (int) (sdiff + 0.5);
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
      if (istar < 1)
	istar = 1;
      if (istar > nstars)
	istar = nstars;
      if (istar == istar1)
	break;
    }
  }
  my_uacstar (istar, &star, fu);
  ra = my_uacra (star.rasec);
  if(my_uacstar (istar-1, &star, fu)) ram=-1;
  else                                ram = my_uacra (star.rasec);
  fclose(fu);
  if(debug) printf(" ra %f rax %f ra-1 %f\n",ra,rax,ram);
  if(ra>rax) {
    if(ram<rax) return(istar);
    else        return(istar-1);
  }
  else       return (istar+1);
}
