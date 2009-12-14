#include "modulos.h"

typedef int gen;

int GAMinFunc(int ndim, double *psol, double *pmin, double *pmax,
	      struct ga_set gs, double (*gafunc)(int, double *))
     
{
  //Genetic algorithm From Charbonneau ApJSS 101, 309
  //Translated into C and some minor changes
  
  //Si la función gafunc necesita información adicional, como puntos x, y, se lo pasaré con variables globales.
  //Esto es diferente a como lo hago para la función que llamo en Amoeba.
  //Por otra parte, pmin y pmax  marcan los límites en los parámetros de búsqueda, de modo que hay que normalizar
  //todas las llamadas a la función gafunc.
  //psol contiene la solución obtenida.

  int ip, k, ig;
  int newtot;
  int ip1, ip2;
  static long idum =-1;
  
  double *pscale, *pzero;
  double *scaledph;

  double **ph;  //Child phenotype. Contains two phenotype.
  double **oldph; //contains phenotype population at a given generation
                  //ph[i] contains individual i 
                  //ph[i][j] contains dimension j for individual i. It differs from the Fortran implementation
  double **newph; //Contains temporary array for the new population
  double *fitness;//fitness[i] contains fitness function value for phenotype i
  int *gen1, *gen2;//Child phenotypes encoding
  unsigned int *ifit, *jfit;//ifit[i] contains the position in array oldph for the ith phenotype. 
                            //jfit[j] contains the osrting index for phenotype j.
  //Matrix allocation
  gen1=vector_i(gs.ncod);
  gen2=vector_i(gs.ncod);
  
  ifit=vector_i(gs.npop);
  jfit=vector_i(gs.npop);
  
  fitness=vector_d(gs.npop);
  oldph=matrix_d(gs.npop, ndim);
  newph=matrix_d(gs.npop, ndim);
  ph=matrix_d(2, ndim);

  pscale=vector_d(ndim);
  pzero=vector_d(ndim);
  scaledph=vector_d(ndim);

  //Random number initialization
  if(idum==-1) idum=-(long)time(NULL)/2;
  //Compute scaling 
  for(k=0;k<ndim;k++) {
    pscale[k]=1;
    pzero[k]=0;
  }
  //Compute initial random population

  for(ip=0;ip<gs.npop;ip++) {
    for(k=0;k<ndim;k++) {
      oldph[ip][k]=ran2(&idum);
      scaledph[k]=pzero[k]+pscale[k]*oldph[ip][k];
    }
    fitness[ip]=gafunc(ndim, scaledph);
  }
  //Rank initial population
  ga_rnkpop(gs.npop, fitness, ifit, jfit);

  //Main generation loop
  for(ig=0;ig<gs.ngeneration;ig++) {
    newtot=0;
    for(ip=0;ip<gs.npop/2;ip++) {
      ga_select(gs.npop, jfit, gs.fdif, &ip1);
      do {
	ga_select(gs.npop, jfit, gs.fdif, &ip2);
      } while(ip2==ip1);
      //Encode parents
      ga_encode(ndim, gs.ncod, oldph[ip1], gen1);
      ga_encode(ndim, gs.ncod, oldph[ip2], gen2);
      //Breed
      ga_cross (ndim, gs.ncod, gs.pcross, gen1, gen2);
      ga_mutate(ndim, gs.ncod, gs.pmut, gen1);
      ga_mutate(ndim, gs.ncod, gs.pmut, gen2);
      //Decode offspring
      ga_decode(ndim, gs.ncod, ph[0], gen1); 
      ga_decode(ndim, gs.ncod, ph[1], gen2); 
      if(gs.irep==1) ga_genrep(ndim, ip, ph, newph);
      else {printf(" irep!=1 not implemented yet\n");exit(1);}
      //End of main population loop
    }
    if(gs.irep==1)  ga_newpop(gafunc, pscale, pzero, gs.ielite, ndim, gs.npop, oldph, newph, ifit, jfit, fitness, &newtot);
    if(gs.imut==2)  ga_adjmut(gs.npop, fitness, ifit, gs.pmutmin, gs.pmutmax, &gs.pmut);
    if(gs.verbose) {
      
    }
    //End of main generation loop
  }
  //Return best phenotype 
  for(k=0;k<ndim;k++) {
    psol[k]=pzero[k]+pscale[k]*oldph[ifit[gs.npop-1]][k];
  }

  //Freeing matrix
  free(scaledph);
  free(pzero);
  free(pscale);
  
  free_matrix_d(ph, 2, ndim);
  free_matrix_d(oldph, gs.npop, ndim);
  free_matrix_d(newph, gs.npop, ndim);
  free(fitness);

  free(ifit);
  free(jfit);

  free(gen1);
  free(gen2);

  return 0;
}


void ga_setdefault(struct ga_set *gs)
{
  gs->npop=100;
  gs->ngeneration=500;
  gs->ncod=8;
  gs->pcross=0.85;
  gs->imut=2;
  gs->pmut=0.005;
  gs->pmutmin=0.0005;
  gs->pmutmax=0.25;
  gs->fdif=1.;
  gs->irep=1;
  gs->ielite=1;
  gs->verbose=0;
}


void ga_newpop(double (*gafunc)(int, double *), double *pscale, double *pzero,
	       int ielite, unsigned int ndim,
	       unsigned int np, double **oldph, double **newph, unsigned int *ifit, unsigned int *jfit, 
	       double *fitns, unsigned int *nnew)
{
  int i,k;
  double *scaledph;
  scaledph=vector_d(ndim);
  *nnew=np;

  for(k=0;k<ndim;k++)   scaledph[k]=pzero[k]+pscale[k]*newph[0][k];

  if(ielite && gafunc(ndim,scaledph) < fitns[ifit[np-1]] ) {
    for(k=0;k<ndim;k++) {
      newph[0][k]=oldph[ifit[np-1]][k];
    }
    *nnew=*nnew-1;
  }
  for(i=0;i<np;i++) {
    for(k=0;k<ndim;k++) {
      oldph[i][k]=newph[i][k];
      scaledph[k]=pzero[k]+pscale[k]*oldph[i][k];
    }
    fitns[i]=gafunc(ndim,scaledph);
  }
  ga_rnkpop(np,fitns,ifit,jfit);
  free(scaledph);
}


void ga_genrep(unsigned int n , unsigned int ip, double **ph, double **newph)
{
  unsigned int i1, i2, k;
  i1=2*ip;
  i2=i1+1;
  for(k=0;k<n;k++) {
    newph[i1][k]=ph[0][k];
    newph[i2][k]=ph[1][k];
  }
}


void ga_rnkpop(unsigned int n, double  *arrin, unsigned int  *indx, unsigned int *rank)
{
  unsigned int i;
  indexx(n,arrin,indx);

  for(i=0;i<n;i++) {
    rank[indx[i]]=n-i;   //  Esto lo defino empezando en 0. He cambiado la definición de rtfit
  }
}


void ga_select(unsigned int np, unsigned int  *jfit, double fdif, int *idad)
{

  unsigned int np1,i;
  double dice, rtfit;
  static long idum =-1;
  
  if(idum==-1) idum=-(long)time(NULL)/2;
  np1  = np + 1;
  dice = ran2(&idum) *np * np1;
  rtfit=0;
  for(i=0;i<np;i++) {
    rtfit=rtfit+np1+fdif*(np - 1 - 2*jfit[i]);  // Definición cambiada porque jfit va de 0 a N-1
    if(rtfit>=dice) {
      *idad=i;
      break;
    }
  }
}


void ga_adjmut(unsigned int np, double *fitns, unsigned int *ifit, double pmutmin, double pmutmax, double *pmut)
{
  double rdif;
  const double rdiflo=0.05;
  const double rdifhi=0.25;
  const double delta =1.5;
  rdif=fabs(fitns[ifit[np-1]]-fitns[ifit[(np-1)/2]]) / (fitns[ifit[np-1]]+fitns[ifit[(np-1)/2]]);

  if(rdif<=rdiflo)      *pmut=minf(pmutmax, *pmut*delta);
  else if(rdif>=rdifhi) *pmut=maxf(pmutmin, *pmut/delta);
}

void  ga_cross(unsigned int ndim, unsigned int ncod, double pcross, int *gen1, int *gen2)
{

  unsigned int i; 
  int ispl, t;
  static long idum =-1;

  if(idum==-1) idum=-(long)time(NULL)/2;
  if(ran2(&idum)<pcross)  {
    ispl=(int)(ran2(&idum) *ndim *ncod) ;
    if(ispl>ndim*ncod-1) {printf(" ga_cross ERROR: ispl > ndim*ncod-1\n");exit(1);}
    for(i=ispl; i<ndim*ncod-1;i++) {
      t=gen2[i];
      gen2[i]=gen1[i];
      gen1[i]=t;
    }
  }
}

void ga_encode(unsigned int ndim, unsigned int ncod, double *ph, int *gen)
{  
  int i, j, ii;
  double z;
  long int ip;
  
  z=pow(10., ncod);
  ii=0;

  
  for(i=0;i<ndim;i++) {
    ip=(long int)(ph[i]*z);
    for(j=ncod-1;j>=0;j--) {
      gen[ii+j]=(int)(ip-10*floor(ip/10.));
      ip=(int)floor(ip/10.);
    }
    ii=ii+ncod;
  }
}

void ga_decode(unsigned int ndim, unsigned int ncod, double *ph, int *gen)
{
  int i, j, ii;
  long int ip;
  double z;

  z=pow(10.,-(float)(ncod));
  ii=0;
  for(i=0;i<ndim;i++) {
    ip=0;
    for(j=0;j<ncod;j++) {
      ip=10*ip+gen[ii+j];
    }
    ph[i]=ip*z;
    ii=ii+ncod;
  }
}



void ga_mutate(unsigned int ndim, unsigned int ncod, double pmut, int *gen)
{
  unsigned int i;
  static long idum =-1;
  
  if(idum==-1) idum=-(long)time(NULL)/2;
  for(i=0;i<ndim*ncod;i++) {
    if(ran2(&idum)<pmut) gen[i]=(int)(ran2(&idum)*10.);
  }
}

