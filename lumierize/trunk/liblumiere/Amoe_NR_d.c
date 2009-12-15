#include "amoeba.h"
#define DEBUG 0

#define ITMIN 5

int Amoe_NR_d(int ndata, double *xdata, double *ydata,
	      double *p, double *y, int ndim, double ftol,int itmax,
	      double (*amofunc)(int, double *, double *, double *))

{
  int mpts,ilo,ihi,inhi,i,j,mp,np;
  double alpha=1.0, beta=0.5, gamma=2.0;
  double *pr,*prr,*pbar;
  double rtol,ypr,yprr;
  double dum1;
  int iter=0;

  mp=ndim+1;
  np=ndim;

  if( (pr   = malloc(ndim*sizeof(double))) == NULL ||
      (prr  = malloc(ndim*sizeof(double))) == NULL ||
      (pbar = malloc(ndim*sizeof(double))) == NULL) {
       printf("Amoe_NR_d: ERROR. No puedo dimensionar la matrices\n");
       exit(1);
       }

    mpts=ndim+1;
    while(1) {
      ilo=0;
      if(y[0] > y[1]) {
        ihi=0;
        inhi=1;
        }
      else {
        ihi=1;
        inhi=0;
        }
      for(i=0; i<mpts; i++) {
        if(y[i] < y[ilo]) ilo=i;
        if(y[i] > y[ihi]) {
          inhi=ihi;
          ihi=i;
          }
        else {
          if(y[i] > y[inhi] && i != ihi) inhi=i;
          }
        }

      dum1=y[ihi]-y[ilo];
      rtol=2.*fabs(dum1)/(fabs(y[ihi])+fabs(y[ilo]));
      if(DEBUG)        printf(" it %d %f %f   rtol %g ftol %g\n",iter,y[ihi],y[ilo],rtol,ftol);

      if(rtol < ftol && iter > ITMIN) {
	if(DEBUG) printf(" Aqui debo haber salido\n");
	free(pr);
	free(prr);
	free(pbar);
	return(iter);
	}
      if(iter == itmax) {
/* 	printf("AMOEBA Exceeding maximum iterations (%d)\n",iter); */
        free(pr);
        free(prr);
        free(pbar);
	return(0);
	}
      iter++;
      for(j=0; j<ndim; j++) pbar[j]=0;
      for(i=0; i<mpts; i++) 
        if(i != ihi) for(j=0; j<ndim; j++) pbar[j] += p[i+j*mp];


      for(j=0; j<ndim; j++) {
	pbar[j] /= (double)ndim;
        pr[j]=(1.+alpha)*pbar[j]-alpha*p[ihi+j*mp];
	}

      ypr=(*amofunc)(ndata,xdata,ydata,pr);

      if(ypr < y[ilo]) {
	for(j=0; j<ndim; j++) prr[j]=gamma*pr[j]+(1.-gamma)*pbar[j];
        yprr=(*amofunc)(ndata,xdata,ydata,prr);
        if(yprr < y[ilo]) {
	  for(j=0; j<ndim; j++) p[ihi+j*mp]=prr[j];
          y[ihi]=yprr;
	  }
        else {
	  for(j=0; j<ndim; j++) p[ihi+j*mp]=pr[j];
          y[ihi]=ypr;
	  }
	}
      else {
        if(ypr > y[inhi]) {
          if(ypr < y[ihi]) {
	    for(j=0; j<ndim; j++) p[ihi+j*mp]=pr[j];
            y[ihi]=ypr;
	    }
	  for(j=0; j<ndim; j++) prr[j]=beta*p[ihi+j*mp]+(1.-beta)*pbar[j];
          yprr=(*amofunc)(ndata,xdata,ydata,prr); 
          if(yprr < y[ihi]) {
	    for(j=0; j<ndim; j++) p[ihi+j*mp]=prr[j];
            y[ihi]=yprr;
	    }
          else {
	    for(i=0; i<mpts; i++) 
	      if(i != ilo) {
	        for(j=0; j<ndim; j++) {
                  pr[j]=0.5*(p[i+j*mp]+p[ilo+j*mp]);
                  p[i+j*mp]=pr[j];
	 	 }
                y[i]=(*amofunc)(ndata,xdata,ydata,pr); 
	      }
            }
	  }
        else  {
	  for(j=1; j<ndim; j++) p[ihi+j*mp]=pr[j];
          y[ihi]=ypr;
	  }
        }
      }
    return(iter);
} 
