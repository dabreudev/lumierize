#include "amoeba.h"

int Amoe_NR(int ndata, float *xdata, float *ydata,
	    float *p, float *y, int ndim, float ftol,int itmax,
	    float (*amofunc)(int, float *, float *, float *))

{
  int mpts,ilo,ihi,inhi,i,j,mp,np;
  float alpha=1.0, beta=0.5, gamma=2.0;
  float *pr,*prr,*pbar;
  float rtol,ypr,yprr;
  float dum1;
  int iter=0;

  mp=ndim+1;
  np=ndim;

  if( (pr   = malloc(ndim*sizeof(float))) == NULL ||
      (prr  = malloc(ndim*sizeof(float))) == NULL ||
      (pbar = malloc(ndim*sizeof(float))) == NULL) {
       printf("Amoe_NR: ERROR. No puedo dimensionar la matrices\n");
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
      if(rtol < ftol) {
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
	pbar[j] /= (float)ndim;
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
