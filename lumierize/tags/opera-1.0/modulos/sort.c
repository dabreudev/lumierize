
#include "modulos.h"

int sort2(n,ra)
int *n;
float *ra;
{
    static  int i__, j, l, ir;
    static float rra;
    int s;
    s=*n;
/*     //printf("hola %d\n",s); */


    /* Parameter adjustments */
    --ra;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    l = *n / 2 + 1;
    ir = *n;
L10:
    if (l > 1) {
	--l;
	rra = ra[l];
    } else {
	rra = ra[ir];
	ra[ir] = ra[1];
	--ir;
	if (ir == 1) {
	    ra[1] = rra;
	    return 0;
	}
    }
    i__ = l;
    j = l + l;
L20:
    if (j <= ir) {
	if (j < ir) {
	    if (ra[j] < ra[j + 1]) {
		++j;
	    }
	}
	if (rra < ra[j]) {
	    ra[i__] = ra[j];
	    i__ = j;
	    j += j;
	} else {
	    j = ir + 1;
	}
	goto L20;
    }
    ra[i__] = rra;
    goto L10;
} /* sort_ */



void Sort(int n,float *x,float *y)
{
  int i,j,k;
  y[0]=x[0];
/*   //printf("Pas %f %f %f %f %f\n",x[0],x[1],x[2],x[3],x[4],x[5]); */
  for(i=1;i<n;i++) {
/*     //    printf(" soritng %d\n",i); */
    for(j=0;j<i;j++) {
/*       //printf("%d %d Pu %f %f %f %f %f\n",i,j,y[0],y[1],y[2],y[3],y[4],y[5]); */
      if(x[i]<y[j]) {
	for(k=i-1;k>j-1;k--)         y[k+1]=y[k];
	y[j]=x[i];
	j=i;
/* 	//printf("Ya %f %f %f %f %f\n",y[0],y[1],y[2],y[3],y[4],y[5]); */
      }
    }
    if(x[i]>y[i-1]) y[i]=x[i];
/*     //printf("Fin %f %f %f %f %f\n",y[0],y[1],y[2],y[3],y[4],y[5]); */
  }
}

void sort(n,ra)
int n;
float ra[];
{
  int l,j,ir,i;
  float rra;
/*   //for(ir=1;ir<n+1;ir++) printf("ZZ %f ",ra[ir]); */
  l=(n >> 1)+1;
/*   //printf("N %d L %d\n",n,l);  */
  ir=n;
  for(;;) {
    if(l>1)
      rra=ra[--l];
    else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if(--ir==1) {
	ra[1]=rra;
	return;
      }
    }
    i=l;
    j=l << 1;
    while(j<=ir) {
      if(j<ir && ra[j]<ra[j+1]) ++j;
      if(rra<ra[j]) {
	ra[i]=ra[j];
	j+=(i=j);
      }
      else j=ir+1;
    }
/*     //for(dum=1;dum<n+1;dum++) printf("ZZ %f ",ra[dum]); */
/*     //printf("\n"); */
  }
  
}


void indexx(unsigned int n, double arr[], unsigned int indx[])
{
  const int nstack=50;
  const int m=7;
  unsigned int i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  double a;
  
  istack=vector_i(nstack);
  for (j=1;j<=n;j++) indx[j -1]=j;
  for (;;) {
    if (ir-l < m) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j -1];
	a=arr[indxt-1];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i -1] -1] <= a) break;
	  indx[i+1 -1]=indx[i -1];
	}
	indx[i+1 -1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack-- -1];
      l=istack[jstack-- -1];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k -1],indx[l+1 -1]);
      if (arr[indx[l+1 -1] -1] > arr[indx[ir -1] -1]) {
	SWAP(indx[l+1 -1],indx[ir -1])
	  }
      if (arr[indx[l -1] -1] > arr[indx[ir -1] -1]) {
	SWAP(indx[l -1],indx[ir -1])
	  }
      if (arr[indx[l+1 -1] -1] > arr[indx[l -1] -1]) {
	SWAP(indx[l+1 -1],indx[l -1])
	  }
      i=l+1;
      j=ir;
      indxt=indx[l -1];
      a=arr[indxt -1];
      for (;;) {
	do i++; while (arr[indx[i -1] -1] < a);
	do j--; while (arr[indx[j -1] -1] > a);
	if (j < i) break;
	SWAP(indx[i -1],indx[j -1])
	  }
      indx[l -1]=indx[j -1];
      indx[j -1]=indxt;
      jstack += 2;
      if (jstack > nstack) {
	printf("NSTACK too small in indexx.\n");
	exit(1);
      }
			
      if (ir-i+1 >= j-l) {
	istack[jstack -1]=ir;
	istack[jstack-1 -1]=i;
	ir=j-1;
      } else {
	istack[jstack -1]=j-1;
	istack[jstack-1 -1]=l;
	l=i;
      }
    }
  }
  free(istack);
  for (j=1;j<=n;j++) indx[j -1] = indx[j -1] -1; //Para que empiece en uno
}
