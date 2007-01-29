#include "modulos.h"

#define DEBUG 0

#define TOL 1.0e-18
static double sqrarg;
static int iminarg1,iminarg2;
static double maxarg1,maxarg2;

#define DMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?  (maxarg1) : (maxarg2))
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?   (iminarg1) : (iminarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)




void svdfit_d( double *x, double *y, double *sig, int n, double *a, int ma, double **cvm, double *chisq,void (*funcs)(double,double [],int)) {


  int j,i;
  double wmax,tmp,thresh,sum,*b,*afunc;
  double **u,**v,*w;


  if((b=malloc(n*sizeof(double))) == NULL) {
    printf("I cannot dimension b of %d elements \n",n);
    exit(1);
  }
  if((afunc=malloc(ma*sizeof(double))) == NULL) {
    printf("I cannot dimension afunc of %d elements \n",ma);
    exit(1);
  }
  u=matrix_d(n,ma);
  v=matrix_d(ma,ma);
  w=vector_d(ma);

  for (i=0;i<n;i++) {
    (*funcs)(x[i],afunc,ma);
    if(DEBUG) printf(" x %f s %f",x[i],sig[i]);
    tmp=1.0/sig[i];
    for (j=0;j<ma;j++) 
      u[i][j]=afunc[j]*tmp;
    if(DEBUG) printf(" a0 %f a1 %f\n",u[i][0],u[i][1]);
    b[i]=y[i]*tmp;
  }
  if(DEBUG) printf("Terino\n");
  svdcmp_d(u,n,ma,w,v);
  if(DEBUG) printf("Paso svdcmp\n");
  if(DEBUG) for(j=0;j<ma;j++) printf("PRI %d w %f\n",j,w[j]);
  wmax=0.0;
  for(j=0;j<ma;j++) 
    if(w[j]>wmax) wmax=w[j];
  thresh=TOL*wmax;
  if(DEBUG) printf(" thres %f\n",thresh);
  for(j=0;j<ma;j++)
    if(w[j]<thresh) w[j]=0.0;
  if(DEBUG) for(j=0;j<ma;j++) printf(" %d w %f\n",j,w[j]);
  if(DEBUG) for(i=0;i<ma;i++) for(j=0;j<ma;j++) printf(" %d %d v %f\n",i,j,v[i][j]);
  if(DEBUG) for(i=0;i<n;i++) for(j=0;j<ma;j++) printf(" %d %d u %f\n",i,j,u[i][j]);
  if(DEBUG) for(i=0;i<n;i++) printf(" %d %d b %f\n",i,j,b[i]);
  svbksb_d(u,w,v,n,ma,b,a);
  *chisq=0.0;
  for(i=0;i<n;i++) {
    (*funcs)(x[i],afunc,ma);
    for(sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
    *chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
  }

  svdvar_d(v,ma,w,cvm);

  free(afunc);
  free(b);
  free(w);
  free_matrix_d(u,n,ma);
  free_matrix_d(v,ma,ma);
  
}


void svdvar_d(double **v,int ma, double w[],double **cvm)
{
  int k,j,i;
  double sum,*wti;

  if((wti=malloc(ma*sizeof(double))) == NULL) {
    printf("I cannot dimension wti of %d elements \n",ma);
    exit(1);
  }

  for(i=0;i< ma ; i++) {
    wti[i]=0.0;
    if (w[i]) wti[i]=1.0/(w[i]*w[i]);
  }
  for(i=0;i< ma ; i++) {
    for(j=0;j<=i ; j++) {
      for(sum=0.0,k=0;k<ma;k++) sum += v[i][k]*v[j][k]*wti[k];
      cvm[j][i]=cvm[i][j]=sum;
    }
  }
  free(wti);
 

}


void svdcmp_d(double **a,int m,int n,double w[],double **v) {

  int flag,i,its,j,jj,k,l=0,nm=0;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  
  if(DEBUG) {
    for(i=1;i<=m;i++) {
      for(k=1;k<=n;k++) {
	printf("EN a  %f ",a[i-1][k-1]);
      }
      printf("\n");
    }
    printf("PASo la primera\n");
  }
  rv1=vector_d(n);
  g=scale=anorm=0.0;
  for(i=1;i<=n;i++ ) {
    l=i+1;
    rv1[i-1]=scale*g;
    g=s=scale=0.0;
    if(i<=m) {
      for(k=i;k<=m;k++) scale += fabs(a[k-1][i-1]);
      if (scale) {
	for(k=i;k<=m;k++) {
	  a[k-1][i-1] /= scale;
	  s += a[k-1][i-1] *a[k-1][i-1];
	}
	f=a[i-1][i-1];
	g=-SIGN(sqrt(s),f);
	h=f*g-s;
	a[i-1][i-1]=f-g;
	for(j=l;j<=n;j++) {
	  for(s=0.0,k=i;k<=m;k++) s += a[k-1][i-1]*a[k-1][j-1];
	  f=s/h;
	  for(k=i;k<=m;k++) a[k-1][j-1] += f*a[k-1][i-1];
	}
	for(k=i;k<=m;k++) a[k-1][i-1] *= scale;
      }
    }
    w[i-1]=scale*g;
    g=s=scale=0.0;
    if(DEBUG) printf(" W[i-1] %f\n",w[i-1]);
    if(i<=m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i-1][k-1]);
      if(DEBUG) printf(" Antes scale \n");
      if (scale) {
	if(DEBUG) printf(" Si scale \n");
	for(k=l;k<=n;k++) {
	  a[i-1][k-1] /= scale;
	  s +=  a[i-1][k-1]* a[i-1][k-1];
	}
	f= a[i-1][l-1];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i-1][l-1]=f-g;
	for(k=l;k<=n;k++) {
	  rv1[k-1]= a[i-1][k-1]/h;
	  if(DEBUG) printf(" %d rv1 %f\n",k-1,rv1[k-1]); 
	}
	for(j=l;j<=m;j++) {
	  for(s=0.0,k=l;k<=n;k++) s+= a[j-1][k-1]* a[i-1][k-1];
	  for(k=l;k<=n;k++)  a[j-1][k-1] += s*rv1[k-1];
	}
	for(k=l;k<=n;k++)  a[i-1][k-1] *=scale;
      }
    }
    anorm=DMAX(anorm,(fabs(w[i-1])+fabs(rv1[i-1])));
  }
  if(DEBUG) {
    for(i=1;i<=m;i++) {
      for(k=1;k<=n;k++) {
	printf("a  %f ",a[i-1][k-1]);
      }
      printf("\n");
    }
    printf("PASo la primera\n");
  }
  
  for(i=n;i>=1;i--) {
    if(DEBUG) printf(" i %d\n",i);
    if (i<n) {
      if (g) {
	for(j=l;j<=n;j++) 
	  v[j-1][i-1]=( a[i-1][j-1]/ a[i-1][l-1])/g;
	for(j=l;j<=n;j++) {
	  for(s=0.0,k=l;k<=n;k++) s+= a[i-1][k-1] * v[k-1][j-1];
	  for(k=l;k<=n;k++)  v[k-1][j-1] += s*v[k-1][i-1];
	}
      }
      for(j=l;j<=n;j++) v[i-1][j-1]=v[j-1][i-1]=0.0;
    }
    v[i-1][i-1]=1.0;
    g=rv1[i-1];
    l=i;
  }
  if(DEBUG) {
    for(i=1;i<=n;i++) {
      for(k=1;k<=n;k++) {
	printf("v  %f ",v[i-1][k-1]);
      }
      printf("\n");
    }
    printf("Paso la segubdo\n");
  }
  for(i=IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i-1];
    for(j=l;j<=n;j++) a[i-1][j-1]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for(s=0.0,k=l;k<=m;k++) s+= a[k-1][i-1]* a[k-1][j-1];
	f=(s/a[i-1][i-1])*g;
	for(k=i;k<=m;k++)  a[k-1][j-1] += f*a[k-1][i-1];
      }
      for(j=i;j<=m;j++) a[j-1][i-1] *=g;
    } else for (j=i;j<=m;j++) a[j-1][i-1] =0.0;
    ++a[i-1][i-1];
  }
  if(DEBUG) {
    for(i=1;i<=n;i++) {
      for(k=1;k<=n;k++) {
	printf("v  %f ",v[i-1][k-1]);
      }
      printf("\n");
    }
    printf("Paso la terce \n");
  }
  if(DEBUG) printf(" ANTES ULTIMO\n");
  for(k=n;k>=1;k--) {
    if(DEBUG) printf(" ESTOY EN k %d \n",k);
    for (its=1;its<=30;its++) {
      if(DEBUG) printf(" Itera %d\n",its);
      flag=1;
      for(l=k;l>=1;l--) {
	nm=l-1;
	if ((double)(fabs(rv1[l-1])+anorm) == anorm) {
	  if(DEBUG) printf(" Se cunple anorm\n");
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm-1])+anorm) == anorm) break;
      }
      if(DEBUG) printf(" Donde estare flag %d\n",flag);
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++ ) {
	  f=s*rv1[i-1];
	  rv1[i-1]=c*rv1[i-1];
	  if ((double)(fabs(f)+anorm) ==anorm) break;
	  g=w[i-1];
	  h=pythag(f,g);
	  if(DEBUG) printf(" f %f g %f pyth %f\n",f,g,h);
	  w[i-1]=h;
	  h=1.0/h;
	  c=g*h;
	  s= -f*h;
	  for(j=1;j<=m;j++ ) {
	    y=a[j-1][nm-1];
	    z=a[j-1][i-1];
	    a[j-1][nm-1]=y*c+z*s;
	    a[j-1][i-1]=z*c-y*s;
	  }
	}
      }
      if(DEBUG) printf(" Estoy intermedio w[k-1] %f\n",w[k-1]);
      z=w[k-1];
      if(l==k) {
	if(z<0.0) {
	  w[k-1]=-z;
	  for(j=1;j<=n;j++) v[j-1][k-1] = -v[j-1][k-1];
	}
	if(DEBUG) printf(" Se da esta condicion z %f\n",z);
	break;
      }
      
      if(DEBUG) printf(" ANtes its\n");
      if(its==30) printf(" Error: no convergence  in 30 svdcmp_d iterations\n");
      if(DEBUG) printf(" l %d\n",l); 
      x=w[l-1];
      nm=k-1;
      if(DEBUG) printf(" nm %d \n",nm);
      y=w[nm-1];
      g=rv1[nm-1];
      h=rv1[k-1];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      if(DEBUG) printf(" ANtes for \n");
      for(j=l;j<=nm;j++) {
	if(DEBUG) printf(" Entro for j %d  (%d)\n",j,nm);
	i=j+1;
	g=rv1[i-1];
	y=w[i-1];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j-1]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g= g*c-x*s;
	h=y*s;
	y *= c;
	if(DEBUG) printf(" Perim for j %d i %d\n",j,i);
	for (jj=1;jj<=n;jj++) {
	  x=v[jj-1][j-1];
	  z=v[jj-1][i-1];
	  v[jj-1][j-1]=x*c+z*s;
	  v[jj-1][i-1]=z*c-x*s;
	}
	if(DEBUG) printf(" Desou forhpy\n");
	z=pythag(f,h);
	w[j-1]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	if(DEBUG) printf(" otro for j %d i %d\n",j,i);
	for( jj=1;jj<=m;jj++) {
	  if(DEBUG) printf(" jj %d n %d\n",jj,n);
	  y=a[jj-1][j-1];
	  z=a[jj-1][i-1];
	  a[jj-1][j-1]=y*c+z*s;
	  a[jj-1][i-1]=z*c-y*s;
	}
	if(DEBUG) printf(" Acabo for\n");
      }
      rv1[l-1]=0.0;
      rv1[k-1]=f;
      if(DEBUG) printf(" LA W %d   %f\n",k-1,x);
      w[k-1]=x;
    }
  }
  free(rv1);
}
			      
  
double pythag_d(double a,double b) 
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa> absb) return(absa*sqrt(1.0+SQR(absb/absa)));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
 



void svbksb_d(double **u,double w[],double **v,int m,int n,double b[],double x[]) {

  int jj,j,i;
  double s,*tmp;

  tmp=vector_d(n);
  
  for(j=0;j<n;j++) {
    s=0.;
    if(w[j]) {
      for (i=0;i<m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
    if(DEBUG) printf(" %d tmp %f\n",j,tmp[j]);
  }
  for(j=0;j<n;j++) {
    s=0.0;
    for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
    if(DEBUG) printf(" s %f tmp %f\n",s,tmp[jj]);
    x[j]=s;
  }
  free(tmp);
}

