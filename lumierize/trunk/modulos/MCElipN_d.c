#include "modulos.h"

#define DEBUG  0
#define DEBUG2 0
#define DEBUG3 0

/* Sirve para ajustar una hyperconica en N dimensiones */
/* La conica debe estar centrada en 0 !!*/
/* La conica se define como:
   x_i x_j C_ij=1
   donde x=(x_0, x_1,   x_N) 
   y C_ij es la matriz, simetrica que define la conica */

/* x[][] es una matriz, en la primera dimension contiene 
   la componente y en la segunda el numero de datos
   x[3][45] sera la componente 3 del punto numero 46 */
/* C[][] es una matriz dimensionada como NxN */

/* Las ecuaciones que permiten resolver el problema por minimos cuadrados son: 
   A_ij_lm C_ij  = b_lm
   donde 
   A_ij_lm= sum_n (x_l x_m sum_i sum_j x_i x_j)
   b_lm   = sum_n x_l x_m  */

/* Este sistema lo resuelvo con gaussj, pero tengo que meter c_ij 
   en una nueva variable lineal cc_i. Lo mismo con a_ij_lm y b_lm 
   que las meto en aa, bb respectivamente */

/* El numero de incognitas en total es N(N+1)/2 */


int  MCElipN_d(int npt, int ndim, double **x, double **C)
{
  int i,j;
  int l,m;
  int n;
  double **aa;
  double **bb;
  double *cc;
  int nvec;

  nvec=ndim*(ndim+1)/2;

  if(DEBUG) printf(" nvec vale %d\n",nvec);

  if(npt < nvec) {
    /*     fprintf(stderr,"MCElip: ERROR. npt < ndim*(ndim+1)/2 (npt=%d, ndim=%d)\n",npt,ndim); */
    return(0);
  }


  
/*  Dimensionado de las matrices */
/*  ---------------------------- */
  aa = matrix_d(nvec,nvec);
  bb = matrix_d(nvec,1);
  cc = vector_d(nvec);
  
  /*  Calculo de los coeficientes del sistema */
  /*  --------------------------------------- */
  
  /* El indice _ij lo transformo en (i+1)i/2+j y el _lm en (l+1)l/2+m 
     El indice _ij_lm lo transformo en i+j+(l+m)*ndim */

  if(DEBUG) {
    printf(" Y ndim vale %d\n",ndim);
    printf(" Por otro npt %d\n",npt);
  }

  for(l=0;l<ndim;l++) {
    for(m=0;m<=l;m++) {
      bb[(l+1)*l/2+m][0]=0.;
      for(n=0;n<npt;n++) {
	bb[(l+1)*l/2+m][0]+=x[l][n]*x[m][n];
	if(DEBUG2) printf("(l+1)l/2+m %d n %d x1 %f x2 %f b %f\n",(l+1)*l/2+m,n,x[l][n],x[m][n],bb[l+m][0]);
      }
      for(i=0;i<ndim;i++) {
	for(j=0;j<=i;j++) {
	  aa[(i+1)*i/2+j][((l+1)*l/2+m)]=0.;
	  for(n=0;n<npt;n++) aa[(i+1)*i/2+j][((l+1)*l/2+m)]+=x[l][n]*x[m][n]*x[i][n]*x[j][n];
	  if(DEBUG3)printf("(i %d j %d l %d m %d idx1 %d idx2 %d aa %f\n",i,j,l,m,(i+1)*i/2+j,(l+1)*l/2+m,aa[(i+1)*i/2+j][((l+1)*l/2+m)]);
	}
      }
    }
  }

  if(DEBUG) {
    for(i=0; i<nvec; i++) {
      printf(" MCELIPN ");
      for(j=0; j<nvec; j++)  printf(" %g + ",aa[i][j]);
      printf(" = %f \n",bb[i][0]);
    }
  }

  if(!gaussj_d(aa , nvec , bb , 1)) return(0) ;
  
  if(DEBUG2) {
    for(i=0; i<nvec ; i++)  printf(" GAUSSJ %d bb %f\n",i,bb[i][0]);
  }
  
  for(i=0;i<ndim;i++) {
    for(j=0;j<=i;j++) {
      C[i][j]=bb[(i+1)*i/2+j][0]/2.;
      C[j][i]=bb[(i+1)*i/2+j][0]/2.;
      if(j==i) C[i][i]*=2;        /* Porque todos los he contado dos veces menos los diagonales */ 
    }
  }

  free_matrix_d(aa,nvec,nvec);
  free_matrix_d(bb,nvec,1);
  free(cc);
  return(1);
}
