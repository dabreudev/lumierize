/* Programa pirateado de pgplot para dibujar imagenes en 3 dimensiones */
/* Para compilar: */
/* 	cc  -I$PGPLOT_DIR -c $s2/Proced/C/modulos/cpgdraw3d.c */
/* 	ARRAY(kx,ny) es la matriz a representar */
/* 
   You must link the resulting object file with the libraries:
   -lf2c -lm   (in that order)
   */

#include "modulos.h"
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)


int fredgo(float *array,int * mn);

struct {
  float deltax, x;
  int step, left, right, it, nx;
  int visble;
} fredcm_;

#define fredcm_1 fredcm_


static float c_b10 = (float)0.;

/* Subroutine */ 
int cpgdraw3d(float *arrayen, int kx, int ny, float size, float angle)
{
  int array_dim1, array_offset, i__1, i__2;
  float r__1, r__2;
  float *array;
  double sin();
  float peak, fmin, fmax, sine;
  int incx, i__, j, ki, kj, mn;
  float dx;
  int mx;
  float height;
  int my;
  float deltav, deltay;


  array=malloc(kx * ny*sizeof(float));
  for(j=0;j<(kx * ny);j++) array[j]=arrayen[j];

  array_dim1 = kx;
  array_offset = array_dim1 + 1;
  array -= array_offset;
  
  /* Function Body */
  mn = kx * ny;
  fredcm_1.nx = kx;
  /*     Check array size: */
  if (fredcm_1.nx < 2 || ny < 2) {
    return 0;
  }
  fmax = array[array_dim1 + 1];
  fmin = fmax;
  i__1 = ny;
  for (j = 1; j <= i__1; ++j) {
    i__2 = fredcm_1.nx;
    for (i__ = 1; i__ <= i__2; ++i__) {
      /* Computing MIN */
      r__1 = array[i__ + j * array_dim1];
      fmin = dmin(r__1,fmin);
      /* Computing MAX */
      r__1 = array[i__ + j * array_dim1];
      fmax = dmax(r__1,fmax);
      /* L10: */
    }
    /* L20: */
  }
  fredcm_1.deltax = size / (fredcm_1.nx + ny);
  sine = sin(angle / (float)58.);
  deltay = fredcm_1.deltax * sine;
  height = size * ((float)1. - dabs(sine));
  deltav = height;
  fmax -= fmin;
  if (fmax < (float)1e-4) {
    fmax = deltav;
  }
  deltav /= fmax;
  mx = fredcm_1.nx + 1;
  my = ny + 1;
  fredcm_1.step = mx;
  
  /* Start PGPLOT buffering. */
  
  cpgbbuf();
  
  /* Work our way down the Y axis, then up the X axis, */
  /* calculating the Y plotter coordinates for each */
  /* column of the plot, doing the hidden-line suppression */
  /* at the same time. */
  
  i__1 = ny;
  for (j = 1; j <= i__1; ++j) {
    kj = my - j;
    ki = 1;
    /*               ( KI,KJ are coordinates of bottom of column) */
    array[ki + kj * array_dim1] = deltay * (ki + kj) + deltav * (array[ki+ kj * array_dim1] - fmin);
  L30:
    peak = array[ki + kj * array_dim1];
  L40:
    ++ki;
    ++kj;
    if (ki > fredcm_1.nx || kj > ny) {
      goto L50;
    }
    array[ki + kj * array_dim1] = deltay * (ki + kj) + deltav * (array[ki+ kj * array_dim1] - fmin);
    if (array[ki + kj * array_dim1] > peak) {
      goto L30;
    }
    if (array[ki + kj * array_dim1] <= peak) {
      array[ki + kj * array_dim1] = -(r__1 = array[ki + kj * array_dim1]
				      , dabs(r__1));
    }
    goto L40;
  L50:
    ;
  }
  
  /* Now to work our way up the X axis */
  
  i__1 = fredcm_1.nx;
  for (i__ = 2; i__ <= i__1; ++i__) {
    ki = i__;
    kj = 1;
    array[ki + kj * array_dim1] = deltay * (ki + kj) + deltav * (array[ki 
								      + kj * array_dim1] - fmin);
  L60:
    peak = array[ki + kj * array_dim1];
  L70:
    ++ki;
    ++kj;
    if (ki > fredcm_1.nx || kj > ny) {
      goto L80;
    }
    array[ki + kj * array_dim1] = deltay * (ki + kj) + deltav * (array[ki 
								      + kj * array_dim1] - fmin);
    if (array[ki + kj * array_dim1] > peak) {
      goto L60;
    }
    if (array[ki + kj * array_dim1] <= peak) {
      array[ki + kj * array_dim1] = -(r__1 = array[ki + kj * array_dim1]
				      , dabs(r__1));
    }
    goto L70;
  L80:
    ;
  }
  
  /* Draw a line along the bottom of the vertical faces */
  
  r__1 = fredcm_1.deltax * (fredcm_1.nx + ny - 2);
  r__2 = deltay * mx;
  cpgmove(r__1, r__2);
  r__1 = fredcm_1.deltax * (ny - 1);
  r__2 = deltay * 2;
  cpgdraw(r__1, r__2);
  r__1 = deltay * my;
  cpgdraw(c_b10, r__1);
  
  /* Array is now ready for plotting.  If a point is */
  /* positive, then it is to be plotted at that Y */
  /* coordinate; if it is negative, then it is */
  /* invisible, but at minus that Y coordinate (the point */
  /* where the line heading towards it disappears has to */
  /* be determined by finding the intersection of it and */
  /* the cresting line). */
  
  /* Plot rows: */
  
  i__1 = ny;
  for (j = 1; j <= i__1; j += 2) {
    kj = my - j;
    dx = fredcm_1.deltax * (j - 2);
    fredcm_1.x = dx + fredcm_1.deltax;
    r__1 = deltay * (kj + 1);
    cpgmove(fredcm_1.x, r__1);
    cpgdraw(fredcm_1.x, array[kj * array_dim1 + 1]);
    fredcm_1.visble = 1;
    i__2 = fredcm_1.nx;
    for (i__ = 2; i__ <= i__2; ++i__) {
      fredcm_1.right = i__ + fredcm_1.nx * (kj - 1);
      fredcm_1.left = fredcm_1.right - 1;
      fredcm_1.it = fredcm_1.right;
      fredcm_1.x = dx + fredcm_1.deltax * i__;
      fredgo(&array[array_offset], &mn);
      /* L90: */
    }
    
    /* Now at far end of row so come back */
    
    --kj;
    if (kj <= 0) {
      goto L170;
    }
    fredcm_1.visble = array[fredcm_1.nx + kj * array_dim1] >= (float)0.;
    dx = fredcm_1.deltax * (fredcm_1.nx + j);
    if (fredcm_1.visble) {
      r__1 = dx - fredcm_1.deltax;
      cpgmove(r__1, array[fredcm_1.nx + kj * array_dim1]);
    }
    fredcm_1.deltax = -fredcm_1.deltax;
    i__2 = fredcm_1.nx;
    for (i__ = 2; i__ <= i__2; ++i__) {
      ki = mx - i__;
      fredcm_1.left = ki + fredcm_1.nx * (kj - 1);
      fredcm_1.right = fredcm_1.left + 1;
      fredcm_1.it = fredcm_1.left;
      fredcm_1.x = dx + fredcm_1.deltax * i__;
      fredgo(&array[array_offset], &mn);
      /* L100: */
    }
    
    fredcm_1.x = dx + fredcm_1.deltax * fredcm_1.nx;
    if (! fredcm_1.visble) {
      cpgmove(fredcm_1.x, array[kj * array_dim1 + 1]);
    }
    r__1 = deltay * (kj + 1);
    cpgdraw(fredcm_1.x, r__1);
    /*               (set DELTAX positive for return trip) */
    fredcm_1.deltax = -fredcm_1.deltax;
    /* L110: */
  }
  
  /* Now do the columns: */
  /* as we fell out of the last DO-loop we do the */
  /* columns in ascending-X order */
  
  incx = 1;
  ki = 1;
  /*               (set DELTAX -ve since scanning R to L) */
L120:
  dx = fredcm_1.deltax * (ki + ny - 1);
  fredcm_1.deltax = -fredcm_1.deltax;
  fredcm_1.x = dx + fredcm_1.deltax;
  cpgmove(fredcm_1.x, array[array_dim1 + 1]);
L130:
  fredcm_1.visble = 1;
  i__1 = ny;
  for (j = 2; j <= i__1; ++j) {
    fredcm_1.left = ki + fredcm_1.nx * (j - 1);
    fredcm_1.right = fredcm_1.left - fredcm_1.nx;
    fredcm_1.it = fredcm_1.left;
    fredcm_1.x = dx + fredcm_1.deltax * j;
    fredgo(&array[array_offset], &mn);
    /* L140: */
  }
  
  /* At far end, increment X and check still inside array */
  
  ki += incx;
  if (ki <= 0 || ki > fredcm_1.nx) {
    goto L180;
  }
  fredcm_1.visble = array[ki + ny * array_dim1] >= (float)0.;
  fredcm_1.deltax = -fredcm_1.deltax;
  dx = fredcm_1.deltax * (ki - 2);
  fredcm_1.x = dx + fredcm_1.deltax;
  if (fredcm_1.visble) {
    cpgmove(fredcm_1.x, array[ki + ny * array_dim1]);
  }
  i__1 = ny;
  for (j = 2; j <= i__1; ++j) {
    kj = my - j;
    fredcm_1.right = ki + fredcm_1.nx * (kj - 1);
    fredcm_1.left = fredcm_1.right + fredcm_1.nx;
    fredcm_1.it = fredcm_1.right;
    fredcm_1.x = dx + fredcm_1.deltax * j;
    fredgo(&array[array_offset], &mn);
    /* L150: */
  }
  
  fredcm_1.x = dx + fredcm_1.deltax * ny;
  if (! fredcm_1.visble) {
    cpgmove(fredcm_1.x, array[ki + array_dim1]);
  }
  if (ki == 1) {
    goto L180;
  }
  r__1 = deltay * (ki + 1);
  cpgdraw(fredcm_1.x, r__1);
  ki += incx;
  if (ki > fredcm_1.nx) {
    goto L180;
  }
  if (ki == 1) {
    goto L120;
  }
L160:
  fredcm_1.deltax = -fredcm_1.deltax;
  dx = fredcm_1.deltax * (1 - ki - ny);
  fredcm_1.x = dx + fredcm_1.deltax;
  r__1 = deltay * (ki + 1);
  cpgmove(fredcm_1.x, r__1);
  cpgdraw(fredcm_1.x, array[ki + array_dim1]);
  goto L130;
  
  /* Do columns backwards because ended rows at far end of X */
  
L170:
  ki = fredcm_1.nx;
  incx = -1;
  dx = fredcm_1.deltax * (ki + ny);
  goto L160;
  
  
L180:
  cpgebuf();
  return 0;
} /* pgdraw3d_ */

/* ----------------------------------------------------------------------- */
/* Subroutine */ int fredgo(float *array,int * mn)
{
  /* System generated locals */
  float r__1;
  
  /* Local variables */
  float y, al, bl, em, ar, xx;
/*   extern  Subroutine  int pgdraw_(), pgmove_(); */
  
  
  
  /* Test visibility */
  
  /* Parameter adjustments */
  --array;
  
  /* Function Body */
  if (array[fredcm_1.it] < (float)0.) {
    goto L80;
  }
  
  /* This point is visible - was last? */
  
  if (fredcm_1.visble) {
    goto L50;
  }
  
  /* No: calculate point where this line vanishes */
  
L10:
  if (fredcm_1.left <= fredcm_1.nx || (fredcm_1.left - 1) % fredcm_1.nx == 
      0 || fredcm_1.right <= fredcm_1.nx || (fredcm_1.right - 1) % 
      fredcm_1.nx == 0) {
    goto L100;
  }
  al = (r__1 = array[fredcm_1.left], dabs(r__1));
  ar = (r__1 = array[fredcm_1.right], dabs(r__1));
  if (array[fredcm_1.left] < (float)0.) {
    goto L70;
  }
  /*               Right-hand point is crested */
L20:
    fredcm_1.right -= fredcm_1.step;
    if (array[fredcm_1.right] < (float)0.) {
      goto L20;
    }
    /*               Left-hand end of cresting line is either */
    /*               RIGHT+NX or RIGHT-1 */
    fredcm_1.left = fredcm_1.right + fredcm_1.nx;
    if (array[fredcm_1.left] < (float)0.) {
      fredcm_1.left = fredcm_1.right - 1;
    }
    
    /*               RIGHT and LEFT index into the endpoints of the */
    /*               cresting line */
L30:
    bl = (r__1 = array[fredcm_1.left], dabs(r__1));
    em = (r__1 = array[fredcm_1.right], dabs(r__1)) - bl;
    xx = em - ar + al;
    if (dabs(xx) < (float)1e-4) {
      goto L60;
    }
    xx = (al - bl) / xx;
L40:
    y = em * xx + bl;
    if (fredcm_1.deltax > (float)0.) {
      xx = (float)1. - xx;
    }
    xx = fredcm_1.x - xx * fredcm_1.deltax;
    if (fredcm_1.visble) {
      goto L90;
    }
    /*               Drawing a line from an invisible point */
    /*               to a visible one */
    cpgmove(xx, y);
    fredcm_1.visble = 1;
L50:
    cpgdraw(fredcm_1.x, array[fredcm_1.it]);
    return 0;
    
L60:
    xx = (float).5;
    goto L40;
    
    /* Left-hand point crested */
    
L70:
    fredcm_1.left -= fredcm_1.step;
    if (array[fredcm_1.left] < (float)0.) {
      goto L70;
    }
    
    /* Right-hand end of cresting line is either LEFT+1 or LEFT-NX */
    
    fredcm_1.right = fredcm_1.left + 1;
    if (array[fredcm_1.right] < (float)0.) {
      fredcm_1.right = fredcm_1.left - fredcm_1.nx;
    }
    goto L30;
    
    /* This point is invisible; if last one was too, then forget it; */
    /* else draw a line towards it */
    
L80:
    if (! fredcm_1.visble) {
      return 0;
    }
    goto L10;
    
L90:
    cpgdraw(xx, y);
L100:
    fredcm_1.visble = 0;
    return 0;
} /* fredgo_ */
