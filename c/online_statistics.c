/* 
 *
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#define NDIM 3
static void online_statistics(unsigned n, double mean[NDIM], double mom[NDIM * NDIM], double val[NDIM]) {
  double omean[NDIM],
         omom[NDIM * NDIM];
  unsigned ix,
           jx;

  if(n == 0){
    return;
  }else if(n == 1){
    for(ix = 0; ix < NDIM; ix++){
      mean[ix] = val[ix];
      for(jx = 0; jx <= ix; jx++){
        mom[ix * NDIM + jx] = 0.;
      }
    }
  }else{
    for(ix = 0; ix < NDIM; ix++) {
      omean[ix] = mean[ix];
      mean[ix] += (val[ix] - mean[ix]) / n;
    }
    for(ix = 0; ix < NDIM; ix++) {
      for(jx = 0; jx <= ix; jx++) {
        mom[ix * NDIM + jx] += (val[ix] - omean[ix]) * (val[jx] - mean[jx]);
      }
    }
    for(ix = 0; ix < NDIM; ix++) {
      for(jx = NDIM - 1; jx > ix; jx--) {
        mom[ix * NDIM + jx] = mom[jx * NDIM + ix];
      }
    }
  }
}
static void gen(double x[NDIM], double cov[NDIM * NDIM]){
  unsigned ix;
  double data[15],
         *xme = data,
         *xmo = data + NDIM,
         *xva = data + NDIM + NDIM * NDIM;
#define BIGN 1000.00
#define URAND (((double)rand() / (double)RAND_MAX - 0.5) * 2.0)
  for(ix = 0; ix < BIGN; ix++){
    xva[0] = 45 + 10. * URAND;
    xva[1] = xva[0];
    xva[2] = -45 + URAND;
    online_statistics(ix + 1, xme, xmo, xva);
  }
  for(ix = 0; ix < NDIM * NDIM; ix++){
    cov[ix] = sqrt(xmo[ix] / ((double)BIGN - 1.));
  }
  for(ix = 0; ix < NDIM; ix++){
    x[ix] = xme[ix];
  }
}
int main(void){
  unsigned ix;
  double x[NDIM], cov[NDIM * NDIM];
  gen(x, cov);
  printf("Mean\n");
  for(ix = 0; ix < NDIM; ix++)
    printf("%10.3e%s", x[ix], (ix + 1) % NDIM == 0 ? "\n" : "\t");
  printf("Covariance\n");
  for(ix = 0; ix < NDIM * NDIM; ix++){
    printf("%10.3e%s", cov[ix], (ix + 1) % NDIM == 0 ? "\n" : "\t");
  }
  return 0;
}
