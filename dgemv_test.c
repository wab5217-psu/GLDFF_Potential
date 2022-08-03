#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>


int main(){

  const MKL_INT n=100;
  const MKL_INT m=50;


  const double alpha=1;
  const double beta=0;
  int sixfour=64;
  int step=1;
  
  double *a;
  a=(double*)mkl_calloc(n*m,sizeof(double),sixfour);
  
  double *b;
  b=(double*)mkl_calloc(n,sizeof(double),sixfour);


  double *res;
  res=(double*)mkl_calloc(m,sizeof(double),sixfour);

  int jr,jc;
  for(jr=0; jr<n; jr++){
    a[jr+n*jr]=1;
    b[jr]=(double)jr;
  }
  
  cblas_dgemv(CblasRowMajor,CblasTrans,m,n,alpha,a,n,b,step,beta,res,step);

  for( jr=0; jr<m; jr++ )fprintf(stderr,"%d %lf\n",jr,res[jr]);

  mkl_free(a);
  mkl_free(b);
  mkl_free(res);
}


      
