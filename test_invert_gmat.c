#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>

#include <mkl.h>

extern int pinv(int m,int n,double* a,double* inv);  


int main(){

  FILE *fp;

  double hold;
  double *G;
  double *GT;
  double *GTG;
  double *INV_ar;
  double *GTG_i;

  int jr,jc,jt,indx,indxT;

  int ngrid=3637;
  int neqn=7000;

  char *vals[5];
  char *token;

  char *line=NULL;
  size_t len=0;

  int stat;
  double alpha=1;
  double beta=0;

  
  G=(double*)mkl_calloc(neqn*ngrid,sizeof(double),64);
  GT=(double*)mkl_calloc(neqn*ngrid,sizeof(double),64);
  GTG=(double*)mkl_calloc(ngrid*ngrid,sizeof(double),64);

  GTG_i=(double*)mkl_calloc(ngrid*ngrid,sizeof(double),64);
  INV_ar=(double*)mkl_calloc(ngrid*ngrid,sizeof(double),64);


  fp=fopen("hold","r");

  while( getline(&line,&len,fp)!= EOF ){
    jt=0;
    token=strtok(line," ");
    while( token != NULL){
      vals[jt]=token;
      token=strtok(NULL," ");      
      jt++;
    }

    if(jt==5){
      sscanf(vals[0],"%d",&jr);
      sscanf(vals[1],"%d",&jc);

      indx=jr*ngrid+jc;
      sscanf(vals[4],"%lf",&hold);
      G[indx]=hold;
    }
  }

  /* for(jr=0; jr<neqn; jr++)for( jc=0; jc<ngrid; jc++){ */
  /*     indx=jr*ngrid+jc; */
  /*     if(G[indx]!=0){ */
  /* 	fprintf(stdout,"%d %d %7.4e\n",jr,jc,G[indx]); */
  /*   } */
  /* } */

  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ngrid,ngrid,neqn,alpha,G,ngrid,G,ngrid,beta,GTG,ngrid);

  stat=pinv((MKL_INT)ngrid,(MKL_INT)ngrid,GTG,INV_ar);


  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ngrid,ngrid,ngrid,alpha,INV_ar,ngrid,GTG,ngrid,beta,GTG_i,ngrid);

  for(jr=0; jr<ngrid; jr++)for(jc=0; jc<ngrid; jc++){
      indx=jr*ngrid+jc;
  if( GTG_i[indx]>1.e-3 )fprintf(stdout,"GTG_i: %d %d %8.4e\n",jr,jc,GTG_i[indx]);
    }
  
  fflush(stdout);

  

}
