#include "stdio.h"
#include "common.h"

int i,j,k;

double fmaxd(double **a){
  double max=a[0][0];
  for(i=0; i<256; i++ ){
    for(j=0; j<256; j++){
      if (a[i][j] > max)
	max = a[i][j];
    }
  }
  return max;
}

double fmind(double **a){
  double min=a[0][0];
  for(i=0; i<256; i++ ){
    for(j=0; j<256; j++){
      if (a[i][j] < min)
	min = a[i][j];
    }
  }
  return min;
}

void wspfile(double **a,int n){
  FILE* f = fopen("strip.txt","w+");
  for(i=0;i<256;i++){
    for(j=0;j<n;j++){
      fprintf(f,"%f ",a[j][i]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
}

double factorial(int k)
{
  double res = 1;
  for(i=1; i<k; i++){
    res = res*(i+1); 
  }
  return res;
}

int norm1d(double *a){
  double sum = 0;
  for(i=0; i<256; i++){
    sum += a[i];
  }
  for(i=0; i<256; i++){
    a[i] /= sum;
  }
  return 1;
}

int norm2d(double **a){
  double max; double min;
  max = fmaxd(a); min = fmind(a);
  for(i=0; i<256; i++){
    for(j=0; j<256; j++){
      a[i][j] = (a[i][j]-min)/(max-min)*2-1;
    }
  }
  return 1;
}

int totsc_fact(double ***a,int n){
  for(i=0; i<256; i++){
    for(j=0; j<256; j++){
      a[n][i][j] = 0;
    }
  }
  for(i=0; i<256; i++){
    for(j=0; j<256; j++){
      for(k=0; k<n; k++){
	a[n][i][j] += a[k][i][j];
      }
    }
  }
  int st = norm2d(a[n]);
  return 1;
}

int subtmean2d(double **a, int n){
  double sum[n];
  for(i=0; i<n; i++) {sum[i]=0;}
  for(i=0; i<n; i++){
    for(j=0; j<256; j++){
      sum[i] += a[i][j]/256;
    }
  }
  for(i=0; i<n; i++){
    for(j=0; j<256; j++){
      a[i][j] -= sum[i];
    }
  }
  return 1;
}

double fsumint(int* a,int n){
  double sum=0;
  for(i=0;i<n;i++){
    sum += a[i];
  }
  return sum;
}

double fsumdub(double* a,int n){
  double sum=0;
  for(i=0;i<n;i++){
    sum += a[i];
  }
  return sum;
}

double mult3sum(double *a, double *b,double *c,int n){
  double sum = 0;
  for(i=0; i<n; i++){
    sum += a[i]*b[i]*c[i];
  }
  return sum;
}


