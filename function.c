#include "stdio.h"
#include "common.h"

double factorial(int k)
{
  double res = 1;
  for(int i=1; i<k; i++){
    res = res*(i+1); 
  }
  return res;
}

int norm(double *a){
  double sum = 0;
  for(int i=0; i<256; i++){
    sum += a[i];
  }
  for(int i=0; i<256; i++){
    a[i] /= sum;
  }
  return 1;
}

int norm(double **a){
  double max = fmax(a);
  double min = fmin(a);
  for(int i=0; i<256; i++){
    for(int j=0; j<256; j++){
      a[i][j] = (a[i][j]-min)/(max-min)*2-1;
    }
  }
  return 1;
}

int totsc_fact(double ***a,int n){
  for(int i=0; i<256; i++){
    for(int j=0; j<256; j++){
      a[n][i][j] = 0;
    }
  }
  for(int i=0; i<256; i++){
    for(int j=0; j<256; j++){
      for(int k=0; k<n; k++){
	a[n][i][j] += a[k][i][j];
      }
    }
  }
  int st = norm(a[n]);
  return 1;
}

int subtmean(double **a, int n){
  double sum[n];
  for(int i=0; i<n; i++) {sum[i]=0;}
  for(int i=0; i<n; i++){
    for(int j=0; j<256; j++){
      sum[i] += a[i][j]/256;
    }
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<256; j++){
      a[i][j] -= sum[i];
    }
  }
  return 1;
}

int subtmean(double ***a){
  double sum =0;
  for(int k=0; k<256; k++){
    sum=0;
    for(int i=0;i<256; i++){
      for(int j=0;j<256; j++){
	sum += a[i][j][k]/(256*256);
      }
    }
    for(int i=0;i<256; i++){
      for(int j=0;j<256; j++){
	a[i][j][k] -=sum;
      }
    }
  }
  return 1;
}

double mult3sum(double *a, double *b,double *c,int n){
  double sum = 0;
  for(int i=0; i<n; i++){
    sum += a[i]*b[i]*c[i];
  }
  return sum;
}

double fmax(double **a){
  double max=a[0][0];
  for(int i=0; i<256; i++ ){
    for(int j=0; j<256; j++){
      if (a[i][j] > max)
	max = a[i][j];
    }
  }
  return max;
}

double fmin(double **a){
  double min=a[0][0];
  for(int i=0; i<256; i++ ){
    for(int j=0; j<256; j++){
      if (a[i][j] < min)
	min = a[i][j];
    }
  }
  return min;
}

