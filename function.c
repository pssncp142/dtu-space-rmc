#include "stdio.h"
#include "stdlib.h"
#include "common.h"

int i,j,k;

double factorial(int k)
{
  double res = 1;
  for(i=1; i<k; i++){
    res = res*(i+1); 
  }
  return res;
}

double mult3sum(double *a,double *b,double *c,int n){
  double sum = 0;
  for(i=0;i<256;i++){
    sum += a[n*256+i]*b[n*256+i]*c[n*256+i];
  }
  return sum;
}

int subtmean(double a[],int n){
  double sum;
  for(j=0; j<n; j++){
    sum = 0;
    for(i=0; i<256; i++){
      sum += a[j*256+i];
    }
    for(i=0; i<256; i++){
      a[j*256+i] -= sum/256;
    }
  }
  return 1;
}

int norm1_1(double a[],int n){
  double max,min;
  for(k=0;k<n;k++){
    max = -1e10; min = +1e10;
    for(i=0;i<256;i++){
      for(j=0;j<256;j++){
	if(a[k*256*256+j*256+i]>max) max = a[k*256*256+j*256+i];
	if(a[k*256*256+j*256+i]<min) min = a[k*256*256+j*256+i];
      }
    }
    for(i=0;i<256;i++){
      for(j=0;j<256;j++){
	a[k*256*256+j*256+i] = (a[k*256*256+j*256+i]-min)/(max-min)*2-1;
      }
    }
  }
  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      a[n*256*256+j*256+i] = 0;
      for(k=0;k<n;k++){
	a[n*256*256+j*256+i] += a[k*256*256+j*256+i];
      }
    }
  }
  max = -1e10; min = +1e10;
  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      if(a[n*256*256+j*256+i]>max) max = a[n*256*256+j*256+i];
      if(a[n*256*256+j*256+i]<min) min = a[n*256*256+j*256+i];
    }
  }
  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      a[n*256*256+j*256+i] = (a[n*256*256+j*256+i]-min)/(max-min)*2-1;
    }
  }
  return 1;
}

int norm0_1(double a[], int n){
  double sum;
  for(j=0; j<n; j++){
    sum = 0;
    for(i=0; i<256; i++){
      sum += a[j*256+i];
    }
    for(i=0; i<256; i++){
      a[j*256+i] /= sum;
    }
  }
  return 1;
}


