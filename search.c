#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "common.h"

#define PI 3.14f

int i,j,k,p,q,r,st,rem;

int main(){

  FILE* f = fopen("strip.txt","w+");
  int sp = n();
  int range=15;
  double sum,asum;
  double max=-1e8;
  double min=1e8;
  double obs[2000];
  double dev[2000]={0};
  double theta[]={0.2,0.3};
  double phi[]={0.2,1.2};
  double nofphot[]={1000.,1000.};
  double offset=0.1;
  double L=20;
  double d=1;
  double n_source=1;

  real(obs,1,theta,phi,offset,L,d,nofphot,10000.,100);
  //mod(obs,theta[0],phi[0],offset,L,d);

  sum = 0;
  for(i=0;i<2000;i++){dev[i]=0;}
  for(i=0;i<sp;i++){
    sum=0;
    for(j=0;j<256;j++){
      sum += obs[i*256+j];
    }
    for(j=0;j<256;j++){
      for(k=-range+j;k<=range+j;k++){
	rem = (k+256)%256;
	if(k>j){
	  dev[i*256+j] += pow(obs[i*256+j]-obs[i*256+rem],2);
	} else if(k<j) {
	  dev[i*256+j] += pow(obs[i*256+j]-obs[i*256+rem],2);
	}
      }
      if(dev[i*256+j]>max) max = dev[i*256+j];
      if(dev[i*256+j]<min) min = dev[i*256+j];
      }
    for(j=0;j<256;j++){
      dev[i*256+j] = (dev[i*256+j]-min)/(max-min);
      dev[sp*256+j] += dev[i*256+j]; 
    }
  }
  max = -1e8;
  min = 1e8;
  for(j=0;j<256;j++){
    if(dev[sp*256+j]>max) max = dev[sp*256+j];
    if(dev[sp*256+j]<min) min = dev[sp*256+j];
  }
  for(j=0;j<256;j++){
    dev[sp*256+j] = (dev[sp*256+j]-min)/(max-min)*2-1;
  }
  for(i=0;i<128;i++){
    sum=0;
    for(j=0;j<256;j++){
      sum += pow(dev[sp*256+j]+sin((j-i)*PI/64),2);
    }
    printf("%d %f\n",i,sum);
  }


  for(j=0;j<256;j++){
    for(i=0;i<sp+1;i++){
      fprintf(f,"%f ",dev[i*256+j]);
    }
    fprintf(f,"\n");
  }
  
  fclose(f);
  system("./plot.py");

  return 1;
}
