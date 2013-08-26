#include "stdio.h"
#include "common.h"
#include "param.h"
#include "math.h"

#define PI 3.14159f

double mean(double*);

int main(){
  
  configure(1);
  offset = 0.0;
  L = 20;
  thick = 0.;
  alpha = -0.*PI;
  beta = 1e-20*PI;
  //sp = 1;
  configure(0);

  FILE* f;
  int i,j,k;
  double model[2000];
  double map[70000];

  //double phi=0.5;
  //double theta=0.2;
  //mod(model,theta,phi,0);
  double theta[] = {0.1,0.1,0.1,0.15,0.05,0.25};
  double phi[] = {1.2,0.6,1.5,1.7,0.2,0.8};
  double nofphot[] = {70,60,70,80,40,70};


  for(k=0;k<50;k++){

  real(model,1,theta,phi,nofphot,1000,100);  
  f=fopen("strip.txt","w+");

    for(j=0;j<256;j++){
      for(i=0;i<sp;i++){
	fprintf(f,"%f ",model[i*256+j]);
      }
      fprintf(f,"\n");
    }
    fclose(f);
    
    corr(map,model,k,0,"asd");
    
    f=fopen("corr.txt","w+");
    
    for(j=0;j<256;j++){
      for(i=0;i<256;i++){
	fprintf(f,"%f ",map[i*256+j]);
      }
      fprintf(f,"\n");
    }
    fclose(f);
    
    printf("%f\n",mean(map));

  }
    //system("./plot.py");

  return 1;
}

double mean(double *map){
  double tot = 0;
  double rms;
  int i,j,k,l;
  int range=5;
  
  for(i=range;i<256-range;i++){
    for(j=range;j<256-range;j++){
       rms = 0;
       for(k=i-range;k<i+range;k++){
	 for(l=j-range;l<j+range;l++){
	   rms += pow(map[k*256+l]-map[i*256+j],2);
	 }
       }
       tot += rms/(pow(range*2+1,2)-1);
     }
  }
  return tot/(256*256);
}
