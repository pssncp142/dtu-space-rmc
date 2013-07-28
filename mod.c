#include "stdio.h"
#include "common.h"
#include "param.h"

#define PI 3.14159f

int main(){
  
  configure(1);
  offset = 0.;
  thick = 0.;
  alpha = -0.75*PI;
  beta = 0.1*PI;
//alpha = 0.5*PI;
  configure(0);

  FILE* f;
  int i,j;
  double model[2000];
  double map[70000];

  double phi=0.4;
  double theta=0.2;
  mod(model,theta,phi);
  /*double theta[] = {0.2,0.1};
  double phi[] = {0.8,0.6};
  double nofphot[] = {10000,1000};
  real(model,1,theta,phi,nofphot,0,100);*/
  
  f=fopen("strip.txt","w+");

  for(j=0;j<256;j++){
    for(i=0;i<sp;i++){
      fprintf(f,"%f ",model[i*256+j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);

  /*corr(map,model,0,"asd");

  f=fopen("corr.txt","w+");
  
  for(j=0;j<256;j++){
    for(i=0;i<256;i++){
      fprintf(f,"%f ",map[i*256+j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);*/

  system("./plot.py");

  return 1;
}
