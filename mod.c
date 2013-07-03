#include "stdio.h"
#include "common.h"
#include "param.h"

int main(){
  
  configure();

  FILE* f;
  int i,j;
  double model[2000];
  double phi=0.4;
  double theta=0.2;

  mod(model,theta,phi);
  
  f=fopen("strip.txt","w+");

  for(j=0;j<256;j++){
    for(i=0;i<sp;i++){
      fprintf(f,"%f ",model[i*256+j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);

  system("./plot.py");

  return 1;
}
