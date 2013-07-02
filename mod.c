#include "stdio.h"
#include "common.h"

int main(){
  
  FILE* f;
  int sp = n();
  int i,j;
  double model[2000];
  double phi=0.4;
  double theta=0.2;
  double offset=0.2;
  double L=20.;
  double d=1.;

  mod(model,theta,phi,offset,L,d);
  
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
