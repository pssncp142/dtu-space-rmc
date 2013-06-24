#include "stdio.h"
#include "stdlib.h"
#include "common.h"

int main(){
  int i,j,st;
  int sp = n();
  double data[500000];
  double theta[5] = {0.3,0.2,0.3,0,0};
  double phi[5] = {1.,0.3,0.7,0,0};
  double nofphot[5] = {1000.,1000.,1000,0,0};

  st = corr(data,3,theta,phi,0.25,50,5,nofphot,0.,50);
  FILE* f = fopen("T.txt","w+");
  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      fprintf(f,"%f ",data[sp*256*256+j*256+i]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  return 0;
}
