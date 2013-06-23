#include "stdio.h"
#include "stdlib.h"
#include "common.h"

int main(){
  int i,j,st;
  int sp = n();
  double data[500000];
  printf("A");
  st = corr(data,0.3,1.2,0.25,50,5,1000.,10.,10);
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
