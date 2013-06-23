#include "stdio.h"
#include "stdlib.h"

#include "common.h"

int main(){

  int i,j,st;
  int sp = n();
  double data[2000];
  st = real(data,0.3,1.2,0.5,50,5,1000.,10.,10);
  FILE* f = fopen("strip.txt","w+");
  for(i=0;i<256;i++){
    for(j=0;j<sp;j++){
      fprintf(f,"%f ",data[j*256+i]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  system("./plot.py");

}
