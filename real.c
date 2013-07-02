#include "stdio.h"
#include "common.h"

int main(){
  
  FILE* f;
  int sp = n();
  int i,j;
  double obs[2000];
  double phi[]={0.,0.,0.};
  double theta[]={0.3,0.2,0.1};
  double nofphot[]={1000.,1000.,1000.};
  double noise=0000.;
  double offset=0.;
  double L=20.;
  double d=1.;
  int n_source=1;
  int turn=100;

  real(obs,n_source,theta,phi,offset,L,d,nofphot,noise,turn);
  
  f=fopen("strip.txt","w+");

  for(j=0;j<256;j++){
    for(i=0;i<sp;i++){
      fprintf(f,"%f ",obs[i*256+j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);

  system("./plot.py");

  return 1;
}
