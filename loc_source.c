#include "stdio.h"
#include "stdlib.h"
#include "common.h"

int main(){
  int i,j,st;
  int sp = n();
  double info[10];
  double theta[5] = {0.3,0.2,0.3,0,0};
  double phi[5] = {1.,0.3,0.7,0,0};
  double nofphot[5] = {1000.,1000.,1000,0,0};

  st = loc_source(info,3,theta,phi,0.,50,5,nofphot,10000.,50);
  return 0;
}
