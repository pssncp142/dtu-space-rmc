#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "common.h"

#define PI 3.14159f

int nofstrip = 6;
double opening = 1./3;

int n(){
  return nofstrip;
}
double op(){
  return opening;
}

double sawtooth(double x, double period)
{
  uint check;
  if(x/(period)<0) {
    check = floor(x/period+600);
  } else {
    check = floor(x/period);
  }
  if (check%6 == 0){
    return -(x-check*period)/period+floor((x-check*period)/period)+1;
  } else if ((check%6 == 4)){
    return -(x-(check+4)*period)/period+floor((x-(check+4)*period)/period)+1;
  } else if ((check%6 == 3)){
    return (x-(check+3)*period)/period-floor((x-(check+3)*period)/period);
  } else if ((check%6 == 5)){
    return (x-(check+5)*period)/period-floor((x-(check+5)*period)/period);
  } else {
    return 0;
  }    
}
