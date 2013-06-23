#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "common.h"

#define PI 3.14159f

int nofstrip = 3;
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
    check = floor(x/period+300);
  } else {
    check = floor(x/period);
  }
  if (check%3 == 0){
    return -(x-check*period)/period+floor((x-check*period)/period)+1;
  } else if (check%3 == 1){
    return 0;
  } else {
    return (x-(check+2)*period)/period-floor((x-(check+2)*period)/period);
  }    
}
