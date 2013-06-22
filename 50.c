#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "common.h"

#define PI 3.14159f

int nofstrip = 2;
double opening = 0.5;

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
    check = floor(x/period+200);
  } else {
    check = floor(x/period);
  }
  if (check%2 == 0){
    return -(x-check*period)/period+floor((x-check*period)/period)+1;
  }  else {
    return (x-(check+1)*period)/period-floor((x-(check+1)*period)/period);
  }    
}
