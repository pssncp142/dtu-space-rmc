#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "common.h"

#define PI 3.14159f

int nofstrip = 4;
double opening = 0.25;

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
    check = floor(x/period+400);
  } else {
    check = floor(x/period);
  }
  if (check%4 == 0){
    return -(x-check*period)/period+floor((x-check*period)/period)+1;
  } else if ((check%4 == 1) || (check%4 ==2)){
    return 0;
  } else {
    return (x-(check+3)*period)/period-floor((x-(check+3)*period)/period);
  }    
}
