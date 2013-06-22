#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "common.h"

#define PI 3.14159f

int nofstrip = 6;
double opening = 1/3;

int n(){
  return nofstrip;
}
double op(){
  return opening;
}

double sawtooth(double x, double period)
{
  uint check;
  if(x/(period*0.5)<0) {
    check = floor(x/(period*0.5)+600);
  } else {
    check = floor(x/(period*0.5));
  }
  if (check%12 == 0 | check%12 == 1){
    return -(x-floor(check*0.5)*period)/period+floor((x-floor(check*0.5)*period)/period)+1;
  } else if (check%12 == 3){
    return (x-(floor(check*0.5)+1)*period)/period-floor((x-(floor(check*0.5)+1)*period)/period)-0.5;
  } else if ((check%12 == 4)){
    return (x-(floor(check*0.5)+2)*period)/period-floor((x-(floor(check*0.5)+2)*period)/period)+0.5;
  } else if ((check%12 == 5)){
    return -(x-(floor(check*0.5)+2)*period)/period+floor((x-(floor(check*0.5)+2)*period)/period)+1.5;
  } else if ((check%12 == 6)){
    return -(x-(floor(check*0.5)+3)*period)/period+floor((x-(floor(check*0.5)+3)*period)/period)+0.5;
  } else if ((check%12 == 10) | (check%12 == 11)){
    return (x-(floor(check*0.5)+5)*period)/period-floor((x-(floor(check*0.5)+5)*period)/period);
  } else {
    return 0;
  } 
}
