#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "common.h"

#define PI 3.14159f

int nofstrip = 5;
double opening = 0.4;

int n(){
  return nofstrip;
}
double op(){
  return opening;
}

double sawtooth(double x, double period){
  uint check;
  if(x/(period*0.5)<0) {
    check = floor(x/(period*0.5)+500);
  } else {
    check = floor(x/(period*0.5));
  }

  if (check%10 == 0 | check%10 == 1){
    return -(x-floor(check*0.5)*period)/period+floor((x-floor(check*0.5)*period)/period)+1;
  } else if (check%10 == 3){
    return (x-(floor(check*0.5)+1)*period)/period-floor((x-(floor(check*0.5)+1)*period)/period)-0.5;
  } else if (check%10 == 4){
    return (x-(floor(check*0.5)+2)*period)/period-floor((x-(floor(check*0.5)+2)*period)/period)+0.5;
  } else if (check%10 == 5){
    return -(x-(floor(check*0.5)+2)*period)/period+floor((x-(floor(check*0.5)+2)*period)/period)+1.5;
  } else if (check%10 == 6){
    return -(x-(floor(check*0.5)+3)*period)/period+floor((x-(floor(check*0.5)+3)*period)/period)+0.5;
  } else if ((check%10 == 8) | (check%10 == 9)){
    return (x-(floor(check*0.5)+4)*period)/period-floor((x-(floor(check*0.5)+4)*period)/period);
  } else {
    return 0;
  } 

}
