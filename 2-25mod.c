/************************************************************
 * Yigit Dallilar 13.06.2013                                *      
 * DTU-Space : 25 perc open grid theoretical modulation     *
 ************************************************************/

#include "math.h"
#include "stdlib.h"
#include "stdio.h"

#define PI 3.14159f

double sawtooth(double,double);

int main() 
{

  double theta = 0.3*PI, phi = 1.2*PI, nofbin = 256, offset = 0.;
  double L = 50, d = 5;

  double *cA = (double*)malloc(nofbin*sizeof(double));
  double *cB = (double*)malloc(nofbin*sizeof(double));
  double *cC = (double*)malloc(nofbin*sizeof(double));
  double *cD = (double*)malloc(nofbin*sizeof(double));
  
  FILE* file = fopen("2-25mod.txt","w+");
  for(int i=0; i<nofbin; i++){
    *(cA+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+offset*PI,PI);
    *(cB+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+1)*PI,PI);
    *(cC+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+2)*PI,PI);
    *(cD+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+3)*PI,PI);
    fprintf(file,"%f %f %f %f\n",*(cA+i),*(cB+i),*(cC+i),*(cD+i));
  }
  fclose(file);

  system("./2-25mod.py");

  return 0;
}

double sawtooth(double x, double period)
{
  uint check;
  if(x/(period)<0) {
    check = floor(x/period);
  } else {
    check = floor(x/period);
  }
  if (check%4 == 0){
    return -(x-check*period)/period+floor((x-check*period)/period)+1;
  } else if (check%4 == 1){
    return 0;
  } else if (check%4 == 2){
    return (x-(check+3)*period)/period-floor((x-(check+3)*period)/period);
  } else {
    return 1;
  }   
}
