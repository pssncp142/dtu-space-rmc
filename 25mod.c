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

  double theta = 0.2*PI, phi = 0.*PI, nofbin = 256, offset = 0.;
  double L = 50, d = 5;

  double *cA = (double*)malloc(nofbin*sizeof(double));
  double *cB = (double*)malloc(nofbin*sizeof(double));
  double *cC = (double*)malloc(nofbin*sizeof(double));
  double *cD = (double*)malloc(nofbin*sizeof(double));
  
  FILE* file = fopen("25mod.txt","w+");
  for(int i=0; i<nofbin; i++){
    *(cA+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+offset,PI);
    *(cB+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+1)*PI,PI);
    *(cC+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset-1)*PI,PI);
    *(cD+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset-1)*PI,PI);
    fprintf(file,"%f %f %f %f\n",*(cA+i),*(cB+i),*(cC+i),*(cD+i));
  }

  return 0;
}

double sawtooth(double x, double period)
{
  uint check;
  if(floor(x/(period))<=0) {
    check = floor(x/period-1);
  } else {
    check = floor(x/period);
  }
  printf("%f %f %d \n",x/period,floor(x/period),check%4);
  if (check%3 == 0){
    return -(x-check*period)/period+floor((x-check*period)/period)+1;
  } else if (check%3 == 1 || check%4 ==1){
    return 0;
  } else {
    return (x-(check+3)*period)/period-floor((x-(check+3)*period)/period);
  }    
}
