/************************************************************
 * Yigit Dallilar 11.06.2013                                *      
 * DTU-Space : 50 perc open grid theoretical modulation     *
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
  
  FILE* file = fopen("50mod.txt","w+");
  for(int i=0; i<nofbin; i++){
    *(cA+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+offset*PI,PI);
    *(cB+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+1)*PI,PI);
    fprintf(file,"%f %f\n",*(cA+i),*(cB+i));
  }
  fclose(file);

  system("./50mod.py");

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
  if (check%2 == 0){
    return -(x-check*period)/period+floor((x-check*period)/period)+1;
  }  else {
    return (x-(check+2)*period)/period-floor((x-(check+2)*period)/period);
  }    
}
