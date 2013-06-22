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

  double theta = 0.25*PI, phi = 1.2*PI, nofbin = 256, offset = 0.5;
  double L = 50, d = 5;

  double *cA = (double*)malloc(nofbin*sizeof(double));
  double *cB = (double*)malloc(nofbin*sizeof(double));
  double *cC = (double*)malloc(nofbin*sizeof(double));
  double *cD = (double*)malloc(nofbin*sizeof(double));
  double *cE = (double*)malloc(nofbin*sizeof(double));
  double *cF = (double*)malloc(nofbin*sizeof(double));
  
  FILE* file = fopen("11-165mod.txt","w+");
  for(int i=0; i<nofbin; i++){
    *(cA+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+offset*PI,PI);
    *(cB+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+1)*PI,PI);
    *(cC+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+2)*PI,PI);
    *(cD+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+3)*PI,PI);
    *(cE+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+4)*PI,PI);
    *(cF+i) = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/nofbin-phi)+(offset+5)*PI,PI);

    fprintf(file,"%f %f %f %f %f %f\n",*(cA+i),*(cB+i),*(cC+i),*(cD+i),*(cE+i),*(cF+i));
  }
  fclose(file);

  system("./11-165mod.py");

  return 0;
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
