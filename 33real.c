/************************************************************
 * Yigit Dallilar 11.06.2013                                *      
 * DTU-Space : Realistic modulation pattern using analytic  *
 * function                                                 *
 ************************************************************/

#include "math.h"
#include "stdlib.h"
#include "stdio.h"

#define PI 3.14159f

double poisson_formula(double, int);
double factorial(int);
double sawtooth(double,double);

int main() 
{
  //initial definition
  double nofphot,nofbins,probint,offset,theta,phi,L,d;

  //initial parameters
  nofphot = 10000; nofbins = 256; probint = 0.05; L = 50; d = 5;  offset = 0.;
  theta = 0.3*PI; phi = 1.2*PI;

  //other calculations and variable definitions
  double PIL_over_d = PI*L/d;
  double var = nofphot/nofbins;
  double prob[101] = {0};
  int dim = 1/probint;
  double randphot[dim];

  int cntr=0; double check;
  for(int i=0; i<100; i++){
    prob[i+1] = prob[i] + poisson_formula(var,i);
    check = (cntr+0.5)*probint;
    if(check < prob[i+1]){
      randphot[cntr] = (check-prob[i])/(prob[i+1]-prob[i]) + i;
      cntr = cntr + 1;
    }
  }
  
  int lastndx = nofbins;
  double countA[lastndx];
  double countB[lastndx];
  double countC[lastndx];
  for(int i=0; i<lastndx; i++){
    double rnd = (float)rand()/RAND_MAX;
    countA[i] = sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+offset*PI,PI)*randphot[(int) floor(rnd/probint)];
    countB[i] = sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+1)*PI,PI)*randphot[(int) floor(rnd/probint)];
    countC[i] = sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+2)*PI,PI)*randphot[(int) floor(rnd/probint)];
  }  
  
  FILE* file;
  file = fopen("33real.txt","w+");
  for(int i=0; i<nofbins; i++){
    fprintf(file,"%f %f %f\n",countA[i],countB[i],countC[i]);
  }  
  fclose(file);

  system("./33real.py");

  return 0;
}

double sawtooth(double x, double period)
{
  uint check;
  if(x/(period)<0) {
    check = floor(x/period-1);
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

double factorial(int k)
{
  double res = 1;
  for(int i=1; i<k; i++){
    res = res*(i+1); 
  }

  return res;
}

double poisson_formula(double var, int k)
{
  return pow(var,k)*exp(-var)/factorial(k);
}
