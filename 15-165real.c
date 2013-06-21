/************************************************************
 * Yigit Dallilar 11.06.2013                                *      
 * DTU-Space : Realistic modulation pattern using analytic  *
 * function                                                 *
 ************************************************************/

#include "math.h"
#include "stdlib.h"
#include "stdio.h"

#define PI 3.14159f

double sawtooth(double,double);
double factorial(int);

int main() 
{
  //initial definition
  double nofphot,nofbins,probint,offset,theta,phi,L,d,noise;

  //initial parameters
  nofphot = 1000.; nofbins = 256; probint = 0.05; L = 50; d = 5;  offset = 0.5;
  theta = 0.3*PI; phi = 1.2*PI; noise = 10.;

  //other calculations and variable definitions
  nofphot /= 3;
  noise *= nofphot;
  double PIL_over_d = PI*L/d;
  double var = nofphot/nofbins;
  double var_n = noise/nofbins;
  int dim = 1/probint;
  double* prob = (double*) malloc(100*sizeof(double));
  double* prob_n = (double*) malloc(100*sizeof(double));
  int* randphot = (int*) malloc(nofbins*sizeof(int));
  int* randphot_n = (int*) malloc(nofbins*sizeof(int));
  int** count = (int**) calloc(6,sizeof(int*));
  for (int i=0; i<6; i++){count[i] = (int*) calloc(nofbins,sizeof(int));}
  double* rate = (double*) malloc(5*sizeof(double));
  double rnd;

  double check; 
  prob[0] = exp(-var); prob_n[0] = exp(-var_n);
  for(int i=0; i<100; i++){
    prob[i+1] = prob[i] + exp(-var)*pow(var,i+1)/factorial(i+1);
    prob_n[i+1] = prob_n[i] + exp(-var_n)*pow(var_n,i+1)/factorial(i+1);
  }
  for(int i=0; i<nofbins; i++){
    rnd = (double)rand()/RAND_MAX;
    for (int j=0; j<100; j++){
      check=prob[j];
      if(rnd < check){
	randphot[i] = j;  
	break;
      }
    }
    rnd = (double)rand()/RAND_MAX;
    for (int j=0; j<100; j++){
      check=prob_n[j];
      if(rnd < check){
	randphot_n[i] = j;  
	break;
      }
    }
  }
  for(int i=0; i<nofbins; i++){
    rate[0] = sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+offset*PI,PI);
    rate[1] = rate[0] + sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+1)*PI,PI);
    rate[2] = rate[1] + sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+2)*PI,PI);
    rate[3] = rate[2] + sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+3)*PI,PI);
    rate[4] = rate[3] + sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+4)*PI,PI);
   
    for(int j=0; j<200; j++){
      if(j==randphot[i]){
	break;
      }
      rnd = (double)rand()/RAND_MAX*2;
      if(rnd < rate[0]){
	count[0][i] += 1;
      } else if (rnd < rate[1]){
	count[1][i] += 1;
      } else if (rnd < rate[2]){
	count[2][i] += 1;
      } else if (rnd < rate[3]){
	count[3][i] += 1;
      } else if (rnd < rate[4]){
	count[4][i] += 1;
      } else {
	count[5][i] += 1;
      }
    }
    for(int j=0; j<200; j++){
      if(j==randphot_n[i]){
	break;
      }
      rnd = (double)rand()/RAND_MAX*6;
      count[(int)rnd][i] += 1;
    }
  }
  
  FILE* file;
  file = fopen("15-165real.txt","w+");
  for(int i=0; i<nofbins; i++){
    fprintf(file,"%d %d %d %d %d %d\n",count[0][i],count[1][i],count[2][i],count[3][i],count[4][i],count[5][i]);
  }  
  fclose(file);

  system("./15-165real.py");

  return 0;
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

double factorial(int k)
{
  double res = 1;
  for(int i=1; i<k; i++){
    res = res*(i+1); 
  }
  return res;
}
