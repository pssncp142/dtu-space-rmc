#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "common.h"

#define PI 3.14159f

int nofstrip = 2;
int opening = 0.5;

int n(){
  return nofstrip;
}


double **real(double theta, double phi, double offset, double L, double d, double nofphot, double noise) 
{
  int i,j,k;
  double nofbins = 256; 
  double probint = 0.05;
  theta = theta*PI; phi = phi*PI; 

  nofphot *= 0.5;
  noise *= nofphot;
  double PIL_over_d = PI*L/d;
  double var = nofphot/nofbins;
  double var_n = noise/nofbins;
  int dim = 1/probint;
  double* prob = (double*) malloc(100*sizeof(double));
  double* prob_n = (double*) malloc(100*sizeof(double));
  int* randphot = (int*) malloc(nofbins*sizeof(int));
  int* randphot_n = (int*) malloc(nofbins*sizeof(int));
  double** count = (double**) calloc(2,sizeof(double*));
  for (int i=0; i<nofstrip; i++){count[i] = (double*) calloc(nofbins,sizeof(double));}
  double* rate = (double*) malloc(sizeof(double));
  double rnd;

  double check; 
  prob[0] = exp(-var); prob_n[0] = exp(-var_n);
  for(i=0; i<100; i++){
    prob[i+1] = prob[i] + exp(-var)*pow(var,i+1)/factorial(i+1);
    prob_n[i+1] = prob_n[i] + exp(-var_n)*pow(var_n,i+1)/factorial(i+1);
  }
  for(i=0; i<nofbins; i++){
    rnd = (double)rand()/RAND_MAX;
    for (int j=0; j<200; j++){
      check=prob[j];
      if(rnd < check){
	randphot[i] = j;  
	break;
      }
    }
    rnd = (double)rand()/RAND_MAX;
    for (j=0; j<200; j++){
      check=prob_n[j];
      if(rnd < check){
	randphot_n[i] = j;  
	break;
      }
    }
  }
  for(i=0; i<nofbins; i++){
    rate[0] = sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+offset*PI,PI);
    for(j=1; j<nofstrip; j++){
      rate[j] = rate[j-1] + sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+j)*PI,PI);
    }
    for(j=0; j<200; j++){
      if(j==randphot[i]){
	break;
      }
      rnd = (double)rand()/RAND_MAX*opening*nofstrip;
      for(k=0; k<nofstrip; k++){
	if(rnd < rate[i]) {
	  ++count[k][i]; break;}
	if(k==nofstrip){
	  ++count[nosftrip][i];}
      }
    }
    for(j=0; j<200; j++){
      if(j==randphot_n[i]){
	break;
      }
      rnd = (double)rand()/RAND_MAX*nofstrip;
      ++count[(int)rnd][i];
    }
  }
  return count;
}

double **mod(double theta, double phi, double offset, double L, double d) 
{
  int i,j;
  theta = theta*PI; phi = phi*PI;

  double **c = (double**)malloc(nofstrip*sizeof(double*));
  for(i=0; i<2; i++){c[i]=(double*)malloc(256*sizeof(double));} 
  for(i=0; i<256; i++){
    for(j=0; j<nofstrip; j++){
      c[j][i] = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/256-phi)+(offset+j)*PI,PI);
    }
  }
  return c;
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
