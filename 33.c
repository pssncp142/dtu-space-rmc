#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "common.h"

#define PI 3.14159f

int nofstrip = 3;
double opening = 1/3;

int n(){
  return nofstrip;
}

double **lsf(double theta,double phi,double offset,double L,double d,double nofphot,double noise)
{
  int i,j,st;
  double **data=(double**)malloc((nofstrip+2)*sizeof(double*));
  for(i=0;i<nofstrip+2;i++){data[i]=(double*)malloc(2*sizeof(double));}
  data[nofstrip][0]=0; data[nofstrip][1]=0;
  double LHS[100]; double RHS[100];
  double **obs=real(theta,phi,offset,L,d,nofphot,noise); 
  double **model=mod(theta,phi,offset,L,d);
 
  for(j=0; j<nofstrip; j++){

    st = norm(model[j]);
    st = norm(obs[j]);
  
    LHS[0] = 1./256; LHS[1]=0; LHS[2]=0; LHS[3]=0; RHS[0] = 0; RHS[1] = 0; 

    for(i=0; i<256; i++){
      LHS[1] += model[j][i]/256;
      LHS[2] += model[j][i]/256;
      LHS[3] += pow(model[j][i],2);
      RHS[0] += obs[j][i]/256;
      RHS[1] += obs[j][i]*model[j][i];
    }

    st =  gaussj_nr(LHS,2,RHS,2);
    data[j][0] = RHS[0];
    data[j][1] = RHS[1];
  }

  for(i=0;i<nofstrip;i++){
    data[nofstrip][0] += data[i][0]/nofstrip;
    data[nofstrip][1] += data[i][1]/nofstrip;
  }

  LHS[0] = (1./256)/nofstrip; LHS[1]=0; LHS[2]=0; LHS[3]=0; RHS[0] = 0; RHS[1] = 0; 

  for(j=0;j<nofstrip;j++){
    for(i=0;i<256;i++){
      LHS[1] += (model[j][i]/256)/(nofstrip*nofstrip);
      LHS[2] += (model[j][i]/256)/(nofstrip*nofstrip);
      LHS[3] += pow(model[j][i],2)/(nofstrip*nofstrip);
      RHS[0] += (obs[j][i]/256)/(nofstrip*nofstrip);
      RHS[1] += (obs[j][i]*model[j][i])/(nofstrip*nofstrip);
    }
  }

  st =  gaussj_nr(LHS,2,RHS,2);
  data[nofstrip+1][0] = RHS[0];
  data[nofstrip+1][1] = RHS[1];
  
  free(obs);
  free(model);
  return data;
}

double **real(double theta, double phi, double offset, double L, double d, double nofphot, double noise) 
{

  phi = PI*phi; theta = theta*PI;
  double nofbins = 256;
  double probint = 0.05;
  nofphot *= 0.5;
  noise *= nofphot;
  double PIL_over_d = PI*L/d;
  double var = nofphot/nofbins;
  double var_n = noise/nofbins;
  int dim = 1/probint;
  double* prob = (double*) malloc(400*sizeof(double));
  double* prob_n = (double*) malloc(400*sizeof(double));
  int* randphot = (int*) malloc(nofbins*sizeof(int));
  int* randphot_n = (int*) malloc(nofbins*sizeof(int));
  double** count = (double**) calloc(nofstrip,sizeof(double*));
  for (int i=0; i<2; i++){count[i] = (double*) calloc(nofbins,sizeof(double));}
  double* rate = (double*) malloc(sizeof(double));
  double rnd;

  double check; 
  prob[0] = exp(-var); prob_n[0] = exp(-var_n);
  for(int i=0; i<400; i++){
    prob[i+1] = prob[i] + exp(-var)*pow(var,i+1)/factorial(i+1);
    prob_n[i+1] = prob_n[i] + exp(-var_n)*pow(var_n,i+1)/factorial(i+1);
  }
  for(int i=0; i<nofbins; i++){
    rnd = (double)rand()/RAND_MAX;
    for (int j=0; j<400; j++){
      check=prob[j];
      if(rnd < check){
	randphot[i] = j;  
	break;
      }
    }
    rnd = (double)rand()/RAND_MAX;
    for (int j=0; j<400; j++){
      check=prob_n[j];
      if(rnd < check){
	randphot_n[i] = j;  
	break;
      }
    }
  }
  for(int i=0; i<nofbins; i++){
    rate[0] = sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+offset*PI,PI);
    for(int j=1; j<nofstrip-1; j++){
      rate[j] = sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+j)*PI,PI);}

    for(int j=0; j<400; j++){
      if(j==randphot[i]){
	break;
      }
      rnd = (double)rand()/RAND_MAX;
      for(int k=0;k<nofstrip-1;k++){
	if(rnd < rate[k]){
	  ++count[k][i]; break;}
	if(k == nofstrip -2 ){
	  ++count[nofstrip-1][i];}
      }
    }
    for(int j=0; j<400; j++){
      if(j==randphot_n[i]){
	break;
      }
      rnd = (double)rand()/RAND_MAX*2;
      count[(int)rnd][i] += 1;
    }
  }
  printf("%f %f %f %f\n",fsum(randphot,256),fsum(randphot_n,256),fsum(count[0],256),fsum(count[1],256));
  return count;
}

double **mod(double theta, double phi, double offset, double L, double d) 
{
  int i,j;
  theta = theta*PI; phi = phi*PI;

  double **c = (double**)malloc(nofstrip*sizeof(double*));
  for(i=0; i<nofstrip; i++){c[i]=(double*)malloc(256*sizeof(double));} 
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
    check = floor(x/period+300);
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
