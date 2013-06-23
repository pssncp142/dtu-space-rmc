#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#include "common.h"

#define PI 3.14159f

int nofstrip;
double opening;

void init(){
  opening = op();
  nofstrip = n();
}


int corr(double data[], double theta,double phi,double offset,double L,double d,double nofphot,double noise,int turn){
  
  
  int i,j,k,l,st;
  init();
  
  double posx,posy;
  double max_angle = PI/3;
  double model[2000],obs[2000],weight[2000];
  st = real(weight,theta,phi,offset,L,d,nofphot,noise,turn); 
  st = real(obs,theta,phi,offset,L,d,nofphot,noise,turn); 
  st = subtmean(obs,nofstrip);
  
  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      posx = (i-256*0.5)*L*tan(max_angle)/(256*0.5)+L*tan(max_angle)/(256*0.5);
      posy = (j-256*0.5)*L*tan(max_angle)/(256*0.5)+L*tan(max_angle)/(256*0.5);
      theta = atan(sqrt((posx*posx+posy*posy)/(L*L)));
      if(posx < 0 && posy > 0){
	phi = atan(posy/posx) + PI;
      }else if(posx < 0 && posy < 0){
	phi = atan(posy/posx) - PI;
      } else {
	phi = atan(posy/posx);
      }
      st = mod(model,theta/PI,phi/PI,offset,L,d);
      st = subtmean(model,nofstrip);
      for(k=0;k<nofstrip;k++){
	data[k*256*256+j*256+i] = mult3sum(model,obs,weight,k);
      }
    }
  }
  
  
  norm1_1(data,nofstrip);  
  
  return 1;
}

int lsf(double data[], double theta,double phi,double offset,double L,double d,double nofphot,double noise,int turn)
{
  init();
  int i,j,k,l,st;

  data[nofstrip]=0;
  double LHS[100]; double RHS[100];
  double obs[2000]; double model[2000];
  st = real(obs,theta,phi,offset,L,d,nofphot,noise,turn); 
  st = mod(model,theta,phi,offset,L,d);

  st = norm0_1(model,nofstrip);
  st = norm0_1(obs,nofstrip);
 
  for(j=0; j<nofstrip; j++){
  
    LHS[0] = 1./256; LHS[1]=0; LHS[2]=0; LHS[3]=0; RHS[0] = 0; RHS[1] = 0; 

    for(i=0; i<256; i++){
      LHS[1] += model[j*256+i]/256;
      LHS[2] += model[j*256+i]/256;
      LHS[3] += pow(model[j*256+i],2);
      RHS[0] += obs[j*256+i]/256;
      RHS[1] += obs[j*256+i]*model[j*256+i];
    }

    st =  gaussj_nr(LHS,2,RHS,2);
    data[j] = RHS[0]/RHS[1];
  }

  LHS[0] = (1./256)/(nofstrip); LHS[1]=0; LHS[2]=0; LHS[3]=0; RHS[0] = 0; RHS[1] = 0; 

  for(i=0;i<256*nofstrip;i++){
    LHS[1] += (model[i]/256)/(nofstrip*nofstrip);
    LHS[2] += (model[i]/256)/(nofstrip*nofstrip);
    LHS[3] += pow(model[i],2)/(nofstrip*nofstrip);
    RHS[0] += (obs[i]/256)/(nofstrip*nofstrip);
    RHS[1] += (obs[i]*model[i])/(nofstrip*nofstrip);
  }

  st =  gaussj_nr(LHS,2,RHS,2);
  data[nofstrip] = RHS[0]/RHS[1];
  
  return 1;
}

int real(double count[], double theta, double phi, double offset, double L, double d, double nofphot, double noise, int turn) 
{
  init();
  int i,j,k,l,st;

  phi = PI*phi; theta = theta*PI;
  for(i=0;i<2000;i++)
    count[i] = 0;
  double nofbins = 256;
  double probint = 0.05;
  nofphot *= opening;
  noise *= nofphot;
  double PIL_over_d = PI*L/d;
  double var = nofphot/nofbins;
  double var_n = noise/nofbins;
  int dim = 1/probint;
  double prob[1000], prob_n[1000], rate[10];
  int randphot[50000], randphot_n[50000];
  double rnd;
  srand(time(NULL));

  double check; 
  prob[0] = exp(-var); prob_n[0] = exp(-var_n);
  for(i=1; i<1000; i++){
    prob[i] = prob[i-1] + exp(-var)*pow(var,i)/factorial(i);
    prob_n[i] = prob_n[i-1] + exp(-var_n)*pow(var_n,i)/factorial(i);
  }
  for(i=0; i<nofbins*turn; i++){
    rnd = (double)rand()/RAND_MAX;
    for (j=0; j<1000; j++){
      check=prob[j];
      if(rnd < check){
	randphot[i] = j;  
	break;
      }
    }
    rnd = (double)rand()/RAND_MAX;
    for (j=0; j<1000; j++){
      check=prob_n[j];
      if(rnd < check){
	randphot_n[i] = j;  
	break;
      }
    }
  }
  for(l=0;l<turn;l++){
    for(i=0; i<nofbins; i++){
      rate[0] = sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+offset*PI,PI); 
      for(j=1; j<nofstrip-1; j++){
	rate[j] = rate[j-1] + sawtooth(PIL_over_d*tan(theta)*cos(i*2*PI/nofbins-phi)+(offset+j)*PI,PI);
      }
    
      for(j=0; j<2000; j++){
	if(j==randphot[l*256+i]){
	  break;
	}
	rnd = (double)rand()/RAND_MAX*nofstrip*opening;
	for(k=0;k<nofstrip-1;k++){
	  if(rnd < rate[k]){
	    ++count[k*256+i]; break;}
	  if(k == nofstrip-2 ){
	    ++count[(nofstrip-1)*256+i];}
	}
      }
      for(j=0; j<2000; j++){
	if(j==randphot_n[i]){
	  break;
	}
	rnd = (double)rand()/RAND_MAX*nofstrip;
	++count[(int)rnd*256+i];
      }
    }
  }
  return 1;
}

int mod(double model[], double theta, double phi, double offset, double L, double d) 
{
  init();
  int i,j,k,l,st;

  theta = theta*PI; phi = phi*PI;

  for(i=0; i<256; i++){
    for(j=0; j<nofstrip; j++){
      model[j*256+i] = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/256-phi)+(offset+j)*PI,PI);
    }
  }
  return 1;
}
