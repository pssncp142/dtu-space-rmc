/*****************************************************************
 * Yigit Dallilar 23.06.2013                                     *
 * Main functions for the rmc simulation and analysis...         *
 *****************************************************************/

/*----------------------------------------------------------------*\
|  NOTE : this file should be compiled with one of the mask        |
|     specific codes.                                              |
|                                                                  |
|  variables needed in calculations :                              |
|  data[] - is the variable for the output                         |
|  n_source - number of sources                                    |
|  *theta - off-axis angle for the source - divided by PI          |
|  *phi - azimuth angle - divided by PI                            |
|  offset - relative slide between upper grid and detector strips  |
|  L - height between detector and the mask                        |
|  d - strip thickness                                             |
|  *nofphot - number of photons to be executed                     |
|  noise - needed for the determination of backgroud photons       |
|      background photons = nofphot*noise                          |
|  turn - number of full rotation                                  |
\-----------------------------------------------------------------*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#include "common.h"

#define PI 3.14159f

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//global variables which are initialized by init function from the specific "mask-name".c files
int nofstrip;
double opening;

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//doing the initalisation for the global variables. It is necessary to call this function
//from the functions below. So that, they know with which mask they are integrated.
void init(){
  opening = op();
  nofstrip = n();
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//this function does the correlation map for the source concerning given values.
//data[k*256*256+j*256+i]
// k --> strip number
// j --> jth grid on the y-axis
// i --> ith grid on the x-axis
int corr(double data[], int n_source, double *theta,double *phi,double offset,double L,double d,double *nofphot,double noise,int turn){  
  
  int i,j,k,l,st;
  init();
  
  double posx,posy;
  double max_angle = PI/3;
  double model[2000],obs[2000],weight[2000];
  st = real(weight,n_source,theta,phi,offset,L,d,nofphot,noise,turn); 
  st = real(obs,n_source,theta,phi,offset,L,d,nofphot,noise,turn); 
  st = subtmean(obs,nofstrip);
  
  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      posx = (i-256*0.5)*L*tan(max_angle)/(256*0.5)+L*tan(max_angle)/(256*0.5);
      posy = (j-256*0.5)*L*tan(max_angle)/(256*0.5)+L*tan(max_angle)/(256*0.5);
      theta[0] = atan(sqrt((posx*posx+posy*posy)/(L*L)));
      if(posx < 0 && posy > 0){
	phi[0] = atan(posy/posx) + PI;
      }else if(posx < 0 && posy < 0){
	phi[0] = atan(posy/posx) - PI;
      } else {
	phi[0] = atan(posy/posx);
      }
      st = mod(model,theta[0]/PI,phi[0]/PI,offset,L,d);
      st = subtmean(model,nofstrip);
      for(k=0;k<nofstrip;k++){
	data[k*256*256+j*256+i] = mult3sum(model,obs,weight,k);
      }
    }
  }
  
  norm1_1(data,nofstrip);  
  
  return 1;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//returns least square fitting result first nofstrip data belongs to each strip nofstrip+1 value is for total fit
int lsf(double data[], int n_source, double *theta,double *phi,double offset,double L,double d,double *nofphot,double noise,int turn)
{
  init();
  int i,j,k,l,st;

  data[nofstrip]=0;
  double LHS[100]; double RHS[100];
  double obs[2000]; double model[2000];
  st = real(obs,1,theta,phi,offset,L,d,nofphot,noise,turn); 
  st = mod(model,theta[0],phi[0],offset,L,d);

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

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//returns result of the monte carlo simulation for every strip
//count[256*j+i]
// j --> strip number
// i --> related phase bin
int real(double count[], int n_source, double* theta, double* phi, double offset, double L, double d, double* nofphot, double noise, int turn) 
{
  init();
  int i,j,k,l,m,st;

  
  for(i=0;i<2000;i++)
    count[i] = 0;
  double nofbins = 256;
  double probint = 0.05;
  noise *= opening;
  double PIL_over_d = PI*L/d;
  double var;
  int dim = 1/probint;
  double prob[1000], rate[10];
  int randphot;
  double rnd;   double check; 
  srand(time(NULL));

  for(m=0;m<n_source;m++){
    var = nofphot[m]*opening/256;
    prob[0] = exp(-var);
    for(l=1; l<1000; l++){
      prob[l] = prob[l-1] + exp(-var)*pow(var,l)/factorial(l);
    }
    printf("%f %f %f\n",nofphot[m],theta[m],phi[m]);
    for(l=0;l<turn;l++){
      for(i=0; i<nofbins; i++){
	rate[0] = sawtooth(PIL_over_d*tan(theta[m]*PI)*cos(i*2*PI/nofbins-phi[m]*PI)+offset*PI,PI); 
	for(j=1; j<nofstrip-1; j++){
	  rate[j] = rate[j-1] + sawtooth(PIL_over_d*tan(theta[m]*PI)*cos(i*2*PI/nofbins-phi[m]*PI)+(offset+j)*PI,PI);
	}
	rnd = (double)rand()/RAND_MAX;
	for (j=0; j<1000; j++){
	  check=prob[j];
	  if(rnd < check){
	    randphot = j;  
	    break;
	  }
	}
	for(j=0; j<2000; j++){
	  if(j==randphot){
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
      }
    }
  }

  var = noise/256;
  prob[0] = exp(-var);
  for(i=1; i<1000; i++){
    prob[i] = prob[i-1] + exp(-var)*pow(var,i)/factorial(i);
  }
  for(l=0;l<turn;l++){
    for(i=0; i<nofbins; i++){
      rnd = (double)rand()/RAND_MAX;
      for (j=0; j<1000; j++){
	check=prob[j];
	if(rnd < check){
	  randphot = j;  
	  break;
	}
      }
      for(j=0; j<2000; j++){
	if(j==randphot){
	  break;
	}
	rnd = (double)rand()/RAND_MAX*nofstrip;
	++count[(int)rnd*256+i];
      }
    }
  }

  return 1;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//returns modulation function
//data[256*j+i]
// j --> strip number
// i --> related phase bin
int mod(double data[], double theta, double phi, double offset, double L, double d) 
{
  init();
  int i,j,k,l,st;

  theta = theta*PI; phi = phi*PI;

  for(i=0; i<256; i++){
    for(j=0; j<nofstrip; j++){
      data[j*256+i] = sawtooth(PI*L/d*tan(theta)*cos(i*2*PI/256-phi)+(offset+j)*PI,PI);
    }
  }
  return 1;
}

/*****************************************************************/
