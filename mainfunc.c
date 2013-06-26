/*******************************************************************
 * Yigit Dallilar 23.06.2013                                       *
 * Main functions for the rmc simulation and analysis...           *
 *******************************************************************/

/*-----------------------------------------------------------------*\
|  NOTE : this file should be compiled with one of the mask         |
|     specific codes.                                               |
|                                                                   |
|  variables needed in calculations :                               |
|                                                                   |
| %  n_source - number of sources                                   |
| %  *theta - off-axis angle for the source - divided by PI         |
| %  *phi - azimuth angle - divided by PI                           |
| %  offset - relative slide between upper grid and detector strips |
| %  L - height between detector and the mask                       |
| %  d - strip thickness                                            |
| %  *nofphot - number of photons to be executed                    |
| %  noise - needed for the determination of backgroud photons      |
|      background photons = nofphot*noise                           |
| %  turn - number of full rotation                                 |
|                                                                   |
|  if the variable is not one of them it should be an output        |
\------------------------------------------------------------------*/

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

double clean(double obs[], double fit[], double sources[], int n_source, double count, double offset, double L, double d){

  int i,j,ndx,st;
  double model[2000];
  double max=-1;
 
  for(j=0;j<n_source;j++){
    st = mod(model,sources[j*2],sources[j*2+1],offset,L,d);
    for(i=0;i<nofstrip*256;i++) obs[i] -= fit[j+1]*model[i]*count/256;
    count -= fit[j+1]*count;
  }

  return count;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//n_source is changed with k !!!!!!
int loc_source(double sources[], double map[], double banned[], int n_source, double L){  
  
  int i,j,k,l,m,n,st;
  init();

  int x_ndx, y_ndx;
  int range=3;
  double sum=0;
  double max, posx=0, posy=0, max_angle=PI/3;
  max = -1;
  
  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      if(banned[j*256+i]==0){
	if(map[j*256+i] > max){
	  max = map[j*256+i];
	  x_ndx = i;
	  y_ndx = j;
	}
      }
    }
  }

  for(i=-range;i<range+1;i++){
    for(j=-range;j<range+1;j++){
      banned[(y_ndx+j)*256+(x_ndx+i)] = 1;
      posx += (((x_ndx+i)-256*0.5)*L*tan(max_angle)/(256*0.5)+0.5*L*tan(max_angle)/(256*0.5))*map[(y_ndx+j)*256+(x_ndx+i)];
      posy += (((y_ndx+j)-256*0.5)*L*tan(max_angle)/(256*0.5)+0.5*L*tan(max_angle)/(256*0.5))*map[(y_ndx+j)*256+(x_ndx+i)];
      sum += map[(y_ndx+j)*256+(x_ndx+i)];
    }
  }

  posx /= sum;
  posy /= sum;

  sources[2*n_source] = atan(sqrt((posx*posx+posy*posy)/(L*L)))/PI;
  if(posx < 0 && posy > 0){
    sources[2*n_source+1] = atan(posy/posx)/PI + 1;
  }else if(posx < 0 && posy < 0){
    sources[2*n_source+1] = atan(posy/posx)/PI - 1;
  } else {
    sources[2*n_source+1] = atan(posy/posx)/PI;
  }

  return n_source+1;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//this function does the correlation map for the source concerning given values.
//data[k*256*256+j*256+i]
// k --> strip number
// j --> jth grid on the y-axis
// i --> ith grid on the x-axis
int corr(double map[], double obs[], double L, double d, double offset){  
  
  int i,j,k,l,st;
  init();
  
  double posx,posy;
  double theta,phi;
  double max_angle = PI/3;
  double data[500000];
  double model[2000],sbtr_obs[2000];
  double half_map = L*tan(max_angle);
  double half_st = 0.5*L*tan(max_angle)/(256*0.5);
  double all_st = L*tan(max_angle)/(256*0.5);
  double L_2 = L*L;
  int grid_2 = 256*256; 
  for(i=0;i<nofstrip*256;i++) sbtr_obs[i] = obs[i];
  st = subtmean(sbtr_obs,nofstrip);

  
  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      posx = i*all_st+half_st-half_map;
      posy = j*all_st+half_st-half_map;
      theta = atan(sqrt((posx*posx+posy*posy)/(L_2)));
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
	data[k*grid_2+j*256+i] = mult3sum(model,sbtr_obs,obs,k);
      }
    }
  }
  
  norm1_1(data,nofstrip);

  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      map[j*256+i] = data[nofstrip*grid_2+j*256+i];
    }
  }
  
  return 1;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//returns least square fitting result first nofstrip data belongs to each strip nofstrip+1 value is for total fit
int lsf(double fit[], double obs[], double sources[], int n_source, double L, double d, double offset)
{
  init();
  int i,j,k,l,st;

  double LHS[700] = {0};
  for(i=0;i<400;i++) fit[i] = 0;
  double model[2000]; double tmp_model[2000]; double tmp_obs[2000];
  double nofstrip_2 = nofstrip*nofstrip;
  
  for(i=0;i<nofstrip*256;i++) tmp_obs[i] = obs[i];
  st = norm0_1(tmp_obs,nofstrip);

  LHS[0] = (1./256)/(nofstrip);
  fit[0] = (1./256)/(nofstrip);
 
  for(j=1;j<n_source+1;j++){
    st = mod(model,sources[(j-1)*2],sources[(j-1)*2+1],offset,L,d);
    st = norm0_1(model,nofstrip);
    for(i=0;i<256*nofstrip;i++){
      fit[j] += (tmp_obs[i]*model[i])/(nofstrip_2);
      LHS[j] += (model[i]/256)/(nofstrip_2);
      LHS[j*(n_source+1)] = LHS[j]; 
      LHS[j*(n_source+1)+j] += pow(model[i],2)/(nofstrip_2);
    }
    for(k=j+1;k<n_source+1;k++){
      st = mod(tmp_model,sources[(k-1)*2],sources[(k-1)*2+1],offset,L,d);
      st = norm0_1(tmp_model,nofstrip);
      for(i=0;i<256*nofstrip;i++){
	LHS[j*(n_source+1)+k] += (tmp_model[i]*model[i])/(nofstrip_2);
	LHS[k*(n_source+1)+j] = LHS[j*(n_source+1)+k];
      }
    }
  }

  st =  gaussj_nr(LHS,n_source+1,fit,n_source+1);
  
  return 1;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//returns result of the monte carlo simulation for every strip
//count[256*j+i]
// j --> strip number
// i --> related phase bin
double real(double obs[], int n_source, double* theta, double* phi, double offset, double L, double d, double* nofphot, double noise, int turn) 
{
  init();
  int i,j,k,l,m,st;
  
  double count = 0;
  for(i=0;i<2000;i++) obs[i] = 0;
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
	      ++obs[k*256+i]; break;}
	    if(k == nofstrip-2 ){
	      ++obs[(nofstrip-1)*256+i];}
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
	++obs[(int)rnd*256+i];
      }
    }
  }

  for(i=0;i<nofstrip*256;i++) count += obs[i];
  
  return count;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//nofstrip*opening

//returns modulation function
//data[256*j+i]
// j --> strip number
// i --> related phase bin
int mod(double model[], double theta, double phi, double offset, double L, double d) 
{
  init();
  int i,j,k,l,st;
  double PIL_over_d = PI*L/d;
  double PI2_over_256 = 2*PI/256;

  for(i=0; i<256; i++){
    for(j=0; j<nofstrip; j++){
      model[j*256+i] = sawtooth(PIL_over_d*tan(theta*PI)*cos(i*PI2_over_256-phi*PI)+(offset+j)*PI,PI)/(nofstrip*opening);
    }
  }
  return 1;
}

/*****************************************************************/
