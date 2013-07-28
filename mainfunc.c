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
| %  *map - correlation map                                         |
| %  *fit - result of lsf calculation 0->background                 | 
| %  *obs - observation data for all strips                         |
| %  *banned - banned locations 1->forbidden 0->available           |
| %  *sources - found source list 2*k-->theta 2*k+1-->phi           |
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
#include "3d_cart_vec.h"
#include "param.h"

#define PI 3.14159f

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/


/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//erase all found sources from the observation data
double clean(double obs[], double fit[], double sources[], int n_source, double count){

  int i,j,ndx,st;
  double model[2000];
  double max=-1;
 
  for(j=0;j<n_source;j++){
    st = mod(model,sources[j*2],sources[j*2+1]);
    for(i=0;i<sp*256;i++) obs[i] -= fit[j+1]*model[i]*count/256;
    count -= fit[j+1]*count;
  }

  return count;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//finds three highest possible sources from the correlation map
//!!!one trick data is outputted as ;
//sources[2*(n_source+k)]   -->theta
//sources[2*(n_source+k)+1] -->phi
//sources[2*(n_source+k)+6] -->x_ndx --> to be able to add the best in the banned map
//sources[2*(n_source+k)+7] -->y_ndx --> to be able to add the best in the banned map
int loc_source(double sources[], double map[], double banned[], int n_source){  
  
  int i,j,k,l,m,n,st;

  int x_ndx[5]={0};
  int y_ndx[5]={0};
  int range=5,search=3,same;
  double max, sum, max_angle=PI/3;
  double posx[5]={0};
  double posy[5]={0}; 

  for(k=0;k<search;k++){

    max = -1;
    sum =  0;
  
    for(i=0;i<256;i++){
      for(j=0;j<256;j++){
	if(banned[j*256+i]==0){
	  same = 0;
	  for(l=0;l<k;l++){
	    if((i<x_ndx[l]+(range+1)) && (i>x_ndx[l]-(range+1)) && (j<y_ndx[l]+(range+1)) && (j>y_ndx[l]-(range+1)))
	      same = 1;
	  }
	  if(!same){
	    if(map[j*256+i] > max){
	      max = map[j*256+i];
	      x_ndx[k] = i;
	      y_ndx[k] = j;
	    }
	  }
	}
      }
    }
    for(i=-1;i<2;i++){
      for(j=-1;j<2;j++){
	posx[k] += (((x_ndx[k]+i)-256*0.5)*L*tan(max_angle)/(256*0.5)+0.5*L*tan(max_angle)/(256*0.5))*map[(y_ndx[k]+j)*256+(x_ndx[k]+i)];
	posy[k] += (((y_ndx[k]+j)-256*0.5)*L*tan(max_angle)/(256*0.5)+0.5*L*tan(max_angle)/(256*0.5))*map[(y_ndx[k]+j)*256+(x_ndx[k]+i)];
	sum += map[(y_ndx[k]+j)*256+(x_ndx[k]+i)];
      }
    }
    posx[k] /= sum;
    posy[k] /= sum;

    sources[2*(n_source+k)] = atan(sqrt((posx[k]*posx[k]+posy[k]*posy[k])/(L*L)))/PI;
    if(posx[k] < 0 && posy[k] > 0){
      sources[2*(n_source+k)+1] = atan(posy[k]/posx[k])/PI + 1;
    }else if(posx[k] < 0 && posy[k] < 0){
      sources[2*(n_source+k)+1] = atan(posy[k]/posx[k])/PI - 1;
    } else {
      sources[2*(n_source+k)+1] = atan(posy[k]/posx[k])/PI;
    }
    sources[2*(n_source+k)+6] = x_ndx[k];
    sources[2*(n_source+k)+7] = y_ndx[k];
  
  }

  return n_source+1;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//this function does the correlation map for the source concerning given values.
//data[k*256*256+j*256+i]
// k --> strip number
// j --> jth grid on the y-axis
// i --> ith grid on the x-axis
int corr(double map[], double obs[], int check, char fname[]){  
  
  int i,j,k,l,st;
  
  FILE* f;
  double posx,posy;
  double theta,phi;
  double max_angle = PI/3;
  double data[500000];
  double model[2000],sbtr_obs[2000];
  double half_map = L*tan(max_angle);
  double half_st = 0.5*L*tan(max_angle)/(256*0.5);
  double st_map = half_st-half_map;
  double all_st = L*tan(max_angle)/(256*0.5);
  double L_2 = L*L;
  int grid_2 = 256*256; 
  char *bin;

  for(i=0;i<sp*256;i++) sbtr_obs[i] = obs[i];
  st = subtmean(sbtr_obs,sp);
  

  if(check == 0){
    f=fopen(fname,"w+");
    printf("Modulation function for all sky map will be written to : %s \n",fname);
    printf("Correlation process :      "); fflush(stdout);
    for(i=0;i<256;i++){
      if(0==i % 25) printf("\b\b\b\b\b%3d%% ",i/25*10); fflush(stdout);
      for(j=0;j<256;j++){
	posx = i*all_st+st_map;
	posy = j*all_st+st_map;
	theta = atan(sqrt((posx*posx+posy*posy)/(L_2)));
	if(posx < 0 && posy > 0){
	  phi = atan(posy/posx) + PI;
	}else if(posx < 0 && posy < 0){
	  phi = atan(posy/posx) - PI;
	} else {
	  phi = atan(posy/posx);
	}
	st = mod(model,theta/PI,phi/PI);
	st = subtmean(model,sp);
	fwrite(model,sizeof(double),256*sp,f);
	for(k=0;k<sp;k++){
	  data[k*grid_2+j*256+i] = mult3sum(model,sbtr_obs,obs,k);
	}
      }
    }
    fclose(f);
  } else {
    f=fopen(fname,"r");
    //printf("Correlation process :      "); fflush(stdout);
    for(i=0;i<256;i++){
      //if(0==i % 25) printf("\b\b\b\b\b%3d%% ",i/25*10); fflush(stdout);
      for(j=0;j<256;j++){
	fread(model,sizeof(double),256*sp,f);
	for(k=0;k<sp;k++){
	  data[k*grid_2+j*256+i] = mult3sum(model,sbtr_obs,obs,k);
	}
      }
    }
    fclose(f);
  }

  
  //printf("\n\n");

  norm1_1(data,sp);

  for(i=0;i<256;i++){
    for(j=0;j<256;j++){
      map[j*256+i] = data[sp*grid_2+j*256+i];
    }
  }
  
  return 1;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

//returns least square fitting result first sp data belongs to each strip sp+1 value is for total fit
int lsf(double fit[], double obs[], double sources[], int n_source)
{

  int i,j,k,l,st;

  double LHS[700] = {0};
  for(i=0;i<400;i++) fit[i] = 0;
  double model[2000]; double tmp_model[2000]; double tmp_obs[2000];
  
  for(i=0;i<sp*256;i++) tmp_obs[i] = obs[i];
  st = norm0_1(tmp_obs,sp);

  LHS[0] = (1./256)/(sp);
  fit[0] = (1./256)/(sp);
 
  for(j=1;j<n_source+1;j++){
    st = mod(model,sources[(j-1)*2],sources[(j-1)*2+1]);
    st = norm0_1(model,sp);
    for(i=0;i<256*sp;i++){
      fit[j] += (tmp_obs[i]*model[i])/(sp_2);
      LHS[j] += (model[i]/256)/(sp_2);
      LHS[j*(n_source+1)] = LHS[j]; 
      LHS[j*(n_source+1)+j] += pow(model[i],2)/(sp_2);
    }
    for(k=j+1;k<n_source+1;k++){
      st = mod(tmp_model,sources[(k-1)*2],sources[(k-1)*2+1]);
      st = norm0_1(tmp_model,sp);
      for(i=0;i<256*sp;i++){
	LHS[j*(n_source+1)+k] += (tmp_model[i]*model[i])/(sp_2);
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
double real(double obs[], int n_source, double* theta, double* phi, double* nofphot, double noise, int turn) 
{
  int i,j,k,l,m,st;
  
  double count = 0;
  for(i=0;i<2000;i++) obs[i] = 0;
  noise *= opening;
  double prob[1000], rate[10], model[2000];
  int randphot;
  double rnd,check,var; 
  double PIL_over_d_tan_theta[50]; for(i=0;i<n_source;i++) PIL_over_d_tan_theta[i] = PI*L*tan(theta[i]*PI)/d;
  double phi_PI[50]; for(i=0;i<n_source;i++) phi_PI[i] = phi[i]*PI;

  srand(time(NULL));

  for(m=0;m<n_source;m++){
    mod(model,theta[m],phi[m]); 
    var = nofphot[m]*opening/256;
    prob[0] = exp(-var);
    for(l=1; l<1000; l++){
      prob[l] = prob[l-1] + exp(-var)*pow(var,l)/factorial(l);
    }
    for(l=0;l<turn;l++){
      for(i=0; i<256; i++){
	rate[0] = model[i]; 
	for(j=1; j<sp; j++){
	  rate[j] = rate[j-1] + model[j*256+i];
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
	  rnd = (double)rand()/RAND_MAX;
	  for(k=0;k<sp;k++){
	    if(rnd < rate[k]){
	      ++obs[k*256+i]; break;}
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
    for(i=0; i<256; i++){
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
	rnd = (double)rand()/RAND_MAX*sp;
	++obs[(int)rnd*256+i];
      }
    }
  }

  for(i=0;i<sp*256;i++) count += obs[i];
  
  return count;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//returns modulation function
//data[256*j+i]
// j --> strip number
// i --> related phase bin
int mod(double model[], double theta, double phi) 
{
  int i,j,k,st;  
  
  double thetap[300], phip[300], frac[300]={0};
  theta *= PI;
  phi *= PI;

  calc_angles(thetap,phip,theta,phi);

  frac_detect(frac,thetap,phip);
  for(i=0;i<256;i++) frac[i] = 1;

  for(i=0; i<256; i++){
    for(j=0; j<sp; j++){
      model[j*256+i] = sawtooth(PIL_over_d*tan(thetap[i])*cos(phip[i])+offset_PI[j],PI)/(st_ob)*frac[i];
    }
  }

  return 1;
}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//calculates the illuminated fraction of the detector for every iteration
//(normalized to 1)

int frac_detect(double frac[], double thetap[], double phip[]){

  int i,j,k;

  double pos[4];
  double t_c_h, t_s_h;
  for(i=0;i<256;i++){
    t_c_h = tan(thetap[i])*cos(phip[i])*height;
    t_s_h = tan(thetap[i])*sin(phip[i])*height;
    pos[0] = mask-t_c_h;
    pos[1] = -mask-t_c_h;
    pos[2] = mask-t_s_h;
    pos[3] = -mask-t_s_h;
    if(pos[0] > det && pos[2] > det && pos[1] < -det && pos[3] < -det)  {frac[i]=1; goto exit;}
    for(j=0;j<2;j++){
      for(k=0;k<2;k++){
	if(pos[j] < det && pos[j] > -det){	  
	  if(pos[2+k] < det && pos[2+k] > -det){
	    if(j==0){
	      if(k==0){
		frac[i] = (det+pos[j])*(det+pos[2+k])/det_2_sq; goto exit;
	      } else {
		frac[i] = (det+pos[j])*(det-pos[2+k])/det_2_sq; goto exit;
	      }
	    } else {
	      if(k==0){
		frac[i] = (det-pos[j])*(det+pos[2+k])/det_2_sq; goto exit;
	      } else {
		frac[i] = (det-pos[j])*(det-pos[2+k])/det_2_sq; goto exit;
	      }
	    }
	  }
	}
      }
    }
    for(j=0;j<2;j++){
      if(pos[j] < det && pos[j] > -det){
	if(pos[2] > det && pos[3] < -det){
	  if(j==0){
	    frac[i] = det_2*(det+pos[j])/det_2_sq; goto exit;
	  } else {
	    frac[i] = det_2*(det-pos[j])/det_2_sq; goto exit;
	  }
	}
      }
      if(pos[2+j] < det && pos[2+j] > -det){
	if(pos[0] > det && pos[1] < -det){
	  if(j==0){
	    frac[i] = det_2*(det+pos[2+j])/det_2_sq; goto exit;
	  } else {
	    frac[i] = det_2*(det-pos[2+j])/det_2_sq; goto exit;
	  }
	}
      }
    }
  exit: ;
  }

  //effective area
  for(i=0;i<256;i++)
    frac[i] *= cos(thetap[i])*(1-thick/d*pow(cos(phip[i]),2));

  return 1;

}

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//returns the coordinates in the reference frame of tilted rotating mask

int calc_angles(double thetap[], double phip[], double theta, double phi){

  int i;
  double tel_axis[3], x_axis[3], source[3];
  double s_beta = sin(beta);
  double s_beta_PI = sin(beta+PI);
  
  source[0] = sin(theta)*cos(phi); source[1] = sin(theta)*sin(phi); source[2] = cos(theta);
  tel_axis[2] = cos(beta);
  x_axis[2] = 0;

  for(i=0;i<256;i++){
    
    tel_axis[0] = s_beta*cos(i*PI2_over_256+alpha_proj);
    tel_axis[1] = s_beta*sin(i*PI2_over_256+alpha_proj);
    x_axis[0] = cos(i*PI2_over_256);
    x_axis[1] = sin(i*PI2_over_256);

    thetap[i] = acos(vec_dotp(tel_axis,source));
    phip[i] = vec_angle(tel_axis,x_axis,source);

  }

  return 1;

}

/*****************************************************************/
