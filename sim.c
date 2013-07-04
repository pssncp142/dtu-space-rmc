#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "time.h"

#define PI 3.14159f
#define RANDOM_SKY 1

#include "common.h"
#include "param.h"

int random_sky(double*,double*,double*,int);

double max_angle = PI/3;

int main(){

  configure();
  
  FILE* f;
  int i,j,k,m,n,ndx,x_ndx,y_ndx,found,st,range=4;
  double max;
  double tmp_s[2];
  double init_obs[2000]={0};
  double obs[2000]={0};
  double map[70000]={0};
  double banned[70000]={0};
  double old_fit[400]={0};
  double fit[400]={0};
  double sources[100]={0};
  double count,totcount;
  double theta[20] = 
    {0.2,0.2,0.3,0.2,0.1,0.15,0.3,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double phi[20] = 
    {0.8,0.2,0.7,-0.2,1.5,1,0.4,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double nofphot[20] = 
    {10000,1000,1200,900,900,1500,700,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double noise = 10000.;
  int n_source = 10;
  int turn = 100;

#if RANDOM_SKY == 1
  random_sky(theta,phi,nofphot,n_source);
#endif

  printf("\n\n");
  printf("DTU SPACE - RMC simulation and source detection simulation\n\n");
  printf(" INPUTS\n---------\n\n");
  printf("  **Configuration Settings\n\n");
  printf("   -Mask Height    : %5.2f\n",L);
  printf("   -Strip Length   : %5.2f\n",d);
  printf("   -Offset         : %5.2f\n",offset);
  printf("   -Rotation       : %d tour\n\n",turn);
  printf("  **Source Settings\n\n");
#if RANDOM_SKY == 1
  printf("   !!Sources are picked randomly");
#endif
  printf("   - # of sources  : %d\n",n_source);
  
  totcount = noise;
  for(i=0;i<n_source;i++) totcount += nofphot[i];

  for(i=0;i<n_source;i++){
    printf("    - %2d --> Theta : %5.2f*PI    Phi : %5.2f*PI    Intensity : %7.2f   Relative : %5.2f perc\n"
	   ,i+1,theta[i],phi[i],nofphot[i],nofphot[i]*100/totcount);
  }
  printf("   - Background noise : %5.2f \n",noise);
  printf("\n\nMonte-Carlo simulation is started...\n");

  totcount = real(obs,n_source,theta,phi,nofphot,noise,turn);
  count = totcount;
  for(i=0;i<2000;i++) {
    init_obs[i] = obs[i];
  }
 
  printf("%5.0f number of photons are created and replaced...\n",count);
  printf("Simulation is completed...\n\n\n");

  printf("Source analysis by removal of the strongest one...\n\n");

  n_source = 0;

  for(k=0;;k++){
   
  again:

    printf("\nIteration %d\n\n",k+1);
    st = corr(map,obs,k);

    n_source = loc_source(sources,map,banned,n_source);

    st = lsf(fit,init_obs,sources,n_source+2);

    printf("Three highest is picked...\n");

    for(i=n_source-1;i<n_source+2;i++){
      printf("  - %2d --> Theta : %5.2f*PI    Phi : %5.2f*PI  \n",i+2-n_source,sources[2*i],sources[2*i+1]);
    }

    max = 0;
    for(i=0;i<3;i++) {
      if(fit[n_source+i] > max) {
	max = fit[n_source+i];
	found = n_source+i-1;
      }
    }

    sources[2*(n_source-1)]   = sources[2*(found)];
    sources[2*(n_source-1)+1] = sources[2*(found)+1];
    for(i=-range;i<range+1;i++){
      for(j=-range;j<range+1;j++){
	banned[((int)sources[2*(found)+7]+j)*256+((int)sources[2*(found)+6]+i)] = 1;
      }
    }

    st = lsf(fit,init_obs,sources,n_source);
 
    if(fit[n_source]<0.02){
      /*printf("Searching limit is reached. Sorting the sources...\n");
      --n_source;
      for(j=0;j<n_source;j++){
	max = 0;
	found = j;
	for(i=j;i<n_source;i++){
	  if(old_fit[i+1]>max){
	    max = old_fit[i+1];
	    found = i;
	  }
	  //printf("%f %d\n ",max*totcount/(opening*turn),found);
	}
	//printf("%d %d\n",found,j);
	if(found!=j){
	  //printf("%d\n",found);
	  tmp_s[0] = sources[2*found];
	  tmp_s[1] = sources[2*found+1];
	  for(i=found-1;i==j;i--){
	    sources[2*(i+1)] = sources[2*i];
	    sources[2*(i+1)+1] = sources[2*i+1];
	  }
	  sources[2*j] = tmp_s[0];
	  sources[2*j+1] = tmp_s[1];
	  n_source = j + 1;
	  ++k;
	  for(i=0;i<70000;i++) banned[i]=0;
	  for(i=0;i<n_source;i++){
	    x_ndx = (int) floor(tan(sources[2*i])*cos(sources[2*i+1])*128/tan(max_angle)+128);
	    y_ndx = (int) floor(tan(sources[2*i])*sin(sources[2*i+1])*128/tan(max_angle)+128);
	    printf("%d %d\n",x_ndx,y_ndx);
	    for(m=-range; m<=range; m++){
	      for(n=-range; n<=range; n++){
		banned[(y_ndx+n)*256+(x_ndx+m)] = 1;
	      }
	    }
	  }
	  for(i=0;i<2000;i++) obs[i] = init_obs[i];
	  count = clean(obs,fit,sources,n_source,count,offset,L,d);
	  goto again;
	}
	}*/	

      printf("Quitting...\n");
      break;
    }

    for(i=0;i<=n_source;i++) old_fit[i] = fit[i];

    printf("\n\nCurrent data : \n\n");

    for(i=0;i<n_source;i++){
      printf("  - %2d --> Theta : %5.3f    Phi : %7.3f    Intensity : %7.2f         Relative Intensity : %5.2f perc\n"
	     ,i+1,sources[2*i],sources[2*i+1],fit[i+1]*totcount/(opening*turn),fit[i+1]*100);
    }

    printf("\n  * Background Intensity : %5.2f                                    Relative Intensity : %5.2f perc\n\n"
	   ,fit[0]*totcount/(opening*turn),fit[0]*100);
            
    for(i=0;i<2000;i++) obs[i] = init_obs[i];
    count = clean(obs,fit,sources,n_source,count);

    printf("Rest counts : %5.2f\n\n",count);
  }  
  return 1;
}

int random_sky(double theta[],double phi[], double nofphot[], int n_sources){
  int i;
  double posx,posy;
  int x_ndx,y_ndx,phot;
  
  srand(time(NULL));
  for(i=0;i<n_sources;i++){
    posx = 2*L*tan(max_angle)*0.9*(((double)rand()/RAND_MAX)-0.5);
    posy = 2*L*tan(max_angle)*0.9*(((double)rand()/RAND_MAX)-0.5);

    theta[i] = atan(sqrt((posx*posx+posy*posy)/(L*L)))/PI;
    if(posx < 0 && posy > 0){
      phi[i] = atan(posy/posx)/PI + 1;
    }else if(posx < 0 && posy < 0){
      phi[i] = atan(posy/posx)/PI - 1;
    } else {
      phi[i] = atan(posy/posx)/PI;
    }
    nofphot[i] = ((double)rand()/RAND_MAX)*1000 + 400;
  }
  return 1; 
}
