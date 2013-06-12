/****************************************************************
 * Yigit Dallilar 12.06.2013                                    *
 * DTU SPACE - cartesian grid correlation                       *
 ****************************************************************/

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define PI 3.14159f

double poisson_formula(double, int);
double factorial(int);
double sawtooth(double,double);

int main()
{
  //printf("12");
  
  double max_angle, nofgrid, L, d,  nofbin, offset;

  max_angle = PI/3; nofgrid = 50; L = 50; d = 5;
  nofbin = 256; offset = 0.25;

  double k = PI*L/d; offset = 2*PI*offset;

  double theta[(int)nofgrid][(int)nofgrid]; double phi[(int)nofgrid][(int)nofgrid];
  for (int i=0; i<nofgrid; i++){
    for (int j=0; j<nofgrid; j++){
      double posx = (i-nofgrid*0.5)*L*tan(max_angle)/(nofgrid*0.5)+L*tan(max_angle)/(nofgrid*0.5);
      double posy = (j-nofgrid*0.5)*L*tan(max_angle)/(nofgrid*0.5)+L*tan(max_angle)/(nofgrid*0.5);
      theta[i][j] = atan(sqrt((posx*posx+posy*posy)/(L*L)));
      if(posx < 0 && posy > 0){
	phi[i][j] = atan(posy/posx) + PI;
      }else if(posx < 0 && posy < 0){
	phi[i][j] = atan(posy/posx) - PI;
      } else {
	phi[i][j] = atan(posy/posx);
      }
    }
  }

  double model[(int)nofgrid][(int)nofgrid][(int)nofbin];
  double sum[(int)nofbin];
  for(int i=0; i<nofbin; i++){
    sum[i] = 0;
  }
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      for(int l=0; l<nofbin; l++){
	model[i][j][l] = sawtooth(k*tan(theta[i][j])*cos(l*2*PI/256.-phi[i][j])+offset,2*PI);
	sum[l] = sum[l] + model[i][j][l];
      }
    }
  }
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      for(int l=0; l<nofbin; l++){
	model[i][j][l] = model[i][j][l] - sum[l]/nofbin;
      }
    }
  }

  FILE* file;
  float real_w[(int)nofbin]; float real_obs[(int)nofbin];
  float dummy; sum[0] = 0;
  file = fopen("realdata.txt","r");
  int ctr = 0;
  for (int i=0; i<nofbin; i++){
    fscanf(file,"%f %f \n",real_w+i,&dummy);
    sum[0] = sum[0] + *(real_w+i);
  }
  fclose(file);
  for (int i=0; i<nofbin; i++){
    real_obs[i] = real_w[i] - sum[0]/nofbin;
  }

  file = fopen("trial.txt","w+");
  double sc_fact[(int)nofgrid][(int)nofgrid];
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      sc_fact[i][j] = 0;
    }
  }
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      sum[0] =0;
      for(int l=0; l<nofbin; l++){
	sc_fact[i][j] = sc_fact[i][j] + real_w[l]*real_obs[l]*model[i][j][l];
      }
      fprintf(file,"%f ",sc_fact[i][j]);
    }
    fprintf(file,"\n");
  }
  fclose(file);

}

double sawtooth(double x, double period)
{
  uint check = floor(x/(2*PI));
  if (check%2 == 0 ){
    return (x-check*2*PI)/period-floor((x-check*2*PI)/period);
  } else {
    return -(x-check*2*PI)/period+floor((x-check*2*PI)/period)+1;
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
