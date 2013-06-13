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
  
  const double nofgrid = 256;
  const double nofbin = 256;

  double max_angle, L, d, offset;

  max_angle = PI/3; L = 50; d = 5; offset = 0.5;

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

  double ***model = (double***)malloc(nofgrid*sizeof(double**));
  for(int i=0; i<nofgrid; i++){
    model[i]=(double**)malloc(nofgrid*sizeof(double*));
    for(int j=0; j<nofbin; j++){
      model[i][j]=(double*)malloc(nofbin*sizeof(double));
    }
  }
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
  float real_w[(int)nofbin][2]; float real_obs[(int)nofbin][2];
  float dummy; sum[0] = 0;
  file = fopen("realdata.txt","r");
  int ctr = 0;
  for (int i=0; i<nofbin; i++){
    fscanf(file,"%f %f \n",*(real_w+i),*(real_w+i)+1);
    sum[0] += **(real_w+i);
    sum[1] += *(*(real_w+i)+1);
  }
  fclose(file);
  for (int i=0; i<nofbin; i++){
    real_obs[i][0] = real_w[i][0] - sum[0]/nofbin;
    real_obs[i][1] = real_w[i][1] - sum[1]/nofbin;
  }

  file = fopen("A.txt","w+");
  double max,min;
  double sc_fact[(int)nofgrid][(int)nofgrid];
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      sc_fact[i][j] = 0;
    }
  }
  
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      for(int l=0; l<nofbin; l++){
	sc_fact[i][j] = sc_fact[i][j] + real_w[l][0]*real_obs[l][0]*model[i][j][l];
      }
      if(i==0 && j==0){
	max = sc_fact[i][j];
	min = sc_fact[i][j];
      } else {
	if(sc_fact[i][j] > max){max = sc_fact[i][j];}
	if(sc_fact[i][j] < min){min = sc_fact[i][j];}
      }
    }
  }
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      fprintf(file,"%f ",2*(sc_fact[i][j]-min)/(max-min)-1);
    }
    fprintf(file,"\n");
  }
  fclose(file);

  file = fopen("B.txt","w+");
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      sc_fact[i][j] = 0;
    }
  }
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      for(int l=0; l<nofbin; l++){
	sc_fact[i][j] = sc_fact[i][j] - real_w[l][1]*real_obs[l][1]*model[i][j][l];
      }
      if(i==0 && j==0){
	max = sc_fact[i][j];
	min = sc_fact[i][j];
      } else {
	if(sc_fact[i][j] > max){max = sc_fact[i][j];}
	if(sc_fact[i][j] < min){min = sc_fact[i][j];}
      }
    }
  }
  for(int i=0; i<nofgrid; i++){
    for(int j=0; j<nofgrid; j++){
      fprintf(file,"%f ",2*(sc_fact[i][j]-min)/(max-min)-1);
    }
    fprintf(file,"\n");
  }
  fclose(file);

}

double sawtooth(double x, double period)
{
  uint check = floor(x/(2*PI));
  if (check%2 == 0 ){
    return (x-(check)*2*PI)/period-floor((x-(check)*2*PI)/period);
  } else {
    return -(x-(check)*2*PI)/period+floor((x-(check)*2*PI)/period)+1;
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
