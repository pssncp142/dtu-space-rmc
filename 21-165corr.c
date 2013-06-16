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
  max_angle = PI/3; L = 50; d = 5; offset = -0.26;
  FILE* file;
  char fname[][10] = {"A.txt","B.txt","C.txt","D.txt","E.txt","F.txt"};
  float real_w[(int)nofbin][6]; float real_obs[(int)nofbin][6];
  double theta[(int)nofgrid][(int)nofgrid]; double phi[(int)nofgrid][(int)nofgrid];
  double sum[(int)nofbin];
  double max,min;
  double sc_fact[(int)nofgrid][(int)nofgrid];
  double ***model = (double***)malloc(nofgrid*sizeof(double**));
  for(int i=0; i<nofgrid; i++){
    model[i]=(double**)malloc(nofgrid*sizeof(double*));
    for(int j=0; j<nofbin; j++){
      model[i][j]=(double*)malloc(nofbin*sizeof(double));
    }
  }

  double k = PI*L/d;

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

  sum[0] = 0; sum[1] = 0; sum[2] = 0; sum[3] = 0; sum[4] = 0; sum[5] = 0;
  file = fopen("21-165real.txt","r");
  int ctr = 0;
  for (int i=0; i<nofbin; i++){
    fscanf(file,"%f %f %f %f %f %f\n",*(real_w+i),*(real_w+i)+1,*(real_w+i)+2,*(real_w+i)+3,*(real_w+i)+4,*(real_w+i)+5);
    sum[0] += **(real_w+i);
    sum[1] += *(*(real_w+i)+1);
    sum[2] += *(*(real_w+i)+2);
    sum[3] += *(*(real_w+i)+3);
    sum[4] += *(*(real_w+i)+4);
    sum[5] += *(*(real_w+i)+5);
  }
  fclose(file);
  for (int i=0; i<nofbin; i++){
    real_obs[i][0] = real_w[i][0] - sum[0]/nofbin;
    real_obs[i][1] = real_w[i][1] - sum[1]/nofbin;
    real_obs[i][2] = real_w[i][2] - sum[2]/nofbin;
    real_obs[i][3] = real_w[i][3] - sum[3]/nofbin;
    real_obs[i][4] = real_w[i][4] - sum[4]/nofbin;
    real_obs[i][5] = real_w[i][5] - sum[5]/nofbin;
  }

  for(int q=0; q<6; q++){
    for(int i=0; i<nofbin; i++){
      sum[i] = 0;
    }
    for(int i=0; i<nofgrid; i++){
      for(int j=0; j<nofgrid; j++){
	for(int l=0; l<nofbin; l++){
	  model[i][j][l] = sawtooth(k*tan(theta[i][j])*cos(l*2*PI/256.-phi[i][j])+(offset+q)*PI,PI);
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
    
    file = fopen(fname[q],"w+");
    for(int i=0; i<nofgrid; i++){
      for(int j=0; j<nofgrid; j++){
	sc_fact[i][j] = 0;
      }
    }
    for(int i=0; i<nofgrid; i++){
      for(int j=0; j<nofgrid; j++){
	for(int l=0; l<nofbin; l++){
	  sc_fact[i][j] = sc_fact[i][j] + real_w[l][q]*real_obs[l][q]*model[i][j][l];
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
}

double sawtooth(double x, double period)
{
  uint check;
  if(x/(period)<0) {
    check = floor(x/period+300);
  } else {
    check = floor(x/period);
  }
  if ((check%6 == 0)){
    return -(x-check*period)/period+floor((x-check*period)/period)+1;
  } else if ((check%6 == 3)){
    return -(x-(check+3)*period)/period+floor((x-(check+3)*period)/period)+1;
  } else if ((check%6 == 1)){
    return 0;
  } else if ((check%6 == 2)){
    return (x-(check+2)*period)/period-floor((x-(check+2)*period)/period);
  } else if ((check%6 == 4)){
    return (x-(check+4)*period)/period-floor((x-(check+4)*period)/period);
  } else {
    return 1;
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
