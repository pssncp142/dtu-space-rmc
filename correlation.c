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
  double max_angle, nofgrid, L;

  max_angle = PI/3; nofgrid = 10; L = 10;

  double theta[(int)nofgrid][(int)nofgrid]; double phi[(int)nofgrid][(int)nofgrid];
  for (int i=0; i<nofgrid; i++){
    for (int j=0; j<nofgrid; j++){
      double posx = (i-nofgrid*0.5)*L*tan(max_angle)/(nofgrid*0.5)+L*tan(max_angle)/nofgrid;
      double posy = (j-nofgrid*0.5)*L*tan(max_angle)/(nofgrid*0.5)+L*tan(max_angle)/nofgrid;
      theta[i][j] = atan(sqrt((posx*posx+posy*posy)/(L*L)));
      if(posx < 0 && posy > 0){
	phi[i][j] = atan(posy/posx) + PI;
	//printf("%f \n",atan(posy/posx) + PI);
      }else if(posx < 0 && posy < 0){
	phi[i][j] = atan(posy/posx) - PI;
      } else {
	phi[i][j] = atan(posy/posx);
      }
      printf("%d %d %f\n",i-5,j-5,theta[i][j]);
    }
  }

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
