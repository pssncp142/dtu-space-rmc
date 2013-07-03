#include "common.h"
#include "param.h"
#include "math.h"

#define PI 3.14159f

//parameters for the simulation
extern double L,d,offset;
extern double det,mask,height;
extern double alpha,beta;
extern double sp,opening;

//some constants to fasten simulation
extern double PI2_over_256;
extern double PIL_over_d;
extern double st_ob;
extern double det_2;
extern double det_2_sq;
extern double offset_PI[6];
extern double sp_2;

int configure(){
  int i;

  L      = 20; 
  d      = 1; 
  offset = 0.;
  det    = 55; 
  mask   = 85; 
  height = 20;
  alpha  = 0.25*PI;
  beta   = 0.2*PI;
  sp     = n();
  opening= op();
  
  PI2_over_256 = 2*PI/256;
  PIL_over_d   = PI*L/d;
  st_ob        = sp*opening;
  det_2        = det;
  det_2_sq     = pow(det_2,2);
  det         *= 0.5;
  mask        *= 0.5;
  sp_2        *= sp;

  for(i=0;i<6;i++) offset_PI[i] = (offset+i)*PI;
}
