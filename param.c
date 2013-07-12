#include "common.h"
#include "param.h"
#include "3d_cart_vec.h"

#include "stdio.h"
#include "math.h"

#define PI 3.14159f

//parameters for the simulation
extern double L,d,offset;
extern double det,mask,height;
extern double alpha,beta,alpha_proj;
extern double opening;
extern double thick;
extern int sp;

//some constants to fasten simulation
extern double PI2_over_256;
extern double PIL_over_d;
extern double st_ob;
extern double det_2;
extern double det_2_sq;
extern double offset_PI[6];
extern double sp_2;

int configure(int a){

  int i;

  if(a){

    L      = 20; 
    d      = 1; 
    offset = 0.25;
    det    = 55; 
    mask   = 85; 
    height = 20;
    thick  = 0.;
    alpha  = 0.*PI;
    beta   = -0.005; //1e-20*PI;
    sp     = n();
    opening= op();
  
    PI2_over_256 = 2*PI/256;
    PIL_over_d   = PI*L/d;
    st_ob        = sp*opening;
    det_2        = det;
    det_2_sq     = pow(det_2,2);
    det         *= 0.5;
    mask        *= 0.5;
    sp_2         = sp*sp;

  }

  for(i=0;i<6;i++) offset_PI[i] = (offset+i)*PI;

  double tel_axis[3], rot_axis[3], sun_axis[3], x_axis[3];
  
  tel_axis[0] = sin(beta)*cos(0); tel_axis[1] = sin(beta)*sin(0); tel_axis[2] = cos(beta);
  rot_axis[0] = 0.; rot_axis[1] = 0.; rot_axis[2] = 1.;
  x_axis[0] = 1.; x_axis[1] = 0; x_axis[2] = 0;

  vec_ort_wc(sun_axis,rot_axis,tel_axis);
  vec_scap(-1,sun_axis);
  vec_rotate(tel_axis,-alpha,sun_axis);
  alpha_proj = vec_angle(rot_axis,sun_axis,x_axis);

  return 0;
}
