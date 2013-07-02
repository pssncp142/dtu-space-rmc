/*****************************************************************\
* Yigit Dallilar 02.07.2013                                       *
*                                                                 *
* C vector algebra library                                        *
\*****************************************************************/

#include "3d_cart_vec.h"

#include "math.h"
#include "stdio.h"

#define PI 3.14159f

int i,j;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//norm of a vector

double vec_norm(double a[]){
  double res=vec_dotp(a,a);
  return sqrt(res);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//unit vector of a

int vec_unit(double a[]){
  double norm = vec_norm(a);
  vec_scap(1/norm,a);
}

int vec_unit_wc(double out[], double a[]){
  double norm = vec_norm(a);
  vec_scap_wc(out,1/norm,a);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//dot product a and b

double vec_dotp(double a[], double b[]){
  double res=0;
  for(i=0;i<3;i++){
    res += a[i]*b[i];
  }
  return res;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//cross product a to b

int vec_cross(double a[], double b[]){
  double cross[3];
  vec_cross_wc(cross,a,b);
  for(i=0;i<3;i++){
    a[i] = cross[i];
  }
}

int vec_cross_wc(double out[], double a[], double b[]){
  out[0] = a[1]*b[2] - a[2]*b[1];
  out[1] = a[2]*b[0] - a[0]*b[2];
  out[2] = a[0]*b[1] - a[1]*b[0];
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//scalar product a and b

int vec_scap(double a, double b[]){
  for(i=0;i<3;i++){
    b[i] = a*b[i];
  }
}

int vec_scap_wc(double out[], double a, double b[]){
  for(i=0;i<3;i++){
    out[i] = a*b[i];
  }
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//orthogonal projection of a into b

int vec_proj(double a[], double b[]){
  double unit[3];
  double norm;

  vec_unit_wc(unit,b);
  norm = vec_dotp(a,unit);
  vec_scap_wc(a,norm,unit);

}

int vec_proj_wc(double out[], double a[], double b[]){
  double unit[3];
  double norm;

  vec_unit_wc(unit,b);
  norm = vec_dotp(a,unit);
  vec_scap_wc(out,norm,unit);

}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//orthogonal vector from b to a

int vec_ort(double a[], double b[]){
  double proj[3];
  vec_proj_wc(proj,a,b);
  vec_subt(a,proj);
}

int vec_ort_wc(double out[], double a[], double b[]){
  double proj[3];
  vec_proj_wc(proj,a,b);
  vec_subt_wc(out,a,proj);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//sum with coefficients x*a+y*b

int vec_csum(double a[], double b[], double x, double y){
  for(i=0;i<3;i++){
    a[i] = x*a[i] + y*b[i];
  }
}

int vec_csum_wc(double out[], double a[], double b[], double x, double y){
  for(i=0;i<3;i++){
    out[i] = x*a[i] + y*b[i];
  }
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//add a to b

int vec_add(double a[], double b[]){
  for(i=0;i<3;i++){
    a[i] += b[i];
  }
}

int vec_add_wc(double out[], double a[], double b[]){
  for(i=0;i<3;i++){
    out[i] = a[i]+b[i];
  }
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//subtract b from a 

int vec_subt(double a[], double b[]){
  for(i=0;i<3;i++){
    a[i] -= b[i];
  }
}

int vec_subt_wc(double out[], double a[], double b[]){
  for(i=0;i<3;i++){
    out[i] = a[i]-b[i];
  }
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//directional angle from b to a projected into axis plane  
//!!counter-clockwise [0,2*PI)

double vec_angle(double axis[], double a[], double b[]){
  double ort_a[3],ort_b[3],cross[3];
  double angle;
  vec_ort_wc(ort_a,a,axis);
  vec_ort_wc(ort_b,b,axis);

  vec_cross_wc(cross,ort_b,ort_a);

  if(vec_dotp(axis,cross) > 0){
    angle = acos(vec_dotp(ort_a,ort_b)/(vec_norm(ort_a)*vec_norm(ort_b)));
  } else {
    angle = 2*PI-acos(vec_dotp(ort_a,ort_b)/(vec_norm(ort_a)*vec_norm(ort_b)));
  }

  return angle;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//rotates a vector around an axis with a given angle (idea : CLHEP)

int vec_rotate(double axis[], double ang, double a[]){
  double out[3];
  vec_rotate_wc(out,axis,ang,a);
  for(i=0;i<3;i++) a[i] = out[i];
}

int vec_rotate_wc(double out[], double axis[], double ang, double a[]){
  double c = cos(ang);
  double s = sin(ang);
  double cmxy = axis[0]*axis[1]*(1-c);
  double cmxz = axis[0]*axis[2]*(1-c);
  double cmyz = axis[1]*axis[2]*(1-c);
  double cmxx = axis[0]*axis[0]*(1-c);
  double cmyy = axis[1]*axis[1]*(1-c);
  double cmzz = axis[2]*axis[2]*(1-c);
  double sx = axis[0]*s;
  double sy = axis[1]*s;
  double sz = axis[2]*s;

  out[0] = (c+cmxx)*a[0] + (cmxy-sz)*a[1] + (cmxz+sy)*a[2];
  out[1] = (cmxy+sz)*a[0] + (c+cmyy)*a[1] + (cmyz-sx)*a[2];
  out[2] = (cmxz-sy)*a[0] + (cmyz+sx)*a[1] + (c+cmzz)*a[2];
}

/*****************************************************************/
