/*****************************************************************\
* Yigit Dallilar 02.07.2013                                       *
*                                                                 *
* C vector algebra library                                        *
\*****************************************************************/

#include "3d_cart_vec.h"

#include "math.h"


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
  vec_subt(a,b);
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

/*****************************************************************/
