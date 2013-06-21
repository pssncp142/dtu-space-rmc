#include "stdio.h"
#include "common.h"

#define PI 3.14159f

int main()
{

  int i,j;
  int nofstrip = n();
  double **tt = lsf(0.3, 1.2, 0., 50, 5, 5000., 9.);
  //FILE* f = fopen("corr.txt","w+");
  for(i=0; i<3; i++){
    //for(j=0; j<256; j++){
    printf("%f %f %f\n",tt[i][0],tt[i][1],tt[i][0]/tt[i][1]);
       }
  //fprintf("\n");
  //}
  //fclose(f);

  return 0;

}
