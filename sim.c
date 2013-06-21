#include "stdio.h"
#include "common.h"

#define PI 3.14159f

int main()
{

  printf("holal");
  int i,j;
  int nofstrip = n();
  double ***model = corr(0.2, 1.2, 0.5, 50, 5, 40000.,0.);
  FILE* f = fopen("corr.txt","w+");
  for(i=0; i<256; i++){
    for(j=0; j<256; j++){
      fprintf(f,"%f \n",model[2][i][j]);
      }
      fprintf(f,"\n");
  }
  fclose(f);

  return 0;

}
