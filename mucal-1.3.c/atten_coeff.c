#include "stdio.h"
#include "mucal.h"

int main(){
  
  FILE* f = fopen("atten_coef.txt","w+");
  int i=0;
  double energy[15],xsec[15],fl_yield[15];
  char err[1000];
  double input[] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  
  for(i=0;i<sizeof(input)/sizeof(input[0]);i++){
    mucal("Si",0,input[i],'c',0,energy,xsec,fl_yield,err);
    fprintf(f,"%8.3f %8.3f \n",input[i],xsec[5]);
    printf("%8.3f %8.3f \n",input[i],xsec[5]);
  }
  fclose(f);

  return 0;

}
