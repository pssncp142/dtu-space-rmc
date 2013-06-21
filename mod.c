#include "stdio.h"
#include "stdlib.h"

#include "common.h"

int main(){

  int sp = n();
  double **data = mod(0.2,1.2,0.5,50,5);
  wspfile(data,sp);
  system("./plot.py");

}
