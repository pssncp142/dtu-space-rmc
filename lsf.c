#include "stdio.h"
#include "stdlib.h"

#include "common.h"

void control(int);
void analysis(char*,int);

int main(){
  
  int sp = n();
  analysis("25lsf.txt",sp);

}

void control(int sp){
  double data[60];
  int i,j,k,l,m,p,q,st;
  
  st = lsf(data,0.2,1.1,0.25,50,5,1000.,25.,150);
  char sp_name[] = "ABCDEF";
  for(i=0;i<sp;i++)
    printf("Strip %c : %f\n",sp_name[i],data[i]);
  printf("Average  : %f\n",data[sp]);
}

void analysis(char fname[],int sp){
  double data[60];
  int i,j,k,l,m,p,q,st;

  double phi[20] = {0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8};
  double theta[20] = {0.05,0.1,0.15,0.2,0.25,0.3};
  double offset[20] = {0.,0.1,0.2,0.3,0.4,0.5};
  int turn[20] = {1,5,20,50,150};
  double nofphot[20] = {500.,1000.,2000.,3000.,4000.};
  double noise[20] = {0.,5.,10.,15.,20.,25.};

  FILE* f = fopen(fname,"w+");

  for(i=0;i<10;i++){
    for(j=0;j<6;j++){
      for(k=0;k<6;k++){
	for(l=0;l<5;l++){
	  for(m=0;m<5;m++){
	    for(p=0;p<6;p++){
	      st = lsf(data,theta[j],phi[i],offset[k],50,5,nofphot[m],
		       noise[p],turn[l]);
	      fprintf(f,"%8.2f %8.1f %8.1f %10d %8.3f %8.3f %8.3f\n",
		      data[sp],noise[p],nofphot[m],turn[l],
		      offset[k],theta[j],phi[i]);
	    }
	  }
	}
      }
    }
  }
  fclose(f);

}
