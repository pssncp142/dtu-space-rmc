#include "stdio.h"

#define PI 3.14159f

#include "common.h"
#include "param.h"

int x[2],y[2],z[2];
int progress(int a, int b, int c);

int main(){

  configure(1);

  FILE* f;
  int i,j,k,l,m,n,o;
  int range = 4;
  char f_corr[] = "cor/Ex1-2-0.bin";
  char f_out[] = "res/Ex2-2-0-000.bin";
  double map[70000] = {0};
  double banned[70000] = {0};
  double obs[2000]={0};
  double init_obs[2000]={0};
  double fit[200]={0};
  double sources[100]={0};
  double theta[3] = {0.15,0.15,0.11};
  double phi[3] = {1.33,0.67,0};
  double nofphot[3] = {200,40,10};
  double noise = 800;
  double off[3] = {1e-20*PI,0.005,-0.005};
  double tcnt,cnt;
  double max;
  int found;
  double keep = nofphot[2];
  int n_source = 2;
  int turn[] = {1,10,100,4000};
  x[0]=0; x[1]=20;
  y[0]=0; y[1]=4;
  z[0]=0; z[1]=2000;

  for(o=x[0];o<x[1];o++){
    nofphot[2] = keep + 4*o;
    f_out[12] = ((int)((int)nofphot[2])/100) + 48;
    f_out[13] = ((int)(((int)nofphot[2])%100)/10) + 48;
    f_out[14] = ((int)(((int)nofphot[2])%10)) + 48;
    configure(0);
    for(i=y[0];i<y[1];i++){
      f_out[10] = 48 + i;
      f = fopen(f_out,"w+");
      for(l=z[0];l<z[1];l++){
	progress(o,i,l);
	tcnt = real(init_obs,n_source,theta,phi,nofphot,noise,turn[i]);
	for(j=0;j<2000;j++) obs[j] = init_obs[j];
	for(j=0;j<4;j++){
	  corr(map,obs,1,f_corr);
	  loc_source(sources,map,banned,j);
	  lsf(fit,init_obs,sources,j+1);
	  for(m=-range;m<range+1;m++){
	    for(n=-range;n<range+1;n++){
	      banned[((int)sources[2*(j)+7]+n)*256+((int)sources[2*(j)+6]+m)]=1;
	    }
	  }
	  if(fit[j+1]<0.025){
	    break;
	  }
	  for(k=0;k<2000;k++) obs[k] = init_obs[k];
	  clean(obs,fit,sources,j+1,tcnt);
	  }
	fwrite(&j,sizeof(int),1,f);
	fwrite(sources,sizeof(double),4,f);
	for(k=0;k<70000;k++) banned[k]=0;
	
	}
      fclose(f);
    }
  }

  printf("\n");

  return 0;
}

int progress(int a, int b, int c){
  
  int i;
  double perc = ((double)a*(z[1]-z[0])*(y[1]-y[0])+b*(z[1]-z[0])+c)*100
    /((x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0])-1);

  printf("\r[");
  
  for(i=0;i<(int)perc;i++) printf(":");
  for(i=0;i<100-(int)perc;i++) printf(" ");
  printf("] %6.3f %%",perc);
  fflush(stdout);

  return 0;
}
