#include "stdio.h"

#define PI 3.14159f

#include "common.h"
#include "param.h"

int progress(int a, int b, int c);

int main(){

  configure(1);

  FILE* f;
  int i,j,k,l,m,n,o;
  int range = 4;
  char f_corr[] = "cor/Ex1-2-0.bin";
  char f_out[] = "res/Ex1-2-0-000.bin";
  double map[70000] = {0};
  double banned[70000] = {0};
  double obs[2000]={0};
  double init_obs[2000]={0};
  double fit[200]={0};
  double sources[100]={0};
  double theta[2] = {0.15,0.15};
  double phi[2] = {1.33,0.67};
  double nofphot[2] = {200,40};
  double noise = 800;
  double off[3] = {1e-20*PI,0.005,-0.005};
  double tcnt,cnt;
  double max;
  int found;
  int n_source = 2;
  int turn = 80;

  for(o=0;o<3;o++){
    beta = off[o];
    f_corr[10] = 48 + o;
    f_out[10] = 48 + o;
    configure(0);
    for(i=0;i<100;i++){
      f_out[12] = ((int)(turn+i)/100) + 48;
      f_out[13] = ((int)((turn+i)%100)/10) + 48;
      f_out[14] = ((int)((turn+i)%10)) + 48;
      f = fopen(f_out,"w+");
      for(l=0;l<2000;l++){
	progress(o,i,l);
	tcnt = real(init_obs,n_source,theta,phi,nofphot,noise,turn+i);
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
	  //printf("%f %f \n ",sources[2*j],sources[2*j+1]);
	  //printf("%f\n",fit[j+1]);
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
  double perc = (a*2000.*100+b*2000+c)*100/(2000*100*3-1);

  printf("\r[");
  
  for(i=0;i<(int)perc;i++) printf(":");
  for(i=0;i<100-(int)perc;i++) printf(" ");
  printf("] %6.3f %%",perc);
  fflush(stdout);

  return 0;
}
