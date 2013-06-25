#include "stdio.h"

#include "common.h"

int main(){
  
  FILE* f;
  int i,j,k,ndx,st;
  double opening = op();
  int nofstrip;
  double output[200];
  double obs[2000]={0};
  double map[70000]={0};
  double fit[400]={0};
  double sources[100]={0};
  double count;
  double theta[20] = 
    {0.3,0.2,0.3,0.2,0.1,0.15,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double phi[20] = 
    {1.2,0.2,0.7,-0.3,1.5,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double nofphot[20] = 
    {1000,1000,1200,700,900,1500,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double noise = 10000.;
  double offset = 0.;
  double L = 50;
  double d = 5;
  int n_source = 6;
  int turn = 50;

  //system("clear");
  printf("\n\n");
  printf("DTU SPACE - RMC simulation and source detection simulation\n\n");
  printf(" INPUTS\n---------\n\n");
  printf("  **Configuration Settings\n\n");
  printf("   -Mask Height    : %5.2f\n",L);
  printf("   -Strip Length   : %5.2f\n",d);
  printf("   -Offset         : %5.2f\n",offset);
  printf("   -Rotation       : %d tour\n\n",turn);
  printf("  **Source Settings\n\n");
  printf("   - # of sources  : %d\n",n_source);
  for(i=0;i<n_source;i++){
    printf("    - %d --> Theta : %5.2f*PI    Phi : %5.2f*PI    Photon per Rotation : %5.2f\n"
	   ,i+1,theta[i],phi[i],nofphot[i]);
  }
  printf("   - Background noise : %5.2f \n",noise);
  printf("\n\nMonte-Carlo simulation is started...\n");

  count = real(obs,n_source,theta,phi,offset,L,d,nofphot,noise,turn);
   
  printf("%5.0f number of photons are created and replaced...\n",count);
  printf("Simulation is completed...\n\n\n");

  printf("Source analysis by removal of the strongest one...\n\n");

  for(k=0;;k++){
   
    printf("\nIteration %d\n\n",k+1);
 
    st = corr(map,obs,L,d,offset);
    n_source = loc_source(sources,map,obs,L);
    if(n_source == -1) {
      printf("More than 5 sources are found. Background dominates. Terminated...\n\n\n");
      break;
    }
    st = lsf(fit,obs,sources,n_source,L,d,offset);

    printf("%d sources are found above 0.8. Possible candidates are :\n",n_source);

    for(i=0;i<n_source;i++){
      printf("  - %d --> Theta : %5.2f*PI    Phi : %5.2f*PI    Intensity : %5.2f\n"
	     ,i+1,sources[2*i],sources[2*i+1],fit[i+1]*count/(opening*turn));
    }
    printf("  - Background Intensity : %5.2f\n\n",fit[0]*count/(opening*turn));
    
    count = clean(obs,fit,sources,n_source,count,offset,L,d);

    output[3*k] = sources[0];
    output[3*k+1] = sources[1];
    output[3*k+2] = sources[2]/(opening*turn);

    printf("Rest counts : %5.2f\n\n",count);
  }

  printf("  OUTPUT\n----------\n\n");
  for(i=0;i<k-1;i++){
    printf("  - %d --> Theta : %5.2f*PI    Phi : %5.2f*PI    Intensity : %5.2f\n"
	   ,i+1,output[3*i],output[3*i+1],output[3*i+2]);
  }

  printf("  - Background Intensity : %5.2f\n\n",count/(opening*turn));
  
  return 1;
}
