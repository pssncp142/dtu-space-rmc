#include "stdio.h"
#include "math.h"
#include "common.h"


double **lsf(double theta,double phi,double offset,double L,double d,double nofphot,double noise)
{

  int i,j,st;
  double data[nofstrip+1][2]; data[nofstrip][0]=0; data[nofstrip][1]=0;
  double LHS[100]; double RHS[100];
  double **obs=real(0.3,1.2,0.,50.,5.,2000.,10.); 
  double **model=mod(0.3,1.2,0.,50.,5.);

  for(j=0; j<nofstrip; j++){

    st = norm(model[j]);
    st = norm(obs[j]);
  
    LHS[0] = 1./256; LHS[1]=0; LHS[2]=0; LHS[3]=0; RHS[0] = 0; RHS[1] = 0; 

    for(i=0; i<256; i++){
      LHS[1] += model[j][i]/256;
      LHS[2] += model[j][i]/256;
      LHS[3] += pow(model[j][i],2);
      RHS[0] += obs[j][i]/256;
      RHS[1] += obs[j][i]*model[j][i];
    }

    st =  gaussj_nr(LHS,2,RHS,2);

    data[j][0] = RHS[0];
    data[j][1] = RHS[1];
  }

  for(i=0;i<nofstrip;i++){
    data[nofstrip][0] += data[i][0]/nofstrip;
    data[nofstrip][1] += data[i][1]/nofstrip;
  }
  
  return data;
}
