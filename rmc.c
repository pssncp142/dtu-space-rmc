/******************************************************************
 * Yigit Dallilar 08.06.2013                                      *
 * DTU Space - Produces normalized count of A strips vs. azimuth  *
 * angle data.                                                    *
 ******************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int main()
{
  //constant pi
  const float pi = 3.1415926f;

  //initialisation parameters
  float phi, theta;            // source parameters.
  float r, strip_grid, T;      // open grid parameters as circle.
  float offset;                // offset value 
  float h;                     // distance between det and grid.
  float alp;                   // starting angle of the grid.
  float grid_length;           // grid length.

  phi = pi*0.3; theta = pi*0.1;
  r = 13.5; strip_grid = 1.; T = 1440.;
  offset = 0.1;
  h = 10.;
  alp = 0.;
  grid_length = 0.1;

  //Result data;
  long totcount = 0;
  long count[(int)T][2];
  for(int i=0; i<T; i++){
    count[i][0] = 0;
    count[i][1] = 0;
  }

  //some constants and some variable definition  
  float length= 2*r/grid_length;             //number of grid
  float length_2 = length*length*0.25;       //constant used in loop
  float r_2 = r*r;                           //constant used in loop
  float r_plus_offset = r+offset;            //constant used in loop
  float disp = h*tan(theta);                 //constant used in loop
  float dispx = -disp*cos(phi);              //constant used in loop
  float dispy = -disp*sin(phi);              //constant used in loop
  float coef = 0.5*grid_length*(1-length);   //constant used in loop
  float posx, posy;                          //x and y positions in upper grid
  float rel_posx;                            //relative x position in turned refrence frame
  float pi_2_over_T = 2*pi/T;                //constant used in loop
  float axis_angle[(int)T][2];               //constant used in loop
  for (int i = 0; i<T; i++){
    axis_angle[i][0] = cos(i*pi_2_over_T-alp);
    axis_angle[i][1] = sin(i*pi_2_over_T-alp);
  }
  float axis_angle_disp[(int)T][2];          //constant used in loop
  for (int i = 0; i<T; i++){
    axis_angle_disp[i][0] = dispx*cos(i*pi_2_over_T-alp);
    axis_angle_disp[i][1] = dispy*sin(i*pi_2_over_T-alp);
  }
  int nofstrip = 2*r/strip_grid;             //number of strips
  int grid[nofstrip];                        //grid opening function
  for (int i = 0; i<nofstrip; i++){
    grid[i] = i % 2;
  }
  int ndx;                                   //value definition

  //iteration in x grids
  for(int i=0 ; i<length ; i++){
    //iteration in y grids
    for(int j=0 ; j<length ; j++){
      //check if it is inside the upper circle
      if (pow(i-length*0.5,2) + pow(j-length*0.5,2) < length_2) {
	//counts total grid 
	totcount = totcount + 1;
	//calculates x and y position
	posx = i*grid_length + coef;
	posy = j*grid_length + coef;
	//turning upper grid
	for(int k=0 ; k<T ; k++){
	  //calculation of x value in new referance frame
	  rel_posx = posx*axis_angle[k][0] + posy*axis_angle[k][1];
	  //checks if passing the upper grid
	  ndx = (int)round(rel_posx+r);
	  if (grid[ndx] == 1){
	    //checks if inside the detector
	    if (pow(posx+dispx,2)+pow(posy+dispy,2) < r_2) {
	      ndx = (int)round(rel_posx+axis_angle_disp[k][0]+axis_angle_disp[k][1]+r_plus_offset);
	      //shadow of the grid choses which strip it is belonged to
	      if (grid[ndx] == 1) {
		count[k][0] = count[k][0] + 1;
	      } else {
	      	count[k][1] = count[k][1] + 1;
	      }
	    }
	  }
	}
      }
    }
  }
  
  //output data
  FILE *file;
  file = fopen("data.txt","w");
  for(int i=0; i<T; i++){
    fprintf(file,"%f %f\n",(float)count[i][0]/(count[i][0]+count[i][1]),i*2/T-alp/pi);
  }
  fclose(file);

  //run python plot script
  system("./plot.py");

  return 0;		 
}

/*********************************************************************/
