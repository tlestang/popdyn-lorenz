#include <iostream>
#include <cmath>
#include <fstream>

int main()
{
  /*Lorenz model parameters*/
  double sigma = 10.0;
  double rho = 28.0;
  double beta = 8./3.;
  /*Simulation parameters*/
  double numberOfIterations;
  double timeStep;
  /*Dynamical qties*/
  double x, y, z;

  /* Open output file*/
  ofstream traj("trajectory_no_noise.datout");
  /*Initialize dynamical qties*/
  x = y = z = 0;
  
  /* START MAIN LOOP OVER TIMESTEPS */

  for(int t=0;t<numberOfIterations;t++)
    {
      /*Displays percentage*/
      if(t%(numberOfIterations/100)==0)
	{
	  k++; cout << k << "%\r"; fflush(stdout);
	}
      /*Compute lorenz model*/
      x = x + timeStep*sigma*(y-x);
      y = y + timeStep*(x*(rho-z) - y);
      z = z + timeStep*(x*y - beta*z);

      /*Write resutl on disk*/
      traj << t*timeStep << " " << x << " " << y << " " << z << endl;
    }

  traj.close();

}