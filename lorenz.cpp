#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
  /*Lorenz model parameters*/
  double sigma = 10.0;
  double rho = 28.0;
  double beta = 8./3.;
  double eps = 0.01; //Amplitude of the noise 
  /*Simulation parameters*/
  int numberOfIterations = 100000;
  double timeStep = 0.001;
  /*Dynamical qties*/
  double x, y, z;

  /* Open output file*/
  ofstream traj("trajectory_no_noise.datout");
  /*Initialize dynamical qties*/
  x = y = z = 1;

  /*Misc*/
  int k=0; // Count for percentage
  int sdt = sqrt(timeStep);
  /* START MAIN LOOP OVER TIMESTEPS */

  for(int t=0;t<numberOfIterations;t++)
    {
      /*Displays percentage*/
      if(t%(numberOfIterations/100)==0)
	{
	  k++; cout << k << "%\r"; fflush(stdout);
	}
      /*Compute lorenz model*/
      x = x + timeStep*sigma*(y-x) + eps*sdt*randN();
      y = y + timeStep*(x*(rho-z) - y) + eps*sdt*randN();
      z = z + timeStep*(x*y - beta*z) + eps*sdt*randN();

      /*Write resutl on disk*/
      traj /*<< t*timeStep << " "*/ << x << " " << y << " " << z << endl;
    }

  traj.close();

}
