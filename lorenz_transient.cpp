#include <iostream>
#include <cmath>
#include <fstream>
#include "randNormal.h"

using namespace std;

void lorenz_transient(double *x0, double totalTime)
{
  /*Lorenz model parameters*/
  double sigma = 10.0;
  double rho = 28.0;
  double beta = 8./3.;
  double eps = 0.001; //Amplitude of the noise 
  /*Simulation parameters*/
  double timeStep = 0.001;
  int numberOfIterations = totalTime/timeStep;
  /*Dynamical qties*/
  double x, y, z;

  /*Misc*/
  int k=0; // Count for percentage
  int sdt = sqrt(timeStep);
  /* First itearation from initial conditions */
      x = x0[0] + timeStep*sigma*(x0[1]-x0[0]) + eps*sdt*randNormal();
      y = x0[1] + timeStep*(x0[0]*(rho-x0[2]) - x0[1]) + eps*sdt*randNormal();
      z = x0[2] + timeStep*(x0[0]*x0[1] - beta*x0[2]) + eps*sdt*randNormal();

  /* START MAIN LOOP OVER TIMESTEPS */
    for(int t=1;t<numberOfIterations;t++)
    {
      /*Compute lorenz model*/
      x = x + timeStep*sigma*(y-x) + eps*sdt*randNormal();
      y = y + timeStep*(x*(rho-z) - y) + eps*sdt*randNormal();
      z = z + timeStep*(x*y - beta*z) + eps*sdt*randNormal();
    }
    /*Update values of coordinates*/
    x0[0] = x;
    x0[1] = y;
    x0[3] = z;
}
