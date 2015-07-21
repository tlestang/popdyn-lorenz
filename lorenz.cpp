#include <iostream>
#include <cmath>
#include <fstream>
#include "randNormal.h"

using namespace std;

double lorenz(double *x0, double *m_x, double alpha, double totalTime)
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
  double s=0.0;
  
  /* First itearation from initial conditions */
      x = x0[0] + timeStep*sigma*(x0[1]-x0[0]) + eps*sdt*randNormal();
      y = x0[1] + timeStep*(x0[0]*(rho-x0[2]) - x0[1]) + eps*sdt*randNormal();
      z = x0[2] + timeStep*(x0[0]*x0[1] - beta*x0[2]) + eps*sdt*randNormal();

  /* START MAIN LOOP OVER TIMESTEPS */
    for(int t=1;t<numberOfIterations;t++)
    {
      /*Displays percentage*/
      /*      if(t%(numberOfIterations/100)==0)
	{
	  k++; cout << k << "%\r"; fflush(stdout);
	  }*/
      /*Compute lorenz model*/
      x = x + timeStep*sigma*(y-x) + eps*sdt*randNormal();
      y = y + timeStep*(x*(rho-z) - y) + eps*sdt*randNormal();
      z = z + timeStep*(x*y - beta*z) + eps*sdt*randNormal();

      /*COmpute integral for weight*/
      s += x;
    }
    /*Update values of coordinates*/
    m_x[0] = x;
    m_x[1] = y;
    m_x[3] = z;
    /*Compute weight and return it*/
    s *= timeStep;
    return exp(alpha*s);
}
