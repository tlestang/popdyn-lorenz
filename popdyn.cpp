#include<iostream>
#include<sstream>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<ctime>
#include "lorenz.h"
using namespace std;

int main()
{  int Nc = 200; // Number of clones
  double T = 1000; // Total simulation time
  double dT = 10; // Cloning timestep
  int K = 100;

  double alpha = -5.0;
  int nbrTimeSteps;
  double **x0; // State of the clones after cloning (initial conditions for next block)
  double **x; // State of the clones after time evolution
  double s[Nc]; //Absolute weight
  int nbrCreatedCopies[Nc];
  double *R_record;
  double mean_err = 0.0; int alphaCount = 0;
  double phi_alpha, phi_theor;
  int NcPrime, deltaN, copyIdx, k;
  int l= 0;
  double **temp_storage;

  // R is the mean weight ; Q is the corresponding variance
  double R; 

  /*Allocate memory for coordinates*/
  x0 = new double*[Nc]; x = new double*[Nc];
  temp_storage = new double*[2*Nc];
  for (int i=0;i<Nc;i++)
    {
      x0[i] = new double[3]; x[i] = new double[3];
    }
  for (int i=0;i<Nc;i++)
    {
      temp_storage[i] = new double[3];
    }
	nbrTimeSteps = floor(T/dT); // Number of cloning steps
	R_record = new double[nbrTimeSteps];
	  
  for(int kreal=0;kreal<K;kreal++) // LOOP on realisations
    {

      cout << "K = " << kreal << "\r"; fflush(stdout);
      
  for(int dum=0;dum<Nc;dum++)
    {
      x0[dum][0]=1.0;
      x0[dum][1]=1.0;
      x0[dum][2]=1.0;
    }      

      l=0;
   
      for(int j=0;j<Nc;j++)
	{
	  lorenz_transient(x0[j], 5.0); //Transient regime
	}

    for(int t=0;t<nbrTimeSteps;t++)
    {
      
      if(t%(nbrTimeSteps/100)==0){l++; cout << l << "%\r"; fflush(stdout);}
      R = 0.0; // Mean weight initilization
      for(int j=0;j<Nc;j++) // Loop on clones
	{
	  s[j] =  lorenz(x0[j], x[j], alpha, dT); // OU simulates the time evolution of 1 
    	                                 //clone durinf dT and return                               
                                         //   the corresponding absolute weight
	R += s[j];
	} // END OF LOOP ON CLONES

      
      R /= Nc;
      NcPrime = 0.0;
      for(int j=0;j<Nc;j++)
	{
	  nbrCreatedCopies[j] = floor(s[j]/R + drand48()); // Nbr of copies that clone j should produce
	  NcPrime += nbrCreatedCopies[j]; //Total number of clones after cloning
	}

      // ------ REAJUSTING NUMBER OF CLONES ----------
      deltaN = NcPrime - Nc;

      	  k=0;
	  for(int j=0;j<Nc;j++)
	    {
	      for(int dum=0;dum<nbrCreatedCopies[j];dum++)
		{
		  temp_storage[k] = x[j];
      		  k++;
		}
	    }
      
      if(deltaN > 0) // if too many copies created we prune
	{
	  for(int i=0;i<deltaN;i++)
	    {
	      copyIdx = floor(NcPrime*drand48());
	      temp_storage[copyIdx] = temp_storage[NcPrime+i];
	      temp_storage[NcPrime+i][0] = 0.0;
	      	      temp_storage[NcPrime+i][1] = 0.0;
		      	      temp_storage[NcPrime+i][2] = 0.0;
		}
	} // END OF CASE DELTA_N > 0
      
      else if(deltaN < 0) // else we add just the right number of clones at random
	    {
	      for(int i=0;i<-deltaN;i++)
		{
		  copyIdx = floor(NcPrime*drand48());
		  temp_storage[NcPrime+i] = temp_storage[copyIdx];
		}
	    } // END OF CASE DELTA N < 0

      //---------------------------------------------------------------------------------

      // CREATE INITIAL CONDITIONS FOR NEXT ITERATION

            for(int j=0;j<Nc;j++)
	{
	  x0[j] = temp_storage[j];
	}
 
      R_record[t] = R; // Store mean weight values for SCGF calculation later on
     
    }// END OF LOOP ON TIMESTEPS
  
  // COMPUTE SCGF
  phi_alpha = 0;
  for(int n=0;n<nbrTimeSteps;n++){phi_alpha += log(R_record[n]);}
  phi_alpha /= T;
  
    } // END OF LOOP ON REALISATIONS K
  
  
    return 0;
}// END OF MAIN



