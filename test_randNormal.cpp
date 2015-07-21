#include <iostream>
#include <fstream>
#include "randNormal.h"

using namespace std;

int main()
{
  double a;
  fstream temp("temp.datout", ios::out);
  for(int i=0;i<10000;i++)
    {
      a = randNormal();
      temp << a << endl;
    }
}  
