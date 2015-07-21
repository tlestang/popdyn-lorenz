#include <cstdlib>
#include <cmath>

double randNormal()
{
  /* Return a random number sampled in N(0.0, 1.0).
     Box-Muller method.
  */

  double x1, x2, w;
  do {
    x1 = 2.0 * (rand () / double (RAND_MAX)) - 1.0;
    x2 = 2.0 * (rand () / double (RAND_MAX)) - 1.0;
    w = x1 * x1 + x2 * x2;
  } while (w >= 1.0);

  w = sqrt (-2.0 * log (w)/w);
  const double y1 = x1 * w;
  const double y2 = x2 * w;

  return y1;
}
