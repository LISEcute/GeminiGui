#ifndef random_
#define random_
#include <cstdlib>
#include <cmath>

/**
 * !\brief Random numbers for a number of distributions
 *
 * Random number generation using the C++ random number function
 */


class CRandom
{
 protected:
  CRandom();
  static  CRandom* fInstance; //!< instance member to make tis a singleton
  bool one; //!< used for Gaus
  double angle; //!< used for Gaus
  double x; //!< parameter
  double  pi; //!< 3.14159
 public:
  static CRandom* instance(); //!< instance member to make this a singleton
  ~CRandom();
  double Rndm();
  double Gaus(double mean,double sigma);
  double expDecayTime(double width);
  double BreitWigner(double mean ,double width);
};


#endif
