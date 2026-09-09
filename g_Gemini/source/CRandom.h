#ifndef random_
#define random_
#include <cstdlib>
#include <cmath>
#include <random>

/**
 * !\brief Random numbers for a number of distributions
 *
 * Random number generation using the C++ random number function
 */


class CRandom
{
 protected:
  CRandom();
  static thread_local CRandom* fInstance; //!< per-thread random-number generator
  bool one; //!< used for Gaus
  double angle; //!< used for Gaus
  double x; //!< parameter
  double  pi; //!< 3.14159
  std::mt19937_64 generator;
  std::uniform_real_distribution<double> uniform;
 public:
  static CRandom* instance(); //!< instance member to make this a singleton
  static void seedCurrentThread(unsigned int seed);
  ~CRandom();
  double Rndm();
  double Gaus(double mean,double sigma);
  double expDecayTime(double width);
  double BreitWigner(double mean ,double width);
};


#endif
