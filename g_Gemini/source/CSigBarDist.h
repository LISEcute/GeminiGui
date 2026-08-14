#include "CSigCharged.h"
#include <QString>
#include <cstdlib>
/**
 *!\brief inverse cross section with barrier distributions
 * 
 * calculates transmission coefficientinverse cross sections 
 * using a simplistic barrier 
 * distribution logic. The final result is the average of three  
 * inverse cross sections.
 *  \f$ T_{l}(Ek) = \frac{T_{l}^{R_{0}-\Delta R}(Ek) + T_{l}^{R_{0}}(Ek) + T_{l}^{R_{0}+\Delta R}(Ek)}{3} \f$ 
 * where \f$ T_{l}^{R_{0}}\f$ is the standard coeff derived using the IWBC and 
 * global optical-model potential. The other two are for when the radius of 
 * the nuclear potential is shifted by \f$ \pm \Delta R \f$. 
 * where \f$  \Delta R = width*\sqrt{temperature} \f$.
 */


class CSigBarDist
{
 private:
  CSigCharged* sigCharged[3]; //!< arrays for standard radii and +- width0
  bool one; //!< if true, no distribution, just standard radius is used
  static double width; //!< width paramter determines shifted radii
  static double const width0; //!< results readin from file for this shift
  double Z; //!< calls to getInverseXsec refer to this residual proton number
  double A; //!< calls to getInverseXsec refer to this residual mass number
  double Zp; //!<proton number of evaporated particle
  double Ap; //!<mass number of evaporated particle
 public:
  CSigBarDist(const QString&, double Zp0, double Ap0 );
  ~CSigBarDist();
  double getInverseXsec(double fEk, double temp);
  static void setBarWidth(double width00);
  static double getBarWidth();
  static void printParameters();
  void prepare(double Z, double A);
  double getBarrier();
};
