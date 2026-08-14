#ifndef sigCharged_
#define sigCharged_

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <QString>
using namespace std;

/**
 *!\brief inverse xsections
 *
 * This class calculates \f$ \sum_{0}^{\infty} (2\ell+1) T_{\ell}(\varepsilon)\f$ which is the inverse xsection divided by \f$ \pi/k^2 \f$
 * 
 */


class CSigCharged
{
 private:
  QString sName; //!< name of input file of coefficients
  double rc0; //!< paramter to calculate radius for Coulomb barrier
  double rc1; //!< paramter to calculate radius for Coulomb barrier
  double rc2; //!< paramter to calculate radius for Coulomb barrier
  double omega0;  //!< paramter for omega
  double omega1;  //!< parameter for omega
  double omega2;  //!< parameter for omega
  double omega3;  //!< parameter for omega
  double rI0; //!< paramter for radius for rotational energy
  double rI1; //!< paramter for radius for rotational energy
  double rI2; //!< paramter for radius for rotational energy
  double aa0; //!< below barrier correction parameter
  double aa1; //!< below barrier correction parameter
  double a0; //!< above barrier correction parameter
  double a1; //!< above barrier correction parameter
  double Zp;  //!< proton number of evaporated particle
  double Ap;  //!< mass number of evaporated particle

  double barrier; //!< Coulomb barrier in MeV
  double InvInertia; //!< inverse of the moment of inertia associated with rotateion
  double omega; //!< amega parameter in MeV
  double a; //!< above barrier correction
  double aa; //!< below barrier correction
  double offset;  //!< offset

  bool neutron; //!<bool to signify neutron calculation
  double n0; //!< neutron parameter
  double n1; //!< neutron parameter
  double n2; //!< neutron parameter

 public:
  CSigCharged(const QString& file, double Zp0, double Ap0); //!< constructor
  void prepare(double Z,double A); //!< prepares for calculations of inverse xsections
  double getInverseXsec(double energy); //!<calculates the inverse xsection 
  double getBarrier(); //!< calculates the barrier

 };
#endif
