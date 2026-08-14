#ifndef _levelDensity
#define _levelDensity

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

/**
 *!\brief returns level density, entropy, temperature
 *
 * This class deals with level densities and related paramters.
 * the level-density parameter is \f$a=\frac{A}{k_{\infty} - \left(k_{\infty} -k_{0} \right) \exp\left( \frac{\kappa}{k_{\infty}-k_{0}}\frac{U}{A}\right)}\f$
 * where \f$ \kappa = a_{\kappa} \exp\left(c_{\kappa} A\right) \f$
 */

class CLevelDensity
{
 private:
  CLevelDensity();
  static CLevelDensity *fInstance; //!< instance member to make this class a singleton
  static double k0; //!< inverse level-density parameter at U=0
  static double kInfinity; //!< inverse level-density parameter at U=infinity
  static double aKappa; //!< mass dependence of kappa
  static double cKappa; //!< mass dependence of kappa
  static double af_an; //!< ratio of sym-fission saddle  to equilibrium level-density para
  static double aimf_an; //!< ratio of asy-fission(imf) saddle to equilibrium level-density para
  static bool normal;
  static double Ucrit0; //!< Ucrit for J= 0, vanishing of pairing
  static double Jcrit; //!< spin for the vanishing of pairing
   double aden; //!< little-density parameter
  double entropy; //!< entropy
  double temp;  //!< temperature
  static double eFade; //!< fade out of shell effects with excitation energy
  static double jFade; //!< fade out of shell effects with spin
  double fU;   //!< thermal excitation energy in MeV
  static double const pi; //!< the mathematical constant \f$\pi\f$

 public:
  static CLevelDensity *instance(); //!< instance member to make this class a singleton
  double getLittleA(int iA,double fU0,double fPairing=0.,double fShell=0.,
                 double fJ=0.,short iFission=0);
  double getLittleA(int iA, short iFission);
  double getU(double fU0,double fPairing, double fShell, double fJ);
  double getAden();
  double getLogLevelDensitySpherical(int iA,double fU0,double fPairing,
	       double fShell,double fJ,double fMinertia,
               short iFission=0);
  double getLogLevelDensitySpherical(int iA,double fU0,double fPairing,
         double fShell);
  double getTemp();
  double getEntropy();
  static void setLittleA(double k00,double aKappa0=0., double cKappa0=0.,
          double kInfinity=12.);
  static void setAfAn(double af_an0);
  static void setAimfAn(double aimf_an0);
  static void setUcrit(double Ucrit00, double Jcrit);
  static double getAKappa();
  static double getCKappa();
  static double getK0();
  static double getKInfinity();
  static double getAfAn();
  static double getAimfAn();
  static void printParameters();
  double getLogLevelDensityScission(int iA, double U, double adenInv=8.);
  double J;
};
#endif
