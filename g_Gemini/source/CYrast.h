#ifndef yrast_
#define yrast_
#include <iostream>
#include <fstream>
#include "CMass.h"
using namespace std;

/**
 *!\brief fission barriers, rotational energies, etc
 *
 * Class to return barriers and rotational energies. It is constructed 
 * from bits of code such as Arnie Sierk Barfit - which gives his Finite-range
 *fission barrier, rotational energies, moments of inertia, surface areas.
 * Also I have used code from the RLDM model. In addition I enclude a 
 * interpolation of conditional fission barriers, using barriers calculated
 * by Sierk
 */

class CYrast
{
 private:

  CYrast();
  static CYrast *fInstance; //!< instance member to make this class a singleton
  static double const pi; //!< 3.14159
  //needed by getYrastRLDM
  static double const x1h[11][6]; //!< number for RLDM
  static double const x2h[11][6]; //!< numbers for RLDM
  static double const x3h[20][10]; //!< number for RLDM
  static double const x1b[11][6]; //!<numbers for RLDM
  static double const x2b[11][6]; //!<numbers for RLDM
  static double const x3b[20][10]; //!<numbers for RLDM
  //needed by Sierk functions
  static double const emncof[4][5]; //!< used in Sierk functions
  static double const elmcof[4][5]; //!<used in Sierk functions
  static double const emxcof[5][7];//!<used in Sierk functions
  static double const elzcof[7][7];//!<used in Sierk functions
  static double const egscof[5][7][5];//!<used in Sierk functions
  static double const aizroc[5][6];//!<used in Sierk functions
  static double const ai70c[5][6];//!<used in Sierk functions
  static double const ai95c[5][6];//!<used in Sierk functions
  static double const aimaxc[5][6];//!<used in Sierk functions
  static double const ai952c[5][6];//!<used in Sierk functions
  static double const aimax2c[5][6];//!<used in Sierk functions
  static double const aimax3c[4][4];//!<used in Sierk functions
  static double const aimax4c[4][4];//!<used in Sierk functions
  static double const bizroc[4][6];//!<used in Sierk functions
  static double const bi70c[4][6];//!<used in Sierk functions
  static double const bi95c[4][6];//!<used in Sierk functions
  static double const bimaxc[4][6];//!<used in Sierk functions
  static double const b[8][5][5];//!<used in Sierk functions
  void lpoly(double,int,double*);

  double A; //!< mass number
  double Z; //!< proton number
  double zz; //!< used in Sierk functions
  double amin; //!< lower limits of application of Sierk routine
  double amax; //!< upper limits of application of Sierk routine
  double pa[7]; //!<used in Sierk routines
  double pz[7];  //!<used in Sierk routines
  //needed by saddlefit
  double c[6][8][2][11][2]; //!< coeff for sadfits
  double cubic(double,double,double,double,double,double);

  static bool first; //!< only write out barrier warning once
  int Narray; //!< number of elements in array of asymmetric barriers
  static double const hbarc; //!< used for asymmetric barriers
  static double const alfinv;//!< used for asymmetric barriers
  static double const srznw;//!< used for asymmetric barriers
  static double const aknw;//!< used for asymmetric barriers
  static double const bb;//!< used for asymmetric barriers
  static double const um;//!< used for asymmetric barriers
  static double const elm;//!< used for asymmetric barriers
  static double const spdlt;//!< used for asymmetric barriers
  static double const asnw;//!< used for asymmetric barriers
  static double const kx[8];//!< used for asymmetric barriers
  static double const ky[6];//!< used for asymmetric barriers
  static double const ka[11];//!< used for asymmetric barriers
  static double const r0; //!< radius parameter
  static double const sep; //!<separation id fm between fragments
  static bool bForceSierk; //!<separation id fm between fragments
  static double addBar; //!<extrapolated Sierk barrier increase by this amount

  double sadArray[300]; //!< array stores the conditional saddle energies
  double sadArrayZA[300]; //!< array stores saddle energies after correction
  CMass * mass; //!< class for mass defects
  static double const deltaJ; //!< used to extend sierk barrier to higher J
  static double const kRotate; //!< constant for rotional energy
  int iZ; //!< proton number
  int iA; //!<mass number
  double fJ; //!< spin
 public:
  static CYrast *instance(); //!< instance member to make this class a singleton
  double Jmax;  //!< max spin where the fission barrier exists
  double getYrast(int,int,double);
  double getYrastModel(int,int,double);
  double getYrastRLDM(int,int,double);
  double getYrastSierk(double);
  double getJmaxSierk(int,int);
  double getBarrierFissionSierk(double);
  double getSymmetricSaddleEnergy(int,int,double);
  double getBarrierFissionRLDM(int,int,double);
  double getBsSierk(double);
  static void forceSierk(bool=1);
  static void printParameters();
  void prepareAsyBarrier(int, int, double);
  void printAsyBarrier();
  double getSaddlePointEnergy(int,int);
  double getSaddlePointEnergy(double);
  double getMomentOfInertiaSierk(double);
  double WignerEnergy(int iZ, int iA);

  double momInertiaMin; //!< minimum saddle-point moment of inertia 
  double momInertiaMid; //!< intermediate saddle-point moment of inertia
  double momInertiaMax; //!< maximum saddle-point moment of inertia
  // double sumPair;
  //double sumShell;
  //double viola(double,double,double,double);
};
#endif
