#include <cmath>
#include "CMass.h"
/**
 *!\brief scission energies
 *
 * calculates information on scission point for saddle-to-scission
 * evaporation and for fission mass distributions
 */


class CScission
{
 protected:
  static double const e2;  //!< \f$e^{2} = 1.44 MeV fm^{-1}\f$
  static double const kRotate; //!< constant for rotational energy 
  static double const slopeViola; //!< for Viola systematics of fission KE
  static double const constViola; //!< for Viola systematics of fission KE

  double mu0; //!< reduced mass
  double Vc0; //!< Coulomb self energy of a symmetric frag
  double R1; //!< radius of  frag1 in fm
  double momInertia1; //!< moment of inertia of  frgament1
  double R2; //!< radius of  frag1 in fm
  double momInertia2; //!< moment of inertia of  frgament1
  double k1; //!< part of the rotational energy
  bool sym; //!< logical symmetic of non symmetric fission
  double Z1;//!< atomic number of lighter fragment after fission
  double Z2;//!< atomic number of heavier fragment
  double A1;//!< mass number of lighter fragment after fission
  double A2;//!< mass number of heavier fragment 

  CMass * mass; //!< gives mass excess 

 public:
  static double const r0; //!< constant R=r0*A^(1/3)
  double sep; //!< separation between the surafces of the 2 spheres
  double sep0; //!< separation determined in init()
  double sep1; //!< separation determined in sigmaFissionSystematics()
  int iA; //!< mass number of system
  int iZ; //!< proton number of system

  double A; //!< double value of iA
  double Z; //!< double value of iZ
  double fJ; //!< spin of system
  double Esymmetric; //!< energy for symmetric mass split
  double ekTot; //!< total fission kinetic energy from getFissionKineticEnergy()
  double Erotate1; //!<rotational energy of fragment1
  double Erotate2; //!<rotational energy of fragment2
  double Epair1; //!< pairing energy of frag1
  double Epair2; //!< pairing energy of frag2
  double Eshell1; //!< shell energy of frag1
  double Eshell2; //!< shell energy of frag2
  double EkCoul; //!< Colomb part of total fisison kinetic energy
  double EkRot; //!< rotational part of total fission kinetic energy

  CScission(int iZ0,int iA0,double fJ0, int iChan);
  CScission();
  void init(int iZ0,int iA0,double fJ0,int ichan,double Z1=0.,double A1=0.);
  double getSep(double EkViola);
  double getScissionEnergy();
  double getScissionEnergy(int iZ1,int iA1);
  double getFissionKineticEnergy(int iZ1, int iA1);
  double sigmaFissionSystematicsScission(int iZ, int iA, double fJ, double fUscission);
  double sigmaFissionSystematicsSaddle(int iZ, int iA, double fJ, double fUscission);
};
