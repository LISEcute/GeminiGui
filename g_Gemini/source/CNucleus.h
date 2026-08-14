// -*- mode: c++ -*-
//
#ifndef nucleus_
#define nucleus_

#include "CMass.h"
#include "CYrast.h"
#include "CLevelDensity.h"
#include "CAngle.h"
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include "CRandom.h"
#include "CAngle.h"
#include "CEvap.h"
#include "CAngleDist.h"
#include "CScission.h"
#include "CWeight.h"
#include "SStoreEvap.h"
#include "CGdr.h"
#include <vector>

#include <QString>
using namespace std;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/**
 *!\brief storage
 *
 * this structure store information of the evaporated particle orbital AM
 */
struct SStoreSub
{
  double gamma; //!< weight factor for the orbital amgular momentum
  short unsigned L; //!< orbital angular momentum
};

typedef vector<SStoreSub> SStoreSubVector;
typedef vector<SStoreSub>::const_iterator SStoreSubIter;

/**
 *!\brief storage
 *
 * this structure info on complex fragment channels
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
struct SStore
{
  double gamma; //!< decay width of channel
  short unsigned iZ; //!< proton number of complex fragment
  short unsigned iA; //!< mass number of complex fragment
};

typedef vector<SStore> SStoreVector;
typedef vector<SStore>::const_iterator SStoreIter;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/**
 *!\brief storage
 *
 * structure stores information of evaporation decay sub channels
 */
struct SStoreChan
{
  double S2; //!< spin of daughter
  double Ek; //!< kinetic energy of evaporated particle (MeV)
  double Ex; //!< excitation energy of daughter (MeV)
  double gamma; //!< partial decay width for this subchannel (MeV)
  double temp; //!< temperature of daughter (MeV)
  double UPA; //!< thermal excitation energy per nucleus of daughter (MeV)
  short unsigned L; //!< orbital angular momentum of evaporated particle
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

/**
 *!\brief Hauser-Feshbach, Bohr-Wheeler, Morretto, formulisms
 *
 * Class CNucleus implements the GEMINI statistical mode code.
 * It follows the decay of a compound nucleus by a sequential series of 
 * binary decays. The decay widths are calculated with Hauser-Feshbash
 * formulism for light particles and Morreto's transition-state formulism
 * for other binary decays
 */

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class CNucleus : public CNuclide, public CWeight
{
public:

    CNucleus();
    CNucleus(int iZ,int iA);
    CNucleus(int iZ,int iA, double fEx, double fJ);
    ~CNucleus();
    static CNucleus* acquire(int iZ,int iA);
    static CNucleus* acquire(int iZ,int iA, double fEx, double fJ);
    static void release(CNucleus* nucleus);

    static CYrast *yrast; //!< gives fission barriers and rotational energies
    static CLevelDensity *levelDensity; //!< gives level densities

 protected:
  double Ecoul; //!< Coulomb barrier (HauserFeshbach)

  bool notStatistical; //!< this does not decay statistically, evap. frag. only

  short unsigned notStatisticalMode;//!< specifies type of nonStatisical decay
  double fPairing; //!< pairing energy
  double fShell; //!<shell correction
  double fU0; //!< thermal excitation energy
  double Erot; //!<yrast energy
  double Jmax; //!< max spin with a fission barrier
  double fMInertia; //!< spherical moment of inertia
  double logLevelDensity; //!< store the log of the level density of the nucleus
  double temp; //!< nuclear temperature
  int fissionZ; //!< proton number of fission fragment
  int fissionA; //!< mass number of fission fragment
  int fissioningZ; //!< proton number of fission parent
  int fissioningA; //!< mass number of fission parent
  int iZ1_IMF_Max; //!< maximum Z for IMF emission

  double fissionU; //!< thermal excitation energy of both fission fragments
  double EdefScission; //!< deformation energy of the scission configuration

  bool saddleToSciss; //!< indicated decay during saddle-to-scission transition
  double timeSinceSaddle; //!< stores the time since the saddle was crossed
  double timeSinceStart; //!< stores the time since the decay began

  void saddleToScission();
  void massAsymmetry(bool);
  bool needSymmetricFission; //!< indicated the Bohr-Wheeler width is needed
  static bool const noSymmetry;//!< true - old gemini with Morreto for all 
  double timeScission; //!< time required to go from saddle to scission
  static double const viscosity_scission; //!< viscosity during saddleTosciss
  static double const viscosity_saddle; //!< viscosity  during saddleTosciss
  static double timeTransient; //!< transient fission delay 
  static double fissionScaleFactor; //!< fission width scaled by this factor
  static double barAdd; //!< adds to Sierk fission barrier
  static unsigned iPoint; //!< pointer to array of stable fragments
  static int iHF; //!< set evaporation mode 
  int HF; //!< evaporation mode chosen for a given decay
  static bool noIMF; //!< no imf emission is considered
  static bool BohrWheeler; //!< no imf emission is considered
  double selectJ(double,double,double,double);

  static short unsigned Zshell; //!< enforce shell effects in evaporation
  CScission scission; //!< gives scission energeis, etc
  static CGdr * GDR ; //!< uder defined GDR line shape
  bool  bStable; //!< indicated this nucleus is particle-stable
  static double const r0; //!< radius const (fm)
  static double const sep; //!< separation between fragments
  static double threshold; //!< used to turn off unlikey evaporations
  void initializeDefaults();
  void initialize(int iZ0, int iA0);
  void initialize(int iZ0, int iA0, double fEx0, double fJ0);
  static vector<CNucleus*> recycled;

  CAngle spin; //!< orientation of the spin axis
  double velocity[3]; //!< velocity vector of nucleus in cm/ns 
  double momentum[3]; //!< momentum vector in MeV/c
  double qVal; //added by MPK 7/30/2015

  static CEvap *evap; //!< stores info on evaporated particles


  CLightP * lightP; //!< points to the light-particle decay mode
  double S2Loop(double Ekvalue);
  double S2Width(double Ekvalue);
  double EkWidth(double ek);
  void getSpin( bool saddle);

  double EkLoop();
  double getSumTl(double,double);
  double getWidthZA(double,short);
  void angleEvap();
  void angleIsotropic();
  void angleGamma();
  double S2Start; //!< Hauser-Feshback spin of daughter
  double UMin; //!< min thermal excitation energy in Hauser-Feshbach
  double EcostMin; //!<the min of the energetic cost of emitting light particles
  static short unsigned const lMaxQuantum; //!< number of l-waves to store angular dist
  static double de;//!< kinetic-energy interval for integrating in Hauser-Feshb
  int lMin; //!< minimum orbital AM for Hauser-Feshbach
  int lMax;//!< maximum orbital AM for Hauser-Feshbach
  double lPlusSMax; //!< max value of l+S of evaporated particle
  double lPlusSMin; //!< min value of l+S of evaporated particle
  double rResidue; //!< radius of daughter
  double rLightP; //!< radius of evaporated particle
  //double fMInertiaOrbit; //!< moment of Inertia for orbital motion
  double S2; //!< spin of daughter
  double EYrast2; //!< rotational energy of daughter
  SStoreEvap * storeEvap; //!< information of evap sub channels
  SStoreSub * storeSub; //!< store info on l distribution

  int iStore; //!< actual number of evap sub channels
  short unsigned EvapZ2; //!< proton number of daughter after evap.
  short unsigned EvapA2; //!< mass number of daughter after evap.
  short unsigned EvapZ1; //!< proton number of evaporated particle
  short unsigned EvapA1; //!< mass number of evaporated particle
  short unsigned EvapL; //!< orbital AM of evaporated particle
  short unsigned EvapMode; //!< ID number of evap channel 
  double EvapEx1; //!< excitation ennergy of evap. particle
  double EvapEx2; //!< excitation energy of daughter after evap.
  double EvapS2; //!< spin of daughter after evap
  double EvapS1; //!< spin of evaporated particle
  double EvapEk; //!< kinetic energy of evaporated particle (MeV)
  double EvapLPlusS; //!< toatl spin plus orbital AM of evaporated particle

  static CAngleDist angleDist; //!< selects angular distributions of decays

  double GammaEx;//!< excitation energy after gamma emission
  double GammaJ; //!< spin after gamma emission
  int GammaL; //!< gamma type E1=1, E2 = 2

  static double const gammaInhibition[3]; 
//!<scaling of gamma width from Weisskopf value
  static double const wue[3]; //!<coeff for Weisskopf units (gamma decay)
  void binaryDecay();
  void exciteScission(double,double,bool sym=1);
  double asyFissionWidth();
  double asyFissionWidthZA();
  double asyFissionWidthBW();
  void force8Be();
  void force5Li();
  void force5He();
  void force9B();

  double evaporationWidthSS();
  double gammaWidth();
  double gammaWidthE1GDR();
  double gammaWidthMultipole(int);
  double hauserFeshbach(int);
  double weiskopf( bool saddle);
  void  asyFissionDivide();
  void recursiveDecay();

  CNucleus*daughterLight;  //!< pointer to the lighter of the decay products
  CNucleus*daughterHeavy;  //!< pointer to the heavier of the decay products
  CNucleus*parent; //!< pointer to the parent nucleus

  static int const Nproducts; 
  //!< total number of possible  decay products from all decays


  static vector<CNucleus *> allProducts;
  //!< array of pointer to all decay products (stable or intermediate) 

  static vector<CNucleus *> stableProducts;
  //!< array of pointers to all stable decay products for all CN decays


  bool bResidue; //!< true if decay produced an evaporation residue
  bool bSymmetricFission; //!< true if decay resulted in symmetric fission
  bool bAsymmetricFission; //!< true if decay resulted in asymmetric fission
  int multPostLight; //!< number of post-fission neutrons for lighter ff
  int multPostHeavy; //!< number of post-fission neutrons for heavier ff
  int multPreSaddle; //!< number of pre-scission neutrons emitted
  int multSaddleToScission; 
  //!< number of neutrons emitted between saddle and scission


  static double const kRotate; //!< constant to calculated rotational energy
  void split(CAngle);

  double sigma2; //!< variance of fission mass distribution
  double symSaddlePoint;//!< symmetric saddle point energy
  static double sumGammaEnergy; //!< store the energy emitted in gamma rays

 static  vector <double> GammaRayEnergy; //!< store each gamma ray energy
  static int  nGammaRays;  //!< number of emitted gamma rays

  static bool  GDRParam; //!< if true, the standard formula for GDR decay width is used, if false the parametrized version


 public:


  bool abortEvent; //!< abort the event
  double evaporationWidth();
  double BohrWheelerWidth();
  double LestoneFissionWidth();
  double LestoneCorrection(double Usaddle, double momInertiaEff,short iAfAn);
 // static CRandom ran; //!< pointer to random number generator
  static double const pi; //!< 3.14159
  static double const EkFraction; // !< calculates the Ek spectra down to this
                                 // fraction of the maximum
  //functions
  double getSumGammaEnergy();
  int   getnGammaRays(); //get number of emitted gamma rays
  double getGammaRayEnergy(int number); // get gamma ray energy
  double getTime();
  QString printHtml();
  void setNewIsotope(int iZ0, int iA0, double fEx0, double fJ0); 

  void setCompoundNucleus(double fEx0,double fJ0);

  void setSpinAxis(CAngle angle);
  void setSpinAxisDegrees(CAngle angle);
  void setVelocityPolar(double =0.,double=0.,double=0.);
  void setVelocityCartesian(double vx=0.,double vy=0.,double vz=0.);

  void reset();
  static void resetGlobal();

  void print(const char *str=nullptr);
  void printStableProducts();
  void printAllProducts();
  void vCMofAllProducts();
  void energyConservation();


  CNucleus* getProducts(int=-1);
  CNucleus* getParent();
  CNucleus* getLightDaughter();
  CNucleus* getHeavyDaughter();
  CNucleus* getCompoundNucleus();


  int getNumberOfProducts();
  int getZmaxEvap();


  void excite(double,double);
  void excite(double);

  double getTheta();
  double getThetaDegrees();
  CAngle getAngle();
  CAngle getAngleDegrees();
  double getKE();
  double getVelocity();
  double getMomentum();
  double* getVelocityVector();
  double* getMomentumVector();

  static void setTimeTransient(double time);
  static void setFissionScaleFactor(double factor);
  static void setBarWidth(double width);
  static void setDeltaE(double de0);
  static void setThreshold(double threshold0);
  static void setAddToFisBarrier(double barAdd0);
  static void setNoIMF();
  static void setYesIMF();
  static void setLestone();
  static void setBohrWheeler();
  static void setSolution(int isol);
  static void setEvapMode(int iHF0=2);
static void setUserGDR(bool mode = true);

  static double getTimeTransient();
  static double getFissionScaleFactor();
  static double getBarWidth();
  static double getDeltaE();
  static double getThreshold();
  static double getAddToFisBarrier();

  void decay();
  bool isAsymmetricFission();
  bool isSymmetricFission();
  bool isNotStatistical();
  bool isSaddleToScission();
  bool isResidue();
  int getMultPost();
  int getMultPre();
  int getMultPostLight();
  int getMultPostHeavy();
  int getMultPreSaddle();
  int getMultSaddleToScission();
  double getFissionTimeSymmetric(double & timeScission);
  double getFissionTimeAsymmetric();
  double getDecayWidth();
  double getLogLevelDensity();
  int origin; //!< specifies the origin of the fragment, prefission, post , etc
  int origin2; //!< specifies the origin of the fragment, prefission, post , etc
  void printParameters();

  //*******ROOT********
    //ClassDef(CNucleus,1)  //Gemini Nucleus
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#endif 
