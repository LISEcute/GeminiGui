#include "CMass.h"
#include "CNuclide.h"

//folowing is needed in ROOT version
//#include "Rtypes.h"

/**
 *!\brief fusion xsection, Bass model, excitation energy in fusion
 *
 * class to determine excitation energy and 
 *the critical angular momemtum in fusion.
 * The fusion xsection \f$ \sigma(\ell) = \pi\lambda^2 \sum \frac{2\ell +1}
{1+\exp\left(\frac{\ell-\ell_0}{\Delta_{\ell}}\right)}\f$
* where \f$\lambda\f$ is really lambdabar
* \f$ \Delta_{\ell}\f$ is the diffuseness
*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class Residual
{
public:
    Residual(){};

    int Z;
    int A;
    int count;
    int index() {return 100*Z+A;}
    string name;
};


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW



class CFus
{
public:

    CFus(double plb0,double dif0);
    CFus(int iZprojectile, int iAprojectile,
         int iZtarget, int iAtarget, double ELab, double dif0);
    void init(double plb0,double dif0);

CNuclide p,t,c; //Oleg

 protected:
  int iZp; //!< projectile proton number
  int iAp; //!<projectle mass number
  int iZt; //!<target proton number
  int iAt; //!<target mass number
  double fElab; //!< lab energy in MeV
  double R12; //!< sum of radii
  double U; //!< reduced mass
  double A; //!< const for Coulomb potential 
  double B; //!< const for centrifugal potential
  double C; //!< nuclear potential constant
  static double const D; //!< Bass potential parameter
  static double const E; //!< Bass potential parameter
  static double const G; //!< Bass potential parameter
  static double const H; //!< Bass potential parameter
  double E1; //!< critical energy 1 in Bass Model
  double E2; //!< critical energy 2 in Bass Model

  double MAX; //!< maximum L for fusion barrier

  double CL1; //!< angular momentum assocaited with E1
  double CL2; //!< angular momentum associated with E2
  double W[300]; //!< fusion barrier for each L

  double F(double R,double AL);
  double FF(double R,double AL);
  double  FFF(double R,double AL);
  double  FFFF(double R, double AL);


 public:

  double plb; //!< pi-lambdabar-squared in mb
  double dif; //!<diffuseness
  double Ecm; //!<reaction center of mass energy in MeV
  double vcm; //!<Compound nucleus velocity in cm/ns
  double vbeam; //!<beam velocity in cm/ns
  double Ex; //!< excitation energy
  int iZcn; //!< compound nucleus atomic number
  int iAcn;//!< compound nucleus mass number
  double getL0(double xsection);
  double getBassL(); 
  double getBassXsec(); 
  double Qval; //added by MPK 7/30/2015
  //following is needed in ROOT version
  //ClassDef(CFus,1); //Gemini CFus
};
