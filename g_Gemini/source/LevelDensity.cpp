#include "CLevelDensity.h"
#include <QString>
#include <QFile>
#include <QTextStream>
//#include <QDebug>
#include <cmath>
#include <iostream>

CLevelDensity* CLevelDensity::fInstance = 0;

bool  CLevelDensity::normal = true;

double const CLevelDensity::pi = acos(-1.);
double CLevelDensity::k0 = 7.3;
double CLevelDensity::kInfinity = 12.;
double CLevelDensity::aKappa = 0.00517;
double CLevelDensity::cKappa = .0345;
double CLevelDensity::af_an = 1.036;
double CLevelDensity::aimf_an = 1.02;
double CLevelDensity::eFade = 18.52;
double CLevelDensity::jFade = 50.;
double CLevelDensity::Ucrit0 = 9.;
double CLevelDensity::Jcrit = 14;
//double  CLevelDensity::Ucrit0 = 0.;
//double  CLevelDensity::Jcrit = -1.;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/**
 * Constructor
 */
CLevelDensity::CLevelDensity()
{
    //constructor read in in level density parameter

    QString fName(":tbl/gemini.inp");
    QFile iff(fName);
    if(!iff.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        std::cout << "Gemini++: unable to open gemini.inp" << std::endl;
        return;
    }

    QTextStream ifFile(&iff);

    QString line;
    QString text1;
    QString text2;
    QString text3;
    QString text4;
    //getline(ifFile,line);
    line = ifFile.readLine();

    ifFile >> text1 >> k0 >> text2 >> aKappa >> text3 >> cKappa >> text4 >> kInfinity;

    ifFile >> text1 >> eFade >> text2 >> jFade >> text3;
    ifFile >> text1 >> Ucrit0 >> text2 >> Jcrit >> text3;

    iff.close();
    ifFile.flush();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CLevelDensity* CLevelDensity::instance() // mod-TU
{
    if (fInstance == 0) {
        fInstance = new CLevelDensity;
    }
    return fInstance;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*********************************************************************
/**
   * Returns the backshifted excitation energy energy to be used in the
   * Fermi gas formula for the level Density
     \param fU0 is the thermal excitation energy in MeV
     \param fPairing is the pairing correction in MeV
     \param fShell is the shell correction in MeV
     \param fJ is the angular momentum in units of hbar
  */
double CLevelDensity::getU(double fU0, double fPairing, double fShell, double fJ)
{

    fU = fU0;
    if (fU <= 0.)
    {
        fU = 0.;
        return fU;
    }
    //simple fade out of pairing
    double shiftP;
    double Ucrit = 0.;
    if (fJ < Jcrit) Ucrit = Ucrit0 * (1.-fJ/Jcrit) * (1.-fJ/Jcrit);
    if (fU > Ucrit) shiftP = fPairing;
    else shiftP = fPairing * (1. - (1.-fU/Ucrit) * (1.-fU/Ucrit));

    fU += shiftP;
    if (fU <= 0.)
    {
        fU = 0.;
        return fU;
    }
    //fU += fShell*(1.0-exp(-fU/eFade-fJ/jFade));
    double shiftS = fShell * tanh(fU/eFade + fJ/jFade);

    fU += shiftS;

    if (fU < 0.) fU = 0.;

    J = fJ;
    return fU;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*********************************************
/**
 * Returns the level-density parameter in units of MeV-1, the getU function
 * must already have been called to use this version
\param iA is the mass number
\param fJ is the angular momentum in units of hbar
\param iFission is a short indicating we are detail with a saddle-point shape
 */

double CLevelDensity::getLittleA(int iA, short iFission/*=0*/)
{
    //calculates the level density parameter
    //iA is nucleus mass number
    //fU is thermal excitation energy
    //fPairing is the pairing energy
    //fShell is the shell correction to the mass

    double fA = (double)iA;
    double kappa = 0.;
    double daden_dU;

    if ((normal && fU/fA < 3.) || aKappa == 0.)
    {
        if (k0 == kInfinity)
        {
            aden = fA/k0;
            daden_dU = 0.;
        }
        else
        {
            if (aKappa > 0.) kappa = aKappa*exp(cKappa*fA);
            //kappa = 1.5+.1143*J;
            double expTerm = exp(-kappa*fU/fA/(kInfinity-k0));
            aden = fA/(kInfinity - (kInfinity-k0)*expTerm);
            daden_dU = -((aden/fA)*(aden/fA)) * kappa * expTerm;
        }
        switch(iFission)
        {
        case 1:
            aden *= af_an;
            daden_dU *= af_an;
            break;
        case 2:
            aden *= aimf_an;
            daden_dU *= aimf_an;
            break;
        }
    }
    else
    {
        double r;
        if (iFission == 1) r = 1.0696;
        else if (iFission == 2) r = 1.05;
        else                    r = 1.0;

        if (aKappa > 0.) kappa = aKappa*exp(cKappa*fA);
        double fofr = ( kInfinity - ( kInfinity - k0 ) * r ) / (k0*r);
        double expTerm = exp(-fofr*kappa*fU/fA/(kInfinity-k0));
        aden = fA/(kInfinity - (kInfinity-k0)*r*expTerm);
        daden_dU = -((aden/fA)*(aden/fA)) * kappa * expTerm * fofr;
    }

    entropy = 2.*sqrt(aden*fU);
    temp = sqrt(fU/aden)/(1.+fU/aden*daden_dU);

    return aden;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*********************************************
/**
 * Returns the level-density parameter in units of MeV-1
\param iA is the mass number
\param fU0 is the thermal excitation energy in MeV
\param fPairing is the pairing correction in MeV
\param fShell is the shell correction in MeV
\param fJ is the angular momentum in units of hbar
\param iFission is a short indicating we are detail with a saddle-point shape
 */

double CLevelDensity::getLittleA(int iA, double fU0, double fPairing/*=0.*/,
                                double fShell/*=0.*/, double fJ/*=0.*/, short iFission/*=0*/)
{
    //calculates the level density parameter
    //iA is nucleus mass number
    //fU is thermal excitation energy
    //fPairing is the pairing energy
    //fShell is the shell correction to the mass

    //if (getU(fU0, fPairing,fShell,fJ) <= 0.) return 0.;
    getU(fU0, fPairing,fShell,fJ);
    return getLittleA(iA,iFission);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//************************************************
/**
 * Returns the natural logrithm of the spin dependent level density in MeV-1
 * if fJ < 0, then littleA is calculated taking into account the spin reduction
 * of paring, but subsequentally the level density is calculated for
 * zero spin. This is done when using Weisskopf formalism for evaporation.
 \param iA is the mass number
 \param fU0 is the thermal excitation energy in MeV
 \param fPairing is the pairing energy in MeV
 \param fShell is the shell correction in MeV
 \param fJ is the spin in hbar
 \param fMinertia is the moment of inertia in nucleon mass* fm^2
\param iFission indicates its for the saddle-point configuration
 */

//spin dependent Fermi_gas level density
double CLevelDensity::getLogLevelDensitySpherical
    (int iA, double fU0, double fPairing,
     double fShell, double fJ, double fMinertia, short iFission/*=0*/)
{
    //calculates the level density
    //iA is nucleus mass number
    //fU is thermal excitation energy
    //fPairing is the pairing energy
    //fShell is the shell correction to the mass

    if (getLittleA(iA,fU0,fPairing,fShell,fabs(fJ),iFission) == 0.) return 0.;
    if (fU <=0.) return 0.;
    if (fJ < 0.) fJ = 0.;
    double sigma =  fMinertia*temp/40.848;

    // pow(fU,1.25), pow(sigma,1.5), pow(aden,0.25) cannot be expanded safely
    double preExp = (2.*fJ+1.)/(1.+pow(fU,(double)(1.25))*pow(sigma,(double)1.5))
                   /pow(aden,(double)0.25)/24./sqrt(2.);
    return entropy + log(preExp);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//******************************************
/**
 * Returns the temperature in MeV
 * getLittleA must be called first
 */
double CLevelDensity::getTemp()
{
    return temp;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*************************************
/**
   * Returns the entropy.
   * getLittleA must be called first
   */
double CLevelDensity::getEntropy()
{
    return entropy;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*************************************
/**
* Returns the level-density parameter in MeV-1
 */
double CLevelDensity::getAden()
{
    return aden;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//************************************************
/**
   * Returns the spin independent Fermi-gas level density
   \param iA is the mass number
   \param fU0 is the thermal excitation energy in MeV
   \param fPairing is the pairing correction in MeV
   \param fShell is the shell correction in MeV
   */
double CLevelDensity::getLogLevelDensitySpherical
    (int iA, double fU0, double fPairing, double fShell)
{
    //calculates the level density
    //iA is nucleus mass number
    //fU is thermal excitation energy
    //fPairing is the pairing energy
    //fShell is the shell correction to the mass

    if (getLittleA(iA,fU0,fPairing,fShell)== 0.) return 0.;
    if (fU <=0.) return 0.;

    // pow(aden,0.25), pow(fU+temp,1.25) cannot be expanded safely
    double preExp = sqrt(pi)/12./pow(aden,(double)0.25)/
                   (1.+pow(fU+temp,(double)1.25));
    return entropy + log(preExp);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//***************************************************
/**
 * Sets the constants for  level-density parameter
 * where \f$a=\frac{A}{k_{\infty} - \left(k_{\infty} -k_{0} \right) \exp\left( \frac{\kappa}{k_{\infty}-k_{0}}\frac{U}{A}\right)}\f$
 * where \f$ \kappa = a_{\kappa} \exp\left(c_{\kappa} A\right) \f$.
 * Note setLittleA(8.) is equivalent to \f$a=A/8\f$
\param k00 is \f$k_{0}\f$
\param aKappa0 is \f$a_{\kappa}\f$
\param cKappa0 is \f$c_{kappa}\f$
\param kInfinity0 is \f$k_{\infty}\f$
*/
void CLevelDensity::setLittleA(double k00, double aKappa0/*=0.*/,
                               double cKappa0 /*=0.*/, double kInfinity0 /*=12.*/)
{
    k0 = k00;
    aKappa = aKappa0;
    cKappa = cKappa0;
    kInfinity = kInfinity0;
    if (aKappa == 0.) kInfinity = k0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//****************************************************
/**
 * returns the inverse level-density parameter at zero excitation energy
 */
double CLevelDensity::getK0()
{
    return k0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//************************************************************
/**
 * returns the inverse level-density parameter at infinite excitation
 * energy
 */
double CLevelDensity::getKInfinity()
{
    return kInfinity;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//****************************************************
/**
 * returns one of the coefficients used to calculate kappa
 * \f$ \kappa = a_{\kappa} \exp\left(c_{\kappa} A\right) \f$
 */
double CLevelDensity::getAKappa()
{
    return aKappa;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//****************************************************
/**
 * returns the other coefficients used to calculate kappa
 * \f$ \kappa = a_{\kappa} \exp\left(c_{\kappa} A\right) \f$

 */
double CLevelDensity::getCKappa()
{
    return cKappa;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*****************************************************
/**
 * returns the ration of saddle-point to equilibrium level-density parameter
 * for symmetyric fission
 */
double CLevelDensity::getAfAn()
{
    return af_an;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*****************************************************
/**
 * returns the ration of saddle-point to equilibrium level-density parameter
 * for asymmetric fission
 */
double CLevelDensity::getAimfAn()
{
    return aimf_an;
}
//*****************************************************
/**
  * set the ration of level-density paramters at the saddle-point to
  * equilibrium shape for symmetric fission
  \param af_an0 is the level-density parameter ratio for saddle-point to equilibirum deformation
  */
void CLevelDensity::setAfAn(double af_an0)
{
    af_an = af_an0;
}
//*****************************************************
/**
  * set the ration of level-density paramters at the saddle-point to
  * equilibrium shape for asymmetric fisison , ie. imf emission
  \param aimf_an0 is the level-density parameter ratio for saddle-point to equilibirum deformation
  */
void CLevelDensity::setAimfAn(double aimf_an0)
{
    aimf_an = aimf_an0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//***************************************************
/**
   * writes out the values of the parameters
   */
void CLevelDensity::printParameters()
{
    std::cout << "k0 = " << k0 << " kInfinity= " << kInfinity << std::endl;
    std::cout << "aKappa= " << aKappa << " cKappa= " << cKappa << std::endl;
    std::cout << "af/an= " << af_an << " for symmetric fission" << std::endl;
    std::cout << "aimf/an= " << aimf_an << " for asymmetric fission" << std::endl;
    std::cout << "eFade= " << eFade << " jFade= " << jFade << std::endl;
    std::cout << "Ucrit0= " << Ucrit0 << " Jcrit= " << Jcrit << std::endl;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//***************************************************
/**
   * returns the level density for a scission configuration.
   * This is used for evaporation from saddle to scission and
   * for determining the fission mass asymmetry
\param iA is mass number of fissioning system
\param fU is thermal excitation of fission system
\param adenInv in the inverse level-density paramter in MeV, is \f$a=A/adenInv\f$
  */
double CLevelDensity::getLogLevelDensityScission(int iA, double fU,
                                                double adenInv/*=8.*/ )
{
    aden = (double)iA/adenInv;
    temp = sqrt(fU/aden);
    entropy = 2.*sqrt(aden*fU);

    // pow(aden,0.25), pow(fU+temp,1.25) cannot be expanded safely
    double preExp = sqrt(pi)/12./pow(aden,(double)0.25)/
                   (1.+pow(fU+temp,(double)1.25));
    return entropy + log(preExp);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*****************************************************
/**
  * set the critical thermal excitation where pairing vanishes
  \param Ucrit0 is critical thermal excitaion energy in MeV
  \param Jcrit  is critical angular momentum where pairing vanishes
  */
void CLevelDensity::setUcrit(double Ucrit00, double Jcrit0)
{
    Ucrit0 = Ucrit00;
    Jcrit = Jcrit0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
