#include <iostream>
#include <cmath>

#include "g_Gemini/source/CNucleus.h"


using namespace std;
extern bool _useAME;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int testTheWidth()
{
  bool save_ame = _useAME;
_useAME = false;

  cout << endl << "===================  test Width  ==================" << endl;

  int iZ = 12; 
  int iA = 24;

  CNucleus CN(iZ,iA);
  CN.setEvapMode(1);  //set Hauser feshbach calculation

  double fEx =30;
  double fJL =10;
  double Erot=CN.yrast->getYrast(iZ,iA,fJL);
  double fU0 = fEx - Erot;
  double r0=1.16;
  cout << "spin (test) for yrast = " << fJL << endl;

  double fPairing = CN.mass->getPairing(iZ,iA);
  double fShell = CN.mass->getShellCorrection(iZ,iA);
  double fMInertia =  0.4*(r0*r0)*pow((double)iA,(double)(5./3.));
  double ld = CN.levelDensity->getLogLevelDensitySpherical(
            iA,fU0,fPairing,fShell,fJL,fMInertia);

/*  double men = CN.mass->getExpMass(0,1);
  double mcn = CN.mass->getCalMass(0,1);
  double mep = CN.mass->getExpMass(1,1);
  double mcp = CN.mass->getCalMass(1,1);

  cout << "masses " << men << " " << mcn << " " << mep << " " << mcp << endl;
*/
  cout << "yrast " << Erot << " fPairing " << fPairing << " fShell " << fShell
       << " fMInertia " << fMInertia << " levelDen " << ld << endl;

  double fJ = 0.;  //compound nucleus spin
  cout << "spin = " << fJ << endl;

  for (int i=0;i<10;i++)
    {
      double fEx = 14. + (double)i; //compound nucleus excitation energy
      CN.excite(fEx,fJ);


      double width = CN.getDecayWidth(); //decay width in MeV 
      double logRho = CN.getLogLevelDensity();
      double rho = exp(logRho); //level density in MeV-1

      cout << "Ex = " << fEx << " width = " << width << " MeV, rho = " <<
    rho << " MeV-1, width*rho = " << width*rho << endl;
    }

_useAME = save_ame;
return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
