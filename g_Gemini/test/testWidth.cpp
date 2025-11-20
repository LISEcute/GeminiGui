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

  float fEx =30;
  float fJL =10;
  float Erot=CN.yrast->getYrast(iZ,iA,fJL);
  double fU0 = fEx - Erot;
  float r0=1.16;
  cout << "spin (test) for yrast = " << fJL << endl;

  float fPairing = CN.mass->getPairing(iZ,iA);
  float fShell = CN.mass->getShellCorrection(iZ,iA);
  float fMInertia =  0.4*pow(r0,2)*pow((float)iA,(float)(5./3.));
  float ld = CN.levelDensity->getLogLevelDensitySpherical(
            iA,fU0,fPairing,fShell,fJL,fMInertia);

/*  float men = CN.mass->getExpMass(0,1);
  float mcn = CN.mass->getCalMass(0,1);
  float mep = CN.mass->getExpMass(1,1);
  float mcp = CN.mass->getCalMass(1,1);

  cout << "masses " << men << " " << mcn << " " << mep << " " << mcp << endl;
*/
  cout << "yrast " << Erot << " fPairing " << fPairing << " fShell " << fShell
       << " fMInertia " << fMInertia << " levelDen " << ld << endl;

  float fJ = 0.;  //compound nucleus spin
  cout << "spin = " << fJ << endl;

  for (int i=0;i<10;i++)
    {
      float fEx = 14. + (float)i; //compound nucleus excitation energy
      CN.excite(fEx,fJ);


      float width = CN.getDecayWidth(); //decay width in MeV 
      float logRho = CN.getLogLevelDensity();
      float rho = exp(logRho); //level density in MeV-1

      cout << "Ex = " << fEx << " width = " << width << " MeV, rho = " <<
    rho << " MeV-1, width*rho = " << width*rho << endl;
    }

_useAME = save_ame;
return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
