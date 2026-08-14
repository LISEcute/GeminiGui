#ifndef lightP_
#define lightP_

#include "CNuclide.h"
#include "CTlBarDist.h"
#include "CSigBarDist.h"
#include "SStoreEvap.h"
#include <vector>

/**
 *!\brief deals with a particular evaporation channel
 *
 * Class CLightP contains all the important information
 * concerning a light particle evaporation channel, it has
 * pointers to a cTlArray class to calculate transmission coeff. 
 */

class CLightP : public CNuclide
{
 protected:
  bool mustDelete; //!< bool to remember that we must delete the tLArray
 public:

  CLightP(int iZ,int iA,double fS,string sName);
  CLightP(int,int,double,CTlBarDist*,CSigBarDist*);

   ~CLightP(); //!< destructor

   CTlBarDist *tlArray; //!< gives transmission coefficients
   CSigBarDist *sigBarDist; //!< gives inverse xsections for charged part
   
  double suppress; //!< suppression factor to scale Hauser-Feshbach width
  double fMInertia; //!< moment of inertia of residue
  double fMInertiaOrbit; //!< moment of inerria for orbit of light P
  double separationEnergy; //!< separation energy in MeV
  double fPair; //!< pairing energy of residue
  double fShell; //!< shell correction for residue
  CNuclide residue; //!< contains info on the daughter nucleus  
  bool odd; //!< true for odd mass
  SStoreEvapVector storeEvap; //!<store evaporation info
  double width; //!<decay width in MeV
  double rLight; //!< radius of light particle in fm
};

#endif
