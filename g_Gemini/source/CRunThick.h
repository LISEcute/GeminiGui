#include "CNucleus.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

/**
 *!\brief example of fusion with root histograms
 *
 * class that I use to simulate statistical decay in fusion 
 * reactions where there is a thick target. 
 * It produces a number of histograms, such a mass, charge
 * and evaporation spectra - The output is stored in root file
 */ 


class CRunThick
{

 public:
  CRunThick(int iZ, int iA, double fEx_min,double fEx_max, double l0_min, 
       double l0_Max,double d0, int lmax, double plb,int nBins,
     int numTot,string title0,double vcm=0.,double thetaDetMin=0.,
     double thetaDetMax = 360.);
};
