#include "g_Gemini/source/CNucleus.h"
#include "g_Gemini/source/CFus.h"
// this is an example of using GEMINI CNucleus class to give the
//statistical decay of a compound nucleus form in a heavy ion fusion reaction
// It calculates the fission, residue and IMF probabilities and xsections.
// Also some multiplicties of evaporated particles
extern bool _useAME;


int testFusion()
{
  bool save_ame = _useAME;
_useAME = false;

    cout << endl << "===================  test Fusion ==================" << endl;

  int Zp = 8; // proton number of projectile
  int Ap = 16; // mass number of projectile
  int Zt = 6; // proton number of target
  int At = 12; //mass number of target
  double Elab = 160.; // labe energy in MeV
  double dif = 2; //diffuseness of fusion spin distribution in hbar

  CFus fus(Zp,Ap,Zt,At,Elab,dif);
  //double l0 = fus.getBassL();  // get maximum spin from Bass model
  double l0 = 23;  // get maximum spin from paper


  //construct fusion spin distribution
  // OT 06/30/2023
  int lmax = qMin((int)l0 +5, 1000);  // maxium spin considered


  // OT 06/30/2023
  //double prob[lmax+1];
  double prob[1000];

  double sum = 0.;
  for (int l=0;l<=lmax;l++)
    {
      prob[l] =  (double)(2*l+1);
      if (dif > 0.) prob[l] /= (1.+exp(((double)l-l0)/dif));
      else if ( l > l0) prob[l] = 0.;
      sum += prob[l];
    }
  for (int l=0;l<=lmax;l++)
    {
      prob[l] /= sum;
      if (l > 0) prob[l] += prob[l-1];
    }

  CNucleus CN(fus.iZcn,fus.iAcn); //constructor
  double fEx = fus.Ex; //excitation energy of compound nucleus


  cout << " compound Nucleus is " << CN.getName() << endl;
  cout << " E* = " << fEx << " critical spin = " << l0 << endl;
  cout << " fusion xsection = " << sum*fus.plb << " mb" << endl;
  double xfus = sum*fus.plb;

  CN.setEvapMode(1);  // force a fully Hauser-Feshbach calculation


  // if IMF probablities are small and you are not interested
  // save time and turn them off.
  //CNucleus::setNoIMF();

  // don't like the default evaporation parameters
  //CLevelDensity::setAfAn(1.036);  // change ratio of saddle to ground state
  //level density parameters
  //CLevelDensity::setLittleA(8.); //change level density parameter to A/8.



  //write out statistical model parameters
  CN.printParameters();


  //define and zero events parameters

  double total = 0.;
  double Nfission = 0.;
  double Nimf = 0.;
  double neutPreSad = 0.;
  double neutSaddleToScission = 0.;
  double neutHeavy = 0.;
  double neutLight = 0.;


  double Ares = 0.;
  double Zres = 0.;
  double resTotal = 0.;
  double neutMultEv = 0.;
  double protMultEv = 0.;
  double alpMultEv = 0.;
  double gammaEnergy = 0.;


  double sum3 = 0.;
  double sum4 = 0.;
  double sum5 = 0.;
  double sum6 = 0.;


  int Nevents = 500;
  double Sconst = xfus/(double)Nevents;

  for (int i=0;i<Nevents;i++)
    {
      //choose the spin of the CN from the determed spin distribution
      double ran = CN.ran->Rndm();
      int l = 0;
      for (;;)
        {
          if (ran < prob[l]) break;
          l++;
        }

      //specify the excitation energy and spin
      CN.setCompoundNucleus(fEx,(double)l);



      //if you are interested in IMF emission at low excitation energy
      //then turn IMF weighting on
      CN.setWeightIMF();// turn on enhanced IMF emission


      CN.decay(); //decay the compound nucleus


      if (CN.abortEvent)  //gemini had trouble with this event
        {
          CN.reset();
          continue;
        }


      // number of stable fragments producted in the decay
      int Nfrag = CN.getNumberOfProducts();

      //set pointer to last fragment which is the evaporation residue
      // if isResidue is true, otherwise the heavy fission fragment
      // this should be the heaviest fragment produced
      CNucleus * products = CN.getProducts(Nfrag-1);

      // the weight will be unity unless setWeightIMF is called
      double weight = products->getWeightFactor();



      if (CN.isSymmetricFission())
        {
          Nfission += weight;  //fission event
          products = CN.getProducts(0);  // go to first evaporated particle
          for (int i=0;i<Nfrag-1;i++)
            {
              if (products->iZ == 0 && products->iA == 1) // look for neutrons
                {
                  if (products->isSaddleToScission())
                    neutSaddleToScission += weight;
                  else if (products->origin == 1) cout << "hell " << endl;
                  else if (products->origin == 0) neutPreSad += weight;
                  else if (products->origin == 2) neutLight += weight;
                  else if (products->origin == 3) neutHeavy += weight;

                }
              // go to next particle
              products = CN.getProducts();
            }

        }

      //intermediate mass fragment
      if (CN.isAsymmetricFission()) Nimf += weight;  //imf event
      total += weight;


      //evaporation resiudes
      if (CN.isResidue())
        {
          Ares += products->iA;
          Zres += products->iZ;
          resTotal += weight;
          gammaEnergy += weight*products->getSumGammaEnergy();
        }
      products = CN.getProducts(0);  // go to first evaporated particle
      for (int i=0;i<Nfrag-1;i++)
        {
          if (products->iZ == 0 && products->iA == 1 && CN.isResidue())  //neutrons
            {
              neutMultEv += weight;
            }
          else if (products->iZ == 1 && products->iA == 1 && CN.isResidue()) //protons
            {
              protMultEv += weight;
            }
          else if (products->iZ == 2 && products->iA == 4 && CN.isResidue())//alpha particles
            {
              alpMultEv += weight;
            }

          else if (products->iZ == 3) sum3 += weight;
          else if (products->iZ == 4) sum4 += weight;
          else if (products->iZ == 5) sum5 += weight;
          else if (products->iZ == 6) sum6 += weight;
          // go to next particle
          products = CN.getProducts();
        }





      //reset the compound nucleus for a new decay
      CN.reset();
    }

  cout << endl;
  cout << "fission probability = " << Nfission/total << "  xsection = " <<
          Nfission*Sconst << " mb" <<endl;

  if (Nfission > 0.)
    {
      cout << " neutron multiplicities in fission" << endl;
      cout << "  presaddle = " << neutPreSad/Nfission << endl;
      cout << "  saddle-to-scission  = " << neutSaddleToScission/Nfission << endl;
      cout << "  post-scission light frag = " << neutLight/Nfission << endl;
      cout << "  post-scission heavy frag = " << neutHeavy/Nfission << endl;
      cout << endl;
    }


  cout << endl;

  cout << "residue (no IMF or Symmetric fission emissions)" << endl;
  cout << " residue xsection = " << resTotal*Sconst << " mb" << endl;
  if (resTotal > 0.)
    {
      cout << "  average residue is Z = " << Zres/resTotal << " A = "
       << Ares/resTotal << endl;

      cout << "  average neutron multiplicity = " << neutMultEv/resTotal << endl;
      cout << "  average proton multiplicity = " << protMultEv/resTotal << endl;
      cout << "  average alpha multiplicity = " << alpMultEv/resTotal << endl;
      cout << " average energy in gamma rays = " << gammaEnergy/resTotal <<
              " MeV" << endl;

      cout << endl;
      cout << " intermediate Mass Fragment (Z > " << CN.getZmaxEvap() << ") (IMF) " << endl;
      cout << " IMF prob = " << Nimf/total
           << " xsection = " << Nimf*Sconst << " mb" <<endl;
      cout << " IMF yields" << endl;
      cout << " xsection of Z=3 " << sum3*Sconst << " mb" << endl;
      cout << " xsection of Z=4 " << sum4*Sconst << " mb" << endl;
      cout << " xsection of Z=5 " << sum5*Sconst << " mb" << endl;
      cout << " xsection of Z=6 " << sum6*Sconst << " mb" << endl;

    }

_useAME = save_ame;
 return 1;

}
