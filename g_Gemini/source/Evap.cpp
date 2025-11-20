#include "CEvap.h"

#include <iostream>
#include <QTextStream>
#include <QFile>
CEvap* CEvap::fInstance = 0;

float const r0 = 1.16;
/**
 * Constructor
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CEvap::CEvap() 
{
  QString fileName(":/tbl/evap.inp");
  QFile iff(fileName);
  if(!iff.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      cout << " file " << fileName.toStdString() << " not found" << endl;
      abort();
    }

  QTextStream ifFile(&iff);

  // read in number of evaporation channels
  ifFile >> nLight;
  prob = new float [nLight];
  // create array of CLightP pointers
  lightP = new CLightP * [nLight];

  tlArray = new CTlBarDist * [nLight];
  sigBarDist = new CSigBarDist * [nLight];

  // read in channel information and initialize pointers
  decay = new SDecay [nLight];

  // skip line
  // string line;
  // getline(ifFile,line);
  // getline(ifFile,line);

  ifFile.readLine();
  ifFile.readLine();


  float fJ, fEx, Ek, suppress;
  int iZ, iA;
  QString nameOld("");
  QString name("");
  bool newTl;
  nTl = 0;
  for (int i=0;i<nLight;i++)
    {
      nameOld = name;
      ifFile >> iZ >> iA >> fJ >> fEx >> name >> suppress >> Ek;
      //  if(i == 0) cout << iZ << " " << iA << " " <<fJ << " " <<fEx << " " << name << " " << suppress<<" "  << Ek << endl; // checked ok


      // determine if we need to create new transmission coeff
      newTl = 1;
      if (i > 0)
        {
          if (name == nameOld) newTl = 0;
        }
      if (newTl)
        {
          tlArray[nTl] = new CTlBarDist(name);
          sigBarDist[nTl] = new CSigBarDist(name,(float)iZ,(float)iA);
          nTl++;
        }
      lightP[i] = new CLightP(iZ,iA,fJ,tlArray[nTl-1],sigBarDist[nTl-1]);
      lightP[i]->rLight = pow((float)iA,(float)(1./3.))*r0;
      lightP[i]->fEx = fEx;
      lightP[i]->suppress = suppress;
      // if an excited state add excitation energy to mass excess
      if (fEx > 0.0) lightP[i]->fExpMass += fEx;
      maxZ = iZ;
      decay[i].Ek = Ek;
      if (Ek == 0.) continue;
      // read in decay information
      ifFile >> decay[i].Z1 >> decay[i].A1 >> decay[i].S1 >>
          decay[i].S2 >> decay[i].L >> decay[i].lPlusS1 >> decay[i].gamma;

      // cout <<  "DECAY " << decay[i].Z1 << " "<< decay[i].A1 << " "<< decay[i].S1<< " " <<
      //decay[i].S2 << " "<< decay[i].L << " "<< decay[i].lPlusS1 << " "<< decay[i].gamma << endl; //checked ok


      // here are sum checks to make sure decay information is possible
      if (decay[i].lPlusS1 > decay[i].S1 + (float)decay[i].L)
        {
          cout << "bad LPlusS1 in evap.cpp for mode " << i << endl;
          abort();
        }

      if (decay[i].lPlusS1 < fabs(decay[i].S1 - (float)decay[i].L))
        {
          cout << "bad LPlusS1 in evap.cpp for mode " << i << endl;
          abort();
        }

      if (fJ > decay[i].lPlusS1+decay[i].S2)
        {
          cout << "bad fJ in evap.cpp for mode " << i << endl;
          abort();
        }

      if (fJ < fabs(decay[i].lPlusS1-decay[i].S2))
        {
          cout << "bad fJ in evap.cpp for mode " << i << endl;
          abort();
        }
    }

  //ifFile.clear();
  //ifFile.close();
  ifFile.flush();
  iff.close();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CEvap* CEvap::instance()
{
  if (fInstance == 0) {fInstance = new CEvap;}
  return fInstance;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//**************************************************
/**
 * Destructor
 */
CEvap::~CEvap()
{
  for (int i=0;i<nLight;i++) delete lightP[i];
  for (int i=0;i<nTl;i++)
    {
      delete tlArray[i];
      delete sigBarDist[i];
    }
  delete [] prob;
  delete [] lightP;
  delete [] tlArray;
  delete [] sigBarDist;
  delete [] decay;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
