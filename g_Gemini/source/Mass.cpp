#include "CMass.h"
#include <cmath>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QSqlRecord>
#include <QElapsedTimer>


extern bool _useAME;
extern QSqlDatabase AmeDB;
extern QSqlDatabase Geminidb;
//int const Zmax = 136; // Value was pulled from the chart.cpp file

CMass* CMass::fInstance = 0;  // mod-TU
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CMass::CMass()
{
  chart = CChart::instance();
  fExpMass = new float [chart->iMassDim];
  fAMEMass = new float [chart->iMassDim];
  fCalMass = new float [chart->iMassDim];
  fFRM     = new float [chart->iMassDim];
  fPair = new float [chart->iMassDim];
  fShell   = new float [chart->iMassDim];
  fShell2 = new float [chart->iMassDim];
  //fBeta2 = new float [chart->iMassDim];


  ReadThomasFermiFile();
  ReadFRDMFile();
  ReadAMEDatabase();

  if(AmeDB.isOpen())      AmeDB.close();
  if(Geminidb.isOpen())   Geminidb.close();


  fExpMass[0] = 8.071;
  fExpMass[1] = 7.289;
  fCalMass[0] = 8.071;
  fCalMass[1] = 7.289;


  int iZmin = 4;
  int iZmax = 136;
  for (int iZ =iZmin;iZ<=iZmax;iZ++)
    {
      int iNmin = chart->getAmin(iZ)-iZ;
      int iNmax = chart->getAmax(iZ)-iZ;
      //cout << iZ << " " << iNmin << " " << iNmax << endl;
      for (int iN = iNmin; iN<=iNmax;iN++)
        {
          int index = chart->getIndex(iZ,iZ+iN);
          fPair[index] = getPairing2(iZ,iZ+iN);

          //redefine the shell correction to be difference between
          //experimental mass and (FRM masses plus given pairing corrections)
          fShell2[index] = fExpMass[index] - fPair[index] - fFRM[index];
        }
    }



  for (int iZ =iZmin;iZ<=iZmax;iZ++)
    {
      int iNmin = chart->getAmin(iZ)-iZ;
      int iNmax = chart->getAmax(iZ)-iZ;

      for (int iN = iNmin; iN<=iNmax; iN++)
        {
          int index = chart->getIndex(iZ,iZ+iN);

          float fpairN=0;

          if(iN%2 == 1) fpairN = 0.;
          else if (iN-1 >= iNmin && iN+1 <= iNmax)
            {
              fpairN = -(getShellCorrection2(iZ,iZ+iN-1) + getShellCorrection2(iZ,iZ+iN+1))/2.
                  +      getShellCorrection2(iZ,iZ+iN);
            }
          else if (iN-1 < iNmin) fpairN = -getShellCorrection2(iZ,iZ+iN+1)
              + getShellCorrection2(iZ,iZ+iN);
          else if (iN+1 > iNmax) fpairN = -getShellCorrection2(iZ,iZ+iN-1)
              + getShellCorrection2(iZ,iZ+iN);


          //then proton pairing
          float fpairZ=0;

          if(iZ%2 == 1) fpairZ = 0.;
          else if (iZ+1 <= iZmax) fpairZ =
              -(getShellCorrection2(iZ-1,iZ-1+iN)
                + getShellCorrection2(iZ+1,iZ+iN+1))/2.
              + getShellCorrection2(iZ,iZ+iN);
          else if (iZ+1 > iZmax)  fpairZ =
              -getShellCorrection2(iZ+1,iZ+iN+1)+ getShellCorrection2(iZ,iZ+iN);


          fPair[index] = fpairN + fpairZ + getPairing2(iZ,iZ+iN);

          fShell[index] = fExpMass[index] - fFRM[index]  - fPair[index];

        }
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CMass* CMass::instance() // mod-TU
{
  if (fInstance == 0) {fInstance = new CMass;}
  return fInstance;
}

//**********************************************
/**
 * Destructor
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CMass::~CMass()
{
  delete [] fExpMass;
  delete [] fCalMass;
  delete [] fFRM;
  // delete [] fBeta2;
  delete [] fShell;
  delete [] fShell2;
  delete [] fAMEMass;
  delete [] fPair;  // Oleg ==> error in original GEMINI
}
//********************************************
/**
 * Returns the experimental mass excess
 *
 * If the experimental excess is not known, then the Moller Nix value
 * is returned
 \param iZ is the proton number
 \param iA is the mass number
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
float CMass::getExpMass(int iZ, int iA)
{
  int i = chart->getIndex(iZ,iA);

  if( :: _useAME)  return fAMEMass[i];
  else             return fExpMass[i];
}
//********************************************
/**
 * Returns the calculated mass excess from Moller and Nix
 *
 \param iZ is the proton number
 \param iA is the mass number
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
float CMass::getCalMass(int iZ, int iA)
{

  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  return fCalMass[i];
}
//********************************************
/**
 * Returns the shell correction from Moller and Nix
 *
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getShellCorrection(int iZ, int iA)
{

  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  return fShell[i];

}
//********************************************
/**
 * Returns the shell correction from Moller and Nix
 *
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getShellCorrection2(int iZ, int iA)
{

  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  return fShell2[i];

}
//********************************************
/**
 * Returns the liquid drop mass from moller and Nix
 \param iZ is the protom number
 \param iA is the mass number
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
float CMass::getLDM(int iZ, int iA)
{
  //find location of nuclide in array where the mass is stored
  int i = chart->getIndex(iZ,iA);
  if ( i == -1)  { return -1000;}

  return fFRM[i];
}
//*************************************************
/**
   *
   * Calculates macroscopic finite range model masses of spherical
   * nucleus using formula of Krappe, Nix, and Sierk.
   *
   * Reference- (Phys Rev C20, 992 (1979))
   * modified to use the parameters of Moller + Nix Nucl. Phys. A361(1981)
   * 117. Pairing correction term for odd-odd nuclei
   * is included, as this is the most appropriate ground state for hot nuclei
   * where shell and pairing effects have washed out.
   \param iZ is the proton number
   \param iA is the mass number
   */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
float CMass::getFRM(int iZ, int iA)
{
  return getFRM((float)iZ,(float)iA);
}
//**********************************************************
/**
   *
   * Calculates macroscopic finite range model masses of spherical
   * nucleus using formula of Krappe, Nix, and Sierk.
   *
   * Reference- (Phys Rev C20, 992 (1979))
   * modified to use the parameters of Moller + Nix Nucl. Phys. A361(1981)
   * 117. Pairing correction term for odd-odd nuclei
   * is included, as this is the most appropriate ground state for hot nuclei
   * where shell and pairing effects have washed out.
   \param fZ is the proton number
   \param fA is the mass number
   */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
float CMass::getFRM(float fZ, float fA)
{
  float fN = fA - fZ;
  float fA13 = pow(fA,(float)(1./3.));
  float fA23 = pow(fA13,2);

  // relative neutron excess
  float fI = (fN-fZ)/fA;
  float fFiss = pow(fZ,2)/fA13;

  // neutron-proton terms
  float const fMassN = 8.071431;
  float const fMassH = 7.289034;
  float fEnz = fMassN*fN + fMassH*fZ;

  //Volume energy
  float const fAv = 15.9937;
  float const fKv = 1.927;
  float fEvol = -fAv*(1.-fKv*pow(fI,2))*fA;

  // Surface energy
  float const fa = 0.68;
  float const fR0 = 1.16;
  float const fAs = 21.13;
  float const fKs = 2.3;
  float fX=fa/(fR0*fA13);
  float fact=1.-3.*pow(fX,2)
      + (1.+1./fX)*(2.+3.*fX+3.*pow(fX,2))*exp(-2./fX);
  float fEsurf = fAs*(1.-fKs*pow(fI,2))*fA23*fact;

  //Coulomb energy
  float const e2 = 1.4399764;
  fact = fFiss-0.76361*pow(fZ,(float)(4./3.))/fA13;
  float fECoul = 0.6*e2/fR0*fact;

  //Wigner term
  float const fW = 36.;
  float const ael = 1.433e-5;
  float fEwigner = fW*fabs(fI)-ael*pow(fZ,(float)2.39);

  //correction to Coulomb energy for diffuse surface
  //see Davies & Nix Phys. Rev. C14 (1976) 1977
  float const b = 0.99;
  float ad=0.7071*b;
  fX = ad/(fR0*fA13);
  fact= 1 - 1.875*fX+2.625*pow(fX,3)
      -.75*exp(-2./fX)*(1.+4.5*fX+7.*pow(fX,2)+3.5*pow(fX,3));
  float fEcd = -3.*pow(fZ,2)*e2*pow(ad,2)/pow(fR0*fA13,3)*fact;

  // correction to coulomb energy due to proton form factor
  float const rp = 0.8;
  float akf=1./fR0*pow(7.06858*fZ/fA,1./3.);
  fX = pow(rp*akf,2);
  float fEcpf=-0.125*pow(rp,2)*e2/pow(fR0,3)
      *(3.0208-0.113541667*fX+0.0012624*pow(fX,2))*pow(fZ,2)/fA;

  //A0 term
  float const c0 = 4.4;
  float fEa0=c0;

  //Charge asymmetry term
  float const ca = 0.212;
  float fEca=ca*(pow(fZ,2)-pow(fN,2))/fA;

  // pairing term for odd-odd nuclei
  float deltau = 12./sqrt(fA);
  float deltal = 20./fA;
  float fEpair =deltau - 0.5 * deltal;

  // add all terms
  return fEnz+fEvol+fEsurf+fECoul+fEwigner+fEcd+fEcpf+fEa0+fEca+fEpair;
}
//****************************************
/**
 * Returns the pairing correction to the mass formula.
 * From from Moller Nix is used.
 \param iZ is the proton number
 \param iA is the mass number
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

float CMass::getPairing(int iZ, int iA)
{
  if(iZ==0 || iZ==iA) return 0.;
  return fPair[chart->getIndex(iZ,iA)];
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//****************************************
/**
 * Returns the pairing correction to the mass formula.
 * From from Moller Nix is used.
 \param iZ is the proton number
 \param iA is the mass number
 */
float CMass::getPairing2(int iZ, int iA)
{

  if(iZ==0 || iZ==iA)
    return 0.;

  float fZ = iZ;
  float fA = iA;
  float fN = fA - fZ;
  int   iN = iA - iZ;
  int  ioez = iZ%2;
  int  ioen = iN%2;
  float fPairing;

  if (iN == iZ && ioez ==1)
    {
      fPairing = 4.8/pow(fN,(float)(1./3.))
          + 4.8/pow(fZ,(float)(1./3.)) - 6.6/pow(fA,(float)(2./3.))
          + 30./fA;
    }
  else if (ioez == 1 && ioen == 1)
    {
      fPairing = 4.8/pow(fN,(float)(1./3.))
          + 4.8/pow(fZ,(float)(1./3.)) - 6.6/pow(fA,(float)(2./3.));
    }
  else if (ioez == 1 && ioen == 0)
    {
      fPairing = 4.8/pow(fZ,(float)(1./3.));
    }
  else if (ioez == 0 && ioen == 1)
    {
      fPairing = 4.8/pow(fN,(float)(1./3.));
    }
  else fPairing = 0.;

  //want to redefine odd-odd to have zero paring energy
  fPairing -=   4.8/pow(fN,(float)(1./3.))
      + 4.8/pow(fZ,(float)(1./3.))
      - 6.6/pow(fA,(float)(2./3.));

  if (iN == iZ) fPairing -= 30./fA;


  return fPairing;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// Reads in the mass table from Moller and Nix
 //========================================================================
//========================== Search Function =============================
//========================================================================
void CMass::FRDMFinder(QVector<QVector<QVariant>> &result)
{
  if (Geminidb.isValid() && Geminidb.isOpen())
    {
      QSqlQuery query(Geminidb);
      query.prepare("SELECT A, Z, k, l, m FROM mass");

      if (!query.exec()) { qDebug() << "Search didn't work"; return; }

      while (query.next())
        {
          QVector<QVariant> values;
          for(int k=0; k<5; k++) values.append(query.value(k));

          result.append(values);
        }
    }
  else  { qDebug() << "Gemini DB connection is not valid.";  }

  return;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void CMass::ReadFRDMFile()
{
  QVector<QVector<QVariant>> Array;
  FRDMFinder(Array);

  for(int i=0; i < Array.size(); i++)
    {
      QList temp = Array[i];
      bool check = true;

      int A = temp[0].toInt();
      int Z = temp[1].toInt();
      double f1 = temp[2].toDouble();
      double f2 = temp[3].toDouble();
      double f3 = 0;

      if(temp[4].toString().isEmpty())  check = false;
      else                              f3 = temp[4].toDouble();

  //    int min = chart->getAmin(Z);
  //    int max = chart->getAmax(Z);
      int index = chart->getIndex(Z,A);

      if (index >= 0)
        {
          float fPair = getPairing2(Z,A);

          if (check)        fExpMass[index] = f3;
          else              fExpMass[index] = f2;

          fCalMass[index] = f2;
          fFRM[index] = f2 - f1 - fPair;
          fShell[index] = f1;
        }
    }
//    Geminidb.close();
//    QSqlDatabase::removeDatabase(Geminidb.connectionName());
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// Reads in the mass table from the Thomas Fermi Model of Myers and Swietcki
 //========================================================================
//========================== Search Function =============================
//========================================================================
void CMass::ThomasFinder(QVector<QVector<QVariant>> &result)
{ 

  if (Geminidb.isValid())
    {
      QSqlQuery query(Geminidb);
      query.prepare("SELECT Z, A,  [M_EX], [M_TH], _SHL_, EVOD FROM mass_tf");

      if (!query.exec())  { qDebug() << "Search didn't work"; return;}

      while (query.next())
        {
          QVector<QVariant> values;
          for(int k=0; k<5; k++) values.append(query.value(k));

          result.append(values);
        }
    } //----------------------------------------------
  else    { qDebug() << "Gemini DB connection is not valid.";  }

  return;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void CMass::ReadThomasFermiFile()
{
  QVector<QVector<QVariant>> Array;
  ThomasFinder(Array);

  for(int i=0; i < Array.size(); i++)
    {
      QList temp = Array[i];
      int A = temp[1].toInt();
      int Z = temp[0].toInt();
      float f1 = temp[3].toFloat();
      float f2 = temp[2].toFloat();
      float fshell = temp[4].toFloat();

 //     int min = chart->getAmin(Z);
 //     int max = chart->getAmax(Z);
      int index = chart->getIndex(Z,A);

      float fPairing;

      if (index >= 0)
        {
          fPairing = getPairing2(Z,A);
          fExpMass[index] = f2;
          if (Z < 5) fCalMass[index] = f2;
          else fCalMass[index] = f1;
          fFRM[index] = fCalMass[index] - fshell - fPairing;
          fShell[index] = fshell*2.;
        }
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//========================== Search Function =============================
//========================================================================
void CMass::AMEFinder(QMap<QString, float> &result)
{
  // OT qDebug() << 4;
  QSqlQuery query(AmeDB);
  query.prepare("SELECT A, Z, MASS_EXCES FROM AME2016 ORDER BY [INDEX] ASC");

  if (!query.exec())
    {
      qDebug() << "Search didn't work";
      return;
    }

  while (query.next())
    {
      QString key = query.value(0).toString() + query.value(1).toString();
      float value = query.value(2).toFloat();
      result.insert(key, value);
    }

  return;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void CMass::ReadAMEDatabase()
{
  QMap<QString, float> AME;
  AMEFinder(AME);

  for (int Z=0; Z<120; Z++)
    {
      int min = chart->getAmin(Z);
      int max = chart->getAmax(Z);

      for (int A = min; A <= max; ++A)
        {
          QString key = QString::number(A) + QString::number(Z);
          float value = AME.contains(key) ? AME.value(key) : getCalMass(Z,A);

          int index = chart->getIndex(Z,A);

          fAMEMass[index] = value;
        }
    }

  fAMEMass[0] = 8.071;
  fAMEMass[1] = 7.289;

//  AME.clear();
//  QSqlDatabase::removeDatabase(AmeDB.connectionName());
//  AmeDB.close();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
