#ifndef mass_
#define mass_
#include "CChart.h"
#include <iomanip>

#include <QMessageBox>
#include <QDir>

using namespace std;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

/**
 *!\brief mass excesses, pairing energy, etc
 *
 * Class associated with returning quanties associated with
 * the mass formula
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class CMass
{
protected:
  CMass();                  //!< constructor
  static  CMass* fInstance; //!< instance member to make tis a singleton

  float * fExpMass;  //!<experimental mass array
  float * fAMEMass; //!< AME mass array MPK
  float * fCalMass;  //!<experimental mass array
  float * fFRM;      //!<finite range mass array
  float * fPair; //!< pairing correction
  float * fShell;    //!< shell correction
  float * fShell2;    //!< shell correction

public:
  // mod-TU CMass();
  ~CMass();
  CChart *chart; //!< contains the considered region of the chart of nuclides
  static CMass* instance(); //!< instance member to make this a singleton

  float getExpMass(int iZ, int iA);
  float getCalMass(int iZ,int iA);
  float getShellCorrection(int iZ, int iA);
  float getShellCorrection2(int iZ, int iA);

  float getFRM(float fZ, float fA);
  float getFRM(int iZ ,int iA);
  float getLDM(int iZ ,int iA);

  float getPairing(int iZ,int iA);
  float getPairing2(int iZ,int iA);
  bool useAME;

private:
  void ReadFRDMFile();
  void ReadThomasFermiFile();
  void ReadAMEDatabase();

  void AMEFinder(QMap<QString, float> &result);
  void ThomasFinder(QVector<QVector<QVariant>> &result);
  void FRDMFinder(QVector<QVector<QVariant>> &result);
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#endif
