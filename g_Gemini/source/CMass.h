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

  double * fExpMass;  //!<experimental mass array
  double * fAMEMass; //!< AME mass array MPK
  double * fCalMass;  //!<experimental mass array
  double * fFRM;      //!<finite range mass array
  double * fPair; //!< pairing correction
  double * fShell;    //!< shell correction
  double * fShell2;    //!< shell correction

public:
  // mod-TU CMass();
  ~CMass();
  CChart *chart; //!< contains the considered region of the chart of nuclides
  static CMass* instance(); //!< instance member to make this a singleton

  double getExpMass(int iZ, int iA);
  double getCalMass(int iZ,int iA);
  double getShellCorrection(int iZ, int iA);
  double getShellCorrection2(int iZ, int iA);

  double getFRM(double fZ, double fA);
  double getFRM(int iZ ,int iA);
  double getLDM(int iZ ,int iA);

  double getPairing(int iZ,int iA);
  double getPairing2(int iZ,int iA);
  bool useAME;

private:
  void ReadFRDMFile();
  void ReadThomasFermiFile();
  void ReadAMEDatabase();

  void AMEFinder(QMap<QString, double> &result);
  void ThomasFinder(QVector<QVector<QVariant>> &result);
  void FRDMFinder(QVector<QVector<QVariant>> &result);
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#endif
