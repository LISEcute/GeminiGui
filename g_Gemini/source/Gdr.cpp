#include "CGdr.h"
#include <QTextStream>
#include <QFile>
#include <QString>

CGdr* CGdr::fInstance = 0; 


/**
 * Constructor
 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CGdr::CGdr()
{
  QString fileName(":/tbl/GDR.inp");
  QFile iff(fileName);
  if(!iff.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      cout << "GDR file " << fileName.toStdString() << " not found" << endl;
      abort();
    }

  QTextStream ifFile(&iff);

  //skip line
  ifFile.readLine();

  //read in paramters of Lorentzian's , up to five are allowed
  double a,b,c;
  N = 0;
  for (int i=0;i<5;i++)
    {
      ifFile >> a >> b >> c;
      if (ifFile.atEnd()) break;
      if (a== 0.) break;
      lineShape[N].strength =  a;
      lineShape[N].energy = b;
      lineShape[N].gamma = c;
      N++;
    }

  if (N == 0)
    {
      cout << "no user-defined GDR" << endl;
      abort();
    }
}
//******************************************************
/**
   * returns the GDR line shape as a sum of up to 5 user-defined
   * Lorentzians
   \param e is the gamma energy in MeV
   */

float CGdr::getLineShape(float e)
{
  float out = 0.;
  for (int i=0;i<N;i++)
    {
      out += lineShape[i].gamma*lineShape[i].strength*pow(e,4)/
          (pow(pow(e,2)-pow(lineShape[i].energy,2),2)+pow(lineShape[i].gamma*e,2));
    }
  return out;
}
//***************************
CGdr* CGdr::instance()
{
  if (fInstance == 0) {
      fInstance = new CGdr;
    }
  return fInstance;
}

