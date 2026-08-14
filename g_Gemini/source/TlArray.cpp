#include "CTlArray.h"
#include <QTextStream>
#include <QFile>
/**
 * Constructor
\param sName0 is the name of the *.tl file in directory /tl where
the parameters specifying the transmission coefficents are contained
*/

CTlArray::CTlArray(const QString& sName0)
{
  sName = ":/tl/"+sName0+".tl";

  QFile iff(sName);
  if(!iff.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      cout << "file " << sName.toStdString() << " not found in CTlArray" << std::endl;
      abort();

    }

  QTextStream ifFile(&iff);


  double fx[7];
  // string line;
  //getline(ifFile,line);
  ifFile.readLine();
  ifFile >> shift;
  iZMin = 2;
  int kk;
  for (int i=1;i<=120;i++)
    {
      
      ifFile >> kk >> fx[0] >> fx[1] >> fx[2]
          >> fx[3] >> fx[4] >> fx[5] >> fx[6];

      //see if data is ok
      if (isnan(fx[0]))
        {
          iZMin = kk;
          break;
        }
      if (fx[0] == 0. && fx[1] == 0. && fx[2] == 0. && fx[3] == 0.)
        {
          iZMin = kk;
          break;
        }
      
      //   if (ifFile.bad() || ifFile.eof())
      if(ifFile.atEnd())
        {
          cout << "trouble reading file " << sName.toStdString() << " in CTlArray" << std::endl;
        }
      for (int j=0;j<7;j++)zcoef.Tl[kk].coef[j] = fx[j];

    }
  ifFile.flush();
  iff.close();
}

//************************************************************
/**
   *  prepares for a series of calls with the same Z value
   \param iZ is the proton number of the daughter
   */

void CTlArray::prepare(int iZ)
{
  trans = &zcoef.Tl[iZ];
}


//***********************************************************
/**
 * The transmission coeff are parameterized as \f$\frac{1}{1+exp(x)}\f$
 * this function returns the value of x
 \param iL is the orbital angular momentum of the evaporated particle
 \param fEk is the kinetic energy of the evaporated particle
 */

double CTlArray::getTermInExp(int iL, double fEk)
{
  double fL = (double)(iL+1);
  double fC1 = trans->coef[6];
  double fC2 = trans->coef[0] + fL*(trans->coef[1] + fL*trans->coef[2]);
  double fC3 = trans->coef[3] + fL*(trans->coef[4] + fL*trans->coef[5]);

  fEk += shift;

  // the interpolation can have a maximum, must find
  double emax; // energy of maximum or minimum
  double d2e; //second derivative at stationary point
  
  if (fC1 == 0.0)
    {
      emax = 0.0;
      d2e = 0.0;
    }
  else
    {
      emax = fC3/fC1;
      if  (emax > 0.0) //real stationary points
        {
          emax = sqrt(emax);
          d2e = 2.0*fC3/(emax*emax*emax);
        }
      else
        {
          emax = 0.0;//no real maximum
          d2e = 0.0;
        }
    }
  double fX;
  if (fEk < emax && d2e < 0.0)
    {
      //fit has a maximum at positive ek, set tl=0 below this energy
      fX = 100.;
    }
  else if (fEk > emax && d2e > 0.0)
    {
      //fit has minimum at positive ek, set eetl at it minimum
      //value above this energy
      fX = fC1*emax + fC2 + fC3/emax;
    }
  else
    {
      fX  = fC1*fEk + fC2 + fC3/fEk;
    }

  return fX;
}
//***************************************************
/**
 * Returns the transmission coefficient
 \param iL is the orbital angular momentum of the evaporated particle
 \param fEk is the kinetic energy of the evaporated particle
*/

double CTlArray::getTl(int iL, double fEk)
{
  double fX = getTermInExp(iL,fEk);
  return 1.0/(1.0+exp(fX));
}
//***********************************
/**
 * returns the quantity
 * \f$S=\sum_{\ell=0}^{\infty} (2\ell+1)T_{\ell}(\varepsilon)\f$
 * which is related to the inverse cross section by
 * \f$S=\frac{\sigma_{inv}}{\pi\lambda^{2}}\f$
 \param fEk is the kinetic energy of the evaporated particle
 */
double CTlArray::getInverseXsec(double fEk)
{
  double tot = 0.;
  double xmax = 0.;
  int iL = 0;
  for(;;)
    {
      double x = (double)(2*iL+1)*getTl(iL,fEk);
      xmax = max(x,xmax);
      tot += x;
      if (x < xmax*.01) break;
      iL++;
      if (iL > 20) break;
    }
  return tot;
}





