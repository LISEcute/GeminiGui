#include "CSigCharged.h"
#include <QTextStream>
#include <QFile>
#include <QDebug>

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CSigCharged::CSigCharged(const QString& sName0, double Zp0, double Ap0)
{

  Zp = Zp0;
  Ap = Ap0;
  if (sName0.compare("neutron")==0)
    {
      neutron = 1;
      return;
    }
  neutron = 0;
  sName = ":/tl/" + sName0 + ".inv";

  QFile iff(sName);
  if(!iff.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      cout << "file " << sName.toStdString() << " not found in CSigCharged" << std::endl;
      abort();
    }

  QTextStream ifFile(&iff);

  ifFile >> rc0 >> rc1 >> rc2;
  ifFile >> rI0 >> rI1 >> rI2;
  ifFile >> omega0 >> omega1 >> omega2 >> omega3;
  ifFile >> a0;
  ifFile >> aa0;

  ifFile.flush();
  iff.close();
}
//******************************************
void CSigCharged::prepare(double Z, double A)
{
  if (neutron)
    {
     n0 = 10.386*(1.-exp((-Z/24.50)));
     n1 = 1.227+0.03444*Z;
     n2 = 2.;
     barrier = .5;
     return;
    }

  double A13 = pow(A,(double)(1./3.));
  double rc = rc0/A13 + rc1 + rc2*A13;
  barrier = Zp*Z*1.44/rc;

  double rI = rc0/A13 + rI1 + rI2*A13;
  double mu = A*Ap/(A+Ap);
  InvInertia = 2.*mu*(rI*rI)/41.563;
  omega = omega0 + omega1*A + omega2*exp(-A/omega3);

  a = a0;
  aa = aa0;

  offset = barrier/2.;
}
//*****************************************
double CSigCharged::getInverseXsec(double energy)
{

  if (neutron) return  (n0+ n1*energy)*(1.-exp(-energy/n2));
  double ddelta = sqrt((offset*offset)+(energy-barrier*energy-barrier)) - offset;
  double delta;

  if (energy > barrier) delta = energy - barrier - a*ddelta;   // ???  error??
  else delta = energy - barrier - aa*ddelta;


  double out = InvInertia*(delta + omega*log(1.+exp(-delta/omega)));
  if (out < 0.) out = 0.;
  return out;
}
//******************************************

double CSigCharged::getBarrier()
{
  return  barrier;
}
