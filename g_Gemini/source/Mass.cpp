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

CMass* CMass::fInstance = 0;  // mod-TU

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CMass::CMass()
{
    chart = CChart::instance();
    fExpMass = new double [chart->iMassDim];
    fAMEMass = new double [chart->iMassDim];
    fCalMass = new double [chart->iMassDim];
    fFRM     = new double [chart->iMassDim];
    fPair    = new double [chart->iMassDim];
    fShell   = new double [chart->iMassDim];
    fShell2  = new double [chart->iMassDim];
    //fBeta2 = new double [chart->iMassDim];

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
    for (int iZ = iZmin; iZ <= iZmax; iZ++)
    {
        int iNmin = chart->getAmin(iZ) - iZ;
        int iNmax = chart->getAmax(iZ) - iZ;

        for (int iN = iNmin; iN <= iNmax; iN++)
        {
            int index = chart->getIndex(iZ, iZ + iN);
            fPair[index] = getPairing2(iZ, iZ + iN);

            //redefine the shell correction to be difference between
            //experimental mass and (FRM masses plus given pairing corrections)
            fShell2[index] = fExpMass[index] - fPair[index] - fFRM[index];
        }
    }

    for (int iZ = iZmin; iZ <= iZmax; iZ++)
    {
        int iNmin = chart->getAmin(iZ) - iZ;
        int iNmax = chart->getAmax(iZ) - iZ;

        for (int iN = iNmin; iN <= iNmax; iN++)
        {
            int index = chart->getIndex(iZ, iZ + iN);

            double fpairN = 0;

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
            double fpairZ = 0;

            if(iZ%2 == 1) fpairZ = 0.;
            else if (iZ+1 <= iZmax) fpairZ =
                    -(getShellCorrection2(iZ-1,iZ-1+iN)
                      + getShellCorrection2(iZ+1,iZ+iN+1))/2.
                    + getShellCorrection2(iZ,iZ+iN);
            else if (iZ+1 > iZmax)  fpairZ =
                    -getShellCorrection2(iZ+1,iZ+iN+1)+ getShellCorrection2(iZ,iZ+iN);

            fPair[index] = fpairN + fpairZ + getPairing2(iZ, iZ + iN);
            fShell[index] = fExpMass[index] - fFRM[index] - fPair[index];
        }
    }
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
CMass* CMass::instance() // mod-TU
{
    if (fInstance == 0) { fInstance = new CMass; }
    return fInstance;
}

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

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getExpMass(int iZ, int iA)
{
    int i = chart->getIndex(iZ,iA);

    if( :: _useAME)  return fAMEMass[i];
    else             return fExpMass[i];
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getCalMass(int iZ, int iA)
{
    int i = chart->getIndex(iZ,iA);
    return fCalMass[i];
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getShellCorrection(int iZ, int iA)
{
    int i = chart->getIndex(iZ,iA);
    return fShell[i];
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getShellCorrection2(int iZ, int iA)
{
    int i = chart->getIndex(iZ,iA);
    return fShell2[i];
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getLDM(int iZ, int iA)
{
    int i = chart->getIndex(iZ,iA);
    if ( i == -1)  { return -1000; }

    return fFRM[i];
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getFRM(int iZ, int iA)
{
    return getFRM((double)iZ,(double)iA);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getFRM(double fZ, double fA)
{
    double fN = fA - fZ;
    double fA13 = pow(fA,(double)(1./3.));
    double fA23 = fA13 * fA13;                 // pow(fA13,2)

    // relative neutron excess
    double fI = (fN-fZ)/fA;
    double fFiss = (fZ*fZ) / fA13;             // pow(fZ,2)/fA13

    // neutron-proton terms
    double const fMassN = 8.071431;
    double const fMassH = 7.289034;
    double fEnz = fMassN*fN + fMassH*fZ;

    //Volume energy
    double const fAv = 15.9937;
    double const fKv = 1.927;
    double fEvol = -fAv*(1.-fKv*(fI*fI))*fA;   // pow(fI,2)

    // Surface energy
    double const fa = 0.68;
    double const fR0 = 1.16;
    double const fAs = 21.13;
    double const fKs = 2.3;
    double fX = fa/(fR0*fA13);
    double fact = 1.-3.*(fX*fX)                // pow(fX,2)
                 + (1.+1./fX)*(2.+3.*fX+3.*(fX*fX))*exp(-2./fX); // pow(fX,2)
    double fEsurf = fAs*(1.-fKs*(fI*fI))*fA23*fact;      // pow(fI,2)

    //Coulomb energy
    double const e2 = 1.4399764;
    fact = fFiss - 0.76361*pow(fZ,(double)(4./3.))/fA13; // fractional power: keep pow
    double fECoul = 0.6*e2/fR0*fact;

    //Wigner term
    double const fW = 36.;
    double const ael = 1.433e-5;
    double fEwigner = fW*fabs(fI) - ael*pow(fZ,(double)2.39); // non-integer: keep pow

    //correction to Coulomb energy for diffuse surface
    //see Davies & Nix Phys. Rev. C14 (1976) 1977
    double const b = 0.99;
    double ad = 0.7071*b;
    fX = ad/(fR0*fA13);
    fact = 1 - 1.875*fX + 2.625*(fX*fX*fX)   // pow(fX,3)
           - .75*exp(-2./fX)*(1.+4.5*fX+7.*(fX*fX)+3.5*(fX*fX*fX)); // pow(fX,2), pow(fX,3)
    double fEcd = -3.*(fZ*fZ)*e2*(ad*ad) /    // pow(fZ,2), pow(ad,2)
                 ((fR0*fA13)*(fR0*fA13)*(fR0*fA13))   // pow(fR0*fA13,3)
                 * fact;

    // correction to coulomb energy due to proton form factor
    double const rp = 0.8;
    double akf = 1./fR0*pow(7.06858*fZ/fA,1./3.); // fractional: keep pow
    fX = (rp*akf)*(rp*akf);                      // pow(rp*akf,2)
    double fEcpf = -0.125*(rp*rp)*e2/(fR0*fR0*fR0) // pow(rp,2), pow(fR0,3)
                  *(3.0208-0.113541667*fX+0.0012624*(fX*fX)) // pow(fX,2)
                  * (fZ*fZ)/fA;                              // pow(fZ,2)

    //A0 term
    double const c0 = 4.4;
    double fEa0 = c0;

    //Charge asymmetry term
    double const ca = 0.212;
    double fEca = ca*((fZ*fZ) - (fN*fN))/fA;       // pow(fZ,2)-pow(fN,2)

    // pairing term for odd-odd nuclei
    double deltau = 12./sqrt(fA);
    double deltal = 20./fA;
    double fEpair = deltau - 0.5 * deltal;

    // add all terms
    return fEnz+fEvol+fEsurf+fECoul+fEwigner+fEcd+fEcpf+fEa0+fEca+fEpair;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getPairing(int iZ, int iA)
{
    if(iZ==0 || iZ==iA) return 0.;
    return fPair[chart->getIndex(iZ,iA)];
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CMass::getPairing2(int iZ, int iA)
{
    if(iZ==0 || iZ==iA)
        return 0.;

    double fZ = iZ;
    double fA = iA;
    double fN = fA - fZ;
    int   iN = iA - iZ;
    int  ioez = iZ%2;
    int  ioen = iN%2;
    double fPairing;

    // fractional powers: keep pow(...)
    if (iN == iZ && ioez ==1)
    {
        fPairing = 4.8/pow(fN,(double)(1./3.))
                   + 4.8/pow(fZ,(double)(1./3.)) - 6.6/pow(fA,(double)(2./3.))
                   + 30./fA;
    }
    else if (ioez == 1 && ioen == 1)
    {
        fPairing = 4.8/pow(fN,(double)(1./3.))
                   + 4.8/pow(fZ,(double)(1./3.)) - 6.6/pow(fA,(double)(2./3.));
    }
    else if (ioez == 1 && ioen == 0)
    {
        fPairing = 4.8/pow(fZ,(double)(1./3.));
    }
    else if (ioez == 0 && ioen == 1)
    {
        fPairing = 4.8/pow(fN,(double)(1./3.));
    }
    else fPairing = 0.;

    //want to redefine odd-odd to have zero paring energy
    fPairing -=   4.8/pow(fN,(double)(1./3.))
                + 4.8/pow(fZ,(double)(1./3.))
                - 6.6/pow(fA,(double)(2./3.));

    if (iN == iZ) fPairing -= 30./fA;

    return fPairing;
}

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

        int index = chart->getIndex(Z,A);

        if (index >= 0)
        {
            double fPair = getPairing2(Z,A);

            if (check)        fExpMass[index] = f3;
            else              fExpMass[index] = f2;

            fCalMass[index] = f2;
            fFRM[index] = f2 - f1 - fPair;
            fShell[index] = f1;
        }
    }
}

//========================================================================
//========================== Search Function =============================
//========================================================================
void CMass::ThomasFinder(QVector<QVector<QVariant>> &result)
{
    if (Geminidb.isValid())
    {
        QSqlQuery query(Geminidb);
        query.prepare("SELECT Z, A,  [M_EX], [M_TH], _SHL_, EVOD FROM mass_tf");

        if (!query.exec())  { qDebug() << "Search didn't work"; return; }

        while (query.next())
        {
            QVector<QVariant> values;
            for(int k=0; k<5; k++) values.append(query.value(k));
            result.append(values);
        }
    }
    else    { qDebug() << "Gemini DB connection is not valid."; }

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
        double f1 = temp[3].toDouble();
        double f2 = temp[2].toDouble();
        double fshell = temp[4].toDouble();

        int index = chart->getIndex(Z,A);

        double fPairing;

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

//========================================================================
void CMass::AMEFinder(QMap<QString, double> &result)
{
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
        double value = query.value(2).toDouble();
        result.insert(key, value);
    }

    return;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void CMass::ReadAMEDatabase()
{
    QMap<QString, double> AME;
    AMEFinder(AME);

    for (int Z=0; Z<120; Z++)
    {
        int min = chart->getAmin(Z);
        int max = chart->getAmax(Z);

        for (int A = min; A <= max; ++A)
        {
            QString key = QString::number(A) + QString::number(Z);
            double value = AME.contains(key) ? AME.value(key) : getCalMass(Z,A);

            int index = chart->getIndex(Z,A);
            fAMEMass[index] = value;
        }
    }

    fAMEMass[0] = 8.071;
    fAMEMass[1] = 7.289;
}
