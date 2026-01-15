#include "CFus.h"
#include <cmath>

// the following is needed in the ROOT version
// ClassImp(CFus)

float const CFus::E = 3.3f;
float const CFus::H = 0.65f;
float const CFus::D = 0.03f;
float const CFus::G = 0.0061f;

/**
 * Constructor
\param plb0 is pi-lambda-squared in mb
\param dif0 is the diffuseness in hbar
*/
CFus::CFus(float plb0, float dif0)
{
    plb = plb0;
    dif = dif0;
}

/**
 * alternative constructor
\param iZprojectile is the projectile atomic number
\param iAprojectile is projectile mass number
\param iZtarget is the target atomic number
\param iAtarget is target mass number
\param Elab is lab. energy of projectile in MeV
\param dif0 is the diffuseness in hbar
 */
CFus::CFus(int iZprojectile, int iAprojectile,
           int iZtarget, int iAtarget, float Elab, float dif0)
{
    iZp = iZprojectile;
    iAp = iAprojectile;
    iZt = iZtarget;
    iAt = iAtarget;
    fElab = Elab;

    // mod-TU CMass mass;
    CMass *mass = CMass::instance(); // mod-TU
    float massProjectile = mass->getExpMass(iZprojectile, iAprojectile); // mod-TU
    float massTarget     = mass->getExpMass(iZtarget, iAtarget);         // mod-TU

    iZcn = iZprojectile + iZtarget;
    iAcn = iAprojectile + iAtarget;

    p = CNuclide(iZp, iAp);
    t = CNuclide(iZt, iAt);
    c = CNuclide(iZcn, iAcn);

    float massCN = mass->getExpMass(iZcn, iAcn); // mod-TU

    float Q = massProjectile + massTarget - massCN;
    Qval = Q;

    Ecm   = (float)iAtarget / (float)(iAtarget + iAprojectile) * Elab;
    vbeam = std::sqrt(2.f * Elab / (float)iAprojectile) * 0.9784f;
    vcm   = (float)iAprojectile * vbeam / (float)(iAtarget + iAprojectile);
    Ex    = Ecm + Q;

    float Ared = (float)(iAtarget * iAprojectile) / (float)(iAtarget + iAprojectile);
    plb = 656.80f / Ared / Ecm;
    dif = dif0;
}

//*********************************************
/**
 * reinitialization
\param plb0 is pi-lambda-squared in mb
\param dif0 is the diffuseness in hbar
 */
void CFus::init(float plb0, float dif0)
{
    plb = plb0;
    dif = dif0;
}

//*******************************
/**
 * returns the critical amgular momentum
\param xsection is the fusion cross section in mb
*/
float CFus::getL0(float xsection)
{
    // first get the value for dif=0
    float l00 = std::sqrt(xsection / plb) - 1.f;

    float l0 = l00;
    float dl0;

    for (;;)
    {
        int il = 0;
        float xsec = 0.f;
        float dxsec = 0.f;
        float maxExtra = 0.f;

        for (;;)
        {
            float fl = (float)il;

            float expo  = std::exp((fl - l0) / dif);
            float trans = 1.f / (1.f + expo);

            // pow(trans,2) -> trans*trans
            float trans2 = trans * trans;
            float DtransDl0 = trans2 * expo / dif;

            float extra = trans * (2.f * fl + 1.f) * plb;
            if (maxExtra < extra) maxExtra = extra;

            xsec  += extra;
            dxsec += DtransDl0 * (2.f * fl + 1.f) * plb;

            if (extra < maxExtra / 1000.f) break;
            il++;
        }

        float difference = xsec - xsection;
        if (std::fabs(difference) < 0.1f) break;

        dl0 = -difference / dxsec;
        l0 += dl0;
    }

    return l0;
}

//**************************************************************
/**
 * returns the Bass-model potential in MeV
 * \param R radius in fm
 * \param AL is orbital angular momentum
 */
float CFus::F(float R, float AL)
{
    float S  = R - R12;
    float R2 = R * R;

    float denom = D * std::exp(S / E) + G * std::exp(S / H);
    return A / R + B * AL * (AL + 1.f) / R2 - C / (denom + G); // NOTE: keeps original structure? (see below)
}

//*************************************************
/**
 * returns the derivative Bass-model potential in MeV withrespect to R
 * \param R radius in fm
 * \param AL is orbital angular momentum
 */
float CFus::FF(float R, float AL)
{
    float S  = R - R12;
    float R2 = R * R;
    float R3 = R2 * R;

    float x = -A / R2 - 2.f * B * AL * (AL + 1.f) / R3;

    float e1 = std::exp(S / E);
    float e2 = std::exp(S / H);
    float denom  = D * e1 + G * e2;
    float denom2 = denom * denom;

    float num = (D / E) * e1 + (G / H) * e2;
    return x + C * num / denom2;
}

//********************************************************
/**
 * returns the 2nd derivative of the Bass-model potential
 * in MeV withrespect to R
 * \param R radius in fm
 * \param AL is orbital angular momentum
 */
float CFus::FFF(float R, float AL)
{
    float S  = R - R12;
    float R2 = R * R;
    float R3 = R2 * R;
    float R4 = R2 * R2;

    float x = 2.f * A / R3 + 6.f * B * AL * (AL + 1.f) / R4;

    float e1 = std::exp(S / E);
    float e2 = std::exp(S / H);

    float Y  = D * e1 + G * e2;
    float Y2 = Y * Y;
    float Y3 = Y2 * Y;

    float Z  = (D / E) * e1 + (G / H) * e2;
    float Z2 = Z * Z;

    // GGG = -2*C*pow(Z,2)/pow(Y,3)
    float GGG = -2.f * C * Z2 / Y3;

    // HHH = C*(D/E/E*exp(S/E)+G/H/H*exp(S/H))/pow(Y,2)
    float XX = (D / (E * E)) * e1 + (G / (H * H)) * e2;
    float HHH = C * XX / Y2;

    return x + GGG + HHH;
}

//**********************************************************
/**
 * returns the 3nd derivative of the Bass-model potential
 * in MeV withrespect to R
 * \param R radius in fm
 * \param AL is orbital angular momentum
 */
float CFus::FFFF(float R, float AL)
{
    float S  = R - R12;
    float R2 = R * R;
    float R3 = R2 * R;
    float R4 = R2 * R2;
    float R5 = R4 * R;

    float X = -6.f * A / R4 - 24.f * B * AL * (AL + 1.f) / R5;

    float e1 = std::exp(S / E);
    float e2 = std::exp(S / H);

    float Y  = D * e1 + G * e2;
    float Y2 = Y * Y;
    float Y3 = Y2 * Y;
    float Y4 = Y2 * Y2;

    float Z  = (D / E) * e1 + (G / H) * e2;
    float Z2 = Z * Z;
    float Z3 = Z2 * Z;

    float XX = (D / (E * E)) * e1 + (G / (H * H)) * e2;

    // YY = D/pow(E,3)*exp(S/E)+G/pow(H,3)*exp(S/H)
    float YY = (D / (E * E * E)) * e1 + (G / (H * H * H)) * e2;

    // ZZ = -6*C*pow(Z,3)/pow(Y,4) + 6*C*Z*XX/pow(Y,3) - C*YY/pow(Y,2)
    float ZZ = -6.f * C * Z3 / Y4
               +  6.f * C * Z * XX / Y3
               -  C * YY / Y2;

    return X - ZZ;
}

//*******************************************************
/**
 * returns the critical orbital angular momentum in the Bass model
 * [Nucl Phys A231 (1974) 45 ] USING THE 1977
 * BASS NUCLEAR POTENTIAL [Phys Rev Letts 39 (1977) 265 ]
 */
float CFus::getBassL()
{
    const float AF = 0.75f;

    float AP3 = std::pow((float)iAp, (float)(1.f / 3.f));
    float AT3 = std::pow((float)iAt, (float)(1.f / 3.f));

    float AP = (float)iAp;
    float AT = (float)iAt;

    A = 1.4399f * (float)(iZp * iZt);
    B = 20.90f * (AP + AT) / (AP * AT);

    float AL, R, DR, DDB = 0.f;

    float RP = 1.16f * AP3 - 1.27f / AP3;
    float RT = 1.16f * AT3 - 1.27f / AT3;

    int MIN = 0;
    R12 = RP + RT;
    U = AP * AT / (AP + AT);

    C = RP * RT / R12;

    if (FF(R12, 0.f) < 0.f)
    {
        CL1 = 0.f;
        return 0.f;
    }
    else
    {
        float R12_2 = R12 * R12;
        float R12_3 = R12_2 * R12;

        CL1 = std::sqrt(0.25f + FF(R12, 0.f) * R12_3 / (2.f * B)) - 0.5f;

        if (FFF(R12, CL1) > 0.f)
        {
            // *********** DETERMINATION OF CRITICAL L AT WHICH POCKET VANISHES***
            float R = R12;
            MIN = (int)CL1;

            float AJ = 0.f;
            float FFB = 0.f;

            for (int J = MIN; J <= 200; J++)
            {
                AJ = (float)J;
                for (;;)
                {
                    float DR = FFF(R, AJ) / FFFF(R, AJ);
                    R = R - DR;
                    if (std::fabs(DR) <= 0.00001f) break;
                }
                if (FF(R, AJ) < 0.f) break;
                FFB = FF(R, AJ);
            }
            CL1 = AJ - 1.f + FFB / (FFB - FF(R, AJ));
        }

        // ********* CALCULATION OF FUSION BARRIERS FOR EACH L VALUE**
        MAX = (int)(CL1 + 1.f);
        R = R12 + 2.f;

        for (int J = 1; J <= MAX; J++)
        {
            AL = (float)J - 1.f;
            for (;;)
            {
                DR = FF(R, AL) / FFF(R, AL);
                R = R - DR;
                if (std::fabs(DR) <= 0.00001f) break;
            }
            W[J] = F(R, AL);
        }

        // ********** DETERMINATION OF THE CRITICAL ANGULAR MOMENTUM L1 **
        if (MIN > 0)
        {
            for (int I = MIN; I <= MAX - 1; I++)
            {
                float AJ = (float)I;
                float DD = W[I + 1] - F(R12, AJ);
                if (DD < 0.f)
                {
                    CL1 = AJ - 1.f + DDB / (DDB - DD);
                    break;
                }
                DDB = DD;
            }
        }

        CL2 = CL1 / AF;
        E1 = F(R12, CL1);
        E2 = F(R12, CL2);
    }

    if (Ecm < W[1]) return 0.f;          // below the barrier
    else if (Ecm > E2) return CL2;
    else if (Ecm > E1) return std::sqrt(0.25f + (Ecm - F(R12, 0.f)) * R12 * R12 / B) - 0.5f;
    else
    {
        int i;
        for (i = 1; i < MAX; i++)
        {
            if (W[i] > Ecm) break;
        }
        float DEL = W[i] - W[i - 1];
        return (float)(i - 2) + (Ecm - W[i - 1]) / DEL;
    }
}

//*********************************************************************
/**
 * returns the Bass model fusion cross section in mb
 */
float CFus::getBassXsec()
{
    float L = getBassL();
    return plb * (L * L + 2.f * L);
}
