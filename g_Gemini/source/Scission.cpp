#include "CScission.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;

namespace {
inline double sqr(double x)   { return x * x; }
inline double cube(double x)  { return x * x * x; }
inline double pow4(double x)  { double x2 = x * x; return x2 * x2; }
inline double pow5(double x)  { return x * x * x * x * x; }
inline double pow6(double x)  { double x2 = x * x; return x2 * x2 * x2; }
inline double cbrt_f(double x){ return static_cast<double>(std::cbrt(static_cast<double>(x))); }
} // namespace

double const CScission::slopeViola = .1189f;
double const CScission::constViola = 7.3f;
double const CScission::e2 = 1.44f;
double const CScission::r0 = 1.2f;
double const CScission::kRotate = 41.563f;

/**
 * simple constructor
 */
CScission::CScission()
{
    mass = CMass::instance();//mass singleton
}

/**
 * constructor
/param iZ0 is proton number of fissioning nucleus
/param iA0 is mass number of fission nucleus
/param fJ0 is spin of fissioning nucleus
/param iChan =1 for imf, =2 symmetric fission
 */
CScission::CScission(int iZ0, int iA0, double fJ0, int iChan)
{
    init(iZ0,iA0,fJ0,iChan);
}


//*******************************************************
/**
 * initialized the calls for a given nucleus
  /param iZ0 is the proton number of the nucleus
  /param iA0 is the mass number of the nucleus
  /param fJ0 is the spin of the nucleus
  /param iChan =1 for imf, =2 symmetric fission
 */
void CScission::init(int iZ0, int iA0, double fJ0,int iChan,
                     double Z10/*=0.*/,double A10/*0.*/)
{
    iA = iA0;
    iZ = iZ0;
    fJ = fJ0;
    A = (double)iA;
    Z = (double)iZ;

    if (A10 == 0)
    {
        sym = 1;
        A1 = A/2.f;
        Z1 = Z/2.f;
    }
    else
    {
        sym = 0;
        Z1 = Z10;
        A1 = A10;
    }

    A2 = (double)iA - A1;
    Z2 = (double)iZ - Z1;

    const double A1_13 = cbrt_f(A1);
    const double A2_13 = cbrt_f(A2);

    R1 = A1_13 * r0;
    momInertia1 = 0.4f * A1 * sqr(R1);

    R2 = A2_13 * r0;
    momInertia2 = 0.4f * A2 * sqr(R2);

    Vc0 = Z1*Z2*e2;
    mu0 = A1*A2/(A1+A2);
    k1 = (kRotate * 0.5f) * mu0 * sqr(fJ);

    double Z2A13;
    if (sym == 1)
    {
        const double A_13 = cbrt_f((double)iA);
        Z2A13 = sqr((double)iZ) / A_13;
    }
    else
    {
        Z2A13 = 6.349f * Z1 * Z2 / (A1_13 + A2_13);
    }

    double ekViola = 0.f;
    if (iChan == 1)
    {
        ekViola = slopeViola*Z2A13 + constViola;
        ekViola += 0.004f*sqr(fJ0);   // angular momentum Viola energy
    }
    else
    {
        //alternatibve from Rusanov et al
        if (Z2A13 < 900.f) ekViola = 0.131f*Z2A13;
        else ekViola = 0.104f*Z2A13 + 24.3f;
    }

    //in a two sphere approx, go find separation energy which gives desired
    //fission kinetic energy
    sep = getSep(ekViola);
    sep0 = sep;
}
//******************************************************

/**
 * approximates the scission configuration by two separated spheres.
 * returns the separation between the surface of the sphere which
 * is consistent with input fission total kinetic energy
 /param ekViola is the total fission kinetic energy
 */
double CScission::getSep(double ekViola)
{
    double R = R1 + R2 + 4.f;
    for(;;)
    {
        const double R2v = R * R;
        double momInertiaTot = momInertia1 + momInertia2 + mu0 * R2v;

        double ek = Vc0 / R + k1 * R2v / sqr(momInertiaTot);

        if (std::fabs(ek - ekViola) < 0.1f) break;

        // dek/dR
        // -Vc0/R^2 - 4*mu0*k1*R^3/(I^3) + 2*R/(I^2)
        double dek =
            -Vc0 / R2v
            - 4.f * mu0 * k1 * (R * R2v) / (momInertiaTot * momInertiaTot * momInertiaTot)
            + 2.f * R / sqr(momInertiaTot);

        double dR = -(ek - ekViola) / dek;

        if (R + dR < 0.f) R = 0.9f * R;
        else R += dR;
    }
    return R - R1 - R2;
}
//***************************************************
/**
 * returns the symmetric fission scission energy
 * from a two sphere approximation using the separation sep
 * previously determined from init() to sigmaFissionSystematics
 */
double CScission::getScissionEnergy()
{
    double R = sep + R1 + R2;
    double momInertiaTot = momInertia1 + momInertia2 + mu0 * sqr(R);
    double EE = Vc0 / R + kRotate * 0.5f / momInertiaTot * sqr(fJ);
    double mass1 = mass->getFRM(Z1,A1);
    double mass2 = mass->getFRM(Z2,A2);
    return EE + mass1 + mass2;
}
//***************************************************
/**
 * returns the asymmetric fission scission energy
 * from a two sphere approximation using the separation sep
 * previously determined from init() to sigmaFissionSystematics()
 /param iZ1 is the proton number of one of the fission fragments
 /param iA1 is the mass number of one of the fission fragments
 */
double CScission::getScissionEnergy(int iZ1, int iA1)
{
    A1 = (double)iA1;
    Z1 = (double)iZ1;
    A2 = A - A1;
    int iA2 = iA - iA1;
    Z2 = Z - Z1;
    int iZ2 = iZ - Z1;

    const double A1_13 = cbrt_f(A1);
    const double A2_13 = cbrt_f(A2);

    R1 = r0 * A1_13;
    R2 = r0 * A2_13;

    double R = sep + R1 + R2;
    momInertia1 = 0.4f * A1 * sqr(R1);
    momInertia2 = 0.4f * A2 * sqr(R2);

    double mu = A1*A2/(A1+A2);
    double momInertiaTot = momInertia1 + momInertia2 + mu * sqr(R);

    double EE = Z1*Z2*e2/R + kRotate * 0.5f / momInertiaTot * sqr(fJ);

    double massLD1 = mass->getFRM(iZ1,iA1);
    double massLD2 = mass->getFRM(iZ2,iA2);

    /*  double massLD1 = mass->getLDM(iZ1,iA1);
      double massLD2 = mass->getLDM(iZ2,iA2);*/

    return  massLD1 + massLD2 + EE;
}
//****************************************
/**
 * returns the total fission kinetic energy from a two sphere approximation
 * using the separation sep determined from either
 * init() to sigmaFissionSystematics(), whichever was called last
 /param iZ1 is the proton number of one of the fragments
 /param iZ2 is the mass number of one of the fragments
 */
double CScission::getFissionKineticEnergy(int iZ1, int iA1)
{
    double A1 = (double)iA1;
    double A2 = A - A1;
    double Z1 = (double)iZ1;
    double Z2 = Z - Z1;

    double R1 = r0 * cbrt_f(A1);
    double R2 = r0 * cbrt_f(A2);

    double R = sep + R1 + R2;

    double momInertia1 = 0.4f * A1 * sqr(R1);
    double momInertia2 = 0.4f * A2 * sqr(R2);

    double mu = A1*A2/(A1+A2);
    double momInertiaOrbit = mu * sqr(R);
    double momInertiaTot = momInertia1 + momInertia2 + momInertiaOrbit;

    double fl = fJ * momInertiaOrbit / momInertiaTot;

    EkCoul =  Z1*Z2*e2 / R;
    EkRot  =  kRotate * 0.5f / momInertiaOrbit * sqr(fl);
    ekTot  =  EkCoul + EkRot;

    Erotate1 = sqr(fJ * momInertia1 / momInertiaTot) * (kRotate * 0.5f) / momInertia1;
    Erotate2 = sqr(fJ * momInertia2 / momInertiaTot) * (kRotate * 0.5f) / momInertia2;

    return ekTot;
}
//**************************************************************
/**
 * estimates the standard deviation of the fission mass distributions
 * from the systematics of Rusanov et al. Physics of the Atomic Nucleus 60
 * (1997) 683 assuming a scission-point logic.
 * subsequentally, it approximates the scission configuration as two separated
 * spheres, where the separation is adjusted to reproduce the mass distribution
\param iZ0 is the proton number
\param iA0 is the mass number
\param fJ is the angular momentum
\param fUScission is the thermal excitation energy at the scission-point in MeV
*/
double CScission::sigmaFissionSystematicsScission(int iZ0, int iA0, double fJ,
                                                 double fUScission)
{
    if (fUScission < 0.f || fUScission > 2000.f)
    {
        cout << "fUScission= " << fUScission << " sigmaFissionSystematics" << endl;
        abort();
    }

    iZ = iZ0;
    iA = iA0;
    A = (double)iA;
    Z = (double)iZ;

    const double A3 = cbrt_f(A);

    // on page 684, the temp is determined with a level-density parameter 0.093A
    //we must do the same to be consistent
    double temp = std::sqrt(fUScission / (0.093f * A));

    double Z2A = sqr(Z) / A;

    //find stiffness from Fig8c
    double d2Vdeta2;
    if (Z2A < 23.49f) d2Vdeta2 = 2.105f;
    else if (Z2A < 30.f) d2Vdeta2 = 1.923f*Z2A - 43.08f;
    else if (Z2A < 33.9f) d2Vdeta2 = 3.643f*Z2A - 94.224f;
    else d2Vdeta2 = -1.3144f*Z2A + 73.45f;

    //use equation 1 to get the variance from the stiffness and temp
    double sigma2 = sqr(A) * temp / (16.f * d2Vdeta2);

    //correction for angular momentum from eq 17 and 18.
    double d2sdl2;
    if (Z2A > 32.7f)
        d2sdl2 = -0.1310f*temp - 0.05147f*Z2A + 0.000766f*sqr(Z2A)
                 + 0.00289f*temp*Z2A + .970f;
    else if (Z2A > 31.f)
        d2sdl2 = .2873f*temp + 0.03687f*Z2A
                 - 0.00974f*temp*Z2A - 1.1143f;
    else
        d2sdl2 = 0.0111f*Z2A - .334f;

    if(d2sdl2 < 0.f) d2sdl2 = 0.f;

    double correction = d2sdl2 * sqr(fJ) * 0.5f;

    //I find that use of the full correction - overestimartes the width
    // so I have scaled it
    // correction*= .75;

    sigma2 += correction;

    //now we determine d2VdA2
    double d2VdA2 = temp / sigma2;

    //this has a component of Coloumb energy of each fragment
    double alpha = Z/A;

    const double cbrt2 = cbrt_f(2.f);
    const double pow2_13 = cbrt2;          // 2^(1/3)
    const double pow2_23 = cbrt2*cbrt2;    // 2^(2/3)

    double Ec0 = 0.7053f*sqr(alpha)*pow2_13*20.f/9.f/A3;

    // from the surface energy
    // pow(A,4/3) = A * A^(1/3) = A * A3
    double Es = -8.f*pow2_13*17.9439f
               * (1.f - 1.7826f*sqr((A - 2.f*Z)/A))
               / 9.f / (A * A3);

    //find unaccounted
    d2VdA2 -= Ec0 + Es;

    const double pow2_53 = pow5(cbrt2); // 2^(5/3) = (2^(1/3))^5

    double fact1 = pow2_53 * A3 * e2 * r0 * sqr(alpha) / 9.f;
    double rr    = pow2_23 * A3 * r0;
    double fact2 = 2.f * e2 * sqr(alpha);

    const double r0_2 = r0*r0;
    const double r0_3 = r0_2*r0;
    const double r0_4 = r0_2*r0_2;
    const double r0_5 = r0_4*r0;
    const double r0_6 = r0_3*r0_3;

    const double A2v = A*A;
    const double A3v = A2v*A;
    const double A4v = A2v*A2v;
    const double A5v = A4v*A;

    const double AA3 = A*A3;          // A * A^(1/3) = A^(4/3)
    const double AA3_2 = AA3*AA3;     // (A*A3)^2

    double fact3 = (88.8945f*A*AA3 + 120.f*A2v + 67.5318f*AA3_2 +
                   46.3521f*A3v*A3 - 76.32f*A4v) * r0_4;

    double fact4 = (384.f*A + 453.572f*A*sqr(A3) + 246.365f*A2v*A3 + 64.8f*A3v) * r0_3;

    double fact5 = (635.f*sqr(A3) + 609.562f*A*A3 + 208.8f*A2v) * r0_2;

    double fact6 = (482.57f*A3 + 288.f*A) * r0;
    double fact7 = 144.f;

    double fact8 =
        A5v * r0_6 *
        (144.f + 544.286f*sqr(A3) + 1131.5f*A*A3 + 1411.2f*A2v +
         1167.49f*AA3_2 + 579.465f*A3v*A3 + 158.184f*A4v);

    double fact9 =
        A4v * sqr(A3) * r0_5 *
        (1088.57f + 3428.79f*sqr(A3) + 5702.4f*A*A3 + 5334.f*A2v +
         2941.9f*AA3_2 + 730.08f*A3v*A3);

    double fact10 = pow4(A*r0) * A3 *
                   (3428.79f + 8640.f*sqr(A3) + 6720.42f*A2v + 1853.28f*AA3_2);

    double fact11 =
        A4v * r0_3 *
        (760.f + 10885.7f*sqr(A3) + 9052.f*A*A3 + 2822.4f*A2v);

    double fact12 =
        A3v * sqr(A3) * r0_2 *
        (5442.86f + 6857.57f*sqr(A3) + 2851.2f*A*A3);

    double fact13 = r0 * (2743.03f*A3v*A3 + 1728.f*A4v);
    double fact14 = 576.f*A3v;

    double fact15 = fJ*(fJ+1.f)*kRotate*64.f;

    double s = 1.f;
    int tries = 0;

    for(;;)
    {
        double sp = s + rr;

        double y = fact1 / sqr(sp) - fact2 / sp;

        //contribution from spin
        double nom = fact3 + fact4*s + fact5*s*s + fact6*cube(s) + fact7*pow4(s);
        double denom = fact8 + fact9*s + fact10*s*s + fact11*cube(s) + fact12*pow4(s)
                      + fact13*pow5(s) + fact14*pow6(s);

        double extra = fact15 * nom / denom;
        y += extra;

        double dy = -2.f*fact1/(sp*sp*sp) + fact2/sqr(sp);

        //contribution from spin
        double dnom = fact4 + 2.f*fact5*s + 3.f*fact6*s*s + 4.f*fact7*cube(s);
        double dden = fact9 + 2.f*fact10*s + 3.f*fact11*s*s + 4.f*fact12*cube(s)
                     + 5.f*fact13*pow4(s) + 6.f*fact14*pow5(s);

        extra = fact15 * dnom / denom - fact15 * nom / sqr(denom) * dden;
        dy += extra;

        double delta = y - d2VdA2;
        double deltaS = -delta/dy;

        if (std::fabs(delta) < .0001f) break;

        s += deltaS;

        if (s < 0.f) s = 1.e-3f;
        else if (s > 50.f) s = 50.f;

        tries++;
        if (tries == 10 || std::isnan(s))
        {
            s = 5.f;
            break;
        }
    }

    sep = s;
    sep1 = sep;
    Esymmetric = getScissionEnergy();

    return sigma2;
}
//**************************************************************
/**
 * estimates the standard deviation of the fission mass distributions
 * from the systematics of Rusanov et al. Physics of the Atomic Nucleus 60
 * (1997) 683 assuming a saddle-point logic.
 * subsequentally, it approximates the scission configuration as two separated
 * spheres, where the separation is adjusted to reproduce the mass distribution
\param iZ0 is the proton number
\param iA0 is the mass number
\param fJ is the angular momentum
\param fUScission is the thermal excitation energy at the scission-point in MeV
*/
double CScission::sigmaFissionSystematicsSaddle(int iZ0, int iA0, double fJ,
                                               double fUScission)
{
    if (fUScission < 0.f || fUScission > 2000.f)
    {
        cout << "fUScission= " << fUScission << " sigmaFissionSystematics" << endl;
        abort();
    }

    iZ = iZ0;
    iA = iA0;
    A = (double)iA;
    Z = (double)iZ;

    const double A3 = cbrt_f(A);

    // on page 684, the temp is determined with a level-density parameter 0.093A
    //we must do the same to be consistent
    double temp = std::sqrt(fUScission / (0.093f * A));

    double Z2A = sqr(Z) / A;

    //find stiffness from Fig8c
    double d2Vdeta2;
    if (Z2A < 23.49f) d2Vdeta2 = 2.105f;
    else if (Z2A < 31.57f) d2Vdeta2 = 1.923f*Z2A - 43.08f;
    else if (Z2A < 34.2f) d2Vdeta2 = 3.19f*Z2A - 83.06f;
    else d2Vdeta2 = -1.7287f*Z2A + 85.42f;

    //use equation 1 to get the variance from the stiffness and temp
    double sigma2 = sqr(A) * temp / (16.f * d2Vdeta2);

    //correction for angular momentum from eq 17 and 18.
    double d2sdl2;
    if (Z2A > 32.7f)
        d2sdl2 = -0.1310f*temp - 0.05147f*Z2A + 0.000766f*sqr(Z2A)
                 + 0.00289f*temp*Z2A + .970f;
    else if (Z2A > 31.f)
        d2sdl2 = .2873f*temp + 0.03687f*Z2A
                 - 0.00974f*temp*Z2A - 1.1143f;
    else
        d2sdl2 = 0.0111f*Z2A - .334f;

    if(d2sdl2 < 0.f) d2sdl2 = 0.f;

    double correction = d2sdl2 * sqr(fJ) * 0.5f;

    //I find that use of the full correction - overestimartes the width
    // so I have scaled it
    // correction*= .75;

    sigma2 += correction;

    //now we determine d2VdA2
    double d2VdA2 = temp / sigma2;

    //this has a component of Coloumb energy of each fragment
    double alpha = Z/A;

    const double cbrt2 = cbrt_f(2.f);
    const double pow2_13 = cbrt2;
    const double pow2_23 = cbrt2*cbrt2;

    double Ec0 = 0.7053f*sqr(alpha)*pow2_13*20.f/9.f/A3;

    // from the surface energy
    double Es = -8.f*pow2_13*17.9439f
               * (1.f - 1.7826f*sqr((A - 2.f*Z)/A))
               / 9.f / (A * A3);

    //find unaccounted
    d2VdA2 -= Ec0 + Es;

    const double pow2_53 = pow5(cbrt2);

    double fact1 = pow2_53 * A3 * e2 * r0 * sqr(alpha) / 9.f;
    double rr    = pow2_23 * A3 * r0;
    double fact2 = 2.f * e2 * sqr(alpha);

    const double r0_2 = r0*r0;
    const double r0_3 = r0_2*r0;
    const double r0_4 = r0_2*r0_2;
    const double r0_5 = r0_4*r0;
    const double r0_6 = r0_3*r0_3;

    const double A2v = A*A;
    const double A3v = A2v*A;
    const double A4v = A2v*A2v;
    const double A5v = A4v*A;

    const double AA3 = A*A3;
    const double AA3_2 = AA3*AA3;

    double fact3 = (88.8945f*A*AA3 + 120.f*A2v + 67.5318f*AA3_2 +
                   46.3521f*A3v*A3 - 76.32f*A4v) * r0_4;

    double fact4 = (384.f*A + 453.572f*A*sqr(A3) + 246.365f*A2v*A3 + 64.8f*A3v) * r0_3;

    double fact5 = (635.f*sqr(A3) + 609.562f*A*A3 + 208.8f*A2v) * r0_2;

    double fact6 = (482.57f*A3 + 288.f*A) * r0;
    double fact7 = 144.f;

    double fact8 =
        A5v * r0_6 *
        (144.f + 544.286f*sqr(A3) + 1131.5f*A*A3 + 1411.2f*A2v +
         1167.49f*AA3_2 + 579.465f*A3v*A3 + 158.184f*A4v);

    double fact9 =
        A4v * sqr(A3) * r0_5 *
        (1088.57f + 3428.79f*sqr(A3) + 5702.4f*A*A3 + 5334.f*A2v +
         2941.9f*AA3_2 + 730.08f*A3v*A3);

    double fact10 = pow4(A*r0) * A3 *
                   (3428.79f + 8640.f*sqr(A3) + 6720.42f*A2v + 1853.28f*AA3_2);

    double fact11 =
        A4v * r0_3 *
        (760.f + 10885.7f*sqr(A3) + 9052.f*A*A3 + 2822.4f*A2v);

    double fact12 =
        A3v * sqr(A3) * r0_2 *
        (5442.86f + 6857.57f*sqr(A3) + 2851.2f*A*A3);

    double fact13 = r0 * (2743.03f*A3v*A3 + 1728.f*A4v);
    double fact14 = 576.f*A3v;

    double fact15 = fJ*(fJ+1.f)*kRotate*64.f;

    double s = 1.f;
    int tries = 0;

    for(;;)
    {
        double sp = s + rr;

        double y = fact1 / sqr(sp) - fact2 / sp;

        //contribution from spin
        double nom = fact3 + fact4*s + fact5*s*s + fact6*cube(s) + fact7*pow4(s);
        double denom = fact8 + fact9*s + fact10*s*s + fact11*cube(s) + fact12*pow4(s)
                      + fact13*pow5(s) + fact14*pow6(s);

        double extra = fact15 * nom / denom;
        y += extra;

        double dy = -2.f*fact1/(sp*sp*sp) + fact2/sqr(sp);

        //contribution from spin
        double dnom = fact4 + 2.f*fact5*s + 3.f*fact6*s*s + 4.f*fact7*cube(s);
        double dden = fact9 + 2.f*fact10*s + 3.f*fact11*s*s + 4.f*fact12*cube(s)
                     + 5.f*fact13*pow4(s) + 6.f*fact14*pow5(s);

        extra = fact15 * dnom / denom - fact15 * nom / sqr(denom) * dden;
        dy += extra;

        double delta = y - d2VdA2;
        double deltaS = -delta/dy;

        if (std::fabs(delta) < .0001f) break;

        s += deltaS;

        if (s < 0.f) s = 1.e-3f;
        else if (s > 50.f) s = 50.f;

        tries++;
        if (tries == 10 || std::isnan(s))
        {
            s = 5.f;
            break;
        }
    }

    sep = s;
    sep1 = sep;
    Esymmetric = getScissionEnergy();

    return sigma2;
}
