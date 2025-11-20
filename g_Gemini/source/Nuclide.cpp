#include "CNuclide.h"

// mod-TU CMass CNuclide::mass;

const char*  CNuclide::name[111]={"n","H","He","Li","Be","B",
                              "C","N","O","F","Ne",
			      "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
			      "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
			      "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
                              "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In",
                              "Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr",
			      "Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
			       "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir",
                              "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
                              "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am",
				  "Cm","Bk","Cf","Es","Fm","Md","No","Lr",
                                  "Rf","Db","Sg","Bh","Hs","Mt","Ds"};

//*******************************************
/**
 * Constructor specifies the isotope
\param iZ0 is the proton number
\param iA0 is the mass number
*/
CNuclide::CNuclide(int iZ0, int iA0)
{
  //constructor
  mass = CMass::instance(); // mod-TU
  ran = CRandom::instance();
  init(iZ0,iA0);
}

CNuclide::CNuclide()  // mod-TU
{
    mass = CMass::instance(); // mod-TU
    ran = CRandom::instance();
}

//*******************************************************
/**
 * Initializes the isotope 
 *
 * Can be used to change the isotope
 \param iZ0 is the proton number
 \param iA0 is the mass number
*/
void CNuclide::init(int iZ0, int iA0)
{
  //initialized the charge and mass
  iZ = iZ0;
  iA = iA0;
  iN = iA - iZ;
 fExpMass = mass->getExpMass(iZ,iA);   // mod-TU   mass now pointer

string sGCN;

       if (iZ == 0 && iA == 1) { sGCN = strChemName = "n";}
  else if (iZ ==1 && iA == 1)  { sGCN = strChemName = "p";}
  else if (iZ == 1 && iA == 2) { sGCN = strChemName = "d";}
  else if (iZ == 1 && iA == 3) { sGCN = strChemName = "t";}
  else if (iZ>110) {
       ostringstream outstring;
       outstring << iA << "X-" << iZ;
       sGCN = strChemName = outstring.str();
      }
  else
    {
       ostringstream outstring;
       outstring << iA;
       strChemName = outstring.str() + string(name[iZ]);
       sGCN = "<sup>" + outstring.str() + " </sup>" +  string(name[iZ]);
       }
 // cout << "fExpMass = " << fExpMass << std::endl;

strName = strChemName;
sGCN = "<font size=+1><font color=\"darkBlue\">" + sGCN + "</font></font>";

strGChemName = QString(sGCN.c_str());
}
//**********************************************************
/**
 * Returns the excess mass of the nuclide
 */
float CNuclide::getExcessMass()
{
  fExpMass = mass->getExpMass(iZ,iA); // no this line in new original source
  return fExpMass;
}
//**********************************************************
/**
 * alternative constructor
 */
CNuclide::CNuclide(int iZ0, int iA0, string strName0)
{
    mass = CMass::instance(); // mod-TU
    ran = CRandom::instance();
    strName = strName0;

    init(iZ0,iA0);
   // cout << iZ0 << " " << iA0 << " " << fExpMass << " " << strChemName << std::endl;
}
//**********************************************************
/**
 * Returns the chemical name of the isotope as a character string
 */
const char * CNuclide::getSymbol()
{
  return strChemName.c_str();
}
//********************************************
/**
 * Returns the chemical name of the isotope as a string
 */
string CNuclide::getName()  {  return strChemName;}
QString CNuclide::getGName() {  return strGChemName;}
