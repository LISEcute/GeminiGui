// -*- mode: c++ -*- 
//
#ifndef angle_
#define angle_

#include <cmath>

//******ROOT*********
  //#include "Rtypes.h"

/**
 *!\brief polar angles
 *
 * Class to deal with polar angles
 */

class CAngle
{
protected:

 public:
  static double const pi; //!< 3.14159
  double theta; //!< polar angle in radians
  double phi; //!< azimuth angle in radians
  CAngle(double,double);
  CAngle(){};
  static CAngle transform(CAngle angle1,CAngle angle2);

  //*****ROOT************
    //ClassDef(CAngle,1) //Gemini CAngle
};


#endif
