// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "CLPTemplateMorph/interface/CLPBin.h"

using namespace std;

ClassImp (CLPBin)

CLPBin::CLPBin() : m_xUpper (0), m_xWidth (0), m_y(0), m_sameAsLastBin (false)
{
}

CLPBin::CLPBin (double xLower, double xWidth, double yValue, bool sameAs) : 
   m_xUpper (xLower + xWidth), m_xWidth (xWidth), m_y(yValue), 
   m_sameAsLastBin (sameAs)
{
}

// friends
ostream& operator<< (ostream& o_stream, const CLPBin &rhs)
{
   o_stream << Form ("[(%7.3f, %7.3f), %9.6f]", rhs.xLower(), rhs.xUpper(),
                     rhs.yValue());
   return o_stream;
} 
