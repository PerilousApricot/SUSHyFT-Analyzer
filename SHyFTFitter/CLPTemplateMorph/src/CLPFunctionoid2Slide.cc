// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include "CLPTemplateMorph/interface/CLPFunctionoid2Slide.h"

using namespace std;

CLPFunctionoid2Slide::CLPFunctionoid2Slide(double percentage) :
m_percentage (percentage)
{
}

double
CLPFunctionoid2Slide::operator() (double alpha, double beta) const
{
   return alpha + m_percentage * (beta - alpha);
}

// friends
ostream& operator<< (ostream& o_stream, const CLPFunctionoid2Slide &rhs)
{
   return o_stream;
} 
