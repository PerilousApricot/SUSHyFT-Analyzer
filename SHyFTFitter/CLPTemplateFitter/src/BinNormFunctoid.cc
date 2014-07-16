// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "CLPTemplateFitter/interface/BinNormFunctoid.h"

using namespace std;

ClassImp (BinNormFunctoid)

BinNormFunctoid::BinNormFunctoid () : m_funcPtr (0)
{
}

BinNormFunctoid::BinNormFunctoid (FuncPtr ptr) : m_funcPtr (ptr)
{
}

double
BinNormFunctoid::operator() (int tempIndex, int binIndex,
                             double value, const double *paramArray) const
{
   assert (m_funcPtr);
   return (*m_funcPtr) (tempIndex, binIndex, value, paramArray);
}

void
BinNormFunctoid::dumpToOstream (ostream& o_stream) const
{
   o_stream << "old" << endl;
   o_stream << (void*) m_funcPtr << endl;
}

// friends
ostream& operator<< (ostream& o_stream, const BinNormFunctoid &rhs)
{
   rhs.dumpToOstream (o_stream);
   return o_stream;
} 
