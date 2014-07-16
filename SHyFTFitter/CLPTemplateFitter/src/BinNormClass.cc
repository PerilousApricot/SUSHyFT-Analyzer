// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "CLPTemplateFitter/interface/BinNormClass.h"

using namespace std;

const double BinNormClass::kMinValue = 1e-10;


BinNormClass::BinNormClass () : 
   m_value (0.), m_paramIndex (kNegativeNumber), 
   m_startStopStep(),
   m_verbose (false),
   m_functoidPtr (0)
{
}

BinNormClass::BinNormClass (double value) : 
   m_value (value), m_paramIndex (kNegativeNumber),  
   m_startStopStep(),
   m_verbose (false),
   m_functoidPtr (0)
{
   double start = -1 * fabs (value);
   double stop  = 10 * fabs (value);
   double step  = fabs (value) / 1e2;
   m_startStopStep.setStartStopStep (start, stop, step);
}

BinNormClass::BinNormClass (double value, 
                            double start, double stop, double step) : 
   m_value (value), m_paramIndex (kNegativeNumber),  
   m_startStopStep(start, stop, step),
   m_verbose (false),
   m_functoidPtr (0)
{
}

BinNormClass::BinNormClass (const ISet &tempSet, const ISet &binSet, 
                            double value) : 
   m_binSet (binSet), 
   m_value (value), m_paramIndex (kNegativeNumber),  
   m_startStopStep(),
   m_verbose (false),
   m_functoidPtr (0)
{
   double start = -1 * fabs (value);
   double stop  = 10 * fabs (value);
   double step  = fabs (value) / 1e2;
   m_startStopStep.setStartStopStep (start, stop, step);
   setTempSet (tempSet);
}

BinNormClass::BinNormClass (const ISet &tempSet, const ISet &binSet, 
                            double value, 
                            double start, double stop, double step) : 
   m_binSet (binSet), 
   m_value (value), m_paramIndex (kNegativeNumber),  
   m_startStopStep(start, stop, step),
   m_verbose (false),
   m_functoidPtr (0)
{
   setTempSet (tempSet);
}

int
BinNormClass::index (int tempIndex, int binIndex) const
{
   // if either set does NOT contain their respective index, return
   // kNegativeNumber
   if ( (m_templateSet.end() == m_templateSet.find (tempIndex) ) ||
        (m_binSet.end()      == m_binSet.find      (binIndex)  ) )
   {
      return kNegativeNumber;
   } else {
      return m_paramIndex;
   }
}

double 
BinNormClass::factor (int tempIndex, int binIndex, 
                      double value, const double *paramArray) const
{
   // Should I be looking at this BinNorm object for this bin and this
   // template?
   if (index (tempIndex, binIndex) < 0)
   {
      // nothing for us here.  Return 1. since this is multiplicitive
      return 1.;
   }
   // Is there a single functoid pointer for all of the templates?
   if (m_functoidPtr)
   {
      return m_functoidPtr->operator() (tempIndex, binIndex, 
                                        value, paramArray);
   }
   BinNormFunctoid::IFMapConstIter iter =
      m_functoidMap.find (tempIndex);
   if (m_functoidMap.end() == iter)
   {
      return value;
   }
   const BinNormFunctoid *functoidPtr = iter->second;
   return functoidPtr->operator() (tempIndex, binIndex, value, paramArray);
}

void 
BinNormClass::insertRangeIntoSet (BinNormClass::ISet &theSet,
                                  int lower, int upper)
{
   for (int loop = lower; loop <= upper; ++loop)
   {
      theSet.insert (loop);
   }
}

void
BinNormClass::outputISet (const ISet &theset, ostream& o_stream)
{
   for (BinNormClass::ISetConstIter iter = theset.begin();
        theset.end() != iter;
        ++iter)
   {
      o_stream << " " << setw(2) << *iter;
   } // for iter
}

void 
BinNormClass::setTempSet (const ISet &tempSet) 
{ 
   m_templateSet = tempSet; 
}

void 
BinNormClass::setFunctoid (BinNormFunctoid &functoid)
{
   m_functoidPtr = &functoid;
}

void 
BinNormClass::setFunctoid (int tempIndex, BinNormFunctoid &functoid)
{
   m_functoidMap[tempIndex] = &functoid;
}

// friends
ostream& operator<< (ostream& o_stream, const BinNormClass &rhs)
{
   o_stream << "template set:";
   BinNormClass::outputISet (rhs.m_templateSet, o_stream);
   o_stream << " bin set:";
   BinNormClass::outputISet (rhs.m_binSet, o_stream);
   o_stream << " : " << rhs.m_paramIndex
            << " : " << rhs.m_startStopStep;
   return o_stream;
} 
