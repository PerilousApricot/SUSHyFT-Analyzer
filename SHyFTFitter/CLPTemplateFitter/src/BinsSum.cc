// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "CLPTemplateFitter/interface/CLPTemplateFitter.h"
#include "CLPTemplateFitter/interface/BinsSum.h"
#include "CLPTemplateFitter/interface/BinsSumCont.h"
#include "CLPTemplateMorph/interface/dumpSTL.icc"
#include "OptionUtils/interface/dout.h"

using namespace std;

BinsSum::BinsSum (std::string name, BinsSumCont *contPtr) :
   m_name (name), m_contPtr (contPtr)
{
}

double
BinsSum::value() const
{
   double total = 0.;
   assert (m_lowerBinVec.size() == m_upperBinVec.size());
   int size = (int) m_lowerBinVec.size();
   for (IVecConstIter tempIter = m_templateIndiciesVec.begin();
        m_templateIndiciesVec.end() != tempIter;
        ++tempIter)
   {
      TH1F *histPtr = m_contPtr->getTemplate( *tempIter );
      for (int loop = 0; loop < size; ++loop)
      {
         total += histPtr->Integral( m_lowerBinVec.at(loop),
                                     m_upperBinVec.at(loop) );
      }
   } // for tempIter
   return total;
}

void 
BinsSum::calcMeanAndRMS (double &mean, double &rms) const
{
   calcMeanAndRMS (mean, rms, m_trialsVec);
}

int
BinsSum::addSelfToContainer()
{
   return m_contPtr->add (*this);
}

void 
BinsSum::calcMeanAndRMS (double &mean, double &rms, const DVec &xvec)
{
   double sumX = 0.;
   double sumX2 = 0.;
   int num = 0;
   for (DVecConstIter iter = xvec.begin();
        xvec.end() != iter;
        ++iter)
   {
      ++num;
      sumX  += *iter;
      sumX2 += *iter * *iter;
   }
   if (! num)
   {
      // we haven't summed over anything
      mean = -999.;
      rms  = -999.;
      return;
   }
   mean = sumX / num;
   double aveX2 = sumX2 / num;
   rms = sqrt (aveX2 - mean * mean);
}


// friends
ostream& operator<< (ostream& o_stream, const BinsSum &rhs)
{
   o_stream << "name " << rhs.m_name << endl;
   // PLEASE come back and fix this!!!
   // dumpNamedSTLendl ((const string)"templates", rhs.m_templateIndiciesVec, o_stream);
   // dumpNamedSTLendl ((const string)"lower    ", rhs.m_lowerBinVec, o_stream);
   // dumpNamedSTLendl ((const string)"upper    ", rhs.m_upperBinVec, o_stream);
   return o_stream;
} 
