// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "CLPTemplateFitter/interface/CLPTemplateFitter.h"
#include "CLPTemplateFitter/interface/BinsSumCont.h"
#include "CLPTemplateMorph/interface/dumpSTL.icc"
#include "OptionUtils/interface/dout.h"

using namespace std;

BinsSumCont::BinsSumCont()
{
}

void
BinsSumCont::cleanHistMap()
{
   for (ITH1FPtrMapIter iter = m_histMap.begin();
        m_histMap.end() != iter;
        ++iter)
   {
      delete iter->second;
   } // for iter
   m_histMap.clear();
}

void 
BinsSumCont::setParamVec (BinsSum::DVec &paramVec)
{
   cleanHistMap();
   m_paramVec = paramVec;
}

int 
BinsSumCont::add (const BinsSum &obj)
{
   int size = (int) m_binsSumVec.size();
   m_binsSumVec.push_back( obj );
   return size;
}

void  
BinsSumCont::recordTrials()
{
   for (BinsSum::VecIter iter = m_binsSumVec.begin();
        m_binsSumVec.end() != iter;
        ++iter)
   {
      iter->recordTrial();
   } // for iter
}

void  
BinsSumCont::reset()
{
   for (BinsSum::VecIter iter = m_binsSumVec.begin();
        m_binsSumVec.end() != iter;
        ++iter)
   {
      iter->reset();
   } // for iter
}

void  
BinsSumCont::getValues (DVec &valuesVec)
{
   valuesVec.clear();
   for (BinsSum::VecIter iter = m_binsSumVec.begin();
        m_binsSumVec.end() != iter;
        ++iter)
   {
      valuesVec.push_back( iter->value() );
   } // for iter
}

void  
BinsSumCont::fillMeansAndRMSs (DVec &meanVec, DVec &rmsVec)
{
   meanVec.clear();
   rmsVec.clear();
   for (BinsSum::VecIter iter = m_binsSumVec.begin();
        m_binsSumVec.end() != iter;
        ++iter)
   {
      double mean, rms;
      iter->calcMeanAndRMS (mean, rms);
      meanVec.push_back (mean);
      rmsVec.push_back  (rms);
   } // for iter
}

void  
BinsSumCont::getAllNames (SVec &nameVec)
{   
   nameVec.clear();
   for (BinsSum::VecIter iter = m_binsSumVec.begin();
        m_binsSumVec.end() != iter;
        ++iter)
   {
      nameVec.push_back( iter->name() );
   } // for iter
}

BinsSum & 
BinsSumCont::createNew (const std::string &name)
{
   BinsSum object (name, this);
   unsigned int index = m_binsSumVec.size();
   m_binsSumVec.push_back( object );
   return m_binsSumVec.at( index );
}

TH1F* 
BinsSumCont::getTemplate (int index) const
{
   if (index < 0) index = kAllTemplates;
   // Do we have it already?
   ITH1FPtrMapConstIter iter = m_histMap.find (index);
   if (m_histMap.end() != iter)
   {
      return iter->second;
   }
   // If we're here, we need to get it
   TH1F *histPtr;
   if (kAllTemplates == index)
   {
      histPtr = m_fitterPtr->updatedHistogram (m_paramVec);
   } else {
      histPtr = m_fitterPtr->updatedHistogram (index, m_paramVec);
   }
   m_histMap[index] = histPtr;
   return histPtr;
}

BinsSum &
BinsSumCont::operator() (int index)
{
   assert (index >= 0 && index < (int) m_binsSumVec.size());
   return m_binsSumVec.at (index);
}


// friends
ostream& operator<< (ostream& o_stream, const BinsSumCont &rhs)
{
   for (BinsSum::VecConstIter iter = rhs.m_binsSumVec.begin();
        rhs.m_binsSumVec.end() != iter;
        ++iter)
   {
      o_stream << *iter << endl;
   }
   return o_stream;
} 
