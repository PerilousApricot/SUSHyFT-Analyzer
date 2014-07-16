// -*- C++ -*-
#include <iostream>
#include <cassert>

#include "PolyNoid.h"
#include "OptionUtils/interface/dout.h"
#include "MRFitterNamespace.h"
#include "CLPTemplateMorph/interface/dumpSTL.icc"

using namespace std;

PolyNoid::IVec   PolyNoid::sm_groupBinBoundariesVec;
int              PolyNoid::sm_lastBin = 1;
int              PolyNoid::sm_lastGroupBin = 0;
int              PolyNoid::sm_groupBinUpper = 0;

PolyNoid::PolyNoid() :
   m_name (""),   m_mean (0.),   m_width (0.), 
   m_lower (0.),  m_upper (0.),  m_step (0.)
{
}

PolyNoid::PolyNoid (std::string &name, double mean, double width,
                    double lower, double upper, double step) :
   m_name (name), m_mean (mean), m_width (width), 
   m_lower (lower),  m_upper (upper),  m_step (step)
{
   if (m_lower == m_upper)
   {
      m_lower = -10. * m_mean;
      m_upper =  10. * m_mean;
      m_step  = 0.001;
      if (0 == lower)
      {
         m_lower = -10.;
         m_upper =  10.;
      } // if no lower value
   } // if 
}

double
PolyNoid::operator() (int tempIndex, int binIndex,
                      double value, const double *paramArray) const
{
   int groupBin = bin2groupBin (binIndex);
   IIPairIMapConstIter iter = 
      m_groupTemplateParamMap.find( make_pair( groupBin, tempIndex ) );
   assert (m_groupTemplateParamMap.end() != iter);
   return polynomial ( value, m_parametersVec.at( iter->second ) );
}

//////////////////////////////
// Regular Member Functions //
//////////////////////////////

void 
PolyNoid::addPolynomial (const IVec &groupVec,
                         const IVec &templateEnumVec,
                         const DVec &parameters)
{
   int index = (int) m_parametersVec.size();
   m_parametersVec.push_back( parameters );
   for (IVecConstIter groupIter = groupVec.begin();
        groupVec.end() != groupIter;
        ++groupIter)
   {
      for (IVecConstIter templateIter = templateEnumVec.begin();
           templateEnumVec.end() != templateIter;
           ++templateIter)
      {
         m_temporaryMap[ make_pair( *groupIter, 
                                    *templateIter ) ] = index;
      } // for templateIter
   } // for groupIter
}

void  
PolyNoid::calculateGroupTemplateParamMap (const IIMap &template2indexMap)
{
   const IIMapConstIter kEndIter = template2indexMap.end();
   for (IIPairIMapConstIter iter = m_temporaryMap.begin();
        m_temporaryMap.end() != iter;
        ++iter)
   {
      int groupBin    = iter->first.first;
      int templateEnum = iter->first.second;
      IIMapConstIter findIter = template2indexMap.find (templateEnum);
      if (kEndIter == findIter)
      {
         // not found.  not good.
         cout << "groupBin     " << groupBin << endl;
         cout << "templateEnum " << templateEnum << endl;
         assert (0);
      } // not found
      m_groupTemplateParamMap[ make_pair( groupBin,
                                          findIter->second ) ] = iter->second;
   } // for iter
}

int 
PolyNoid::addSelfToFitter (CLPTemplateFitter &fitter)
{
   BinNormClass::ISet templateSet, binSet;
   for (IIPairIMapConstIter iter = m_groupTemplateParamMap.begin();
        m_groupTemplateParamMap.end() != iter;
        ++iter)
   {
      // just for clarification
      int groupBin = iter->first.first;
      int tempIndex = iter->first.second;
      templateSet.insert (tempIndex);
      int lower, upper;
      lowerUpperBoundaries (lower, upper, groupBin);
      BinNormClass::insertRangeIntoSet (binSet, lower, upper);
   } // for iter
   // dumpSTL (binSet);
   // cout << endl;     
   // Make a bin norm
   BinNormClass tempNorm (templateSet,
                          binSet,
                          mean(),
                          lower(),
                          upper(),
                          step());
   // Add it to the fitter
   int index = fitter.addBinNorm (name(), tempNorm,
                                  mean(), width());
   // Now add the PolyNoiid functoid
   fitter.setBinNormFunctoid (index, *this);
   return index;
}


/////////////////////////////
// Static Member Functions //
/////////////////////////////


void 
PolyNoid::setGroupBinWidths (const IVec &widthVec)
{
   sm_groupBinBoundariesVec.clear();
   sm_groupBinBoundariesVec.push_back( 0 );
   int total = 0;
   for (IVecConstIter iter = widthVec.begin();
        widthVec.end() != iter;
        ++iter)
   {
      total += *iter;
      sm_groupBinBoundariesVec.push_back( total );
   } // for iter
   sm_groupBinUpper = (int) sm_groupBinBoundariesVec.size() - 1;
}


double
PolyNoid::polynomial (double x, const DVec &params)
{
   Double_t xToTheN = 1.;
   Double_t total = 0.;
   unsigned int size = params.size();
   const DVecConstIter kEnd = params.end();
   for (DVecConstIter iter = params.begin();
        kEnd != iter;
        ++iter)
   {
      total += xToTheN * (*iter);
      xToTheN *= x;
   } // for iter
   return total;
}

int 
PolyNoid::bin2groupBin (int bin)
{
   if (sm_lastBin == bin)
   {
      return sm_lastGroupBin;
   }
   if (bin > sm_lastBin)
   {
      for (int index = sm_lastGroupBin;
           index <= sm_groupBinUpper;
           ++index)
      {
         if (sm_groupBinBoundariesVec.at (index)     <  bin &&
             sm_groupBinBoundariesVec.at (index + 1) >=   bin)
         {
            // got it
            sm_lastGroupBin = index;
            sm_lastBin = bin;
            return index;
         }
      } // for index
      assert (0);
   } else {
      for (int index = sm_lastGroupBin;
           index >= 0;
           --index)
      {
         if (sm_groupBinBoundariesVec.at (index)     <  bin &&
             sm_groupBinBoundariesVec.at (index + 1) >=   bin)
         {
            // got it
            sm_lastGroupBin = index;
            sm_lastBin = bin;
            return index;
         }
      } // for index
      assert (0);
   } // else (smaller bin)
   return -1;
}

void 
PolyNoid::lowerUpperBoundaries (int &lower, int &upper, int groupBin)
{
   // sm_groupBinUpper should already be set before this is called
   assert (sm_groupBinUpper);
   assert (groupBin >= 0 && groupBin < sm_groupBinBoundariesVec.size());
   lower = sm_groupBinBoundariesVec.at( groupBin ) + 1;
   upper = sm_groupBinBoundariesVec.at( groupBin + 1 );
}

const PolyNoid::IVec &
PolyNoid::groupBinBoundariesVec() 
{ 
   // sm_groupBinUpper should already be set before this is called
   assert (sm_groupBinUpper);
   return sm_groupBinBoundariesVec; 
}

int
PolyNoid::lastBin()
{
   unsigned int size = sm_groupBinBoundariesVec.size();
   if (! size)
   {
      return -1;
   }
   return sm_groupBinBoundariesVec.at (size - 1);
}


void
PolyNoid::dumpToOstream (ostream& o_stream) const
{
   o_stream << name() << " : " 
            << Form ("%.2f +- %.2f, [%.2f, %.2f] %.2f",
                     mean(), width(), 
                     lower(), upper(), step())
            << endl;
   for (PolyNoid::IIPairIMapConstIter iter = 
           m_groupTemplateParamMap.begin();
        m_groupTemplateParamMap.end() != iter;
        ++iter)
   {
      o_stream << "   [" << iter->first.first << ", " << iter->first.second
               << "] : ";
      dumpSTL (m_parametersVec.at (iter->second), o_stream);
      o_stream << endl;               
   } // for iter
}

ostream& operator<< (ostream& o_stream, const PolyNoid &rhs)
{
   rhs.dumpToOstream (o_stream);
   return o_stream;
}
