// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "TCanvas.h"

#include "CLPTemplateMorph/interface/CLPTemplateMorph.h"
#include "CLPTemplateMorph/interface/dout.h"
#include "CLPTemplateMorph/interface/dumpSTL.icc"

using namespace std;

ClassImp (CLPTemplateMorph)

CLPTemplateMorph::CLPTemplateMorph() : TH1F(), 
   m_dontRemovePoints (true), m_smoothDeltas (false), m_outputAll (true),
   m_isMorphContainer (false), m_noUnderOverflow (false)
{
}

CLPTemplateMorph::CLPTemplateMorph(const char* name, const char *title) :
   TH1F(name, title, 10, 0, 1), 
   m_dontRemovePoints (true), m_smoothDeltas (false), m_outputAll (true),
   m_isMorphContainer (false), m_noUnderOverflow (false)
{
}

CLPTemplateMorph::CLPTemplateMorph(const char* name, const char *title,
                        Int_t numBins, Axis_t xlow, Axis_t xhigh) : 
   TH1F(name, title, numBins, xlow, xhigh), 
   m_dontRemovePoints (true), m_smoothDeltas (false), m_outputAll (true),
   m_isMorphContainer (false), m_noUnderOverflow (false)
{
}

void
CLPTemplateMorph::addDefaultTemplate (TH1F *histPtr)
{
   if (isMorphContainer())
   {
      cerr << "CLPTemplateMorph::addDefaultTemplate() "
           <<"should not be called on container morphs.  Aborting." << endl;
      return;
   }
   if (m_defaultBC.size())
   {
      // this isn't the first time...
      return;
   }
   m_defaultBC.integrateTH1F (histPtr);
   recalculate();
}

int
CLPTemplateMorph::addMorphVariable (const string &name, double center)
{
   if (isMorphContainer())
   {
      cerr << "CLPTemplateMorph::addMorphVariable() "
           <<"should not be called on container morphs.  Aborting." << endl;
      return -1;
   }
   int retval = (int) m_paramVec.size();
   m_paramVec.push_back( center ); // start this off with no morph by default
   m_centerVec.push_back( center );
   m_nameVec.push_back( name );
   m_indexVec.push_back( -1 );
   m_deltaMorphMapVec.push_back( CLPBinCont::DMap() );
   return retval;
}

void
CLPTemplateMorph::addMorphTemplate (int index, double value, TH1F *histPtr)
{
   if (isMorphContainer())
   {
      cerr << "CLPTemplateMorph::addMorphTemplate() "
           <<"should not be called on container morphs.  Aborting." << endl;
      return;
   }
   if (! m_defaultBC.size())
   {
      // you need to have defined a default template first
      return;
   }
   // make sure we have a valid index
   assert (index >= 0 && index <= size());
   // make sure this is not too close to the center value
   assert (! CLPBinCont::areTheSame (value, m_centerVec.at (index)) );
   CLPBinCont::DMap &dmap = m_deltaMorphMapVec.at (index);
   CLPBinCont cont (histPtr);
   CLPBinCont delta;
   delta.horizontalSubtract (m_defaultBC, cont);
   if (m_smoothDeltas)
   {
      delta.smoothBeginningOfDelta();
   }
   dmap[value] = delta;
}

void
CLPTemplateMorph::addMorphTemplate (const string &name, double value, 
                                    TH1F *histPtr)
{
   int index = morphIndex (name);
   if (index >= 0)
   {
      addMorphTemplate (index, value, histPtr);
      return;
   }
   // if we're still here, then we messed up
   cerr << "CLPTemplateMorph::addMorphTemplate() error: '" 
        << name << "' not found.  Aborting." << endl;
   assert (0);
}

void
CLPTemplateMorph::setParamVec (const DVec &paramVec)
{
   assert (paramVec.size() >= m_paramVec.size());
   // it's exactly the right size
   if (paramVec.size() == m_paramVec.size())
   {
      m_paramVec = paramVec;
   } else {
      // Copy the first N entries and skip the rest
      for (int loop = 0; loop < size(); ++loop) 
      {
         m_paramVec[loop] = paramVec.at (loop);
      } // for loop
   }
   if (isMorphContainer())
   {
      for (PtrVecIter iter = m_ptrVec.begin();
           m_ptrVec.end() != iter;
           ++iter)
      {
         (*iter)->setParamVec (m_paramVec);
      } // for iter
   } // if morph container
   recalculate ();
}

void
CLPTemplateMorph::setParamVecFromArray (const double *arrayPtr)
{
   int numMorph = size();
   for (int loop = 0; loop < numMorph; ++loop)
   {
      int index = m_indexVec.at (loop);
      if (index < 0)
      {
         // no match
         continue;
      }
      m_paramVec.at (loop) = arrayPtr [index];
   } // for loop
   if (isMorphContainer())
   {
      for (PtrVecIter iter = m_ptrVec.begin();
           m_ptrVec.end() != iter;
           ++iter)
      {
         (*iter)->setParamVec (m_paramVec);
      } // for iter
   }
   recalculate ();
}

int
CLPTemplateMorph::morphIndex (const string &name) const
{
   // Since we aren't loading strings very often, I'm not worried
   // about this taking a lot of time.  If I were, I'd use a map
   // instead.
   int numMorphs = size();
   for (int index = 0; index < numMorphs; ++index)
   {
      if (m_nameVec.at(index) == name)
      {
         return index;
      }
   } // for loop
   // if we're still here, we didn't find it
   return -1;
}

void
CLPTemplateMorph::setInternalIndicies (const SVec &nameVec)
{   
   if (isMorphContainer())
   {
      for (PtrVecIter iter = m_ptrVec.begin();
           m_ptrVec.end() != iter;
           ++iter)
      {
         (*iter)->setInternalIndicies (nameVec);
      } // for iter
      // we don't want to do this for this copy, so we're done
   } // if is a morph container

   // set all of them to -1 to be safe
   for (IVecIter iter = m_indexVec.begin();
        m_indexVec.end() != iter;
        ++iter)
   {
      *iter = -1;
   }
   int size = (int) nameVec.size();
   for (int which = 0; which < size; ++which)
   {
      int index = morphIndex ( nameVec.at (which) );
      if (index >= 0)
      {
         m_indexVec.at (index) = which;
      } // if found
   } // for which
}

void
CLPTemplateMorph::recalculate ()
{
   if (isMorphContainer())
   {
      for (PtrVecIter iter = m_ptrVec.begin();
           m_ptrVec.end() != iter;
           ++iter)
      {
         (*iter)->recalculate();
      } // for iter
      _fillMorphContainer();
      return;
   } // if is a morph container

   // O.k. If we're here, we're not a morph container, so lets get
   // down to work.
   CLPBinCont::setDontRemovePointsDefault( dontRemovePoints() );
   // start off with your default container
   CLPBinCont cont (m_defaultBC);
   int numMorph = size();
   for (int loop = 0; loop < numMorph; ++loop)
   {
      // are there any entries?
      const CLPBinCont::DMap &dmap = m_deltaMorphMapVec.at (loop);
      unsigned int dmapsize = dmap.size();
      if (! dmapsize)
      {
         // Nothing to see here folks.  Let's go on to the next morph
         continue;
      }
      double value = m_paramVec.at (loop);
      double center = m_centerVec.at (loop);
      // is there only 1 entry
      if (1 == dmapsize)
      {
         // There's only one other template, so
         CLPBinCont::DMapConstIter iter = dmap.begin();
         const CLPBinCont &delta = iter->second;
         double fraction = ( value - center ) / 
            (iter->first - center);
         CLPBinCont tempCont;
         tempCont.horizontalAdd (cont, delta, fraction);
         cont = tempCont;
         continue;
      }
      // if we're here, then we have at least two entries.  Our next
      // job is to find the two best consecuative iterators.  At the
      // end of the following loop, both iteraters MUST be valid.
      CLPBinCont::DMapConstIter leftIter  = dmap.begin();
      CLPBinCont::DMapConstIter rightIter = leftIter; ++rightIter;
      bool keepLooping = true;
      if (value < leftIter->first)
      {
         // don't bother looping
         keepLooping = false;
      }
      while (keepLooping)
      {
         // are we there yet?
         if (leftIter->first <= value &&
             rightIter->first > value)
         {
            // Perfect
            break;
         }
         ++leftIter;
         ++rightIter;
         if (dmap.end() == rightIter)
         {
            // too far
            --leftIter;
            --rightIter;
            break;
         } // if past the end
      } // while keepLooping;

      double left  = leftIter->first;
      double right = rightIter->first;
      //cout << "left " << left << " right " << right << endl;

      // There are three cases where we want to use the center
      // 1) L < C < R
      // 2) V and C < L
      // 3) V and C > R

      // Case 1 - Using Center
      if (left < center && center < right)
      {
         // Using center case 1:
         // if V < C, we use LC, else we use CR
         if (value <= center)
         {
            // LC - fraction >= 0
            double fraction = (center - value) / (center - left);
            CLPBinCont tempCont;
            const CLPBinCont &delta = leftIter->second;
            tempCont.horizontalAdd (cont, delta, fraction);         
            cont = tempCont;
            continue;            
         } else {
            // CR - fraction > 0
            double fraction = (value - center) / (right - center);
            CLPBinCont tempCont;
            const CLPBinCont &delta = rightIter->second;
            tempCont.horizontalAdd (cont, delta, fraction);         
            cont = tempCont;
            continue;            
         }
      } // center in between left and center

      // Case 2 - Using Center
      if (value <= left && center < left)
      {
         // fraction can be positive or negative
         double fraction = (value - center) / (left - center);
         CLPBinCont tempCont;
         const CLPBinCont &delta = leftIter->second;
         tempCont.horizontalAdd (cont, delta, fraction);         
         cont = tempCont;
         continue;            
      }

      // Case 3 - Using Center
      if (value >= right && center > right)
      {
         // fraction can be positive or negative
         double fraction = (center - value) / (center - right);
         CLPBinCont tempCont;
         const CLPBinCont &delta = rightIter->second;
         tempCont.horizontalAdd (cont, delta, fraction);         
         cont = tempCont;
         continue;            
      }

      // If we're here, then we want to use the right and left iterators
      assert (leftIter->first != rightIter->first);
      CLPBinCont delta;
      double fraction = (value - left) / (right - left);
      delta.horizontalDeltaMorph (leftIter->second, rightIter->second, 
                                  fraction);
      // are we above or below the center
      // if (value > center)
      // {
      //    double fraction = (value - left) / (right - left);
      //    delta.horizontalDeltaMorph (leftIter->second, rightIter->second, 
      //                                fraction);
      // } else {
      //    double fraction = (value - right) / (left - right);
      //    delta.horizontalDeltaMorph (rightIter->second, leftIter->second, 
      //                                fraction);
      // }
      CLPBinCont tempCont;
      tempCont.horizontalAdd (cont, delta);
      _printTGraph (cont); _printTGraph (delta);
      cont = tempCont;      

   } // for loop
   cont.fillTH1F (this, m_noUnderOverflow);
}

void
CLPTemplateMorph::allTGraphs (TGraphPtrVec &tgraphPtrVec) const
{
   tgraphPtrVec.clear();
   TGraph *ptr = m_defaultBC.tgraphPtr();
   ptr->SetTitle ("Default");
   tgraphPtrVec.push_back (ptr);
   for (int loop = 0; loop < size(); ++loop)
   {
      const CLPBinCont::DMap &tempMap = m_deltaMorphMapVec.at (loop);
      const string &name = m_nameVec.at (loop);
      for (CLPBinCont::DMapConstIter iter = tempMap.begin();
           tempMap.end() != iter;
           ++iter)
      {
         ptr = iter->second.tgraphPtr();
         ptr->SetTitle (Form ("%s_%f", name.c_str(), iter->first) );
         tgraphPtrVec.push_back (ptr);         
      } // for iter
   } // for loop
}

void
CLPTemplateMorph::takeMorphVariable  (const CLPTemplateMorph &rhs,
                                      const string &name)
{
   if (m_isMorphContainer)
   {
      assert (! m_isMorphContainer );
   }
   int index = rhs.morphIndex (name);
   if (-1 == index )
   {
      cerr << "rhs does not have '" << name << "' defined." << endl;
      assert (0);
   }
   if (-1 != morphIndex (name) )
   {
      cerr << "This object does have '" << name << "' defined." << endl;
      assert (0);
   }
   m_nameVec.push_back (name);
   m_deltaMorphMapVec.push_back( rhs.m_deltaMorphMapVec.at ( index ) );
   m_paramVec.        push_back( rhs.m_paramVec.at         ( index ) );
   m_centerVec.       push_back( rhs.m_centerVec.at        ( index ) );
   m_indexVec.        push_back( rhs.m_indexVec.at         ( index ) );
}

void
CLPTemplateMorph::addMorph (CLPTemplateMorph *morphPtr)
{
   m_isMorphContainer = true;
   // We need for the container morph (this one) to have the same
   // number of variables as all of the morphs it contains.  First,
   // lets make sure that the names that we have are in the same
   // order.
   int minsize = std::min (size(), morphPtr->size());
   for (int loop = 0; loop < minsize; ++loop)
   {
      if (m_nameVec.at(loop) != morphPtr->m_nameVec.at(loop))
      {
         cerr << "CLPTemplateMorph::addMorph Error: "
              << "Variable Names do not match: " << loop
              << " : " << m_nameVec.at(loop) << " != " 
              << morphPtr->m_nameVec.at(loop) << endl;
         assert (0);
      } // if strings don't match
   } // for loop;
   if (size() < morphPtr->size())
   {
      // this is the first one
      m_paramVec = morphPtr->m_paramVec;
      m_nameVec  = morphPtr->m_nameVec;
      m_indexVec = morphPtr->m_indexVec;
   }
   m_ptrVec.push_back (morphPtr);
   int totalBins = 0;
   for (PtrVecIter iter = m_ptrVec.begin();
        m_ptrVec.end() != iter;
        ++iter)
   {
      totalBins += (*iter)->GetNbinsX();
   } // for iter
   SetBins (totalBins, 0.5, totalBins + 0.5);
   _fillMorphContainer();
}

TH1F*
CLPTemplateMorph::th1fCopy() const
{
   TString name = GetName();
   name += "_";
   int bins = GetNbinsX();
   double lower = GetXaxis()->GetXmin();
   double upper = GetXaxis()->GetXmax();
   TH1F *retval = new TH1F (name, name, bins, lower, upper);
   for (int bin = 0; bin <= bins + 1; ++bin)
   {
      retval->SetBinContent (bin, GetBinContent (bin));
   }
   return retval;
}

void 
CLPTemplateMorph::setParamVecs (PtrVec &morphPtrVec, const DVec &paramVec)
{
   for (PtrVec::iterator iter = morphPtrVec.begin();
        morphPtrVec.end() != iter;
        ++iter)
   {
      (*iter)->setParamVec (paramVec);
   } // for iter
}

void
CLPTemplateMorph::_printTGraph (const CLPBinCont &cont) const
{
   static int times = 0;
   return;
   TCanvas *c1 = new TCanvas;
   TGraph *tgPtr = cont.tgraphPtr();
   tgPtr->Draw("ALP");
   c1->Print ( Form ("bincont_%02d.eps", times++) );
   delete tgPtr;
   delete c1;
}

void
CLPTemplateMorph::_fillMorphContainer()
{
   int binOffset = 0;
   Reset();
   for (PtrVecIter iter = m_ptrVec.begin();
        m_ptrVec.end() != iter;
        ++iter)
   {
      int numBins = (*iter)->GetNbinsX();
      for (int bin = 1; bin <= numBins; ++bin)
      {
         SetBinContent (binOffset + bin, (*iter)->GetBinContent (bin));
      } // for bin
      binOffset += numBins;
   } // for iter
   assert (GetNbinsX() == binOffset);
}

// friends
ostream& operator<< (ostream& o_stream, const CLPTemplateMorph &rhs)
{
   int size = rhs.size();
   for (int loop = 0; loop < size; ++loop)
   {
      o_stream << setw(10) << rhs.m_nameVec.at (loop)
               << " : "    << rhs.m_paramVec.at (loop)
               << endl;
      if (rhs.m_outputAll)
      {
         const CLPBinCont::DMap &deltaMap = rhs.m_deltaMorphMapVec.at (loop);
         for (CLPBinCont::DMapConstIter iter = deltaMap.begin();
              deltaMap.end() != iter;
              ++iter)
         {
            cout << "  " << rhs.m_nameVec.at (loop) << " " 
                 << iter->first << endl << iter->second << endl << endl;
         } // for iter
      } // if outputAll
   } // for loop
   return o_stream;
} 
