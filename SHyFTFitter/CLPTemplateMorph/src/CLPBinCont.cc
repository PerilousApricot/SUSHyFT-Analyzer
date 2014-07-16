// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "CLPTemplateMorph/interface/CLPBinCont.h"
#include "CLPTemplateMorph/interface/CLPFunctionoid2Slide.h"
#include "CLPTemplateMorph/interface/dout.h"

using namespace std;

ClassImp (CLPBinCont)

// static member data
bool CLPBinCont::sm_dontRemovePointsDefault = true;


CLPBinCont::CLPBinCont() : 
   m_vec(), m_integral(0), m_noOverflow (false), 
   m_debug (false)
{
   m_dontRemovePoints = sm_dontRemovePointsDefault;
}

CLPBinCont::CLPBinCont (TH1F *histPtr) : 
   m_vec(), m_integral(0), m_noOverflow (false), 
   m_debug (false)
{
   m_dontRemovePoints = sm_dontRemovePointsDefault;
   integrateTH1F (histPtr);
}

void
CLPBinCont::integrateTH1F (TH1F *histPtr)
{
   clear();
   int numBins = histPtr->GetNbinsX();
   double integral = histPtr->Integral (0, numBins + 1);
   if (! integral)
   {
      // this histogram is empty, so stop
      return;
   }
   m_integral = integral;
   double total = 0;
   // before we start, put an entry of 0 below the underflow bin
   double low   = histPtr->GetBinLowEdge (0);
   double width = histPtr->GetBinWidth (0);
   addBin (low - width, width, 0, true);
   for (int bin = 0; bin <= numBins + 1; ++bin)
   {
      low     = histPtr->GetBinLowEdge (bin);
      width   = histPtr->GetBinWidth (bin);
      double content = histPtr->GetBinContent (bin);
      total += content;
      bool isZero = (0 == content);
      addBin (low, width, total / integral, isZero);
   } // for bin
   // after we've started, lets just make sure it ends the way we
   // think it should
   addBin (low + width, width, 1, true);
}

void
CLPBinCont::addBin (double low, double width, double total, bool isZero)
{
   bool sameAsLast = isZero;
   // first, is this a zero change bin?
   if (isZero)
   {
      // are there at least two bins?
      if (m_vec.size() >= 2)
      {
         // is the previous bin a 'sameAsLastBin?'
         CLPBin::VecIter iter = m_vec.end();
         --iter;
         if (iter->isSameAsLastBin())
         {
            // Aha!  So, lets just modify the X values and NOT add
            // another bin.  The lower edge will be the same and the
            // new width is the old width + 'width'.
            double lower = iter->xLower();
            width += iter->xWidth();
            iter->setXLowerWidth (lower, width);
            return;
         } // if same as last bin
      } // if size >= 2
      else if (0 == m_vec.size())
      {
         // if this is the first bin, then don't mark it same as the
         // last.
         sameAsLast = false;
      } // if first entry
   } // if this is a zero bin
   m_vec.push_back( CLPBin (low, width, total, sameAsLast) );
}

void
CLPBinCont::addBin (double upper, double total)
{
   m_vec.push_back( CLPBin (upper, 0., total, false) );
}

void
CLPBinCont::addBinCarefully (double low, double width, double total)
{
   if (m_dontRemovePoints)
   {
      // don't do this
      return;
   }
   // If this function is being called, it is because there are
   // ordering issues with the previous bin.
   assert (! m_vec.empty());
   CLPBin::VecIter iter = m_vec.end() - 1;
   double prevY = -1;
   do
   {
      // can we get away with modifying this bin?
      double prevWidth = low - iter->xLower();
      if (prevWidth >= 0)
      {
         // good, we can just update this one.
         iter->setXLowerWidth ( iter->xLower(), prevWidth );
         prevY = iter->yValue();
         break;
      }
      // If we're here, then we need to delete this one and try again.
      // 'erase' returns the element AFTER the iterator erased.
      iter = m_vec.erase (iter) - 1;
   } while (! m_vec.empty());
   addBin (low, width, total, (total == prevY));
}


void
CLPBinCont::fillTH1F (TH1F *histPtr, bool noUnderOverflow) const
{
   histPtr->Reset();
   if (m_vec.size() < 2)
   {
      // Nothing to do here, folks.
      return;
   }
   ///////////////////////////////////////////////
   // Find the pair of iterators that encompass //
   // the lower edge of the histogram.          //
   ///////////////////////////////////////////////
   CLPBin::VecConstIter lowerIter = m_vec.begin();
   CLPBin::VecConstIter upperIter = lowerIter + 1;   
   double startX = histPtr->GetBinLowEdge(1);
   double startY = 0;
   if (! moveIterators (lowerIter, upperIter, startX, startY))
   {
      cerr << "CLPBinCont::fillTH1F() error: Could not find lower edge. "
           << "Aborting." << endl;
      return;
   }
   // set the underflow bin
   histPtr->SetBinContent (0, startY * integral());
   int numBins = histPtr->GetNbinsX();
   double endY = 0;
   for (int bin = 1; bin <= numBins; ++bin)
   {
      double endX = histPtr->GetBinLowEdge (bin) + histPtr->GetBinWidth (bin);
      if (! moveIterators (lowerIter, upperIter, endX, endY))
      {
         // we've reached the end, so lets stop filling this histogram
         return;
      }
      histPtr->SetBinContent (bin, (endY - startY) * integral());
      startY = endY;
   } // for bin
   // if we're still here, set the overflow bin
   histPtr->SetBinContent (numBins + 1, (1 - endY) * integral());

   if (noUnderOverflow)
   {
      // add under and overflow to first and last bin
      histPtr->SetBinContent ( 1,
                               histPtr->GetBinContent (0) +
                               histPtr->GetBinContent (1) );
      histPtr->SetBinContent ( numBins,
                               histPtr->GetBinContent (numBins) +
                               histPtr->GetBinContent (numBins + 1) );
      // Reset under and overflow bins
      histPtr->SetBinContent (0, 0.);
      histPtr->SetBinContent (numBins + 1, 0.);
   }
   
}

TGraph*
CLPBinCont::tgraphPtr() const
{
   vector< double > xVec, yVec;
   for (CLPBin::VecConstIter iter = m_vec.begin();
        m_vec.end() != iter;
        ++iter)
   {
      xVec.push_back( iter->xUpper() );
      yVec.push_back( iter->yValue() );
   } // for iter
   int size = (int)xVec.size();
   TGraph *retval = new TGraph (size, &(xVec[0]), &(yVec[0]));
   return retval;
}

bool
CLPBinCont::moveIterators (CLPBin::VecConstIter &lowerIter,
                          CLPBin::VecConstIter &upperIter,
                          double xValue, double &yValue) const
{
   bool found = false;
   double inBetween = -1;
   do
   {
      inBetween = (xValue - lowerIter->xUpper()) / upperIter->xWidth();
      if (inBetween >= 0 && inBetween < 1)
      {
         found = true;
         break;
      }
      // if we're still here, then we need to increment the iterators
      ++lowerIter;
      ++upperIter;
      if (m_vec.end() == upperIter)
      {
         return false;
      }
   } while (! found);
   yValue = lowerIter->yValue() + 
      inBetween * (upperIter->yValue() - lowerIter->yValue());
   return true;
}

void
CLPBinCont::horizontalMorph (const CLPBinCont &alpha, const CLPBinCont &beta, 
                             double param)
{
   clear();
   // set the integral first.
   m_integral = percentInBetween (alpha.m_integral, beta.m_integral, param);
   assert (alpha.m_vec.size() >= 3 && beta.m_vec.size() >= 3);
   // the first bin is the 0 bin which is always there
   CLPBin::VecConstIter alphaLowerIter = alpha.m_vec.begin();
   CLPBin::VecConstIter betaLowerIter  = beta.m_vec.begin();
   m_vec.push_back (*alphaLowerIter);
   double alphaXLower = alphaLowerIter->xUpper() + alphaLowerIter->xWidth();
   double betaXLower  = betaLowerIter->xUpper() + betaLowerIter->xWidth();;
   if (alphaXLower != betaXLower)
   {
      // now check percentage difference
      double percent = 2 * (alphaXLower - betaXLower) / 
         (fabs(alphaXLower) + fabs(betaXLower));
      assert (percent < 0.0001);
   }
   ++alphaLowerIter;
   ++betaLowerIter;
   CLPBin::VecConstIter alphaUpperIter = alphaLowerIter + 1;
   CLPBin::VecConstIter betaUpperIter  = betaLowerIter  + 1;
   // we are assuming that the two CDFs are starting at the same lower
   // value   
   double lower = percentInBetween (alphaLowerIter->xLower(), 
                                    betaLowerIter->xLower(),
                                    param);   
   double upper = percentInBetween (alphaLowerIter->xUpper(), 
                                    betaLowerIter->xUpper(),
                                    param);   
   double prevY = percentInBetween (alphaLowerIter->yValue(), 
                                    betaLowerIter->yValue(),
                                    param);
   addBin (lower, upper - lower, prevY, areTheSame (0., prevY));
   double oldUpper = upper;
   while (alpha.m_vec.end() != alphaUpperIter &&
          beta.m_vec.end()  != betaUpperIter)
   {
      // are the two iterators pointing at the same Y value?
      if ( areTheSame ( alphaUpperIter->yValue(), betaUpperIter->yValue() ) )
      {
         // are either alpha or beta the same for the next iterator
         bool nextSame = false;
         CLPBin::VecConstIter alphaNextIter = alphaUpperIter + 1;
         CLPBin::VecConstIter betaNextIter = betaUpperIter + 1;
         if (alpha.m_vec.end() != alphaNextIter)
         {
            if ( areTheSame (alphaUpperIter->yValue(), 
                             alphaNextIter->yValue() ) )
            {
               nextSame = true;
               ++alphaUpperIter; ++alphaLowerIter;               
            } // if this and next iters are the same
         } // not the last alpha iter
         if (beta.m_vec.end() != betaNextIter)
         {
            if ( areTheSame (betaUpperIter->yValue(), 
                             betaNextIter->yValue() ) )
            {
               nextSame = true;
               ++betaUpperIter; ++betaLowerIter;               
            } // if this and next iters are the same
         } // not the last beta iter
         if (nextSame)
         {
            // let's go again
            continue;
         }
         double upper = percentInBetween (alphaUpperIter->xUpper(), 
                                          betaUpperIter->xUpper(),
                                          param);
         lower = oldUpper;
         double width = upper - lower;
         double newY = std::min( alphaUpperIter->yValue(), 
                                 betaUpperIter->yValue() );
         // make sure this is not a zero width bin
         if (width >=0)
         {
            addBin (lower, width, newY, areTheSame (newY, prevY));
            oldUpper = upper;
            prevY = newY;
         }
         ++alphaUpperIter; ++alphaLowerIter;
         ++betaUpperIter;  ++betaLowerIter;
         continue;
      }
      // So there not equal.  Who has the smaller Y value
      double alphaXUpper, betaXUpper, newY;
      if (alphaUpperIter->yValue() < betaUpperIter->yValue())
      {
         // alpha is the small guy
 
         // first check to make sure that the upper and lower iters
         // don't have the same Y value
         if (areTheSame (alphaUpperIter->yValue(), alphaLowerIter->yValue()))
         {
            // They are the same
            ++alphaUpperIter; ++alphaLowerIter;
            continue;
         }
         newY = alphaUpperIter->yValue();
         alphaXUpper = alphaUpperIter->xUpper();
         double frac = (newY - betaLowerIter->yValue()) /
            (betaUpperIter->yValue() - betaLowerIter->yValue());
         betaXUpper = percentInBetween (betaLowerIter->xUpper(),
                                        betaUpperIter->xUpper(),
                                        frac);
         // now increment alpha for next time
         ++alphaUpperIter; ++alphaLowerIter;        
      } else {
         // beta is the small guy
         // first check to make sure that the upper and lower iters
         // don't have the same Y value
         if (areTheSame (betaUpperIter->yValue(), betaLowerIter->yValue()))
         {
            // They are the same
            ++betaUpperIter; ++betaLowerIter;
            continue;
         }
         newY = alphaUpperIter->yValue();
         betaXUpper = betaUpperIter->xUpper();
         double frac = (newY - alphaLowerIter->yValue()) /
            (alphaUpperIter->yValue() - alphaLowerIter->yValue());
         alphaXUpper = percentInBetween (alphaLowerIter->xUpper(),
                                         alphaUpperIter->xUpper(),
                                         frac);
         ++betaUpperIter; ++betaLowerIter;
      }
      lower = oldUpper;
      double upper = percentInBetween (alphaXUpper, betaXUpper, param);
      double width = upper - lower;
      // make sure this is not a zero width bin
      if (width >= 0)
      {
         addBin (lower, width, newY, areTheSame (newY, prevY));
         oldUpper = upper;
         prevY = newY;
      }
      // cout << endl;
   } // while iters are valid
   // finish it off.
   lower = oldUpper;
   upper = percentInBetween (alphaLowerIter->xUpper(), 
                             betaLowerIter->xUpper(),
                             param);
   double yVal = percentInBetween (alphaLowerIter->yValue(), 
                                   betaLowerIter->yValue(),
                                   param);
   addBin (lower, upper - lower, yVal, areTheSame (yVal, prevY));
}

void
CLPBinCont::horizontalSubtract (const CLPBinCont &alpha, 
                                const CLPBinCont &beta)
{
   static CLPFunctionoid2 subtract;
   horizontalFunction (alpha, beta, subtract);
}

void 
CLPBinCont::horizontalDeltaMorph (const CLPBinCont &alpha,
                                  const CLPBinCont &beta,
                                  double percentage)
{
   CLPFunctionoid2Slide slide (percentage);
   horizontalFunction (alpha, beta, slide);
}


void
CLPBinCont::horizontalFunction (const CLPBinCont &alpha, 
                                const CLPBinCont &beta,
                                const CLPFunctionoid2 &functionoid)
{
   clear();
   m_integral = beta.m_integral;
   assert (alpha.m_vec.size() >= 2 && beta.m_vec.size() >= 2);
   CLPBin::VecConstIter alphaLowerIter = alpha.m_vec.begin();
   CLPBin::VecConstIter betaLowerIter  = beta.m_vec.begin();
   CLPBin::VecConstIter alphaUpperIter = alphaLowerIter + 1;
   CLPBin::VecConstIter betaUpperIter  = betaLowerIter  + 1;
   // before we add anything, they should both be at 0;
   addBin (functionoid (alphaLowerIter->xLower(), 
                        betaLowerIter->xLower()), 
           0.);
   while (alpha.m_vec.end() != alphaUpperIter &&
          beta.m_vec.end()  != betaUpperIter)
   {
      // are the two iterators pointing at the same Y value?
      if ( areTheSame( alphaUpperIter->yValue(), betaUpperIter->yValue() ) )
      {
         addBin (functionoid (alphaLowerIter->xUpper(),
                              betaLowerIter->xUpper()),
                 std::min (alphaUpperIter->yValue(),
                           betaUpperIter->yValue()) );
         ++alphaUpperIter; ++alphaLowerIter;
         ++betaUpperIter;  ++betaLowerIter;
         continue;
      } // both the same
      // So there not equal.  Who has the smaller Y value
      if (alphaUpperIter->yValue() < betaUpperIter->yValue())
      {
         // alpha is the small guy
         double frac = (alphaUpperIter->yValue() - betaLowerIter->yValue()) /
            (betaUpperIter->yValue() - betaLowerIter->yValue());
         double betaXUpper = percentInBetween (betaLowerIter->xUpper(),
                                               betaUpperIter->xUpper(),
                                               frac);
         addBin (functionoid (alphaUpperIter->xUpper(),
                              betaXUpper),                              
                 alphaUpperIter->yValue());
         // now increment alpha for next time
         ++alphaUpperIter; ++alphaLowerIter;
      } else {
         // beta is the small guy
         double frac = (betaUpperIter->yValue() - alphaLowerIter->yValue()) /
            (alphaUpperIter->yValue() - alphaLowerIter->yValue());
         double alphaXUpper = percentInBetween (alphaLowerIter->xUpper(),
                                                alphaUpperIter->xUpper(),
                                                frac);
         addBin (functionoid (alphaXUpper,
                              betaUpperIter->xUpper()), 
                 betaUpperIter->yValue());
         // now increment beta for next time
         ++betaUpperIter; ++betaLowerIter;
      }
      // finish it off
   } // while not at the end
}

void
CLPBinCont::horizontalAdd (const CLPBinCont &alpha, 
                           const CLPBinCont &delta,
                           double param)
{
   clear();
   m_integral = percentInBetween (alpha.m_integral,
                                  delta.m_integral,
                                  param);
   assert (alpha.m_vec.size() >= 3 && delta.m_vec.size() >= 3);
   // the first bin is the 0 bin which is always there.
   CLPBin::VecConstIter alphaLowerIter = alpha.m_vec.begin();
   CLPBin::VecConstIter deltaLowerIter = delta.m_vec.begin();
   m_vec.push_back (*alphaLowerIter);
   double xLower = alphaLowerIter->xUpper();
   ++alphaLowerIter;
   ++deltaLowerIter;
   CLPBin::VecConstIter alphaUpperIter = alphaLowerIter + 1;
   CLPBin::VecConstIter deltaUpperIter = deltaLowerIter + 1;
   // add the first point
   double prevY = 0;
   while (alpha.m_vec.end() != alphaUpperIter &&
          delta.m_vec.end() != deltaUpperIter)
   {
      // are the two iterators pointing at the same Y value?
      if ( areTheSame (alphaUpperIter->yValue(), deltaUpperIter->yValue() ) )
      {
         // are either alpha or delta the same for the next iterator
         bool nextSame = false;
         CLPBin::VecConstIter alphaNextIter = alphaUpperIter + 1;
         CLPBin::VecConstIter deltaNextIter = deltaUpperIter + 1;
         if (alpha.m_vec.end() != alphaNextIter)
         {
            if ( areTheSame (alphaUpperIter->yValue(), 
                             alphaNextIter->yValue() ) )
            {
               nextSame = true;
               ++alphaUpperIter; ++alphaLowerIter;               
            } // if this and next iters are the same
         } // not the last alpha iter
         if (delta.m_vec.end() != deltaNextIter)
         {
            if ( areTheSame (deltaUpperIter->yValue(), 
                             deltaNextIter->yValue() ) )
            {
               nextSame = true;
               ++deltaUpperIter; ++deltaLowerIter;               
            } // if this and next iters are the same
         } // not the last delta iter
         if (nextSame)
         {
            // let's go again
            continue;
         }
         double xUpper = alphaUpperIter->xUpper() + 
            param * deltaUpperIter->xUpper();         
         double newY = std::min( alphaUpperIter->yValue(),
                                 deltaUpperIter->yValue() );
         // Make sure we are not adding a negative width bin.  If we
         // are, skip it.
         if (xLower <= xUpper)
         {
            addBin (xLower, xUpper - xLower, newY, areTheSame (newY, prevY));
            xLower = xUpper;
            prevY = newY;
         } else {
            addBinCarefully (xUpper, xLower - xUpper, newY);
         }
         ++alphaUpperIter; ++alphaLowerIter;
         ++deltaUpperIter; ++deltaLowerIter;
         continue;
      } // both the same
      // So there not equal.  Who has the smaller Y value
      if (alphaUpperIter->yValue() < deltaUpperIter->yValue())
      {
         // alpha is the small guy

         // first check to make sure that the upper and lower iters
         // don't have the same Y value
         if (areTheSame (alphaUpperIter->yValue(), alphaLowerIter->yValue()))
         {
            // They are the same
            ++alphaUpperIter; ++alphaLowerIter;
            continue;
         }
         double frac = (alphaUpperIter->yValue() - deltaLowerIter->yValue()) /
            (deltaUpperIter->yValue() - deltaLowerIter->yValue());
         double deltaXUpper = percentInBetween (deltaLowerIter->xUpper(),
                                                deltaUpperIter->xUpper(),
                                                frac);
         double xUpper = alphaUpperIter->xUpper() + param * deltaXUpper;
         double newY = alphaUpperIter->yValue();
         // Make sure we are not adding a negative width bin.  If we
         // are, skip it.
         if (xLower <= xUpper)
         {
            addBin (xLower, xUpper - xLower, newY, areTheSame (newY, prevY));
            xLower = xUpper;
            prevY = newY;
         }
         // now increment alpha for next time
         ++alphaUpperIter; ++alphaLowerIter;
      } else {
         // delta is the small guy
         // first check to make sure that the upper and lower iters
         // don't have the same Y value
         if (areTheSame (deltaUpperIter->yValue(), deltaLowerIter->yValue()))
         {
            // They are the same
            ++deltaUpperIter; ++deltaLowerIter;
            continue;
         }
         double frac = (deltaUpperIter->yValue() - alphaLowerIter->yValue()) /
            (alphaUpperIter->yValue() - alphaLowerIter->yValue());
         double alphaXUpper = percentInBetween (alphaLowerIter->xUpper(),
                                                alphaUpperIter->xUpper(),
                                                frac);
         double xUpper = alphaXUpper + param * deltaUpperIter->xUpper();
         double newY = deltaUpperIter->yValue();
         // Make sure we are not adding a negative width bin.  If we
         // are, skip it.
         if (xLower <= xUpper)
         {
            addBin (xLower, xUpper - xLower, newY, areTheSame (newY, prevY));
            xLower = xUpper;
            prevY = newY;
         } else {
            addBinCarefully (xUpper, xLower - xUpper, newY);
         }
         // now increment delta for next time
         ++deltaUpperIter; ++deltaLowerIter;
      }
   } // while not at the end
   // // finish it off
   // double xUpper = alphaLowerIter->xUpper() + 
   //    param *  deltaLowerIter->xUpper();
   // addBin (xLower, xUpper - xLower, 1, areTheSame (1., prevY));
}

void
CLPBinCont::smoothBeginningOfDelta (double percentStep)
{
   double kStart = 0.10;
   double minVal = 0;
   double maxVal = 0;
   bool first = true;
   bool foundIter = false;
   CLPBin::VecIter leftIter (m_vec.begin());
   for (CLPBin::VecIter iter = m_vec.begin();
        m_vec.end() != iter;
        ++iter)
   {
      if (! foundIter && iter->yValue() > kStart)
      {
         leftIter = iter;
         foundIter = true;
      }
      if (first)
      {
         minVal = maxVal = iter->xUpper();
         first = false;
         continue;
      }
      minVal = std::min (minVal, iter->xUpper());
      maxVal = std::max (maxVal, iter->xUpper());
   } // for iter
   dout << "min " << minVal << " max " << maxVal << endl;
   double step = percentStep * (maxVal - minVal);
   dout << "step " << step << endl;
   CLPBin::VecIter rightIter = leftIter + 1;
   assert (m_vec.end() != rightIter);
   while ( (m_vec.begin() != rightIter) )
   {
      if (0 == rightIter->yValue())
      {
         // if we're here, we want to delete everything from the
         // beginning up to, but not including, the rightIter.
         dout << "found 0" << endl;
         m_vec.erase( m_vec.begin(), rightIter );
         return;
      }
      if ( fabs( leftIter->xUpper() - rightIter->xUpper() ) > step)
      {
         dout << "step" << endl;
         // if we're here, we want to
         // 1) fix leftIter such that:
         //   a) it has a y value of 0, and
         //   b) it has an x value that is linear from the next two points, and
         // 2) delete everything from the beginning up to, but not
         //    including, leftIter.
         CLPBin::VecIter nextRightIter = rightIter + 1;
         assert (m_vec.end() != nextRightIter);
         if ( ( rightIter->xUpper() == nextRightIter->xUpper() ) || 
              ( rightIter->yValue() >= nextRightIter->yValue() ))
         {
            // we have a vertical line OR no lever arm
            leftIter->setXUpperWidth( rightIter->xUpper(), 0. );
            leftIter->setYvalue (0.);
         } else {
            double deltaX =  nextRightIter->xUpper() - rightIter->xUpper();
            double deltaY =  nextRightIter->yValue() - rightIter->yValue();
            double numUnits = rightIter->yValue() / deltaY;
            leftIter->setXUpperWidth( rightIter->xUpper() - deltaX * numUnits,
                                      0.);
            leftIter->setYvalue (0.);            
         }
         m_vec.erase( m_vec.begin(), leftIter);
         return;
      } // if step too big
      --leftIter;
      --rightIter;
   } // while     
   // If we're here, then there's nothing to do
}

bool
CLPBinCont::areTheSame (double alpha, double beta)
{
   const double kEpsilon = 1e-7;
   return (fabs (alpha - beta) < kEpsilon);
}

// friends
ostream& operator<< (ostream& o_stream, const CLPBinCont &rhs)
{
   int bin = 0;
   for (CLPBin::VecConstIter iter = rhs.m_vec.begin();
        rhs.m_vec.end() != iter;
        ++iter)
   {
      o_stream << setw(4) << bin++ << ") " << *iter << endl;
   }
   o_stream << "Integral: " << rhs.integral() << endl;
   return o_stream;
} 
