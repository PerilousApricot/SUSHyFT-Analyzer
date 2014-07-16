// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "CLPTemplateFitter/interface/CLPTrapezoid.h"
#include "CLPTemplateFitter/interface/dout.h"
#include "TGraph.h"
#include "CLPTemplateMorph/interface/dumpSTL.icc"

using namespace std;

ClassImp (CLPTrapezoid)

CLPTrapezoid::CLPTrapezoid (double x, double y) : 
m_lowerX (x), m_upperX (x), m_lowerY (y), m_upperY (y), 
   m_totalIntegral (0.), m_position (-1)
{
}

bool
CLPTrapezoid::operator> (const CLPTrapezoid &rhs) const
{
   return maxY() > rhs.maxY();
}

bool
CLPTrapezoid::operator< (const CLPTrapezoid &rhs) const
{
   return maxY() < rhs.maxY();
}

double
CLPTrapezoid::integral() const
{
   return 0.5 * (m_upperY + m_lowerY) * (m_upperX - m_lowerX);
}

void
CLPTrapezoid::scale (double factor)
{
   m_lowerY *= factor;
   m_upperY *= factor;
   m_totalIntegral *= factor;
}

int
CLPTrapezoid::positionWRT (double xVal) const
{
   if (xVal < lowerX())
   {
      // we are to the right, so postive
      return kToTheRight;
   }
   if (xVal > upperX())
   {
      // we are to the left, so negative
      return kToTheLeft;
   }
   // If we're still here, the porrage must be just right...
   return kInside;
}

bool
CLPTrapezoid::hasVerticalOverlapWith (const CLPTrapezoid &rhs) const
{
   double myMax  = maxY();
   double myMin  = minY();
   double rhsMax = rhs.maxY();
   double rhsMin = rhs.maxY();
   if ((myMax > rhsMax && myMin < rhsMax) ||
       (myMax > rhsMin && myMin < rhsMin) )
   {
      // we've got overlap
      return true;
   } else {
      // no overlap
      return false;
   }
}

double
CLPTrapezoid::absSlope() const
{
   return fabs (m_upperY - m_lowerY) / (m_upperX - m_lowerX);
}

// static member functions
void
CLPTrapezoid::limits (Vec &curve,
                      double &upperLimit, double &lowerLimit, 
                      double percent)
{
   upperLimit = 0;
   lowerLimit = upperLimit + 1;
   assert (percent > 0 && percent < 1);

   // for (VecConstIter iter = curve.begin();
   //      curve.end() != iter;
   //      ++iter)
   // {
   //    cout << Form ("%.4f, %.8f", iter->upperX(), 
   //                  iter->upperY() / 7.06654e-24)
   //         << endl;
   // }

   ////////////////////////////
   // Make Proper Trapezoids //
   ////////////////////////////
   VecIter leftIter = curve.begin();
   double lowestX = leftIter->lowerX();
   for (VecIter rightIter = leftIter + 1;
        curve.end() != rightIter;
        ++leftIter, ++rightIter)
   {
      leftIter->setUpperX( rightIter->upperX() );
      leftIter->setUpperY( rightIter->upperY() );
   } // for iter
   // get rid of the last (zero width) trapezoid
   double biggestX = leftIter->upperX();
   curve.erase (leftIter);

   ////////////////////////////////////////
   // Sort Trapezoids from Tall to Short //
   ////////////////////////////////////////
   //sort (curve.rbegin(), curve.rend());
   sort ( curve.begin(), curve.end(), std::greater< CLPTrapezoid >() );
   // integrate this puppy up
   double total = 0;
   int position = 0;
   for (VecIter iter = curve.begin();
        curve.end() != iter;
        ++iter)
   {
      total += iter->integral();
      iter->setTotalIntegral (total);
      iter->setPosition (position++);
   } // for iter
   assert (total);

   /////////////////////////////
   // Rescale and Get Maximum //
   /////////////////////////////
   double maxX = -1;
   double max = 0;
   double factor = 1. / total;
   total = 0.;
   for (VecIter iter = curve.begin();
        curve.end() != iter;
        ++iter)
   {
      total += iter->integral();
      iter->scale (factor);
      if (iter->maxY() > max)
      {
         // take maxX to be the middle
         maxX = (iter->upperX() + iter->lowerX()) / 2;
         max = iter->maxY();
      } // if new maximum
   } // for iter

   ///////////////////////////////////
   // Scan from Short to Tall       //
   // Is this a Single Sided Limit? //
   ///////////////////////////////////
   bool seenLeft  = false;
   bool seenRight = false;
   // I am not changing the vector 'curve' anymore.  Since I want to
   // move from back to front, I will define and end and begin that
   // work for me.
   const VecConstIter kRevBegin = curve.end() - 1;
   const VecConstIter kRevEnd   = curve.begin() - 1;
   VecConstIter iter;
   for (iter = kRevBegin;
        kRevEnd != iter;
        --iter)
   {
      int position = iter->positionWRT (maxX);
      // is this far enough?
      if (iter->totalIntegral() < percent)
      {
         // we're done
         break;
      }
      if (kToTheRight == position)
      {
         seenRight = true;
      } else if (kToTheLeft == position)
      {
         seenLeft = true;
      }
   } // while iter
   assert (kRevEnd != iter);

   //////////////////////
   // One Sided Limits //
   //////////////////////
   // Is this an upper or lower limit
   if ( (seenRight && ! seenLeft) ||
        (seenLeft  && ! seenRight) )
   {
      bool isUpperLimit = seenRight && ! seenLeft;
      double areaNeeded = percent - iter->totalIntegral();
      ++iter;
      double slope = ( iter->upperY() - iter->lowerY() ) /
         ( iter->upperY() - iter->lowerY() );
      double deltaX;
      if (slope)
      {
         double b0m;
         if (isUpperLimit)
         {
            // upper limit
            b0m = iter->lowerY() / slope;
         } else {
            slope *= -1;
            b0m = iter->upperY() / slope;
         }
         if (slope > 0)
         {
            deltaX = - b0m + sqrt (b0m * b0m + 2 * areaNeeded / slope);
         } else {
            deltaX = - b0m - sqrt (b0m * b0m + 2 * areaNeeded / slope);
         } // else if slope < 0         
      } // if slope
      else
      {
         deltaX = areaNeeded / iter->lowerY();
      }
      if (isUpperLimit)
      {
         // Upper Limit
         lowerLimit = lowestX - 1;
         upperLimit = iter->lowerX() + deltaX;
      } else {
         // Lower Limit
         upperLimit = biggestX + 1;
         lowerLimit = iter->upperX() - deltaX;
      }
      return;
   }

   //////////////////////////////////////////
   // If we're still here, then we have... //
   // Two Sided Limits                     //
   //////////////////////////////////////////
   int sameAs, other;
   if (kToTheRight == iter->positionWRT (maxX))
   {
      sameAs = kToTheRight;
      other  = kToTheLeft;
   } else {
      sameAs = kToTheLeft;
      other  = kToTheRight;
   }

   // First, lets look for a previous iterator so we know we are too
   // small.
   VecConstIter otherIter = iter - 1;
   while (curve.end() != otherIter)
   {
      if (otherIter->positionWRT (maxX) == other)
      {
         // done
         break;
      } // if
      --otherIter;
   } // while
   assert (kRevEnd != otherIter);
   double minArea, maxArea;
   minMaxArea (minArea, maxArea, iter, otherIter);

   int count = 0;
   //cout << curve << endl;
   while (maxArea < percent || minArea > percent)
   {
      if (0) cout << "iter  " << *iter << endl 
                  << "other " << *otherIter << endl 
                  << "min " << minArea << " max " << maxArea << endl << endl;
      if (++count > 10)
      {
         // We shoudln't be here anymore.  Put out an error message
         // and return with the initial values so that it's clear we
         // failed.

         return;
      }
      // which way do we need to go
      VecConstIter nextIter;
      if (maxArea < percent)
      {
         // we don't have enough yet.  Which one is smaller?
         if (*iter < *otherIter)
         {
            nextIter = iter + 1;
         } else {
            nextIter = otherIter + 1;
         }
      } else {
         // we have too much.  Which one is bigger
         if (*iter > *otherIter)
         {
            nextIter = iter - 1;
         } else {
            nextIter = otherIter - 1;
         }
      } // else if minArea > percent
      // Which one should we replace?
      if (nextIter->positionWRT (maxX) == other)
      {
         // replace the other iter
         otherIter = nextIter;
      } else {
         iter = nextIter;
      }
      // Don't forget to recalculate everything.
      minMaxArea (minArea, maxArea, iter, otherIter);
   } // while
   if (0) cout << "iter  " << *iter << endl 
               << "other " << *otherIter << endl 
               << "min " << minArea << " max " << maxArea << endl << endl;
  
   // We now have the right two iterators and the right range for the
   // area.  Let's get the final answer.
   VecConstIter biggerIter  = iter;
   VecConstIter smallerIter = otherIter;
   if (*smallerIter > *biggerIter)
   {
      // Ooops.  Got 'em backwards
      biggerIter  = otherIter;
      smallerIter = iter;
   }
   // If these iterators are NOT next to each other, then we need to
   // grab the area of the intervening iterators.
   VecConstIter middleIter = biggerIter + 1;
   double additionalArea = 0;
   assert (smallerIter > biggerIter);
   while (middleIter != smallerIter)
   {
      // cout << "  Adding " << *middleIter << endl;
      additionalArea += middleIter->integral();
      ++middleIter;
   }
   // cout << "addArea " << additionalArea << endl;
   double prevTotal = biggerIter->totalIntegral() - biggerIter->integral();
   double deltaY    = biggerIter->maxY() - smallerIter->maxY();
   double deltaX    = deltaY / biggerIter->absSlope();
   minArea = prevTotal + additionalArea 
      + deltaX * smallerIter->maxY() 
      + 0.5 * deltaX * deltaY;
   double neededArea = percent - minArea;
   if (0) cout << " pT " << prevTotal 
               << " dY " << deltaY
               << " mA " << minArea 
               << " nA " << neededArea << endl;

   double biggerSlope  = biggerIter->absSlope();
   double smallerSlope = smallerIter->absSlope();
   if (! biggerSlope || ! smallerSlope)
   {
      // At least one of the two slopes is 0
      return;
   }
   assert (biggerSlope > 0 && smallerSlope > 0);

   // We want to start with the smaller maximum height, h.  If I did
   // my math correctly, then we want 
   // deltaH = h - sqrt (h^2 - 2 * mEff * neededArea), where
   // 1 / mEff = 1 / m1 + 1 / m2.  Given deltaH, we have
   // deltaX1 = deltaH / m1, and deltaX2 = deltaH / m2.
   //
   // Here goes:
   double height = smallerIter->maxY();
   double mEff = 1. / ((1. / biggerSlope) + (1. / smallerSlope));
   double deltaH = height - sqrt (height * height - 2 * mEff * neededArea);
   double deltaXbigger  = deltaH / biggerSlope;
   double deltaXsmaller = deltaH / smallerSlope;   

   // Now we need to figure out which one is on the right
   if (0) cout << " smLX " << smallerIter->lowerX() 
               << " smUX " << smallerIter->upperX()
               << " biLX " << biggerIter->lowerX() 
               << " biUX " << biggerIter->upperX()
               << " dXs  " << deltaXsmaller
               << " dXb  " << deltaXbigger << endl
               << "added Area " 
               << 0.5 * deltaXsmaller * deltaH + 0.5 * deltaXbigger * deltaH 
               << " + " << height * (deltaXsmaller + deltaXbigger) << endl;

   // Don't forget to add the smidgen of deltaXbigger that is needed
   // to match the heights.
   double additionalDeltaXbigger = 
      (biggerIter->maxY() - smallerIter->maxY()) / biggerSlope;
   if (0) cout << "aDXb " << additionalDeltaXbigger << endl;
   deltaXbigger += additionalDeltaXbigger;

   if (kToTheLeft == biggerIter->positionWRT (maxX))
   {
      // upper limit is the 'lower' part to the right + its deltaX
      upperLimit = smallerIter->lowerX() + deltaXsmaller;
      // lower limit is the 'upper' part to the left - its deltaX
      lowerLimit = biggerIter->upperX()  - deltaXbigger;
   } else {
      // upper limit is the 'lower' part to the right + its deltaX
      upperLimit = biggerIter->lowerX()  + deltaXbigger;
      // lower limit is the 'upper' part to the left - its deltaX
      lowerLimit = smallerIter->upperX() - deltaXsmaller;
   }
}

bool
CLPTrapezoid::minMaxArea (double &minArea, double &maxArea,
                          const VecConstIter &alpha,
                          const VecConstIter &beta)
{

   // I want alpha to be bigger than beta
   
   VecConstIter smallerIter;
   VecConstIter biggerIter;
   if (*beta > *alpha)
   {
      biggerIter  = beta;
      smallerIter = alpha;
   } else {
      biggerIter  = alpha;
      smallerIter = beta;
   }
   if (! biggerIter->hasVerticalOverlapWith (*smallerIter))
   {
      // Nothing we can do.  So return the total area of the taller
      // one (the smaller total area) for both values and be done with
      // it.
      maxArea = minArea = biggerIter->totalIntegral();
      return false;
   }
   // If these iterators are NOT next to each other, then we need to
   // grab the area of the intervening iterators.
   VecConstIter middleIter = biggerIter + 1;
   double additionalArea = 0;
   assert (smallerIter > biggerIter);
   while (middleIter != smallerIter)
   {
      additionalArea += middleIter->integral();
      ++middleIter;
   }
   double prevTotal = biggerIter->totalIntegral() - biggerIter->integral(); 
   // The smallest area that is allowed is the previous total + the
   // area included on the bigger piece until it hits the maximum
   // height of the smaller piece.  Note that if we're here, we know
   // that the bigger piece can not have 0 slope.
   // The additional area = 0.5 * deltaX * deltaY, where deltaX =
   // deltaY / slope.
   double deltaY = biggerIter->maxY() - smallerIter->maxY();
   double deltaX = deltaY / biggerIter->absSlope();
   minArea = prevTotal + additionalArea 
      + deltaX * smallerIter->maxY() 
      + 0.5 * deltaX * deltaY;
   //minArea = prevTotal + 0.5 * deltaY * deltaY / biggerIter->absSlope();
   // dout << "pT " << prevTotal << " dY " << deltaY 
   //      << " m " << biggerIter->absSlope() << " minA " << minArea << endl;
   // the maximum area allowed is the previous total + the areas of
   // the two trapezoids minus the small sliver of one the trapezoids
   // where the height is below the minimum height of the other.
   maxArea = biggerIter->totalIntegral() + smallerIter->integral()
      + additionalArea;
   // now correct it
   double slope;
   double alphaMin = biggerIter->minY();
   double betaMin  = smallerIter->minY();
   double min;
   deltaY = alphaMin - betaMin;
   if (deltaY < 0)
   {
      // alphaMin is smaller
      deltaY *= -1;
      slope = biggerIter->absSlope();
      min = alphaMin;
   } else {
      slope = smallerIter->absSlope();
      min = betaMin;
   }
   deltaX = deltaY / slope;
   maxArea -= deltaX * min + 0.5 * deltaY * deltaX;
   return true;
}

double
CLPTrapezoid::findYintersect (const Vec &curve, 
                              double signedStep, 
                              double yDelta, double epsilon)
{
   DVec xVec, yVec;
   assert (curve.size());
   // get TGraph and xMin and yMin
   double xMin = 0, yMin = 0;
   double first = true;
   for (CLPTrapezoid::VecConstIter iter = curve.begin();
        curve.end() != iter;
        ++iter)
   {
      double xVal( iter->lowerX() );
      double yVal( iter->lowerY() );
      xVec.push_back( xVal );
      yVec.push_back( yVal );
      if (first)
      {
         xMin = xVal;
         yMin = yVal;
         first = false;
      } else {
         if (yVal < yMin)
         {
            xMin = xVal;
            yMin = yVal;
         }
      }      
   } // for iter   
   // dumpNamedSTLendl ("xVec", xVec);
   // dumpNamedSTLendl ("yVec", yVec);
   TGraph graph (xVec.size(), &(xVec[0]), &(yVec[0]));
   double yValue = yDelta + yMin;
   // cout << "xMin " << xMin << " yValue " << yValue << endl;
   double upper, lower;
   bool positive = false;
   if (signedStep > 0)
   {
      positive = true;
      upper = xMin + signedStep;
      lower = xMin;
   } else {
      lower = xMin + signedStep;
      upper = xMin;
   }
   double value = graph.Eval (xMin + signedStep, 0, "S") - yValue;
   // cout << "value " << value << " " << yValue << " " << xMin << " " 
   //     << signedStep << endl;
   int count = 0;
   bool contained = false;
   while (fabs (value) > epsilon  || ! contained)
   {
      if (++count >= 40)
      {
         cerr << "breaking" << endl;
         break;
      }
      double upperVal  = graph.Eval (upper,  0, "S") - yValue;
      double lowerVal  = graph.Eval (lower,  0, "S") - yValue;
      if (upperVal * lowerVal < 0)
      {
         contained = true;
      } 
      else
      {
         // we're not bounded yet
         if (positive)
         {
            lower = upper;
            upper = upper + signedStep;
            value = graph.Eval (upper, 0, "S") - yValue;
            // cout << "pos " << value << endl;
         } else {
            upper = lower;
            lower = lower + signedStep;
            value = graph.Eval (lower, 0, "S") - yValue;
            // cout << "neg " << value << endl;
         }
         continue;
      }
      // we're bounded.  Pick a point in the middle
      double middle = (upper + lower) / 2;
      double middleVal = value = graph.Eval (middle,  0, "S") - yValue;
      if (upperVal * middleVal > 0)
      {
         // both are same sign
         upper = middle;
         // cout << "up " << value << endl;
      } else {
         lower = middle;
         // cout << "lo " << value << endl;
      }
   }
   return (upper + lower) / 2;
}

// friends
ostream& operator<< (ostream& o_stream, const CLPTrapezoid &rhs)
{
   o_stream << "[ (" 
            << Form ("%.4f", rhs.lowerX() ) << ", " 
            << Form ("%11.8f", rhs.lowerY() ) << ") , (" 
            << Form ("%.4f", rhs.upperX() ) << ", " 
            << Form ("%11.8f", rhs.upperY() ) << ") ] "
            << Form ("%.8f", rhs.integral() ) << ", "
            << Form ("%.8f", rhs.totalIntegral() );
   if (rhs.lowerY() > rhs.upperY())
   {
      o_stream << "   right " << rhs.position();
   } else {
      o_stream << "   left " << rhs.position();
   }
   return o_stream;
} 

ostream& operator<< (ostream& o_stream, const CLPTrapezoid::Vec &rhs)
{
   int count = 0;
   for (CLPTrapezoid::VecConstIter iter = rhs.begin();
        rhs.end() != iter;
        ++iter)
   {
      o_stream << setw(4) << count++ << ") " << *iter << endl;
   } // for iter   
   return o_stream;
} 
