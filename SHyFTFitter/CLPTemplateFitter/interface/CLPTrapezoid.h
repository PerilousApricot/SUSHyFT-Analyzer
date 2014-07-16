// -*- C++ -*-

#if !defined(CLPTrapezoid_HH)
#define CLPTrapezoid_HH

#include <algorithm>
#include <vector>

#include "TSystem.h"

class CLPTrapezoid : public TObject
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////
      
      typedef std::vector< CLPTrapezoid >  Vec;
      typedef Vec::iterator                VecIter;
      typedef Vec::const_iterator          VecConstIter;
      // typedef Vec::const_reverse_iterator  VecConstRevIter; // bug in GCC
      typedef Vec::reverse_iterator        VecConstRevIter;
      typedef std::vector< double >        DVec;

      enum
      {
         kToTheRight =  1,
         kToTheLeft  = -1,
         kInside     =  0,
      };

      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const CLPTrapezoid &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      CLPTrapezoid (double x, double y);

      bool operator> (const CLPTrapezoid &rhs) const;
      bool operator< (const CLPTrapezoid &rhs) const;

       ////////////////
      // One Liners //
      ////////////////
      double lowerX() const { return m_lowerX; }
      double upperX() const { return m_upperX; }
      double lowerY() const { return m_lowerY; }
      double upperY() const { return m_upperY; }
      double totalIntegral() const { return m_totalIntegral; }
      double minY()   const { return std::min (m_lowerY, m_upperY); }
      double maxY()   const { return std::max (m_lowerY, m_upperY); }
      int position()  const { return m_position; }

      void setLowerX (double lowerX) { m_lowerX = lowerX; }
      void setUpperX (double upperX) { m_upperX = upperX; }
      void setLowerY (double lowerY) { m_lowerY = lowerY; }
      void setUpperY (double upperY) { m_upperY = upperY; }
      void setTotalIntegral (double totalIntegral) 
      { m_totalIntegral = totalIntegral; }
      void setPosition (int position) { m_position = position; }

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      // returns the area of this trapezoid
      double integral() const;

      // scales all y values and integral numbers by 'factor'
      void scale (double factor);

      // returns a negative number if this is left of xVal, a positive
      // number if this is right of xVal, and 0 if this centered on
      // xVal.
      int positionWRT (double xVal) const;

      // returns true if there is an overlap in the Y values
      bool hasVerticalOverlapWith (const CLPTrapezoid &rhs) const;

      // returns the absolute value of the slope
      double absSlope() const;

      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////

      static void limits (Vec &curve,
                          double &upperLimit, double &lowerLimit, 
                          double percent = 0.95);

      // returns the minimum and maximum areas that are allowed if we
      // use these two end points.  Returns 'true' if the two
      // iterators have an overlap..
      static bool minMaxArea (double &minArea, double &maxArea,
                              const VecConstIter &alpha,
                              const VecConstIter &beta);


      // A very rough and inefficient method for finding an Y intersect
      static double findYintersect (const Vec &curve, double signedStep, 
                                    double yDelta, double epsilon = 0.00001);


  private:
      ///////////////////////
      // Private Constants //
      ///////////////////////

      /////////////////////////
      // Private Member Data //
      /////////////////////////
      double m_lowerX;
      double m_upperX;
      double m_lowerY;
      double m_upperY;
      double m_integral;
      double m_totalIntegral;
      int    m_position; // for debugging purposes

      ClassDef (CLPTrapezoid, 1) // CLPTrapezoid Class
};


////////////////////////
// ostream Definition //
////////////////////////
extern ostream& operator<< (ostream& o_stream, const CLPTrapezoid::Vec &right);

#endif // CLPTrapezoid_HH
