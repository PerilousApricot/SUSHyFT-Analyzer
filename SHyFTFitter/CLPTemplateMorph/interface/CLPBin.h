// -*- C++ -*-

#if !defined(CLPBin_HH)
#define CLPBin_HH

#include "TSystem.h"

#include <vector>

class CLPBin : public TObject
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////
      
      typedef std::vector< CLPBin >    Vec;
      typedef Vec::iterator            VecIter;
      typedef Vec::const_iterator      VecConstIter;

      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const CLPBin &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      CLPBin();
      CLPBin (double xLower, double xWidth, double yValue, bool sameAs = false);

       ////////////////
      // One Liners //
      ////////////////
      void setXLowerWidth (double xlower, double width)
      { m_xWidth = width; m_xUpper = xlower + width; }
      void setXUpperWidth (double xupper, double width)
      { m_xWidth = width; m_xUpper = xupper; }
      void setYvalue (double yvalue) { m_y = yvalue; }
      void setSameAsLastBin (bool sameAs = true) 
      { m_sameAsLastBin = sameAs; }
      

      double xUpper()  const { return m_xUpper; }
      double xMiddle() const { return m_xUpper - m_xWidth / 2.; }
      double xLower()  const { return m_xUpper - m_xWidth; }
      double xWidth()  const { return m_xWidth; }
      double yValue()  const { return m_y; }
      bool   isSameAsLastBin() const { return m_sameAsLastBin; }

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////

  private:
      ///////////////////////
      // Private Constants //
      ///////////////////////

      /////////////////////////
      // Private Member Data //
      /////////////////////////
      double m_xUpper;
      double m_xWidth;
      double m_y;
      bool   m_sameAsLastBin;

      ClassDef (CLPBin, 1) // CLPBin Class
};

#endif // CLPBin_HH
