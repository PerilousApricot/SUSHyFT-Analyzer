// -*- C++ -*-

#if !defined(CLPFunctionoid2Slide_HH)
#define CLPFunctionoid2Slide_HH

#include "TSystem.h"

#include "CLPTemplateMorph/interface/CLPFunctionoid2.h"

class CLPFunctionoid2Slide : public CLPFunctionoid2
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////

      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const CLPFunctionoid2Slide &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      CLPFunctionoid2Slide (double percentage);

       ////////////////
      // One Liners //
      ////////////////

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////
      virtual double operator() (double alpha, double beta) const;

      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////

  protected:
      /////////////////////////
      // Protected Constants //
      /////////////////////////

      ///////////////////////////
      // Protected Member Data //
      ///////////////////////////

  private:
      ///////////////////////
      // Private Constants //
      ///////////////////////

      /////////////////////////
      // Private Member Data //
      /////////////////////////
      double m_percentage;

};

#endif // CLPFunctionoid2Slide_HH
