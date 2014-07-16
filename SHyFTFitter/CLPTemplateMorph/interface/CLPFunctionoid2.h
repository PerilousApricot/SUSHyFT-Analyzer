// -*- C++ -*-

#if !defined(CLPFunctionoid2_HH)
#define CLPFunctionoid2_HH

#include "TSystem.h"

class CLPFunctionoid2
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
                                       const CLPFunctionoid2 &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      CLPFunctionoid2();

      virtual ~CLPFunctionoid2() {};

       ////////////////
      // One Liners //
      ////////////////

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      virtual double operator() (double alpha, double beta) const
      { return beta - alpha; }

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

};

#endif // CLPFunctionoid2_HH
