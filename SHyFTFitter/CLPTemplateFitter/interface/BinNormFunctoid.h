// -*- C++ -*-

#if !defined(BinNormFunctoid_HH)
#define BinNormFunctoid_HH

#include <map>

#include "TSystem.h"

class BinNormFunctoid : public TObject
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////

      typedef double (*FuncPtr) (int tempIndex, int binIndex,
                                 double value, const double *paramArray);
      typedef std::map< int, BinNormFunctoid* >  IFMap;
      typedef IFMap::const_iterator              IFMapConstIter;


      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const BinNormFunctoid &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      BinNormFunctoid ();
      BinNormFunctoid (FuncPtr ptr);

       ////////////////
      // One Liners //
      ////////////////

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      // function used by fitter
      virtual double operator() (int tempIndex, int binIndex,
                                 double value, const double *paramArray) const;

      // print out information about this class
      virtual void dumpToOstream (ostream& o_stream) const;

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

      FuncPtr m_funcPtr; //!

      ClassDef (BinNormFunctoid, 1) // BinNormFunctoid Class


};

#endif // BinNormFunctoid_HH
