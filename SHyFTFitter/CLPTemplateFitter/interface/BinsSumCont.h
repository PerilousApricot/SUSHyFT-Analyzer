// -*- C++ -*-

#if !defined(BinsSumCont_H)
#define BinsSumCont_H

#include "TH1F.h"

#include "CLPTemplateFitter/interface/BinsSum.h"

class CLPTemplateFitter;

class BinsSumCont
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////
      // Containers
      typedef std::map< int, TH1F* >      ITH1FPtrMap;
      typedef std::vector< std::string >  SVec;
      typedef std::vector< double >       DVec;
      // Iterators
      typedef ITH1FPtrMap::iterator       ITH1FPtrMapIter;
      typedef ITH1FPtrMap::const_iterator ITH1FPtrMapConstIter;
      typedef DVec::const_iterator        DVecConstIter;

      enum { kAllTemplates = -1 };

      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const BinsSumCont &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      BinsSumCont();

       ////////////////
      // One Liners //
      ////////////////

      void setFitter (CLPTemplateFitter *fitterPtr)
      { m_fitterPtr = fitterPtr; }

      int size() const { return (int) m_binsSumVec.size(); }

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      // loops over map, deleting any histograms then clearing map.
      // Called automatically when 'setParamVec()' is called.
      void cleanHistMap();

      // set the parameters used
      void setParamVec (DVec &paramVec);

      // add a BinsSum to self
      int add (const BinsSum &obj);

      // loop over all binsSum object and record trials
      void recordTrials();

      // reset all containted binsSums
      void reset();

      // loop over all binsSum and get their current values
      void getValues (DVec &valuesVec);

      // loop over all binsSum and fill means and RMSs
      void fillMeansAndRMSs (DVec &meanVec, DVec &rmsVec);

      // fills an SVec with all names
      void getAllNames (SVec &nameVec);

      // gets from fitter if necessary and returns TH1F*
      TH1F* getTemplate (int index) const;

      // creates a new object with name, puts it in the vector, and
      // returns a reference to object.  Reference is only guaranteed
      // to be good until 'createNew()' or 'add()' is called.
      BinsSum &createNew (const std::string &name);

      // access to individual binsSum
      BinsSum &operator() (int index);

  private:
      /////////////////////////
      // Private Member Data //
      /////////////////////////

      CLPTemplateFitter  *m_fitterPtr;
      BinsSum::DVec       m_paramVec;
      BinsSum::Vec        m_binsSumVec;
      mutable ITH1FPtrMap m_histMap;

};


#endif // BinsSumCont_H
