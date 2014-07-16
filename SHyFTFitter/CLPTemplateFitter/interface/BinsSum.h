// -*- C++ -*-

#if !defined(BinsSum_H)
#define BinsSum_H

#include <vector>
#include <map>
#include <string>

class BinsSumCont;

class BinsSum
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////
      // Containers
      typedef std::vector< BinsSum >      Vec;
      typedef std::vector< int >          IVec;
      typedef std::vector< double >       DVec;
      // Iterators
      typedef Vec::iterator               VecIter;
      typedef Vec::const_iterator         VecConstIter;
      typedef IVec::const_iterator        IVecConstIter;
      typedef DVec::const_iterator        DVecConstIter;

      /////////////
      // friends //
      /////////////

      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const BinsSum &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      BinsSum (std::string name, BinsSumCont *containerPtr);

      ////////////////
      // One Liners //
      ////////////////

      // return name
      const std::string &name() const { return m_name; }

      // add a template index
      void addTemplateIndex (int index) 
      { m_templateIndiciesVec.push_back (index); }

      // add a bin range
      void addBinRange (int lower, int upper)
      { m_lowerBinVec.push_back (lower); m_upperBinVec.push_back (upper); }

      // reset counting
      void reset() { m_trialsVec.clear(); }

      // calculate value and record as trial
      void recordTrial() { m_trialsVec.push_back( value() ); }

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      // calculate and return value given this paramVec
      double value() const;

      // calculate mean and RMS from trials
      void calcMeanAndRMS (double &mean, double &rms) const;

      // adds self to container.  Returns index
      int addSelfToContainer();

      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////

      static void calcMeanAndRMS (double &mean, double &rms, const DVec &xvec);

  private:
      /////////////////////////
      // Private Member Data //
      /////////////////////////
      std::string  m_name;
      IVec         m_templateIndiciesVec;
      IVec         m_lowerBinVec;
      IVec         m_upperBinVec;
      DVec         m_trialsVec;

      // useful pointers
      BinsSumCont *m_contPtr;

};


#endif // BinsSum_H
