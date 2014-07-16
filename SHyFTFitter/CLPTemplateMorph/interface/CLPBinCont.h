// -*- C++ -*-

#if !defined(CLPBinCont_HH)
#define CLPBinCont_HH

#include <map>

#include "TSystem.h"
#include "TH1F.h"
#include "TGraph.h"

#include "CLPTemplateMorph/interface/CLPBin.h"
#include "CLPTemplateMorph/interface/CLPFunctionoid2.h"

class CLPBinCont : public TObject
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////

      typedef std::map< double, CLPBinCont > DMap;
      typedef DMap::const_iterator           DMapConstIter;
      typedef std::vector< DMap >            DMapVec;
      typedef DMapVec::const_iterator        DMapVecConstIter;

      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const CLPBinCont &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      CLPBinCont();
      CLPBinCont (TH1F *histPtr);

       ////////////////
      // One Liners //
      ////////////////

      // clears vector
      void clear() { m_vec.clear(); }
      
      // integral
      double integral() const { return m_integral; }

      // returns the size
      int size() const { return (int) m_vec.size(); }

      // sets and returns the value of dontRemovePoints
      bool dontRemovePoints() const { return m_dontRemovePoints; }
      void setDontRemovePoints (bool dontRemovePoints) 
      { m_dontRemovePoints = dontRemovePoints; }
      static void setDontRemovePointsDefault (bool dontRemoveDefault)
      { sm_dontRemovePointsDefault = dontRemoveDefault; }

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      // fills itself with the integral of histPtr
      void integrateTH1F (TH1F *histPtr);

      // Adds a bin to the vector.  If this bin is zero, will check
      // whether to add a new bin or modify the last bin.
      void addBin (double low, double width, double total, bool isZero);
      void addBin (double upper, double total);

      // Adds a bin to the vector MODIFYING the previous bins so that
      // we ensure that our bins only grow in X.
      void addBinCarefully (double low, double width, double total);

      // differentiates self to fill histPtr
      void fillTH1F (TH1F *histPtr, bool noUnderOverFlow = false) const;

      // returns a new TGraph pointer.  User is responsible for the memeory
      TGraph *tgraphPtr() const;

      // Increments the two iterators until the given x value is
      // included in the range of the the upperX() values.  Returns
      // true if both iterators are still valid.
      bool moveIterators (CLPBin::VecConstIter &lowerIter,
                          CLPBin::VecConstIter &upperIter,
                          double xValue, double &yValue) const;

      // sets this object as the horizontal morph between two other
      // BinConts.  'param' = 0 means all alpha, 'param' = 1 means all
      // beta.
      void horizontalMorph (const CLPBinCont &alpha, const CLPBinCont &beta, 
                            double param);

      // sets this object as the difference between two BinCounts.
      void horizontalSubtract (const CLPBinCont &alpha, 
                               const CLPBinCont &beta);

      void horizontalDeltaMorph (const CLPBinCont &alpha,
                                 const CLPBinCont &beta,
                                 double percentage);

      void horizontalFunction (const CLPBinCont &alpha, 
                               const CLPBinCont &beta,
                               const CLPFunctionoid2 &functionoid);

      // sets this object as the sum between a default BinCount and a
      // subtracted one.
      void horizontalAdd (const CLPBinCont &alpha, 
                          const CLPBinCont &delta,
                          double percentage = 1.);

      // gets rid of the first few points to avoid any bumps
      void smoothBeginningOfDelta (double percentStep = 0.05);


      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////

      // returns the number 'param' percent in between alpha and beta
      static double percentInBetween (double alpha, double beta, double param)
      { return alpha + param * (beta - alpha); }

      // returns true if the two numbers are with in kEpsilon of each other
      static bool areTheSame (double alpha, double beta);

  private:
      ///////////////////////
      // Private Constants //
      ///////////////////////

      /////////////////////////
      // Private Member Data //
      /////////////////////////

      CLPBin::Vec m_vec;
      double      m_integral;
      bool        m_noOverflow;
      bool        m_dontRemovePoints;
      bool        m_debug;

      static bool sm_dontRemovePointsDefault;

      ClassDef (CLPBinCont, 1) // CLPBinCont Class
};

#endif // CLPBinCont_HH
