// -*- C++ -*-

#if !defined(CLPTemplateMorph_HH)
#define CLPTemplateMorph_HH

#include <vector>
#include <string>

#include "TSystem.h"
#include "TH1F.h"

#include "CLPTemplateMorph/interface/CLPBinCont.h"

class CLPTemplateMorph : public TH1F
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////

      // Vectors
      typedef std::vector< double >            DVec;
      typedef std::vector< int >               IVec;
      typedef std::vector< std::string >       SVec;
      typedef std::vector< TGraph* >           TGraphPtrVec;
      typedef std::vector< CLPTemplateMorph* > PtrVec;
      // Iterators
      typedef DVec::const_iterator             DVecConstIter;
      typedef IVec::iterator                   IVecIter;
      typedef IVec::const_iterator             IVecConstIter;
      typedef SVec::const_iterator             SVecConstIter;
      typedef TGraphPtrVec::const_iterator     TGraphPtrVecConstIter;
      typedef PtrVec::iterator                 PtrVecIter;
      typedef PtrVec::const_iterator           PtrVecConstIter;

      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const CLPTemplateMorph &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      CLPTemplateMorph ();
      CLPTemplateMorph (const char* name, const char *title);
      CLPTemplateMorph (const char* name, const char *title,
                        Int_t numBins, Axis_t xlow, Axis_t xhigh);

      ////////////////
      // One Liners //
      ////////////////
      int size() const { return (int) m_paramVec.size(); }

      // sets and returns the value of dontRemovePoints
      bool dontRemovePoints() const { return m_dontRemovePoints; }
      void setDontRemovePoints (bool dontRemovePoints) 
      { m_dontRemovePoints = dontRemovePoints; }

      // if set true, will smooth the beginning of the delta CDFs
      bool smoothDeltas() const { return m_smoothDeltas; }
      void setSmoothDeltas (bool smoothDeltas) 
      { m_smoothDeltas = smoothDeltas; }

      // turns on full output
      void setOutputAll (bool outputAll) { m_outputAll = outputAll; }

      // returns true if this is a container of other CLPTemplateMorphs
      bool isMorphContainer() const { return m_isMorphContainer; }

      // setting and checking noUnderOverflow
      bool noUnderOverflow() const { return m_noUnderOverflow; }
      void setNoUnderOverflow (bool noUnderOverflow = true) 
      { m_noUnderOverflow = noUnderOverflow; }



      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////
      // adds default template.  Can only be done once.  Must be done
      // BEFORE adding any morphing variables or templates.
      void addDefaultTemplate (TH1F *histPtr);

      // adds a new morphing variable.  Must provide name.  'center'
      // describes where the default template should be used.  Returns
      // index.
      int addMorphVariable (const std::string &name, double center = 0.);

      // adds a new template and variable  to a morphing variable.
      void addMorphTemplate (const std::string &name, double value, 
                             TH1F *histPtr);
      void addMorphTemplate (int index, double value, TH1F *histPtr);

      // set the paramVec and recalculates the histogram.  
      void setParamVec (const DVec &paramVec);

      // this function works as the above, except that you must have
      // already called 'setInternalIndicies()'.  The user is taking
      // responsibility that arrayPtr is at least as large as the
      // highest index set.
      void setParamVecFromArray (const double *arrayPtr);

      // returns the index of a name.  -1 if not found
      int morphIndex (const std::string &name) const;

      // sets the internal index on the position of a name in a vector
      // of strings.
      void setInternalIndicies (const SVec &nameVec);

      // take the current parameter values and recalculates the histogram
      void recalculate();

      // fills tgraphPtrVec with all CLPBinConts
      void allTGraphs (TGraphPtrVec &tgraphPtrVec) const;

      // take the CDF for a given variable from another
      // CLPTemplateMorph.  Does not work with Morph containers and
      // 'name' must exist in rhs, but not yet in this morph.
      void takeMorphVariable (const CLPTemplateMorph &rhs,
                              const std::string &name);

      // Adds a Morph to the ptrVec.  This should only be done after
      // the morph being passed in is fully configured.  
      // NOTE: As this only keeps pointers, changes to the morphs that
      // are added will affect this morph.
      void addMorph (CLPTemplateMorph *morphPtr);
      void addMorph (CLPTemplateMorph &morph)
      { addMorph( &morph ); }

      // makes a TH1F* copy of itself.  User is responsible for memory.
      TH1F *th1fCopy () const;

      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////
      static void setParamVecs (PtrVec &morphPtrVec, const DVec &paramVec);

  private:

      //////////////////////////////
      // Private Member Functions //
      //////////////////////////////
      // For debugging purposes
      void _printTGraph (const CLPBinCont &cont) const;

      // fills the TH1F of a morph container from the m_ptrVec
      // histograms.
      void _fillMorphContainer();

      ///////////////////////
      // Private Constants //
      ///////////////////////

      /////////////////////////
      // Private Member Data //
      /////////////////////////
      
      bool       m_dontRemovePoints;
      bool       m_smoothDeltas;
      bool       m_outputAll;
      bool       m_isMorphContainer;
      bool       m_noUnderOverflow;
      DVec       m_paramVec;
      DVec       m_centerVec;
      SVec       m_nameVec;
      IVec       m_indexVec;
      PtrVec     m_ptrVec;
      CLPBinCont m_defaultBC;
      CLPBinCont::DMapVec m_deltaMorphMapVec;
      

      ClassDef (CLPTemplateMorph, 1) // CLPTemplateMorph Class
};



#endif // CLPTemplateMorph_HH
