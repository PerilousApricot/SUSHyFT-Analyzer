// -*- C++ -*-

#if !defined(CLPTemplateFitter_HH)
#define CLPTemplateFitter_HH

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include "TMinuit.h"
#include "TH1F.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "CLPTemplateFitter/interface/BinNormClass.h"
#include "CLPTemplateFitter/interface/CLPTrapezoid.h"
#include "CLPTemplateFitter/interface/StartStopStep.h"
#include "CLPTemplateMorph/interface/CLPTemplateMorph.h"
#include "CLPTemplateFitter/interface/BinsSumCont.h"

class CLPTemplateFitter
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////
      // Vectors
      typedef std::vector< TH1F* >            TH1FPtrVec;
      typedef std::vector< int >              IVec;
      typedef std::vector< double >           DVec;
      typedef std::vector< std::string >      SVec;
      typedef std::vector< DVec >             DVecVec;
      // Maps
      typedef std::map< std::string, int >    SIMap;
      typedef std::map< std::string, double > SDMap;
      // DPairs
      typedef std::pair< double, double >     DPair;
      typedef std::set< DPair >               DPairSet;
      // Iterators
      typedef IVec::const_iterator            IVecConstIter;
      typedef DVec::const_iterator            DVecConstIter;
      typedef SVec::const_iterator            SVecConstIter;
      typedef SIMap::const_iterator           SIMapConstIter;
      typedef SDMap::const_iterator           SDMapConstIter;
      typedef DPairSet::const_iterator        DPairSetConstIter;

      enum
      {
         kOutputFitResults  = 1 << 0,
         kOutputCovarMatrix = 1 << 1,
         kOutputBinNorms    = 1 << 2,
         kOutputFunctoids   = 1 << 3,
      };

      /////////////
      // friends //
      /////////////

      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const CLPTemplateFitter &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////

      CLPTemplateFitter();

      ///////////////////
      // The Big Three //
      ///////////////////

      // Copy constructor (NOT COMPLETE)
      CLPTemplateFitter (const CLPTemplateFitter& rhs);

      // Destructor
      ~CLPTemplateFitter();

      // Operator=  (NOT COMPLETE)
      CLPTemplateFitter &operator= 
      (const CLPTemplateFitter& rhs);

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////
      
      // adds data template.  Returns the old TH1F pointer (if any)
      TH1F* addData (TH1F *dataHist);
      const TH1F* getData () const { return m_dataHPtr; }
      void  copyDataFromOtherFitter (const CLPTemplateFitter &otherFitter);
      
      // adds a MC template.  Returns template index
      int addTemplate (const std::string& name, TH1F *tempHist,
                       double norm, double mean = 0., double sigma = 0.,
                       double start = 0., double stop = 0., double step = 0.);

      // adds a bin-by-bin norm.  Returns bin norm's index
      int addBinNorm (const std::string &name, const BinNormClass &binNorm,
                      double mean = 1., double sigma = 0.);

      // adds a functoid to a binNorm (either one binNorm for all
      // templates or a binNorm for only one template.
      void setBinNormFunctoid (const std::string &name, 
                               BinNormFunctoid &functoid);
      void setBinNormFunctoid (const std::string &name, int tempIndex,
                               BinNormFunctoid &functoid);
      void setBinNormFunctoid (int binNormIndex, BinNormFunctoid &functoid);
      void setBinNormFunctoid (int binNormIndex, int tempIndex,
                               BinNormFunctoid &functoid);

      // tell fitter to fix/release a parameter in Minuit
      bool fixParameter     (const std::string &name);
      bool releaseParameter (const std::string &name);

      // Tells you whether a parameter is fixed or not
      bool isParameterFixed (const std::string &name) const;
      bool isParameterFixed (int paramIndex) const;

      // tells fitter to update value of a parameter
      bool setParameter (const std::string &name, double value);

      // sets a constraint mean and sigma.  Negative sigma means
      // unconstrained.
      void setConstraint (const std::string &name, double mean, double sigma);

      // get constraint for a given parameter
      void getConstraint (const std::string &name, double &mean, double &sigma);
      void getConstraint (int paramIndex, double &mean, double &sigma);

      // sets up and runs Minuit minimization to fit templates to data
      void fit();

      // fits letting everything float except the variable given is
      // fixed to the value given.  Returns -2 * ln likelihood.
      double fitEverythingBut (const std::string &name, double value);

      // scans a variable in the given range.  Return value is a
      // DPairSet.
      void scanVariable (CLPTrapezoid::Vec &retval, const std::string &name,
                         double lower, double upper, int numPoint = 1000);

      // returns the value of the log likelihood given the values of
      // the parameter stored in 'parameterVec'.
      double logLikelihoodValue (const DVec &parameterVec);
      double logLikelihoodValue ()
      { return logLikelihoodValue (m_fitVec); }

      // reset Minuit parameters to default starting values
      void resetMinuitParameters();

      // returns a string name of a given parameter
      const std::string &paramName (int param) const;

      // returns a TH1F* to the new histogram.  User is responsible
      // for the memory.  If no name is passed in, then a histogram of
      // the sum of templates will be returned.
      // NOTE: 
      // 1) These member functions are NOT const because in order
      // to process the request, we may need to update the Morphs
      // internally.
      // 2) 'name' is passed by copy instead of by reference because
      // it is modified internally.
      TH1F* updatedHistogram ();
      TH1F* updatedHistogram (std::string name);
      TH1F* updatedHistogram (int tempIndex);
      TH1F* updatedHistogram (const DVec &paramVec);
      TH1F* updatedHistogram (std::string name, const DVec &paramVec);
      TH1F* updatedHistogram (int tempIndex, const DVec &paramVec);

      // returns a TH1F* that is a Fit/Data ratio.  Errors are based
      // off of sqrt(N) in each data bin.  Bins with 0 data events are
      // ignored.  User is responsible for the memory
      TH1F* ratioHistogram () const;

      // returns a TH1F* that is the residual of the fit (data -
      // fit). Bins with 0 data events are ignored.  User is
      // responsible for the memory.
      TH1F* residualHistogram () const;

      // outputs fit result to STDOUT
      void outputFitResults() const;

      // outputs the covariance matrix to STDOUT
      void outputCovarianceMatrix (bool nonZeroOnly = true) const;

      // outputs the covariance matrix to STDOUT
      void outputCorrelationMatrix (bool nonZeroOnly = true) const;

      // removes bins from fitting consideration
      void removeBinFromFit  (int bin);
      void removeBinsFromFit (const BinNormClass::ISet &binset);

      // return the total bin contents.  If no DVec is provided, then
      // the results from the fit are returned.
      double getTotalBinContent (int binNormIndex) const;
      double getTotalBinContent (int binNormIndex, const DVec &paramVec) const;

      // returns the value or error of a parameter after the fit
      double getValue (const std::string &name) const;
      double getError (const std::string &name) const;
      double getPosError (const std::string &name) const;
      double getNegError (const std::string &name) const;

      // returns the covariance matrix as well as fills the SVec as
      // the order of the variables in the covariance matrix.  By
      // default, the function only returns the subset of variables
      // that are not fixed.
      TMatrixD getCovarianceMatrix (SVec &paramNamesVec, 
                                    bool nonZeroOnly = true) const;
      TMatrixD getCovarianceMatrix (IVec &paramIndiciesVec, 
                                    bool nonZeroOnly = true) const;

      // calculates and stores internally square root of covariance
      // matrix as well as indicies of non-zero elements
      void storeSqrtMatrix ();

      // after storeSqrtMatrix has been called, generates correlated
      // random numbers.  Any fixed parameters will not be modified.
      void generateRandomParams (DVec &paramVec) const;

      // returns paramVec with the current values of all parameters,
      // symetric and positive and negeative errors..
      void fillParamVec (DVec &paramVec) const;
      void fillErrorVec (DVec &paramVec) const;
      void fillPosErrorVec (DVec &paramVec) const;
      void fillNegErrorVec (DVec &paramVec) const;

      // Returns SIMap which contains map between parameter name and
      // indices used.
      SIMap getNameIndexMap() const;

      // fills namesVec with all of the names of the parameters in order
      void fillNameVec (SVec &nameVec) const;

      // Given a string, it returns the index of the parameter (-1 if
      // not found).  This function is expected to be somewhat
      // expensive and should be called sparingly.
      int parameterIndex (const std::string &name) const;

      // initialize Minuit
      void _initializeMinuit();

      /////////////////////////////////////
      // Inline Regular Member Functions //
      /////////////////////////////////////

      // When set true, tells 'updatedHistogram()' and
      // 'ratioHistogram()' to use all bins, regardless of whether the
      // bin is used in the fit.
      void setUseAllBinsInPlots (bool useAll = true)
      { m_useAllBinsInPlots = useAll; }
      bool usingAllBinsInPlots () const { return m_useAllBinsInPlots; }

      // sets minuit's verbose level.  Only works if used BEFORE first
      // fit.
      void setMinuitVerboseLevel (int level) { m_minuitVerboseLevel = level; }

      // return the number of parameters used
      int size() const { return m_numTemplates + m_numBinNorms; }
      int numTemplates() const { return m_numTemplates; }
      int numBinNorms()  const { return m_numBinNorms; }

      // setting and checking verbose levels
      int verbose() const           { return m_verbose; }
      int verbose (int flag) const  { return m_verbose & flag; }
      void setVerbose (int verbose) { m_verbose = verbose; }

      // set and check 'doMinos' variable
      bool doMinos() const { return m_doMinos; }
      void setDoMinos (bool doMinos) { m_doMinos = doMinos; }

      // set and value for minExpectedBinContent
      double minExpectedBinContent() const { return m_minExpectedBinContent; }
      void setMinExpectedBinContent (double minExpectedBinContent) 
      { m_minExpectedBinContent = minExpectedBinContent; }

      // returns reference to binsSumCont.  Note that m_binsSumCont is
      // mutable, so it can return a non-const ref for a const
      // function call.
      BinsSumCont &binsSumCont() const { return m_binsSumCont; }


      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////

      // function that Minuit minimizes
      static void minimizeFcn (int &npar, double *gin, double &retval, 
                               double *parameterArray, int iflag);

      // returns log of the probability that we see 'observed' when we
      // expected 'expected'.
      static double logPoisson (double observed, double expected);

      // Return the index for a given name.  Returns -1 if name is not
      // found.  User must provide nameIndexMap (from 'getNameIndexMap').
      static int indexOfName (const std::string &name, 
                              const SIMap &nameIndexMap);

      // returns square root of a matrix
      static void sqrtMatrix (TMatrixD &retval, const TMatrixD &mat);

      // Given square root of covariance matrix and a vector of means,
      // generates correlated random numbers
      static void gaussCorrRand (TVectorD &retVec, const TMatrixD &sqrtCovMat,
                                 const TVectorD &meanVec);


  private:

      //////////////////////////////
      // Private Member Functions //
      //////////////////////////////

      // returns multiplicative factor that fitter has assigned to a
      // bin of a template.
      double _getBinFactor (int tempIndex, int binNormIndex,
                            double *arrayAddress = 0) const;
      double _getBinFactor (int tempIndex, int binNormIndex,
                            const DVec &paramVec) const;
      
      // returns bin contents for a template
      double _getBinContent (int tempIndex, int binNormIndex,
                             double *arrayAddress = 0) const;
      double _getBinContent (int tempIndex, int binNormIndex,
                             const DVec &paramVec) const;
      
      // returns the total bin contents for all templates
      double _getTotalBinContent (int binNormIndex,
                                  double *arrayAddress = 0) const;
      double _getTotalBinContent (int binNormIndex,
                                  const DVec &paramVec) const;

      // returns multiplicative factor by bin
      double _byBin (int tempIndex, int binNormIndex,
                     double *arrayAddress = 0) const;
      double _byBin (int tempIndex, int binNormIndex,
                     const DVec &paramVec) const;

      // returns the index of a parameter given its name.  Return -1
      // if the parameter name doesn't find any matches
      int _paramIndex (const std::string &name, 
                       int &tempIndex, int &binNormIndex) const;
      int _paramIndex (const std::string &name) const
      { int temp, binNorm; return _paramIndex (name, temp, binNorm); }

      // returns the probability of the data fluctuating to the
      // current templates
      double _totalLogProb () const;

      // destroys Minuit's variables
      void _destroyMinuit();

      // number of events in the bin range
      double _eventsInRange (TH1F* histPtr);

      // validates bin index
      void _validateBinIndex (const std::string &function, 
                              int binNormIndex) const;

      // validates template index
      void _validateTemplateIndex (const std::string &function, 
                                   int tempIndex) const;

      // validates minuit is initialized
      void _validateMinuitInit (const std::string &function) const;

      // validates that a histogram is o.k.  Also sets informations
      // about histogram sizes.
      void _validateHistogram (const std::string &function, 
                               TH1F *histPtr);

      // returns the ISet of parameters to use (all if nonZeroOnly is
      // false, just the non-fixed parameter indicies if true)
      void _fillParamISet (BinNormClass::ISet &paramSet, 
                           bool nonZeroOnly) const;

      void _setAllMorphingParameters (const double *parameterArray);

      ///////////////////////////////////////
      // One Line Private Member Functions //
      ///////////////////////////////////////

      // set the address of Minuit parameters array
      void _setArrayAddress (double *arrayAddress)
      { m_arrayAddress = arrayAddress; }

      // returns true if bin is in set of bins being fit
      bool _isBinFit (int bin) const
      { return m_binsSet.end() != m_binsSet.find (bin); }

      /////////////////////////
      // Private Member Data //
      /////////////////////////
      TH1F               *m_dataHPtr;
      TMinuit            *m_minuitPtr;
      mutable double     *m_arrayAddress;
      int                 m_numTemplates;
      int                 m_numBinNorms;
      int                 m_numBins;
      int                 m_numMinuitParams;
      int                 m_minuitVerboseLevel;
      int                 m_debugLevel;
      int                 m_verbose;
      bool                m_useAllBinsInPlots;
      bool                m_doMinos;
      double              m_minExpectedBinContent;
      // below here, all private member data work fine with their
      // default values
      TH1FPtrVec          m_templateHPtrVec;
      DVec                m_normVec;
      DVec                m_fitVec;
      DVec                m_errorVec;
      DVec                m_posErrorVec;
      DVec                m_negErrorVec;
      DVec                m_constraintMeanVec;
      DVec                m_constraintSigmaVec;
      StartStopStep::Vec  m_startStopStepVec;
      DVec                m_binNormConstraintMeanVec;
      DVec                m_binNormConstraintSigmaVec;      
      SVec                m_templateNameVec;
      SVec                m_binNormNameVec;
      DVecVec             m_covarMatrixDVV;
      SIMap               m_templateIndexMap;
      SIMap               m_binNormIndexMap;
      BinNormClass::Vec   m_binNormVec;
      BinNormClass::ISet  m_fixedParameterSet;
      BinNormClass::ISet  m_fixedBinNormSet;
      BinNormClass::ISet  m_binsSet;
      BinNormClass::ISet  m_morphSet;
      mutable BinsSumCont m_binsSumCont;

      // used for generating correlated random numbers
      TMatrixD           m_sqrtCovarMatrix;
      IVec               m_nonZeroElementsVec;

      ////////////////////////////////
      // Private Static Member Data //
      ////////////////////////////////
      static CLPTemplateFitter *sm_tfPtr;
      static const std::string  kEmptyString;
};

#endif // CLPTemplateFitter_HH
