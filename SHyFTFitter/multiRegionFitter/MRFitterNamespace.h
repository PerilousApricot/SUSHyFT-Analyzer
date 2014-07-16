// -*- C++ -*-

#if !defined(MIIIFitterNamespace_HH)
#define MIIIFitterNamespace_HH

#include <string>
#include <vector>
#include <map>
#include <set>

#include "TFile.h"
#include "TH1F.h"

#include "CLPTemplateFitter/interface/CLPTemplateFitter.h"
#include "OptionUtils/interface/CommandLineParser.h"
#include "PolyNoid.h"


namespace mrf
{
   /////////////////////////////////////////////////
   // /////////////////////////////////////////// //
   // // Enums, Typedefs and Constants, Oh My! // //
   // /////////////////////////////////////////// //
   /////////////////////////////////////////////////

   struct BlurStruct
   {
         int    m_index;
         double m_mean;
         double m_sigma;
   };

   struct ScanStruct
   {
         std::string m_name;
         double      m_lower;
         double      m_upper;
         int         m_numSteps;
         double      m_lowerError;
         double      m_upperError;
   };

   // Containers
   typedef std::vector< std::vector< int > >    TwoDimIVec;
   typedef std::vector< TH1F* >                 TH1FPtrVec;
   typedef std::vector< std::string >           SVec;
   typedef std::vector< int >                   IVec;
   typedef std::vector< double >                DVec;
   typedef std::pair< std::string, TH1F* >      STH1FPtrPair;
   typedef std::pair< std::string, double >     SDPair;
   typedef std::vector< STH1FPtrPair >          STH1FPtrPairVec;
   typedef std::vector< SDPair >                SDPairVec;
   typedef std::map< std::string, SDPairVec >   SSDPairVecMap;
   typedef std::vector< BlurStruct >            BlurStructVec;
   typedef std::vector< ScanStruct >            ScanStructVec;
   typedef std::map< std::string,               
                     BlurStruct >               BlurStructMap;
   typedef std::map< std::string, TH1F* >       STH1FPtrMap;
   typedef std::map< std::string, double >      SDMap;
   typedef std::map< std::string, int >         SIMap;
   typedef std::map< std::string, std::string > SSMap;
   typedef std::map< std::string,               
                     std::pair < double,        
                                 double > >     SPDMap;
   typedef std::map< int, int >                 IIMap;
   typedef std::set< std::string >              SSet;
   // Iterators                                 
   typedef TH1FPtrVec::const_iterator           TH1FPtrVecConstIter;
   typedef SVec::const_iterator                 SVecConstIter;
   typedef STH1FPtrPairVec::const_iterator      STH1FPtrPairVecConstIter;
   typedef SDPairVec::const_iterator            SDPairVecConstIter;
   typedef SSDPairVecMap::const_iterator        SSDPairVecMapConstIter;
   typedef STH1FPtrMap::iterator                STH1FPtrMapIter;
   typedef STH1FPtrMap::const_iterator          STH1FPtrMapConstIter;
   typedef BlurStructVec::const_iterator        BlurStructVecConstIter;
   typedef BlurStructMap::const_iterator        BlurStructMapConstIter;
   typedef ScanStructVec::iterator              ScanStructVecIter;
   typedef ScanStructVec::const_iterator        ScanStructVecConstIter;
   typedef SDMap::const_iterator                SDMapConstIter;
   typedef SIMap::iterator                      SIMapIter;
   typedef SIMap::const_iterator                SIMapConstIter;
   typedef SSMap::const_iterator                SSMapConstIter;
   typedef SPDMap::const_iterator               SPDMapConstIter;
   typedef IIMap::const_iterator                IIMapConstIter;
   typedef SSet::const_iterator                 SSetConstIter;

   struct TemplateMorphInfo
   {
         std::string    m_defaultName;
         SDMap          m_defaultValuesMap;
         SSDPairVecMap  m_additionalTemplatesVecMap;
   };

   typedef std::map< std::string, TemplateMorphInfo > STemplateMorphInfoMap;
   typedef STemplateMorphInfoMap::iterator       STemplateMorphInfoMapIter;
   typedef STemplateMorphInfoMap::const_iterator STemplateMorphInfoMapConstIter;

   enum
   {
      kNotAssigned = -1,
   };

   // template names (as stored in histogram)
   //const std::string kDataName  = "Data";
   extern std::string ns_dataName;

   /////////////////////
   // /////////////// //
   // // Variables // //
   // /////////////// //
   /////////////////////

   // vector to hold histogram pointers and names

   // constants
   const DVec     kEmptyDVec;
   const SDMap    kEmptySDMap;

   extern optutl::CommandLineParser ns_parser;

   // dummy variable to force initalization of namespace
   extern int ns_dummy;

   ////////////////////////////////
   // ////////////////////////// //
   // // Function Definitions // //
   // ////////////////////////// //
   ////////////////////////////////


   // blurs paramVec according to the blurVec
   void blurParamVec (DVec &paramVec, const BlurStructVec &blurVec);



   // hook up this namespaces OptionUtils
   void _setupOptions();

   // does any necessary calcuation given command line options and
   // updates 'tag' string for output root file.
   void processOptions (std::string &tag);

   void setDataName();

   // Assumes that outputPtr goes from 0.5 to N + 0.5
   void addToHistogram (TH1F *outputPtr, TH1F *inputPtr,  
                        int offset, int rebin, float scale = 1.);



   void _2vecsToMap (const SVec &namesVec, const DVec &valuesVec,
                     SDMap &valueMap);
   void _2vecsToMap (const SVec &namesVec, const IVec &valuesVec,
                     SIMap &valueMap);
   void _2vecsToMap (const SVec &namesVec, const SVec &valuesVec,
                     SSMap &valueMap);

   void _svecToSDMap (const SVec &nameAndValueVec, SDMap &valueMap);
   
   void _makeEverythingPretty();

   int  _initializeMRF();


   // sets error bars based on sqrt of bin contents
   void _addSqrtErrors (TH1F *histPtr);

   //////////////////////////
   // //////////////////// //
   // // MRFitter Class // //
   // //////////////////// //
   //////////////////////////


   class MRFitter
   {

      public:

         // constructor
         MRFitter (const std::string &configFilename, 
                   const std::string &prefix = "");
                      
         // Fit the branching fraction.
         void fit (bool verbose = true);

         // throw PE.  If 'ideal' is true, then it will add the
         // properly normalized templates to make the PE.
         //
         // Note: The fitter is not passed as a constant reference as
         // it may need to modify itself if template morphing is
         // involved.
         void throwPE (const DVec &paramVec,
                       bool ideal = false);

         // saves a canvas of the fitter result.  Note that we are
         // passing a string and not a constant reference because we
         // will change the value of this.
         void saveCanvasResult (std::string outputName = "output",
                                const DVec &paramVec = kEmptyDVec);

         // draw lines on plots separating jet/tag bins
         void drawLines (double maxHeight, bool unitSpacing = false);


         // test fitter. This function will change over time and is
         // used for debugging.
         void testFitter();

         // Convert fitter output paramVec to tree order.  This works
         // for the fitter's paramVec as well as it error vectors.
         // Note that the vectors must already be the right size.
         void convertFitterToTreeOrder (const DVec &paramVec,
                                        DVec &treeParamVec);


         // setup a paramVec with reasonable variables.  If initialize
         // is set to true, the function will make sure the vector is
         // the right size and has reasonable values for all entries.
         // If set to false, only the 'additionalValueMap' values will
         // be used.
         void setupParamVec (DVec &paramVec,
                             bool initialize = false,
                             const SDMap &additionalValueMap = kEmptySDMap);

         // Using fitter and user-supplied inputs, decides which variables
         // will be varied for PE throwing.
         void setupBlurStruct (BlurStructVec &blurVec, 
                               const BlurStructMap &blurMap);
         void setupBlurStruct (BlurStructVec &blurVec, const SVec &names, 
                               const DVec &means, const DVec &sigmas);

         BinsSumCont &binsSumCont() { return m_fitter.binsSumCont(); }


         // Grab the "data" from the fitter that generates the PEs and
         // pull it into this fitter.
         void getDataFromOtherFitter (const MRFitter &otherFitter);
         

         ////////////////////////
         // One line functions //
         ////////////////////////

         const SVec &toFillNames() const
         { return m_toFillNames; }

         int numToFill() const
         { return m_numToFill; }       

         const ScanStructVec &scanStructVec() const
         { return m_scanStructVec; }

         ///////////////////////////////////////
         // Promoting Fitter Member Functions //
         ///////////////////////////////////////

         // At some point, it probably makes sense to have MRFitter
         // derive from CLPTemplateFitter.  In the mean time...

         void fillParamVec (DVec &paramVec) 
         { m_fitter.fillParamVec (paramVec); }

         void fillErrorVec (DVec &errorVec) 
         { m_fitter.fillErrorVec (errorVec); }

         void fillPosErrorVec (DVec &errorVec) 
         { m_fitter.fillPosErrorVec (errorVec); }

         void fillNegErrorVec (DVec &errorVec) 
         { m_fitter.fillNegErrorVec (errorVec); }

         void storeSqrtMatrix()
         { m_fitter.storeSqrtMatrix(); }

         void fillNameVec (SVec &nameVec) const
         { m_fitter.fillNameVec (nameVec); }

         void generateRandomParams (DVec &paramVec) const
         { m_fitter.generateRandomParams (paramVec); }

         void scanVariable (CLPTrapezoid::Vec &retval, 
                            const std::string &name,
                            double lower, double upper, int numPoints = 1000)
         { m_fitter.scanVariable (retval, name, lower, upper, numPoints); }

         // returns the value or error of a parameter after the fit
         double getValue (const std::string &name) const
         { return m_fitter.getValue (name); }
         double getError (const std::string &name) const
         { return m_fitter.getError (name); }
         double getPosError (const std::string &name) const
         { return m_fitter.getPosError (name); }
         double getNegError (const std::string &name) const
         { return m_fitter.getNegError (name); }


         // Given a string, it returns the index of the parameter (-1 if
         // not found).  This function is expected to be somewhat
         // expensive and should be called sparingly.
         int parameterIndex (const std::string &name) const
         { return m_fitter.parameterIndex (name); }

         //////////////////////
         // Access Functions //
         //////////////////////
         const std::string &configFilename() const 
         { return m_configFilename; } 
         const std::string &templateFilename() const 
         { return m_templateFilename; }
         const std::string &prefix() const 
         { return m_prefix; }
         

         int         &integerValue  (std::string key) 
         { return m_parser.integerValue (key); }
         double      &doubleValue   (std::string key) 
         { return m_parser.doubleValue  (key); }
         std::string &stringValue   (std::string key) 
         { return m_parser.stringValue  (key); }
         bool        &boolValue     (std::string key) 
         { return m_parser.boolValue    (key); }
         IVec        &integerVector (std::string key) 
         { return m_parser.integerVector(key); }
         DVec        &doubleVector  (std::string key) 
         { return m_parser.doubleVector (key); }
         SVec        &stringVector  (std::string key) 
         { return m_parser.stringVector (key); }



      private:

         // load configuration file
         bool _loadConfigFile (const std::string &configFilename, 
                              bool mainLoop = true);

         // initialize the fitter according to the mode flag
         void _initializeFitter();



         // Figure out which variables will be saved when PEs are run.
         // Needs to be called after 'initializeFitter()'
         void _decideWhichVariablesToSave();

         void _loadTemplates (const std::string &filename);
         void _hookupFitter ();
         void _createBinsSums ();
   
         // Converts '_'s to ' 's
         void _fixGroupStrings();

         // adds polynoid to internal map
         void _addPolyToMap (const PolyNoid &poly);
         // provides access to these polynoids
         PolyNoid &_poly (const std::string &name);
         // Adds polynoids to the fitter
         void _addPolyMapToFitter();
         // print out morphMap
         void _printMorphMap();

         // rescales defaults
         void _rescaleDefaults();

         //////////////////
         // Data Members //
         //////////////////

         CLPTemplateFitter m_fitter;
         std::string m_configFilename;
         std::string m_prefix;
         std::string m_templateFilename;
         STH1FPtrMap m_histPtrMap;
         SDMap       m_defaultValuesMap;
         SIMap       m_colorMap;
         SDMap       m_constraintMap;
         SIMap       m_nameFitterIndexMap;
         SIMap       m_nameTemplateIDMap;    // name -> enum
         SIMap       m_nameTemplateIndexMap; // name -> fitter index
         SIMap       m_groupIndexMap;
         SIMap       m_groupBinsMap;
         SDMap       m_groupLowerEdgeMap;
         SDMap       m_groupUpperEdgeMap;

         // place to hold binning information
         IVec        m_numBinsVec;
         IVec        m_lowerEdgeBinVec;
         IVec        m_upperEdgeBinVec;
         int         m_totalBins;

         // other variables
         SIMap              m_ignoredNamesMap;
         BinNormClass::ISet m_allBinsSet;
         bool               m_fixedGroupStrings;
         ScanStructVec      m_scanStructVec;

         // names and indicies to be filled for PEs
         SVec        m_toFillNames;
         IVec        m_toFillIndicies;
         int         m_numToFill;
         optutl::CommandLineParser m_parser;

         PolyNoid::Map m_polynoidMap;
         bool          m_canAddToPolyMap;
         IIMap         m_polyTemplate2indexMap;
         SSet          m_usesTemplateMorphing;
         STemplateMorphInfoMap m_morphMap;         
   };

};

#endif // MIIIFitterNamespace_HH

