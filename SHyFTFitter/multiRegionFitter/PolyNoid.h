// -*- C++ -*-

#if !defined(Polynoid_HH)
#define Polynoid_HH

#include <vector>
#include <map>
#include <string>

#include "CLPTemplateFitter/interface/BinNormFunctoid.h"
#include "CLPTemplateFitter/interface/CLPTemplateFitter.h"

class PolyNoid : public BinNormFunctoid
{
   public:

      ////////////////////////
      // Typedefs and Enums //
      ////////////////////////
      typedef std::map< std::string, PolyNoid > Map;
      typedef std::vector< double >             DVec;
      typedef std::vector< int >                IVec;
      typedef std::vector< DVec >               DVecVec;
      typedef std::pair< int, int >             IIPair;
      typedef std::map< IIPair, int>            IIPairIMap;
      typedef std::map< int, int >              IIMap;
      typedef Map::iterator                     MapIter;
      typedef Map::const_iterator               MapConstIter;
      typedef DVec::const_iterator              DVecConstIter;
      typedef IVec::const_iterator              IVecConstIter;
      typedef IIMap::const_iterator             IIMapConstIter;
      typedef IIPairIMap::const_iterator        IIPairIMapConstIter;

      /////////////
      // friends //
      /////////////

      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const PolyNoid &rhs);

      ////////////////////////////////
      // Constructors and Operators //
      ////////////////////////////////
      PolyNoid();
      PolyNoid (std::string &name, double mean, double width,
                double lower = 0., double upper = 0., double step = 0.);

      // operator()
      virtual double operator() (int tempIndex, int binIndex,
                                 double value, 
                                 const double *paramArray = 0) const;

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      //////////////////////
      // Inline Functions //
      //////////////////////

      // access
      const std::string &name()  const { return  m_name;    }
      double mean()  const { return  m_mean;    }
      double width() const { return  m_width;   }
      double lower() const { return  m_lower;   }
      double upper() const { return  m_upper;   }
      double step()  const { return  m_step;    }

      // setting
      void setName  (const std::string &name) { m_name  = name;  }
      void setMean  (double mean)             { m_mean  = mean;  }
      void setWidth (double width)            { m_width = width; }
      void setLower (double lower)            { m_lower = lower; }
      void setUpper (double upper)            { m_upper = upper; }
      void setStep  (double step)             { m_step  = step;  }


      
      // add a polynomial for a given set of group bins and templates
      // (described by enum)
      void addPolynomial (const IVec &groupVec,
                          const IVec &templateEnumVec,
                          const DVec &parameters);
      
      // given information in template2indexMap and m_temporaryMap,
      // calculates m_groupTemplateParamMap
      void calculateGroupTemplateParamMap (const IIMap &template2indexmap);

      // // adds self to polynoidMap.  Can only be done before fitter is
      // // hooked up (will abort if tried afterwards).
      // void addSelfToMap();

      // adds self to a fitter.  Returns index.
      int addSelfToFitter (CLPTemplateFitter &fitter);
      
      // print out information about this class
      virtual void dumpToOstream (ostream& o_stream) const;



      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////
      

      // sets jet bin widths (i.e., a vector of widths: e.g., 20 20 10
      // 10 ...)
      static void setGroupBinWidths (const IVec &widthVec);

      // returns the polynomial evaluated at x.
      static double polynomial (double x, const DVec &params);

      // given a bin, returns jetbin index
      static int bin2groupBin (int bin);

      // returns upper and lower bin boundaries for a given group bin
      static void lowerUpperBoundaries (int &lower, int &upper, int groupBin);

      // returns reference to bin bounaries vector
      static const IVec &groupBinBoundariesVec();

      // returns number of last bin
      static int  lastBin();

   private:

      /////////////////////////
      // Enums and constants //
      /////////////////////////

      enum 
      {
         kError   = -1,
      };      


      /////////////////
      // Member Data //
      /////////////////

      std::string m_name;
      double      m_mean;
      double      m_width;
      double      m_lower;
      double      m_upper;
      double      m_step;

      // Holds all polynomials
      DVecVec        m_parametersVec;
      // holds map between template index stored in fitter and
      // parameters in 'm_paramsVecVec'
      IIPairIMap     m_groupTemplateParamMap;
      // temporarily holds map between enum of template and parameters
      // (used before all fitter templates have been added)
      IIPairIMap     m_temporaryMap;

      ////////////////////////
      // Static Member Data //
      ////////////////////////

      // static int    sm_numJets;
      // static IIMap  sm_template2indexMap;

      // variables needed for bin to group bin conversion
      static IVec   sm_groupBinBoundariesVec;
      static int    sm_lastBin;
      static int    sm_lastGroupBin;
      static int    sm_groupBinUpper;


};

#endif // Polynoid_HH


