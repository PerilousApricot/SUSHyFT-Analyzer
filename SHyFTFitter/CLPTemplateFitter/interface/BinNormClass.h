// -*- C++ -*-

#if !defined(BinNormClass_HH)
#define BinNormClass_HH

#include <iostream>
#include <set>
#include <vector>
#include <map>

#include "CLPTemplateFitter/interface/BinNormFunctoid.h"
#include "CLPTemplateFitter/interface/StartStopStep.h"

class BinNormClass
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////

      static const double kMinValue;

      typedef std::set< int >             ISet;
      typedef std::vector< BinNormClass > Vec;
      typedef ISet::iterator              ISetIter;
      typedef ISet::const_iterator        ISetConstIter;
      typedef Vec::const_iterator         VecConstIter;

      enum {
         kNegativeNumber = -1
      };

      enum {
         kNumParams = 5
      };

      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const BinNormClass &rhs);

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      BinNormClass ();
      BinNormClass (double startingValue);
      BinNormClass (double startingValue, 
                    double start, double stop, double step);
      BinNormClass (const ISet &tempSet, const ISet &binSet, 
                    double startingValue);
      BinNormClass (const ISet &tempSet, const ISet &binSet, 
                    double startingValue, 
                    double start, double stop, double step);

       ////////////////
      // One Liners //
      ////////////////

      // Setting functions
      void setBinSet  (const ISet &binSet)  { m_binSet      = binSet;  }
      void setParamIndex    (int index)     { m_paramIndex  = index;   }
      void setStartingValue (double value)  { m_value       = value;   }
      void setVerbose (bool verbose) const  { m_verbose     = verbose; }

      // access functions
      double startingValue() const { return m_value; }
      int    paramIndex()    const { return m_paramIndex; }

      double start() const { return m_startStopStep.start(); }
      double stop()  const { return m_startStopStep.stop(); }
      double step()  const { return m_startStopStep.step(); }
      void setStart (double start) { m_startStopStep.setStart (start); }
      void setStop  (double stop)  { m_startStopStep.setStop  (stop); }
      void setStep  (double step)  { m_startStopStep.setStep  (step); }
      void setStartStopStep (double start, double stop, double step)
      { m_startStopStep.setStartStopStep (start, stop, step); }

      const StartStopStep &startStopStep() const
      { return m_startStopStep; }

      // returns functoid pointer
      const BinNormFunctoid *functoidPtr () const { return m_functoidPtr; }

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      // returns the index (paramIndex) if template and bin IDs are
      // contained, else kNegativeNumber
      int index (int tempIndex, int binIndex) const;

      // returns the factor by which this bin should be multiplied
      double factor (int tempIndex, int binIndex, 
                     double value, const double *paramArray) const;

      // adds tempSet
      void setTempSet (const ISet &tempSet);

      // set one functoid for all templates
      void setFunctoid (BinNormFunctoid &functoid);    

      // set a functoid for only one template
      void setFunctoid (int tempIndex, BinNormFunctoid &functoid);    

      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////

      // insert a range into a set
      static void insertRangeIntoSet (ISet &theset, int lower, int upper);

      // prints contents of ISet to STDOUT
      static void outputISet (const ISet &theset, 
                              std::ostream& o_stream = std::cout);

      // // returns true if the set containts the element
      // static bool setContains (const ISet &theSet, int element)
      // {  return ( theSet.end() != theSet.find (element) ); }

  private:
      /////////////////////////
      // Private Member Data //
      /////////////////////////
      ISet   m_templateSet;
      ISet   m_binSet;
      double m_value;
      int    m_paramIndex;      
      StartStopStep          m_startStopStep;
      mutable bool           m_verbose;

      BinNormFunctoid       *m_functoidPtr;
      BinNormFunctoid::IFMap m_functoidMap;
};

#endif // BinNormClass_HH
