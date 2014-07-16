// -*- C++ -*-

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLine.h"
#include "TText.h"
#include "THStack.h"
#include "TLegend.h"
#include "TRandom.h"

#include "MRFitterNamespace.h"
#include "OptionUtils/interface/dout.h"
#include "OptionUtils/interface/CommandLineParser.h"
#include "CLPTemplateMorph/interface/dumpSTL.icc"

using namespace std;

//////////////////////////
// Variable Definitions //
//////////////////////////

optutl::CommandLineParser 
mrf::ns_parser ("config.file [options]",
                optutl::CommandLineParser::kNoDefaultArgs);

string mrf::ns_dataName = "Data";
                 
// Initialize the namespace
int mrf::ns_dummy = mrf::_initializeMRF();

//////////////////////////
// Function Definitions //
//////////////////////////

mrf::MRFitter::MRFitter (const std::string &configFilename,
                         const std::string &prefix) : 
   m_configFilename (configFilename), m_prefix (prefix),
   m_totalBins (0), m_fixedGroupStrings (false),
   m_numToFill (0), m_parser (mrf::ns_parser),
   m_canAddToPolyMap (true)
{
   cout << "creating MRFitter object";
   if (m_prefix.length())
   {
      cout << " (" << m_prefix << ")";
   }   
   cout << endl << "loading configuration: " << configFilename << endl;
   _loadConfigFile (configFilename);
   // // before returning, make sure the 'Data' is one of the names of
   // // the samples
   SVec &namesVec = stringVector("names");
   if ( namesVec.end() ==
        find (namesVec.begin(),
              namesVec.end(),
              ns_dataName) )
   {
      namesVec.insert( namesVec.begin(), 1, ns_dataName );
   }
   _initializeFitter();
   // WARNING:: We now want to tell the PolyNoids about our binning.
   // This is somewhat dubious as I may call the static member
   // function multiple times.  This is okay _IFF_ all of my MRFitter
   // objects all have the same binning.  Since this is the case for
   // the foreseeable future, it should be fine, but I thought I
   // should leave a warning here just in case.
   if (boolValue("print"))
   {
      m_parser.printOptionValues();
   }
}

bool
mrf::MRFitter::_loadConfigFile (const std::string &configFilename, 
                                bool mainLoop)
{
   using namespace optutl;
   if (mainLoop)
   {
      // If we're here, then we're the main version of this subroutine
      // call.  Check anything that needs to be checked exactly once
      // before loading configuration files.
      m_ignoredNamesMap.clear();
      for (SVecConstIter iter = 
              stringVector ("ignoreInclude").begin();
           stringVector ("ignoreInclude").end() != iter;
           ++iter)
      {
         m_ignoredNamesMap[ *iter ] = 0;
      } // for iter
   }
   string fullConfigName (configFilename);
   if (stringValue ("inputPrefix").length())
   {
      fullConfigName = stringValue("inputPrefix") + "/" + configFilename;
   }
   ifstream source (fullConfigName.c_str(), ios::in);
   if (! source)
   {
      cerr << "file " << fullConfigName << " could not be opened" << endl;
      exit (-1);            
   }
   string line;
   string polynoidName;
   while (getline (source, line))
   {
      string command, rest;
      CommandLineParser::findCommand (line, command, rest);
      if (! command.length())
      {
         // blank line
         continue;
      }
      if ('#' == command.at(0))
      {
         // comment - ignore
         continue;
      }
      if ("-" == command)
      {
         // function declaration
         SVec words;
         CommandLineParser::split (words, rest);
         if (words.size() < 3)
         {
            dout << "loadConfigFile() error: I don't understand line: "
                 << line << endl
                 << "Not enough words." << endl;
            exit (-1);
         }
         polynoidName = words.at(0);
         double mean  = atof( words.at(1).c_str() );
         double width = atof( words.at(2).c_str() );
         if (words.size() >= 6)
         {
            // The user has given us the bounds.  Let's use them.
            double   lower = atof( words.at(3).c_str() );
            double   upper = atof( words.at(4).c_str() );
            double   step  = atof( words.at(5).c_str() );
            PolyNoid poly (polynoidName, mean, width, lower, upper, step);
            _addPolyToMap (poly);
         } else {
            // Let PolyNoid determine the bounds and step
            PolyNoid poly (polynoidName, mean, width);
            _addPolyToMap (poly);
         }
         continue;
      } // if -
      if ("--" == command)
      {
         // function definition
         if (! polynoidName.length() )
         {
            dout << "loadConfigFile() error: You can't add polynomials with"
                 << " out first defining: " << endl
                 << "a PolyNoid" << endl
                 << line << endl
                 << "Not enough words." << endl;
            exit (-1);
         }
         SVec words;
         CommandLineParser::split (words, rest, ":");
         if (words.size() != 3)
         {
            dout << "loadConfigFile() error: I don't understand line: "
                 << line << endl
                 << "Not enough words." << endl;
            exit (-1);
         } // if too small
         string groupString  = words.at (0);
         string samplesString = words.at (1);
         string polyString    = words.at (2);

         CommandLineParser::removeLeadingAndTrailingSpaces (groupString);
         CommandLineParser::removeLeadingAndTrailingSpaces (samplesString);
         CommandLineParser::removeLeadingAndTrailingSpaces (polyString);   

         // group bins
         SVec groupWords;
         CommandLineParser::split (groupWords, groupString);
         IVec groupVec;
         for (SVecConstIter iter = groupWords.begin();
              groupWords.end() != iter;
              ++iter)
         {
            bool found = false;
            int groupBin = 0;
            for (SVecConstIter groupIter = 
                    stringVector("groupNames").begin();
                 stringVector("groupNames").end() != groupIter;
                 ++groupIter,++groupBin)
            {
               if (*iter == *groupIter)
               {
                  found = true;
                  break;
               }
            }
            if (! found)
            {
               // dumpNamedSTLendl ("groups: ", stringVector("groupNames") );
               dout << "loadConfigFile() error: I don't understand line: "
                    << line << endl
                    << "'" << *iter << "'" << endl;;
               exit (-1);
            }
            groupVec.push_back (groupBin);
         } // for iter
         // samples
         SVec samplesWords;
         IVec samplesVec;
         CommandLineParser::split (samplesWords, samplesString);
         for (SVecConstIter iter = samplesWords.begin();
              samplesWords.end() != iter;
              ++iter)
         {
            SIMapConstIter siIter = m_nameTemplateIDMap.find( *iter );
            if (m_nameTemplateIDMap.end() == siIter)
            {
               int index = m_nameTemplateIDMap[ *iter ] = 
                  (int) m_nameTemplateIDMap.size();
               samplesVec.push_back( index );
            } else {
               samplesVec.push_back( siIter->second );            
            }
         } // for iter
         // polynomial parameters
         SVec polyWords;
         DVec polyVec;
         CommandLineParser::split (polyWords, polyString);
         for (SVecConstIter iter = polyWords.begin();
              polyWords.end() != iter;
              ++iter)
         {
            polyVec.push_back( atof( iter->c_str() ) );
         }         
         _poly (polynoidName).addPolynomial (groupVec,
                                             samplesVec,
                                             polyVec);
         continue;            
      } // if --
      if ("@" == command)
      {
         // include file
         SVec filenames;
         CommandLineParser::split (filenames, rest);
         for (SVecConstIter iter = filenames.begin();
              filenames.end() != iter;
              ++iter)
         {
            // check to make sure I'm not supposed to ignore this file
            SIMapIter ignoreIter = m_ignoredNamesMap.find (*iter);
            if (m_ignoredNamesMap.end() != ignoreIter)
            {
               // ignore this one
               ignoreIter->second += 1;
               cout << "Ignoring '" << ignoreIter->first << "'." << endl;
               // O.k.  Skip it.
               continue;
            }               
            // I currently have no protection against file A load file
            // B and file B loading file A.  It would be easy enough
            // to implement.
            cout << "loading " << *iter << endl;
            if (! _loadConfigFile (*iter, false))
            {
               dout << "Failed to load " << *iter << ".  Aborting."
                    << endl;
               exit (0);
            }
         } // for iter
         continue;
      } // if *
      if ("+" == command)
      {
         // we're setting a variable
         string nospaces;
         for (int position = 0; 
              position < (int) rest.length(); 
              ++position)
         {
            char ch = rest.at (position);
            if ('#' == ch)
            {
               // start of a comment
               break;
            } else if (' ' == ch || '\t' == ch)
            {
               continue;
            }
            nospaces += ch;
         } // for position
         // set the variable from the string UNLESS this variable has
         // already been set (=> ', true')
         if (! m_parser.setVariableFromString (nospaces, true) )
         {
            dout << "Don't understand line" << endl << line << endl
                 << "in options file '" << configFilename << "'.  Aborting."
                 << endl;
            exit(0);
         } // if setting variable failed
         continue;
      } // if +
      dout << "Don't understand line" << endl << line << endl
           << "in options file '" << configFilename << "'.  Aborting."
           << endl;
      exit(0);
   } // while getline
   if (mainLoop)
   {
      // If we're here, then we're the main version of this subroutine
      // call.  Check anything that needs to be checked exactly once
      // after loading configuration files.
      
      // Make sure we ignored all config files that we were asked to:
      bool problem = false;
      for (SIMapConstIter iter = m_ignoredNamesMap.begin();
           m_ignoredNamesMap.end() != iter;
           ++iter)
      {
         if ( ! iter->second)
         {
            problem = true;
            cerr << "Error: Did not ignore '" << iter->first
                 << "' as requested." << endl;
         }
      } // for iter
      if (problem)
      {
         cerr << "Aborting." << endl;
         exit(0);
      } // if problem

      // load any additional configuration files requested now
      for (SVecConstIter iter = stringVector ("includeFiles").begin();
           stringVector ("includeFiles").end() != iter;
           ++iter)
      {
         if (! _loadConfigFile (*iter, false))
         {
            dout << "Failed to load " << *iter << ".  Aborting."
                 << endl;
            exit (0);
         }         
         cout << "loading " << *iter << endl;
      } // for iter
   }

   return true;
}


void
mrf::MRFitter::fit (bool verbose)
{
   m_fitter.resetMinuitParameters();
   m_fitter.fit();
   if (verbose)
   {
      m_fitter.outputFitResults();
   }
   if (boolValue ("showCorrelations"))
   {
      m_fitter.outputCorrelationMatrix();
   }
   if (verbose || boolValue ("showCorrelations"))
   {
      cout << endl;
   }
   // Do we scan any variables?
   if (m_scanStructVec.size())
   {
      for (ScanStructVecIter ssvIter = m_scanStructVec.begin();
           m_scanStructVec.end() != ssvIter;
           ++ssvIter)
      {
         double fitValue = getValue (ssvIter->m_name);
         CLPTrapezoid::Vec scan;
         scanVariable (scan, 
                       ssvIter->m_name, 
                       ssvIter->m_lower + fitValue, 
                       ssvIter->m_upper + fitValue,
                       ssvIter->m_numSteps);
         double step = (ssvIter->m_upper - ssvIter->m_lower) / 
            ssvIter->m_numSteps;
         ssvIter->m_upperError = 
            CLPTrapezoid::findYintersect (scan,  step, 1.) - fitValue;
         ssvIter->m_lowerError = 
            CLPTrapezoid::findYintersect (scan, -step, 1.) - fitValue;
         if (verbose)
         {
            cout << "scanning " << ssvIter->m_name << endl;
            cout << scan << endl << endl;
            cout << fitValue 
                 << " +" << ssvIter->m_upperError 
                 << " "  << ssvIter->m_lowerError
                 << endl;
         }
      } // for ssvIter
   } // if anything to scan
}


void
mrf::MRFitter::throwPE (const DVec &paramVec, bool ideal)
{
   TH1F *dataHPtr = m_histPtrMap[ns_dataName];
   dataHPtr->Reset();
   // we want to get the total template for a given set of parameters:
   TH1F *templateHPtr = m_fitter.updatedHistogram (paramVec);
   assert (templateHPtr);
   if (ideal)
   {
      // Throw an 'ideal PE'
      dataHPtr->Add (templateHPtr);
   } else {
      // generate PEs from throwing random numbers
      int numEvents = gRandom->Poisson( templateHPtr->Integral() );
      // // DEBUG ONLY
      // int numEvents = templateHPtr->Integral();
      for (int loop = 0; loop < numEvents; ++loop)
      {
         dataHPtr->Fill( templateHPtr->GetRandom() );
      } // for loop
   } // else if ideal
   delete templateHPtr;
}


void
mrf::_setupOptions()
{
   using namespace optutl;
   ns_parser.addOption ("templateFile",
                        VariableMapCont::kString,
                        "Name of root files containing templates",
                        "");
   ns_parser.addOption ("intLumi",
                        VariableMapCont::kDouble,
                        "Integrated Luminosity (fb^-1)",
                        1.0);
   ns_parser.addOption ("lowerAxisValue",
                        VariableMapCont::kDouble,
                        "Lower axis value (needed for morphing)",
                        0.0);
   ns_parser.addOption ("upperAxisValue",
                        VariableMapCont::kDouble,
                        "Upper axis value (needed for morphing)",
                        0.0);
   ns_parser.addOption ("minuitVerbose",
                        VariableMapCont::kInteger,
                        "Minuit verbosity level",
                        -1);
   ns_parser.addOption ("doMinos",
                        VariableMapCont::kBool,
                        "Whether to run Minos",
                        false);
   ns_parser.addOption ("showCorrelations",
                        VariableMapCont::kBool,
                        "Show correlations matrix",
                        false);
   ns_parser.addOption ("minPredictionBin",
                        VariableMapCont::kDouble,
                        "Minimum value prediction in a bin can have",
                        1e-4);
   ns_parser.addOption ("dataErrorBars",
                        VariableMapCont::kBool,
                        "(Incorrectly) put error bars on data",
                        false);
   ns_parser.addOption ("useManyTemplates",
                        VariableMapCont::kBool,
                        "Use the many templates (e.g.,'Wbb_bb' instead of 'Wbb'",
                        false);
   ns_parser.addOption ("outputFitterSetup",
                        VariableMapCont::kBool,
                        "Output initial fitter setup",
                        false);
   ns_parser.addOption ("dontPlotData",
                        VariableMapCont::kBool,
                        "Do not plot data points",
                        false);
   ns_parser.addOption ("logPlots",
                        VariableMapCont::kBool,
                        "Use log Y axis for plots",
                        false);
   ns_parser.addOption ("rescaleDefaults",
                        VariableMapCont::kBool,
                        "Rescale the defaults so total prediction is equal to total data ",
                        false);
   ns_parser.addOption ("imageExtention",
                        VariableMapCont::kString,
                        "Extention of graphics files",
                        "png");
   ns_parser.addOption ("plotTitle",
                        VariableMapCont::kString,
                        "Title of plot ",
                        "SHyFT Fit");
   ns_parser.addOption ("saveTemplates",
                        VariableMapCont::kBool,
                        "Whether or not to  store templates (only with savePlots)",
                        false);
   ns_parser.addOption ("fixParameters",
                        VariableMapCont::kStringVector,
                        "Names of parameters to fix");
   ns_parser.addOption ("fixParamVal",
                        VariableMapCont::kStringVector,
                        "Names and value of parameters to fix");
   ns_parser.addOption ("names",
                        VariableMapCont::kStringVector,
                        "Names of templates");
   ns_parser.addOption ("groupNames",
                        VariableMapCont::kStringVector,
                        "Names of groups");
   ns_parser.addOption ("excludeGroups",
                        VariableMapCont::kStringVector,
                        "Name of groups to exclude from fit");
   ns_parser.addOption ("groupStrings",
                        VariableMapCont::kStringVector,
                        "Pretty Names of groups");
   ns_parser.addOption ("endings",
                        VariableMapCont::kStringVector,
                        "Possible endings of templates");
   ns_parser.addOption ("includeFiles",
                        VariableMapCont::kStringVector,
                        "Additional configuration files to include");
   ns_parser.addOption ("ignoreInclude",
                        VariableMapCont::kStringVector,
                        "Names of include files to NOT load");
   ns_parser.addOption ("binsSum",
                        VariableMapCont::kStringVector,
                        "Information for BinsSums Objects");
   ns_parser.addOption ("binsPerGroup",
                        VariableMapCont::kIntegerVector,
                        "Number of bins in each group (per group)");       
   ns_parser.addOption ("constraintNames",
                        VariableMapCont::kStringVector,
                        "constraint Names");   
   ns_parser.addOption ("constraintValues",
                        VariableMapCont::kDoubleVector,
                        "constraint Values");  
   ns_parser.addOption ("defaultNames",
                        VariableMapCont::kStringVector,
                        "default Names");      
   ns_parser.addOption ("defaultValues",
                        VariableMapCont::kDoubleVector,
                        "default Values");     
   ns_parser.addOption ("colorNames",
                        VariableMapCont::kStringVector,
                        "color Names");        
   ns_parser.addOption ("colorValues",
                        VariableMapCont::kIntegerVector,
                        "color Values");       
   ns_parser.addOption ("rangeNames",
                        VariableMapCont::kStringVector,
                        "range Names");        
   ns_parser.addOption ("lowerRange",
                        VariableMapCont::kDoubleVector,
                        "Lower range values");       
   ns_parser.addOption ("upperRange",
                        VariableMapCont::kDoubleVector,
                        "Upper range values");       
   ns_parser.addOption ("stepRange",
                        VariableMapCont::kDoubleVector,
                        "Step range values");       
   ns_parser.addOption ("legendCoordinates",
                        VariableMapCont::kDoubleVector,
                        "legend coordinates");
   ns_parser.addOption ("morphNames",
                        VariableMapCont::kStringVector,
                        "Names of templates that use template morphing");
   ns_parser.addOption ("morphDefaultTemplates",
                        VariableMapCont::kStringVector,
                        "Histogram name for each default template");
   ns_parser.addOption ("morphDefaultValues",
                        VariableMapCont::kStringVector,
                        "Default value for default template for given systematic");
   ns_parser.addOption ("morphAddTemplates",
                        VariableMapCont::kStringVector,
                        "Additional templates:values for particular "
                        "sysemtatic uncertainty");
   ns_parser.addOption ("lowerEdgeVector",
                        VariableMapCont::kStringVector,
                        "Lower edge vector");
   ns_parser.addOption ("upperEdgeVector",
                        VariableMapCont::kStringVector,
                        "Upper edge vector");
   ns_parser.addOption ("lowerEdge",
                        VariableMapCont::kString,
                        "Lower edge");
   ns_parser.addOption ("upperEdge",
                        VariableMapCont::kString,
                        "Upper edge");
   ns_parser.addOption ("ScanVariables",
                        VariableMapCont::kStringVector,
                        "name:lower:upper:steps of variables to scan");
   ns_parser.addOption ("inputPrefix",
                        VariableMapCont::kString,
                        "Prefix to prepend to all input files",
                        "");
   ns_parser.addOption ("dataName",
                        VariableMapCont::kString,
                        "Name of Data Histogram",
                        "Data");
   ns_parser.stringVector ("endings").push_back ("");
}

void
mrf::processOptions (string &tag)
{
}


void
mrf::setDataName ()
{ 
   ns_dataName = ns_parser.stringValue ("dataName");
}


void 
mrf::addToHistogram (TH1F *outputPtr, TH1F *inputPtr,  
                     int offset, int rebin, float scale)
{
   int numBins = inputPtr->GetNbinsX();
   for (int bin = 1; bin <= numBins; ++bin)
   {
      int tobin = (bin - 1) / rebin + 1 + offset;
      outputPtr->Fill( (float) tobin, scale * inputPtr->GetBinContent(bin) );
   } // for bin   
}


void
mrf::MRFitter::_initializeFitter()
{
   // Convert pairs of vectors to maps
   _2vecsToMap (stringVector ("constraintNames"), 
                doubleVector ("constraintValues"), 
                m_constraintMap);
   _2vecsToMap (stringVector ("defaultNames"),    
                doubleVector ("defaultValues"),    
                m_defaultValuesMap);
   _2vecsToMap (stringVector ("colorNames"),      
                integerVector ("colorValues"),      
                m_colorMap);
   // First figure out what is going to use morphing
   for (SVecConstIter iter = stringVector("morphNames").begin();
        stringVector("morphNames").end() != iter;
        ++iter)
   {
      // make sure this sample exists.  Unfortunately, I do not have a
      // map or a set of names, so I'm going to just loop through the
      // vector of names each time.  This certainly isn't the quick
      // way of doing this.
      SVecConstIter result = find (stringVector ("names").begin(),
                                   stringVector ("names").end(),
                                   *iter);
      if ( result == stringVector ("names").end())
      {
         // messed this one up
         cerr << "Error: There is no group '" << *iter << 
            "' to morph.  Aborting." << endl;
         assert (0);
      }
      m_usesTemplateMorphing.insert( *iter );
   } // for iter
   // Make sure stringVector ("groupNames") exist and m_groupIndexMap is setup.
   if (! stringVector ("groupNames").size())
   {
      stringVector ("groupNames").push_back ("");
   }
   int count = 0;
   SSMap lowerEdgeStringMap, upperEdgeStringMap;
   if (integerVector ("binsPerGroup").size())
   {
      // if we have it, we'll set everything up now.  If not, we'll
      // get it from the data.
      _2vecsToMap (stringVector ("groupNames"),
                   integerVector ("binsPerGroup"), 
                   m_groupBinsMap);   
   }
   if (stringVector ("groupNames").size() == 
       stringVector ("upperEdgeVector").size())
   {
      _2vecsToMap (stringVector ("groupNames"),
                   stringVector ("upperEdgeVector"), 
                   upperEdgeStringMap);   
   }
   if (stringVector ("groupNames").size() == 
       stringVector ("lowerEdgeVector").size())
   {
      _2vecsToMap (stringVector ("groupNames"),
                   stringVector ("lowerEdgeVector"), 
                   lowerEdgeStringMap);   
   }
   for (SVecConstIter iter = stringVector ("groupNames").begin();
        stringVector ("groupNames").end() != iter;
        ++iter)
   {
      m_groupIndexMap[ *iter ] = count++;
      SSMapConstIter upperIter = upperEdgeStringMap.find( *iter );
      SSMapConstIter lowerIter = lowerEdgeStringMap.find( *iter );
      if (upperEdgeStringMap.end() != upperIter &&
          lowerEdgeStringMap.end() != lowerIter)
      {
         m_groupUpperEdgeMap[ *iter ] = atof( upperIter->second.c_str() );
         m_groupLowerEdgeMap[ *iter ] = atof( lowerIter->second.c_str() );
      } else {    
         m_groupUpperEdgeMap[ *iter ] = atof( stringValue("upperEdge").c_str() );
         m_groupLowerEdgeMap[ *iter ] = atof( stringValue("lowerEdge").c_str() ); 
      }
      // if the bin map isn't filled, fill it
   }
   _loadTemplates (stringValue ("templateFile"));
   PolyNoid::setGroupBinWidths (integerVector("binsPerGroup"));
   _makeEverythingPretty();
   if (boolValue ("rescaleDefaults"))
   {
      _rescaleDefaults();
   }
   _hookupFitter ();
   _decideWhichVariablesToSave();
   m_fitter.setDoMinos (boolValue ("doMinos"));
   _createBinsSums ();
   // take care of any groups that aren't supposed to be in the fit
   BinNormClass::ISet nofitSet;
   const IVec &binBoundariesVec = PolyNoid::groupBinBoundariesVec();
   for (SVecConstIter iter = 
           stringVector ("excludeGroups").begin();
        stringVector ("excludeGroups").end() != iter;
        ++iter)
   {
      SIMapConstIter groupIter = m_groupIndexMap.find (*iter);
      if (m_groupIndexMap.end() == groupIter)
      {
         cerr << "Error: There is no group '" << *iter << 
            "' to exclude.  Aborting." << endl;
         exit (0);
      }
      int lowerEdge = binBoundariesVec.at( groupIter->second ) + 1;
      int upperEdge = binBoundariesVec.at( groupIter->second + 1 );
      BinNormClass::insertRangeIntoSet (nofitSet, lowerEdge, upperEdge);
   } 
   if (! nofitSet.empty())
   {
      m_fitter.removeBinsFromFit (nofitSet);
   }    
   // setup scan structs if requested
   for (SVecConstIter svIter = stringVector("scanVariables").begin();
        stringVector("scanVariables").end() != svIter;
        ++svIter)
   {
      // split into pieces on ':' and make sure we have the right
      // number of entries.
      SVec pieces;
      optutl::CommandLineParser::split (pieces, *svIter, ":");
      if (pieces.size() != 4)
      {
         dout << "Don't understand '" << *svIter 
              << "' in 'scanVariables': ";
         dumpSTLendl (stringVector("scanVariables"));
         exit (0);
      }
      ScanStruct entry;
      entry.m_name     = pieces.at (0);
      entry.m_lower    = atof( pieces.at(1).c_str() );
      entry.m_upper    = atof( pieces.at(2).c_str() );
      entry.m_numSteps = atoi( pieces.at(3).c_str() );
      m_scanStructVec.push_back( entry );
   } // for iter
   
}


void
mrf::MRFitter::saveCanvasResult  (std::string outputName,
                                  const DVec &paramVec)
{
   const string kExtention = (string) "." + stringValue ("imageExtention");
   // we don't expect the output name to have any extention (either
   // '.root' or '.png', but it's better to be safe than sorry:
   outputName = 
      optutl::CommandLineParser::removeEnding (outputName, 
                                               kExtention);
   outputName = 
      optutl::CommandLineParser::removeEnding (outputName, 
                                               ".root");
   TCanvas canvas ("canvas", "canvas", 800, 600);
   if (boolValue ("logPlots"))
   {
      canvas.SetLogy();
   }
   const int kNumTemplates = m_fitter.numTemplates() + 1; // + 1 for data
   TH1FPtrVec  histPtrVec (kNumTemplates, (TH1F*) 0);
   // get the histograms:
   // For data, we want to clone it since we're going to be deleting
   // all of the histograms at the end.
   int dataIndex = 0;
   TH1F *dataHist = histPtrVec.at(dataIndex) =  
      (TH1F*) m_fitter.getData()->Clone ("newdata");
   // don't try and get data, so start at 1
   for (int histIndex = 1; histIndex < kNumTemplates; ++histIndex)
   {
      const string &name = stringVector ("names").at (histIndex);
      if (paramVec.size())
      {
         // use the vector passed in
         histPtrVec.at (histIndex) = m_fitter.updatedHistogram (name, paramVec);
      } else {
         // use the last fit
         histPtrVec.at (histIndex) = m_fitter.updatedHistogram (name);
      }
      TH1F *histPtr = histPtrVec.at (histIndex);
      histPtr->SetFillColor( m_colorMap [name] ); 
      histPtr->SetLineColor (kBlack);
   } // for histIndex
   THStack stack ("Stack", stringValue ("plotTitle").c_str());
   for (int histIndex = kNumTemplates - 1; histIndex > 0; --histIndex)
   {
      stack.Add( histPtrVec.at (histIndex) );
   }
   double maxHeight = stack.GetMaximum();
   maxHeight *= 1.5;
   if (boolValue ("logPlots"))
   {
      maxHeight *= 20;
   }
   stack.SetMaximum (maxHeight);
   stack.Draw();
   if (boolValue ("dataErrorBars"))
   {
      _addSqrtErrors (dataHist);
   }
   dataHist->SetMarkerStyle (20);
   if (! boolValue ("dontPlotData"))
   {
      dataHist->Draw ("p same");      
   }
   // save templates to root file
   if (boolValue ("saveTemplates"))
   {
      TFile *filePtr = 
         TFile::Open ( (outputName + "_templates.root").c_str(),
                       "recreate");
      assert (filePtr);
      for (TH1FPtrVecConstIter iter = histPtrVec.begin();
           histPtrVec.end() != iter;
           ++iter)
      {
         CLPTemplateMorph *morphPtr = dynamic_cast< CLPTemplateMorph* > (*iter);
         if ( morphPtr )
         {
            // Yup.  This is a morph.
            TH1F *tempHist = morphPtr->th1fCopy();
            tempHist->Write();
            delete tempHist;
         } else {
            (*iter)->Write();
         }
      }
      // save the data too unless we were asked not to
      if (! boolValue ("dontPlotData"))
      {
         dataHist->Write();
      }
      filePtr->Close();
      delete filePtr;
   } // if saving templates
   drawLines (maxHeight);
   if (doubleVector ("legendCoordinates").size() < 4)
   {
      doubleVector ("legendCoordinates").clear();
      doubleVector ("legendCoordinates").push_back (0.3);
      doubleVector ("legendCoordinates").push_back (0.55); 
      doubleVector ("legendCoordinates").push_back (0.40);
      doubleVector ("legendCoordinates").push_back (0.85);
   }
   TLegend legend ( doubleVector ("legendCoordinates").at(0),
                    doubleVector ("legendCoordinates").at(1),
                    doubleVector ("legendCoordinates").at(2),
                    doubleVector ("legendCoordinates").at(3) );
   if (! boolValue ("dontPlotData"))
   {
      legend.AddEntry (dataHist, ns_dataName.c_str(), "p");
   }
   for (int histIndex = 1; histIndex < kNumTemplates; ++histIndex)
   {
      legend.AddEntry (histPtrVec.at (histIndex),
                       stringVector ("names").at (histIndex).c_str(), "f");
   } // for histIndex
   legend.SetFillColor(kWhite);
   legend.Draw();
   canvas.Print( (outputName + kExtention).c_str() );   

   // Make Njet plots
   TH1FPtrVec  histNjetPtrVec (kNumTemplates, (TH1F*) 0);
   int njetBins = (int) stringVector ("groupNames").size();
   if (njetBins < 2)
   {
      // don't bother
      return;
   }
   for (int histIndex = 0; histIndex < kNumTemplates; ++histIndex)
   {
      string histName = stringVector ("names").at (histIndex) + "Njet";
      TH1F* histPtr = histNjetPtrVec.at (histIndex) =
         new TH1F (histName.c_str(), histName.c_str(),
                   njetBins, 0, njetBins);
      if (stringVector ("names").at (histIndex) != ns_dataName)
      {
         // set colors
         histPtr->SetFillColor( m_colorMap [stringVector ("names").at (histIndex)] ); 
         histPtr->SetLineColor (kBlack);
      } else {
         // set marker for data
         histPtr->SetMarkerStyle (20);         
      }
      int bin = 0;
      const TH1F *bigHistPtr = histPtrVec.at (histIndex);
      int numGroups = (int) stringVector ("groupNames").size();
      const IVec &binBoundariesVec = PolyNoid::groupBinBoundariesVec();
      for (int groupBin = 0; groupBin < numGroups; ++groupBin)
      {
         int lower, upper;
         lower = binBoundariesVec.at (groupBin);
         upper = binBoundariesVec.at (groupBin + 1);
         histPtr->SetBinContent (++bin, 
                                 bigHistPtr->Integral (lower,
                                                       upper));
      } // for groupBin
   } // for histIndex
   
   THStack njetStack ("njetStack", (stringValue ("plotTitle") + 
                                    " N Jet Distributions").c_str());
   for (int histIndex = kNumTemplates - 1; histIndex > 0; --histIndex)
   {
      njetStack.Add( histNjetPtrVec.at (histIndex) );
   }
   maxHeight = njetStack.GetMaximum();
   maxHeight *= 1.2;
   njetStack.SetMaximum (maxHeight);
   njetStack.Draw();
   if (boolValue ("dataErrorBars"))
   {
      _addSqrtErrors( histNjetPtrVec.at (dataIndex) );
   }
   if (! boolValue ("dontPlotData"))
   {
      histNjetPtrVec.at (dataIndex)->Draw ("p same");      
   }
   drawLines (maxHeight, true);
   canvas.Print( (outputName + "_njet" + kExtention).c_str() );      
}


void
mrf::MRFitter::drawLines (double maxHeight, bool unitSpacing)
{
   double textFactor = 1.;
   if (boolValue ("logPlots"))
   {
      textFactor /= 20;
   }
   _fixGroupStrings();
   const IVec &binBoundariesVec = PolyNoid::groupBinBoundariesVec();
   assert (binBoundariesVec.size());
   int numGroups = (int) stringVector ("groupNames").size();
   if (numGroups < 2)
   {
      // If there's only one group, then don't bother with lines and
      // text, etc.
      return;
   }
   for (int groupBin = 0; groupBin < numGroups; ++groupBin)
   {
      double lower, upper;
      lower = (double) binBoundariesVec.at (groupBin);
      upper = (double) binBoundariesVec.at (groupBin + 1);
      if (unitSpacing)
      {
         lower = (double) groupBin;
         upper = lower + 1;
      }
      TText *textPtr = new TText ( (upper + lower) / 2.,
                                   0.6 * maxHeight * textFactor,
                                   stringVector ("groupStrings").at(groupBin).c_str() );
      textPtr->SetTextAlign (12);
      textPtr->SetTextAngle (90);
      textPtr->Draw();
      TLine *linePtr = new TLine (upper, 0., upper, maxHeight);
      // if (kNumJetBins != jetBin)
      // {
      linePtr->SetLineStyle (3);
      // }
      linePtr->Draw();
   } // for groupBin
}


void
mrf::MRFitter::testFitter()
{
   // get the top template as it is by default
   DVec paramVec;   
   const string kTtbarName = "Top";
   TCanvas canvas;
   setupParamVec (paramVec, true);
   TH1F *top1hPtr = m_fitter.updatedHistogram (kTtbarName, paramVec);
   // now find the test variable and it it to 1 instead of 0
   SIMapConstIter iter = m_nameFitterIndexMap.find ("test");
   assert (m_nameFitterIndexMap.end() != iter);
   paramVec.at( iter->second ) = 1;
   TH1F *top2hPtr = m_fitter.updatedHistogram (kTtbarName, paramVec);
   top2hPtr->Divide (top1hPtr);
   top2hPtr->SetTitle ("Ratio Stair Plot");
   top2hPtr->Draw();
   top2hPtr->SetMinimum (0.);
   drawLines (top2hPtr->GetMaximum());   
   canvas.Print ("test.png");
}


void
mrf::MRFitter::convertFitterToTreeOrder (const DVec &paramVec,
                                         DVec &treeParamVec)
{
   for (int loop = 0; loop < m_numToFill; ++loop)
   {
      treeParamVec.at( loop ) = paramVec.at( m_toFillIndicies.at( loop ) );
   } // for loop
}

void
mrf::MRFitter::setupParamVec (DVec &paramVec,
                              bool initialize,
                              const SDMap &additionalValueMap)
{
   if (initialize)
   {
      m_fitter.fillParamVec (paramVec);
   }
   for (SDMapConstIter valueIter = additionalValueMap.begin();
        additionalValueMap.end() != valueIter;
        ++valueIter)
   {
      SIMapConstIter indexIter = m_nameFitterIndexMap.find (valueIter->first);
      assert (m_nameFitterIndexMap.end() != indexIter);
      paramVec.at( indexIter->second ) = valueIter->second;
   } // for valueIter
}


void
mrf::MRFitter::setupBlurStruct (BlurStructVec &blurVec, 
                                const BlurStructMap &blurMap)
{
   blurVec.clear();
   
   int size = m_fitter.size();
   for (int paramIndex = 0; paramIndex < size; ++paramIndex)
   {
      if (m_fitter.isParameterFixed (paramIndex))
      {
         // we don't mess with fixed parameters
         continue;
      }
      double mean, sigma;
      m_fitter.getConstraint (paramIndex, mean, sigma);
      if (sigma > 0)
      {
         BlurStruct blur;
         blur.m_mean  = mean;
         blur.m_sigma = sigma;
         blur.m_index = paramIndex;
         blurVec.push_back (blur);
         continue;
      }
      // if we're still here, then this parameter is not fixed or
      // constraint.  Let's see if the user has passed us any useful
      // information?
      const string &name = m_fitter.paramName (paramIndex);
      BlurStructMapConstIter iter = blurMap.find (name);
      if (blurMap.end() != iter)
      {
         BlurStruct blur = iter->second;
         blur.m_index = paramIndex;
         blurVec.push_back (blur);
         continue;
      }
      // if we're here that means that this parameter is not Gaussian
      // constrained and that we haven't explicitly asked for it to be
      // varied in the PEs.  In this case, we just leave it alone and
      // do nothing.
   } // for paramIndex
   // print everything out
   cout << "Throwing:" << endl;
   for (BlurStructVecConstIter iter = blurVec.begin();
        blurVec.end() != iter;
        ++iter)
   {
      const string &name = m_fitter.paramName (iter->m_index);
      cout << Form ("%10s (%2d) %6.3f +- %6.3f",
                    name.c_str(), iter->m_index, iter->m_mean, iter->m_sigma)
           << endl;
   } // for iter
}


void
mrf::MRFitter::setupBlurStruct (BlurStructVec &blurVec, 
                                const SVec &names, 
                                const DVec &means, 
                                const DVec &sigmas)
{
   assert (names.size() == means.size() &&
           names.size() == sigmas.size());
   BlurStructMap blurMap;
   int size = (int) names.size();
   for (int loop = 0; loop < size; ++loop)
   {
      BlurStruct blur;
      blur.m_mean  = means.at(loop);
      blur.m_sigma = sigmas.at(loop);
      blur.m_index = -1; // doesn't really matter
      blurMap[ names.at(loop) ] = blur;
   } // for loop
   setupBlurStruct (blurVec, blurMap);
}


void
mrf::blurParamVec (DVec &paramVec, const BlurStructVec &blurVec)
{
      for (BlurStructVecConstIter iter = blurVec.begin();
        blurVec.end() != iter;
        ++iter)
   {
      paramVec.at( iter->m_index ) =
         gRandom->Gaus( iter->m_mean, iter->m_sigma);
   } // for iter
}

void
mrf::MRFitter::getDataFromOtherFitter (const mrf::MRFitter &otherFitter)
{
   m_fitter.copyDataFromOtherFitter( otherFitter.m_fitter );
}


void
mrf::MRFitter::_loadTemplates (const string &templateFilename)
{
   using namespace optutl;
   cout << "loading templates:     " << templateFilename << endl << endl;
   m_templateFilename = templateFilename;
   string fullTemplateName (templateFilename);
   if (stringValue ("inputPrefix").length())
   {
      fullTemplateName = stringValue("inputPrefix") + "/" + templateFilename;
   }
   TFile *filePtr = TFile::Open( fullTemplateName.c_str() );
   if (! filePtr)
   {
      cerr << "Can not open file: '" << templateFilename << ".'  Aborting.";
      assert (0);
   }
   TH1::AddDirectory (false);
   // First, let's use the 'data' histogram to get the sizes of
   // everything if we haven't already.
   if (! m_totalBins)
   {
      for (SVecConstIter groupIter = 
              stringVector ("groupNames").begin();
           stringVector ("groupNames").end() != groupIter;
           ++groupIter)
      {
         string histName = ns_dataName + *groupIter;
         TH1F *histPtr = (TH1F*) filePtr->Get( histName.c_str() );
         if (! histPtr)
         {
            dout << "data: " << histName << endl;
         }
         assert (histPtr);
         int numBins = histPtr->GetNbinsX();
         m_numBinsVec.push_back( numBins );
         m_lowerEdgeBinVec.push_back( m_totalBins + 1 );
         m_upperEdgeBinVec.push_back( m_totalBins + numBins );
         m_totalBins += numBins;
      } // for groupIter
      BinNormClass::insertRangeIntoSet (m_allBinsSet, 1, 
                                        m_totalBins);
      // if we don't have the integer vector binsPerGroup already
      // filled, let's fill it now.
      if (! integerVector ("binsPerGroup").size() )
      {
         integerVector ("binsPerGroup") = m_numBinsVec;
         _2vecsToMap (stringVector ("groupNames"),
                      integerVector ("binsPerGroup"), 
                      m_groupBinsMap);   
      }
   } // if not totalBins
   ///////////////////////////////////////////
   // Load all non-morphing templates first //
   ///////////////////////////////////////////
   for (SVecConstIter groupIter = 
           stringVector ("groupNames").begin();
        stringVector ("groupNames").end() != groupIter;
        ++groupIter)
   {
      int groupIndex = m_groupIndexMap[ *groupIter ];
      int offset     = m_lowerEdgeBinVec.at(groupIndex) - 2;
      int numBins    = m_numBinsVec.at(groupIndex);
      for (SVecConstIter nameIter = 
              stringVector ("names").begin();
           stringVector ("names").end() != nameIter;
           ++nameIter)
      {
         // if this name is being used in morphing, don't load it now
         if (m_usesTemplateMorphing.end() != 
             m_usesTemplateMorphing.find( *nameIter ))
         {
            continue;
         } // don't do it now
         for (SVecConstIter endingIter = 
                 stringVector ("endings").begin();
              stringVector ("endings").end() != endingIter; 
              ++endingIter)
         {
            string histName = *nameIter + *groupIter + *endingIter;
            TH1F *histPtr = (TH1F*) filePtr->Get( histName.c_str() );
            if (! histPtr)
            {
               if ( ! m_usesTemplateMorphing.size() )
               {
                  cerr << "Warning: '" << histName << "' is missing." << endl;
               }
               continue;
            }
            assert (histPtr);
            assert (numBins == histPtr->GetNbinsX());
            double scale = doubleValue ("intLumi");
            if (ns_dataName == *nameIter)
            {
               scale = 1;
            }
            if (! scale)
            {
               scale = 1.;
            }
            string newHistName = *nameIter  + *endingIter;
            // check to make sure that we have already made this histogram
            TH1F *newHistPtr = 0;
            STH1FPtrMapIter histPtrIter = m_histPtrMap.find (newHistName);
            if (m_histPtrMap.end() == histPtrIter)
            {
               // didn't find it
               newHistPtr = m_histPtrMap [newHistName] =  
                  new TH1F (newHistName.c_str(), 
                            newHistName.c_str(), 
                            m_totalBins, 0, 
                            m_totalBins);
            } else {
               // found it
               newHistPtr = histPtrIter->second;
            }
            addToHistogram (newHistPtr, histPtr, offset, 1, scale);
         } // for endingIter
      } // for nameIter
   } // for groupIter
   //////////////////////////////////////////////
   // Now load all template-morphing templates //
   //////////////////////////////////////////////
   // Before anything, see if we need to bother
   if ( ! m_usesTemplateMorphing.size())
   {
      delete filePtr;
      return;
   }
   // now, initialize STMIMap object
   const TemplateMorphInfo kEmptyTMI;
   for (SSetConstIter utmIter = m_usesTemplateMorphing.begin();
        m_usesTemplateMorphing.end() != utmIter;
        ++utmIter)
   {
      for (SVecConstIter groupIter = 
              stringVector ("groupNames").begin();
           stringVector ("groupNames").end() != groupIter;
           ++groupIter)
      {
         const string fullName = *utmIter + *groupIter;
         m_morphMap[ fullName ] = kEmptyTMI;
         // unless overridden, the default name should be what it
         // would have been if there is no morphing.
         m_morphMap[ fullName ].m_defaultName = fullName;
      } // for groupIter
   } // for utmIter
   // Default templates
   for (SVecConstIter mdtIter = stringVector("morphDefaultTemplates").begin();
        stringVector("morphDefaultTemplates").end() != mdtIter;
        ++mdtIter)
   {
      SVec pieces;
      CommandLineParser::split (pieces, *mdtIter, ":");
      if (pieces.size() != 2)
      {
         dout << "Don't understand '" << *mdtIter 
              << "' in 'morphDefaultTemplates': ";
         dumpSTLendl (stringVector("morphDefaultTemplates"));
         exit (0);
      }
      STemplateMorphInfoMapIter tmiIter = m_morphMap.find( pieces.at(0) );
      if (m_morphMap.end() == tmiIter)
      {
         dout << "Don't understand '" << *mdtIter 
              << "' in 'morphDefaultTemplates': ";
         dumpSTLendl (stringVector("morphDefaultTemplates"));
         exit (0);
      }
      tmiIter->second.m_defaultName = pieces.at(1);
   } // for mdtIter
   for (SVecConstIter mdvIter = stringVector("morphDefaultValues").begin();
        stringVector("morphDefaultValues").end() != mdvIter;
        ++mdvIter)
   {
      SVec pieces;
      CommandLineParser::split (pieces, *mdvIter, ":");
      if (pieces.size() != 3)
      {
         dout << "Don't understand '" << *mdvIter 
              << "' in 'morphDefaultValues': ";
         dumpSTLendl (stringVector("morphDefaultValues"));
         exit (0);
      }
      STemplateMorphInfoMapIter tmiIter = m_morphMap.find( pieces.at(0) );
      if (m_morphMap.end() == tmiIter)
      {
         dout << "Don't understand '" << *mdvIter 
              << "' in 'morphDefaultValues': ";
         dumpSTLendl (stringVector("morphDefaultTemplates"));
         exit (0);
      }
      tmiIter->second.m_defaultValuesMap[ pieces.at(1) ] =
         atof( pieces.at(2).c_str() );
   } // for mdvIter
   for (SVecConstIter matIter = stringVector("morphAddTemplates").begin();
        stringVector("morphAddTemplates").end() != matIter;
        ++matIter)
   {
      SVec pieces;
      CommandLineParser::split (pieces, *matIter, ":");
      if (pieces.size() != 4)
      {
         dout << "Don't understand '" << *matIter 
              << "' in 'morphAddTemplates': ";
         dumpSTLendl (stringVector("morphAddTemplates"));
         exit (0);
      }
      STemplateMorphInfoMapIter tmiIter = m_morphMap.find( pieces.at(0) );
      if (m_morphMap.end() == tmiIter)
      {
         dout << "Don't understand '" << *matIter 
              << "' in 'morphAddTemplates': ";
         dumpSTLendl (stringVector("morphAddTemplates"));
         exit (0);
      }
      tmiIter->second.m_additionalTemplatesVecMap[ pieces.at(2) ].
         push_back( make_pair (pieces.at(1), atof( pieces.at(3).c_str() ) ) );
   } // for matIter

   _printMorphMap();   

   for (SSetConstIter utmIter = m_usesTemplateMorphing.begin();
        m_usesTemplateMorphing.end() != utmIter;
        ++utmIter)
   {
      CLPTemplateMorph::PtrVec morphPtrVec;
      for (SVecConstIter groupIter = 
              stringVector ("groupNames").begin();
           stringVector ("groupNames").end() != groupIter;
           ++groupIter)
      {
         const string fullName = *utmIter + *groupIter;
         STemplateMorphInfoMapConstIter mmIter = m_morphMap.find (fullName);
         if (m_morphMap.end() == mmIter)
         {
            // this shouldn't happen
            assert (0);
         }
         const TemplateMorphInfo &morphInfo = mmIter->second;
         CLPTemplateMorph *morphPtr =
            new CLPTemplateMorph( fullName.c_str(), fullName.c_str(),
                                  m_groupBinsMap[ *groupIter ],
                                  m_groupLowerEdgeMap[ *groupIter ],
                                  m_groupUpperEdgeMap[ *groupIter ] );
         morphPtr->setNoUnderOverflow (true);
         morphPtrVec.push_back( morphPtr );
         // get and adddefault histogram
         TH1F *defaultHistPtr = 
            (TH1F*) filePtr->Get( morphInfo.m_defaultName.c_str() );
         if (! defaultHistPtr )
         {
            dout << "Can not find '" << morphInfo.m_defaultName
                 << "' in template file." << endl;
            exit (0);
         }
         if ( doubleValue ("intLumi") )
         {
            defaultHistPtr->Scale( doubleValue ("intLumi") );
         }
         morphPtr->addDefaultTemplate (defaultHistPtr);
         // now add all appropriate morphing variables
         for (SDMapConstIter dvmIter = morphInfo.m_defaultValuesMap.begin();
              morphInfo.m_defaultValuesMap.end() != dvmIter;
              ++dvmIter)
         {
            morphPtr->addMorphVariable( dvmIter->first, dvmIter->second );
         } // for dvmIter
         // now add all other morph templaltes
         for (SSDPairVecMapConstIter atvIter =
                 morphInfo.m_additionalTemplatesVecMap.begin();
              morphInfo.m_additionalTemplatesVecMap.end() != atvIter;
              ++atvIter)
         {
            const string &variable = atvIter->first;
            const SDPairVec &addInfoVec = atvIter->second;
            for (SDPairVecConstIter aiIter = addInfoVec.begin();
                 addInfoVec.end() != aiIter;
                 ++aiIter)
            {
               const string &histName = aiIter->first;
               double value = aiIter->second;
               TH1F *histPtr = 
                  (TH1F*) filePtr->Get( histName.c_str() );
               if (! histPtr )
               {
                  dout << "Can not find '" << histName
                       << "' in template file." << endl;
                  exit (0);
               }          
               if ( doubleValue ("intLumi") )
               {
                  histPtr->Scale( doubleValue ("intLumi") );
               }
               morphPtr->addMorphTemplate (variable, value, histPtr);
            } // for aiIter
         } // for atvIter
      } // for groupIter
      // now I've got the vector of template morphs, let's add them
      // together into a single morph and add that to the fitter.
      CLPTemplateMorph *morphPtr = 
         new CLPTemplateMorph ( utmIter->c_str(), utmIter->c_str() );
      for (CLPTemplateMorph::PtrVecConstIter mpvIter = morphPtrVec.begin();
           morphPtrVec.end() != mpvIter;
           ++mpvIter)
      {
         morphPtr->addMorph( *mpvIter );
      } // for mpvIter
      m_histPtrMap[ *utmIter ] = morphPtr;
   } // for utmIter

   delete filePtr;
}

void::
mrf::MRFitter::_rescaleDefaults()
{
   double total = 0.;
   double data = 0.;
   double lumiScale = doubleValue ("intLumi");
   if (! lumiScale )
   {
      lumiScale = 1.;
   }
   for (SVecConstIter nameIter = stringVector ("names").begin();
        stringVector ("names").end() != nameIter;
        ++nameIter)
   {
      const string &name = *nameIter;
      if (ns_dataName == name)
      {
         // Don't add data as a template
         data = m_histPtrMap[name]->Integral();
         continue;
      }
      // get  default values
      double defaultValue = 1.;
      SDMapConstIter defaultIter = m_defaultValuesMap.find( name );
      if (m_defaultValuesMap.end() != defaultIter)
      {
         defaultValue = defaultIter->second;
      }
      total += m_histPtrMap[name]->Integral() * defaultValue * lumiScale;
   }
   double factor = data / total;
   for (SVecConstIter nameIter = stringVector ("names").begin();
        stringVector ("names").end() != nameIter;
        ++nameIter)
   {
      const string &name = *nameIter;
      if (ns_dataName == name)
      {
         // Don't add data as a template
         continue;
      }
      // get and reset default values
      double defaultValue = 1.;
      SDMapConstIter defaultIter = m_defaultValuesMap.find( name );
      if (m_defaultValuesMap.end() != defaultIter)
      {
         defaultValue = defaultIter->second;
      }
      m_defaultValuesMap[name] = defaultValue * factor;
   }
}

void
mrf::MRFitter::_printMorphMap()
{
   for (STemplateMorphInfoMapConstIter iter = m_morphMap.begin();
        m_morphMap.end() != iter;
        ++iter)
   {
      cout << iter->first << " : " << iter->second.m_defaultName << endl;
      for (SDMapConstIter dvmIter = iter->second.m_defaultValuesMap.begin();
           iter->second.m_defaultValuesMap.end() != dvmIter;
           ++dvmIter)
      {
         cout << "   " << dvmIter->first << " : " << dvmIter->second;
      } // for dvmIter
      cout << endl;
      for (SSDPairVecMapConstIter advIter = 
              iter->second.m_additionalTemplatesVecMap.begin();
           iter->second.m_additionalTemplatesVecMap.end() != advIter;
           ++advIter)
      {
         cout << "     " << advIter->first << " : ";
         for (SDPairVecConstIter sdIter = advIter->second.begin();
              advIter->second.end() != sdIter;
              ++sdIter)
         {
            cout<< " (" << sdIter->first << ", " << sdIter->second << ")";
         } // for sdIter
         cout << endl;
      } // for advIter
   } // for iter
}


void
mrf::MRFitter::_hookupFitter ()
{
   m_fitter.setVerbose (CLPTemplateFitter::kOutputFunctoids &
                      CLPTemplateFitter::kOutputBinNorms);
   // setup range information
   SDMap lowerRangeMap, upperRangeMap, stepRangeMap;
   _2vecsToMap (stringVector ("rangeNames"), 
                doubleVector ("lowerRange"), lowerRangeMap);
   _2vecsToMap (stringVector ("rangeNames"), 
                doubleVector ("upperRange"), upperRangeMap);
   _2vecsToMap (stringVector ("rangeNames"), 
                doubleVector ("stepRange"),  stepRangeMap);
   // Add templates.
   for (SVecConstIter nameIter = stringVector ("names").begin();
        stringVector ("names").end() != nameIter;
        ++nameIter)
   {
      const string &name = *nameIter;
      if (ns_dataName == name)
      {
         // Don't add data as a template
         continue;
      }
      // get constraints and default values
      double defaultValue =   1.;
      double lowerValue   = -10.;
      double upperValue   =  10.;
      double stepSize     =   0.001;
      SDMapConstIter defaultIter = m_defaultValuesMap.find( name );
      if (m_defaultValuesMap.end() != defaultIter)
      {
         defaultValue = defaultIter->second;
         if (defaultValue)
         {
            lowerValue = -10. * fabs (defaultValue);
            upperValue =  10. * fabs (defaultValue);
         }
         // if defaultValue is 0, leave the +- 10 as the upper and
         // lower values
      }
      // check RangeMaps
      SDMapConstIter rangeIter = lowerRangeMap.find( name );
      if (lowerRangeMap.end() != rangeIter)
      {
         lowerValue = lowerRangeMap [name];
         upperValue = upperRangeMap [name];
         stepSize   = stepRangeMap  [name];
      } // found it
      // Null for width means no constraint
      double constraintMean  = 0.;
      double constraintWidth = 0.;
      SDMapConstIter constraintIter = m_constraintMap.find( name );
      if (m_constraintMap.end() != constraintIter)
      {
         // we've got a constraint.  Since there's no sense in making
         // a constraint mean different than the default value, just
         // use the default value for the constraint.
         constraintMean  = defaultValue;
         constraintWidth = constraintIter->second;
      } // if constraint found
      // Figure out how many templates I have
      STH1FPtrPairVec histVec;
      for (SVecConstIter endingIter = 
              stringVector ("endings").begin();
           stringVector ("endings").end() != endingIter; 
           ++endingIter)
      {
         string newHistName = name  + *endingIter;
         STH1FPtrMapConstIter histPtrIter =
            m_histPtrMap.find (newHistName);
         if (m_histPtrMap.end() != histPtrIter)
         {
            histVec.push_back( make_pair( name, histPtrIter->second ) );
         }
      } // for endingIter
      // Better have something
      int size = (int) histVec.size();
      assert (size);
      // let's do simple first
      if (1 == size)
      {
         const STH1FPtrPair &histPair = histVec.at(0);
         // Only one.  Just add it.
         m_nameTemplateIndexMap[name] =
            m_fitter.addTemplate (histPair.first,
                                histPair.second,
                                defaultValue,
                                constraintMean, 
                                constraintWidth,
                                lowerValue,
                                upperValue,
                                stepSize );             // step size
         continue;
      } // size == 1
      // if we're here, it's because we have to add more than one
      // template.  Let's first make sure that we are self-consistent.
      assert( m_histPtrMap.end() == m_histPtrMap.find (name) );
      // o.k.  Let's add these puppies
      BinNormClass::ISet templateSet;
      for (STH1FPtrPairVecConstIter histIter = histVec.begin();
           histVec.end() != histIter;
           ++histIter)
      {
         int templateIndex = m_fitter.addTemplate (histIter->first,
                                                 histIter->second,
                                                 1.0, // to be fixed at 1
                                                 1.0, 
                                                 0,
                                                 0,
                                                 2,
                                                 0.001 ); 
         m_nameTemplateIndexMap[name] = templateIndex;
         templateSet.insert (templateIndex);
         // tell the fitter to fix this parameter
         m_fitter.fixParameter (histIter->first);
      } // for histIter
      // we've added and fixed all of the templates.  We now need to
      // add a bin norm to float and tie them all together.
      BinNormClass tempNorm (templateSet, 
                             m_allBinsSet,
                             defaultValue,
                             lowerValue,
                             upperValue,
                             stepSize );             // step size
      m_fitter.addBinNorm (name, tempNorm, constraintMean, constraintWidth);
   } // for histIndex
   // Now add data
   m_fitter.addData( m_histPtrMap[ns_dataName] );
   ///////////////////////////////////////
   // Update PolyNoid Template Indicies //
   ///////////////////////////////////////
   // first tell the Polynoid class about the information it needs.
   for (SIMapConstIter ntiIter = m_nameTemplateIndexMap.begin();
        m_nameTemplateIndexMap.end() != ntiIter;
        ++ntiIter)
   {
      // make sure our other map has this one, too
      SIMapConstIter otherIter = m_nameTemplateIDMap.find( ntiIter->first );
      if (m_nameTemplateIDMap.end() == otherIter)
      {
         // we don't use this one so skip it
         continue;
      }
      assert( m_nameTemplateIDMap.end() != otherIter );
      m_polyTemplate2indexMap[ otherIter->second ] = ntiIter->second;
      // PolyNoid::setTemplate2index (otherIter->second,
      //                              ntiIter->second);
   } // for ntiIter
   // loop over all PolyNoids and have them update this info and add
   // them to the fitter.
   _addPolyMapToFitter ();
   // Need to fix any parameters?
   // Fix parameters with new values
   for (SVecConstIter fixIter = stringVector ("fixParamVal").begin();
        stringVector ("fixParamVal").end() != fixIter;
        ++fixIter)
   {
      SVec words;
      optutl::CommandLineParser::split (words, *fixIter, "=");
      if (2 != words.size())
      {
         cerr << "Can not fix parameter value using '"
              << *fixIter << "'.  Aborting." << endl;
         exit(0);
      }
      const string &paramName  =       words.at(0);
      double        paramValue = atof( words.at(1).c_str() );
      assert( m_fitter.setParameter (paramName, paramValue) );
      assert( m_fitter.fixParameter (paramName) );
   } // for fixIter
   // Fix the parameters with values as they are
   for (SVecConstIter fixIter = 
           stringVector ("fixParameters").begin();
        stringVector ("fixParameters").end() != fixIter;
        ++fixIter)
   {
      assert( m_fitter.fixParameter (*fixIter) );
   } // for fixIter
   // output fitter setup?
   if (boolValue ("outputFitterSetup"))
   {
      cout << "fitter" << endl << m_fitter << endl;
   }
   // Other fitter options
   m_fitter.setMinuitVerboseLevel (integerValue ("minuitVerbose"));
   m_fitter.setMinExpectedBinContent (doubleValue ("minPredictionBin"));
   m_fitter._initializeMinuit();   
   m_nameFitterIndexMap = m_fitter.getNameIndexMap();
}

void
mrf::MRFitter::_createBinsSums ()
{
   using namespace optutl;
   assert( stringVector("binsSum").size() % 3 == 0 );
   int size = (int) stringVector("binsSum").size() / 3;
   for (int loop = 0; loop < size; ++loop)
   {
      string name   = stringVector("binsSum").at( 3 * loop     );
      string temps  = stringVector("binsSum").at( 3 * loop + 1 );
      string groups = stringVector("binsSum").at( 3 * loop + 2 );
      BinsSum &binsSum = m_fitter.binsSumCont().createNew (name);
      SVec tempWords;
      CommandLineParser::split (tempWords, temps, ";");
      for (SVecConstIter tempIter = tempWords.begin();
           tempWords.end() != tempIter;
           ++tempIter)
      {
         int tempIndex = m_fitter.parameterIndex( *tempIter );
         assert (tempIndex >= 0);
         binsSum.addTemplateIndex (tempIndex);
      } // for tempIter
      SVec groupWords;
      CommandLineParser::split (groupWords, groups, ";");
      for (SVecConstIter groupIter = groupWords.begin();
           groupWords.end() != groupIter;
           ++groupIter)
      {
         if ("all" == *groupIter)
         {
            binsSum.addBinRange (1, m_totalBins);
            continue;
         }
         // eventually I should do something here.  In the mean time...
         assert (0);
      } // for groupIter
   } // for loop
}

void
mrf::MRFitter::_decideWhichVariablesToSave()
{
   m_toFillNames.clear();
   m_toFillIndicies.clear();
   m_numToFill = 0;
   int size = m_fitter.size();
   for (int loop = 0; loop < size; ++loop)
   {
      if (! m_fitter.isParameterFixed (loop) )
      {
         // we've got a live one, here
         m_toFillNames.push_back( m_fitter.paramName (loop) );
         m_toFillIndicies.push_back( loop );
         ++m_numToFill;
      }
   } // for loop
}


void 
mrf::_makeEverythingPretty()
{
   // Can be used to set styles....
   gROOT->SetStyle("Plain");  
   gStyle->SetOptStat(0000);
}


int
mrf::_initializeMRF()
{
   _setupOptions();
   return 0;
}


void
mrf::_2vecsToMap (const SVec &namesVec, const DVec &valuesVec,
                  SDMap &valueMap)
{
   assert( namesVec.size() == valuesVec.size() );
   int size = (int) namesVec.size();
   for (int loop = 0; loop < size; ++loop)
   {
      valueMap[ namesVec.at (loop) ] = valuesVec.at (loop);
   }
}


void
mrf::_2vecsToMap (const SVec &namesVec, const IVec &valuesVec,
                  SIMap &valueMap)
{
   assert( namesVec.size() == valuesVec.size() );
   int size = (int) namesVec.size();
   for (int loop = 0; loop < size; ++loop)
   {
      valueMap[ namesVec.at (loop) ] = valuesVec.at (loop);
   }
}

void
mrf::_2vecsToMap (const SVec &namesVec, const SVec &valuesVec,
                  SSMap &valueMap)
{
   assert( namesVec.size() == valuesVec.size() );
   int size = (int) namesVec.size();
   for (int loop = 0; loop < size; ++loop)
   {
      valueMap[ namesVec.at (loop) ] = valuesVec.at (loop);
   }
}

void 
mrf::_svecToSDMap (const SVec &nameAndValueVec, SDMap &valueMap)
{
   assert( ! (nameAndValueVec.size() % 2) );
   int size = ( (int) nameAndValueVec.size() ) / 2;
   for (int loop = 0; loop < size; ++loop)
   {
      const string &name = nameAndValueVec.at( 2 * loop );
      double value = atof( nameAndValueVec.at( 2 * loop + 1).c_str() );
      valueMap [name] = value;
   }
}


void
mrf::MRFitter::_fixGroupStrings()
{
   if (m_fixedGroupStrings) return;
   m_fixedGroupStrings = true;
   SVec newStrings;
   int count = 0;
   if (! stringVector ("groupStrings").size() )
   {
      stringVector ("groupStrings") = stringVector("groupNames");
   }
   for (SVecConstIter iter = 
           stringVector ("groupStrings").begin();
        stringVector ("groupStrings").end() != iter;
        ++iter)
   {
      string line = *iter;
      string::size_type where = line.find ("_");
      while (string::npos != where)
      {
         line = line.replace (where, 1, " ");
         where = line.find ("_");
         assert (++count < 1000);
      } // while
      newStrings.push_back (line);
   } // for iter
   stringVector ("groupStrings") = newStrings;
}

void
mrf::MRFitter::_addPolyToMap (const PolyNoid &poly)
{
   if (! m_canAddToPolyMap)
   {
      cerr << "MRFitter::_addPolyToMap() error:"<< endl;
      assert (0);
   }
   if (m_polynoidMap.end() != m_polynoidMap.find( poly.name() ))
   {
      // already an entry there.  No good
      cerr << "MRFitter::_addPolyToMap() error:"
           << poly.name() << "' already exists.  Aborting." << endl;
      assert (0);
   }
   m_polynoidMap[ poly.name() ] = poly; 
}

PolyNoid &
mrf::MRFitter::_poly (const string &name)
{
   PolyNoid::MapIter iter = m_polynoidMap.find( name );
   if (m_polynoidMap.end() == iter )
   {
      // already an entry there.  No good
      cerr << "MRFitter::_poly() error:'" 
           << name << "' doesn't exists.  Aborting." << endl;
      assert (0);
   }
   return iter->second;
}

void
mrf::MRFitter::_addPolyMapToFitter ()
{
   // O.k.  So after we do this, we can't add any more PolyNoids to
   // the map.
   m_canAddToPolyMap = false;
   for (PolyNoid::MapIter iter = m_polynoidMap.begin();
        m_polynoidMap.end() != iter;
        ++iter)
   {
      //cout << "loading " << iter->first << endl;
      iter->second.calculateGroupTemplateParamMap (m_polyTemplate2indexMap);
      iter->second.addSelfToFitter (m_fitter);
   } // for iter

}




void
mrf::_addSqrtErrors (TH1F *histPtr)
{
   int numBins = histPtr->GetNbinsX();
   for (int bin = 1; bin <= numBins; ++bin)
   {
      double contents = histPtr->GetBinContent (bin);
      if (contents > 0)
      {
         histPtr->SetBinError (bin, sqrt (contents));
      }
   } // for bin
}
