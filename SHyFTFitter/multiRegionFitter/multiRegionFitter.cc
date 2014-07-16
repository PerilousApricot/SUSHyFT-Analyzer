#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBenchmark.h"
#include "TRandom3.h"

#include "OptionUtils/interface/CommandLineParser.h"
#include "OptionUtils/interface/dout.h"
#include "CLPTemplateMorph/interface/dumpSTL.icc"
#include "MRFitterNamespace.h"

#include <iostream>
#include <string>
#include <cassert>

using namespace std;

// forward declarations
void localInitialize();

int main(int argc, char** argv)
{
   /////////////////////////////////////////
   // Setup fitters and parameter vectors //
   /////////////////////////////////////////
   localInitialize();     // hook up local options
   mrf::ns_parser.parseArguments (argc, argv, true);
   mrf::setDataName();
   const optutl::CommandLineParser::SVec arguments = 
      mrf::ns_parser.nonOptionArgsVec();

   if ( ! arguments.size() or arguments.size() > 2)
   {
      cerr << "Must provide exactly 1 or 2 configuration filenames." << endl;
      mrf::ns_parser.help();
   }
   string fitterConfigName   ( arguments.at(0) );
   string peFitterConfigName ( arguments.at(0) );
   if (arguments.size() > 1)
   {
      peFitterConfigName = arguments.at(1);
   }
   mrf::MRFitter fitter   ( fitterConfigName );
   mrf::MRFitter peFitter ( peFitterConfigName, "PE");

   gRandom->SetSeed (fitter.integerValue ("randomSeed"));
   mrf::DVec paramVec;
   mrf::DVec peParamVec;   
   BinsSumCont &binsSumCont = fitter.binsSumCont();
   // setup parameter vectors to fitter initial values
   fitter.setupParamVec   (paramVec,   true);
   peFitter.setupParamVec (peParamVec, true);

   ///////////////////////
   // Setup output name //
   ///////////////////////
   string tag;
   mrf::processOptions (tag);
   fitter.stringValue ("output") = 
      optutl::CommandLineParser::removeEnding (fitter.
                                               stringValue ("output"), 
                                               ".root");
   fitter.stringValue ("output") += tag;
   if (fitter.integerValue ("randomSeed"))
   {
      fitter.stringValue ("output") += 
         Form ("_seed%03d", fitter.integerValue ("randomSeed"));
   }
   if (! fitter.doubleValue ("intLumi") )
   {
      // do nothing
   } else if (fitter.doubleValue ("intLumi") < 0.0001 ||
       fitter.doubleValue ("intLumi") > 1e4)
   {
      fitter.stringValue ("output") += 
         Form ("_lum%f", fitter.doubleValue ("intLumi"));
   } else if (fitter.doubleValue ("intLumi") >= 0.0001 &&
              fitter.doubleValue ("intLumi") < 1.0)
   {
      fitter.stringValue ("output") += 
         Form ("_lum%.4f", fitter.doubleValue ("intLumi"));
   } else if (fitter.doubleValue ("intLumi") >= 1.0 &&
              fitter.doubleValue ("intLumi") < 1e3)
   {
      fitter.stringValue ("output") += 
         Form ("_lum%.1f", fitter.doubleValue ("intLumi"));
   } else {
      fitter.stringValue ("output") += 
         Form ("_lum%.0f", fitter.doubleValue ("intLumi"));
   }

   //////////////////////////////////////////////
   // Save templates with specified parameters //
   //////////////////////////////////////////////
   if (fitter.stringVector ("saveCanvasWithParams").size())
   {
      mrf::DVec tempParamVec;
      fitter.setupParamVec   (tempParamVec, true);
      mrf::SDMap valueMap;
      mrf::_svecToSDMap (fitter.stringVector ("saveCanvasWithParams"),
                         valueMap);
      for (mrf::SDMapConstIter iter = valueMap.begin();
           valueMap.end() != iter;
           ++iter)
      {
         int index = fitter.parameterIndex( iter->first );
         if (index < 0)
         {
            // uh oh
            cerr << "Fitter does not have a parameter '" << iter->first
                 << "'." << endl;
            assert (0);
         }
         tempParamVec.at (index) = iter->second;
      } // for iter
      fitter.saveCanvasResult (fitter.stringValue ("output") + "_params",
                                tempParamVec);
   }

   //////////////
   // Fit Data //
   //////////////
   if (fitter.boolValue ("fitData"))
   {
      fitter.fit (true);
      if (fitter.boolValue ("savePlots"))
      {
         fitter.saveCanvasResult (fitter.stringValue ("output"));
      }
      if (binsSumCont.size() && fitter.integerValue ("covarPEs"))
      {
         mrf::DVec resultVec;
         fitter.fillParamVec (resultVec);
         mrf::SVec nameVec;
         binsSumCont.getAllNames (nameVec);
         mrf::DVec centerVec;
         binsSumCont.setParamVec (resultVec);
         binsSumCont.getValues (centerVec);
         fitter.storeSqrtMatrix();
         for (int loop = 0; 
              loop < fitter.integerValue ("covarPEs"); 
              ++loop)
         {
            mrf::DVec randomVec (resultVec);
            fitter.generateRandomParams (randomVec);
            binsSumCont.setParamVec (randomVec);
            binsSumCont.recordTrials();
         } // for loop
         mrf::DVec meanVec, rmsVec;
         binsSumCont.fillMeansAndRMSs (meanVec, rmsVec);
         // dumpNamedSTLendl ("meas   ", meanVec);
         // dumpNamedSTLendl ("center ", centerVec);
         // dumpNamedSTLendl ("rms    ", rmsVec);
      } // if we have something to play with
   } // if fitting data

   /////////////////
   // Test Fitter //
   /////////////////
   if (fitter.boolValue ("runFitterTest"))
   {
      fitter.testFitter();
   }

   // only bother going further if we are going to generate PEs
   if (! fitter.integerValue ("numPEs"))
   {
      return 0;
   }

   ////////////////////////////
   // ////////////////////// //
   // // //////////////// // //
   // // // Throw PEs! // // //
   // // //////////////// // //
   // ////////////////////// //
   ////////////////////////////

   //////////////////////
   // //////////////// //
   // // Setup Tree // //
   // //////////////// //
   //////////////////////
   TFile *outputFilePtr = 0;
   TTree *treePtr = 0;
   string outputName = fitter.stringValue ("output") + ".root";
   cout << "Opening '" << outputName<< "' for output." << endl;
   outputFilePtr = TFile::Open (outputName.c_str(), "RECREATE");
   assert (outputFilePtr);
   // write out names of objects that are saved in vectors
   outputFilePtr->WriteObject (& (fitter.toFillNames()), "names");
   outputFilePtr->WriteObject (& (peFitter.toFillNames()), "PEnames");
   mrf::SVec scanNamesVec;
   int scanVecSize = (int) fitter.scanStructVec().size();
   if (scanVecSize)
   {
      for (mrf::ScanStructVecConstIter ssvIter = fitter.scanStructVec().begin();
           fitter.scanStructVec().end() != ssvIter;
           ++ssvIter)
      {
         scanNamesVec.push_back( ssvIter->m_name );
      } // for ssvIter
   }
   outputFilePtr->WriteObject ( &scanNamesVec, "scanNames");
   treePtr = new TTree ("PEtree", "PEtree");

   ////////////////////
   // Define Vectors //
   ////////////////////
   mrf::DVec *genVecPtr     = new mrf::DVec (fitter.numToFill(), 0.);
   mrf::DVec *measVecPtr    = new mrf::DVec (fitter.numToFill(), 0.);
   mrf::DVec *errVecPtr     = new mrf::DVec (fitter.numToFill(), 0.);
   // minos errors
   mrf::DVec *posErrVecPtr  = 0;
   mrf::DVec *negErrVecPtr  = 0;
   // scan errors
   mrf::DVec *scanPosErrVecPtr = 0;
   mrf::DVec *scanNegErrVecPtr = 0;
   // Bins sums numbers
   mrf::DVec *bscMeasVecPtr = 0;
   mrf::DVec *bscCentVecPtr = 0;
   mrf::DVec *bscGenVecPtr  = 0;
   mrf::DVec *bscErrVecPtr  = 0;

   ////////////////////
   // Setup Branches //
   ////////////////////
   treePtr->Branch ("measVec", "vector<double>", &measVecPtr);
   treePtr->Branch ("genVec",  "vector<double>", &genVecPtr);
   treePtr->Branch ("errVec",  "vector<double>", &errVecPtr);
   if (fitter.boolValue ("doMinos"))
   {
      treePtr->Branch ("posErrVec",  "vector<double>", &posErrVecPtr);
      treePtr->Branch ("negErrVec",  "vector<double>", &negErrVecPtr);
      posErrVecPtr  = new mrf::DVec (fitter.numToFill(), 0.);
      negErrVecPtr  = new mrf::DVec (fitter.numToFill(), 0.);
   }
   if (scanVecSize)
   {
      scanPosErrVecPtr = new mrf::DVec (scanVecSize, 0.);
      scanNegErrVecPtr = new mrf::DVec (scanVecSize, 0.);
      treePtr->Branch ("scanPosErr", "vector<double>", &scanPosErrVecPtr);
      treePtr->Branch ("scanNegErr", "vector<double>", &scanNegErrVecPtr);
   }
   if (binsSumCont.size() && fitter.integerValue ("covarPEs"))
   {
      bscMeasVecPtr = new mrf::DVec;
      bscCentVecPtr = new mrf::DVec;
      bscGenVecPtr  = new mrf::DVec;
      bscErrVecPtr  = new mrf::DVec;
      treePtr->Branch ("measBsc", "vector<double>", &bscMeasVecPtr);
      treePtr->Branch ("centBsc", "vector<double>", &bscMeasVecPtr);
      treePtr->Branch ("genBsc",  "vector<double>", &bscGenVecPtr);
      treePtr->Branch ("errBsc",  "vector<double>", &bscErrVecPtr);
      mrf::SVec bscNamesVec;
      binsSumCont.getAllNames (bscNamesVec);
      outputFilePtr->WriteObject (&bscNamesVec, "bscNames");      
   } // if we are using binsSumContainers

   ////////////////////////
   // Final Preparations //
   ////////////////////////
   int everyN = 1000;
   if (fitter.integerValue ("numPEs") % 20 < everyN && 
       fitter.integerValue ("numPEs") > 20)
   {
      everyN = fitter.integerValue ("numPEs") / 20;
      if (! everyN)
      {
         everyN = 1;
      }
   }
   // setup vector, initializing it
   // dumpNamedSTLendl ("names", peFitter.stringVector ("PEnames"));
   mrf::BlurStructVec blurVec;
   peFitter.setupBlurStruct (blurVec, 
                             peFitter.stringVector ("PEnames"), 
                             peFitter.doubleVector ("PEmeans"), 
                             peFitter.doubleVector ("PEsigmas"));
   mrf::SVec peParamNamesVec;
   peFitter.fillNameVec (peParamNamesVec);

   /////////////////////////
   // Throw, Baby, Throw! //
   /////////////////////////
   for (int loop = 0; loop < peFitter.integerValue ("numPEs"); ++loop)
   {
      if ((loop + 1) % everyN == 0 || fitter.boolValue ("printFitResults"))
      {
         cout << "event " << loop + 1 << endl;
      }
      // only put in new values into vector
      mrf::blurParamVec (peParamVec, blurVec);
      if (fitter.boolValue ("printFitResults"))
      {
         assert (peParamVec.size() == peParamNamesVec.size());
         int size = peParamNamesVec.size();
         cout << "Thrown:" << endl;
         for (int loop = 0; loop < size; ++loop)
         {
            cout << Form ("%-10s: %f",
                          peParamNamesVec.at (loop).c_str(),
                          peParamVec.at(loop))
                 << endl;
         }
      } // if printing results
      // throw using PE fitter
      peFitter.throwPE (peParamVec, peFitter.boolValue ("idealPE"));
      // move PE from PE fitter to fitting fitter
      fitter.getDataFromOtherFitter (peFitter);
      // fit
      fitter.fit (fitter.boolValue ("printFitResults"));
      mrf::DVec fitterValues, fitterErrors;
      fitter.fillParamVec (fitterValues);
      fitter.fillErrorVec (fitterErrors);
      // generated values is from peFitter; results from fitter
      peFitter.convertFitterToTreeOrder (peParamVec,   *genVecPtr);
      fitter.convertFitterToTreeOrder   (fitterValues, *measVecPtr);
      fitter.convertFitterToTreeOrder   (fitterErrors, *errVecPtr);      
      if (fitter.boolValue ("doMinos"))
      {
         mrf::DVec posErrors, negErrors;
         fitter.fillPosErrorVec (posErrors);
         fitter.fillNegErrorVec (negErrors);
         fitter.convertFitterToTreeOrder   (posErrors, *posErrVecPtr);      
         fitter.convertFitterToTreeOrder   (negErrors, *negErrVecPtr);      
      }
      if (scanVecSize)
      {
         for (int loop = 0; loop < scanVecSize; ++loop)
         {
            scanPosErrVecPtr->at(loop) = 
               fitter.scanStructVec().at(loop).m_upperError;
            scanNegErrVecPtr->at(loop) = 
               fitter.scanStructVec().at(loop).m_lowerError;
         } // for loop
      } // if scanning variables
      if (fitter.boolValue ("savePlots"))
      {
         fitter.saveCanvasResult( Form ("%s_pe%d", 
                                        fitter.stringValue("output").c_str(), 
                                        loop) );
         fitter.saveCanvasResult( Form ("%s_pe%d_should", 
                                        fitter.stringValue("output").c_str(), 
                                        loop),
                                  paramVec );
      }
      if (binsSumCont.size() && fitter.integerValue ("covarPEs"))
      {
         binsSumCont.reset();
         fitter.storeSqrtMatrix();
         binsSumCont.setParamVec (fitterValues);
         binsSumCont.getValues (*bscCentVecPtr);
         for (int loop = 0; loop < fitter.integerValue ("covarPEs"); 
              ++loop)
         {
            mrf::DVec randomVec (fitterValues);
            fitter.generateRandomParams (randomVec);
            binsSumCont.setParamVec (randomVec);
            binsSumCont.recordTrials();
         }
         binsSumCont.fillMeansAndRMSs (*bscMeasVecPtr, *bscErrVecPtr);
         binsSumCont.setParamVec (paramVec);
         binsSumCont.getValues (*bscGenVecPtr);
         if (fitter.boolValue ("printFitResults"))
         {
            // dumpNamedSTLendl ("meas   ", *bscMeasVecPtr);
            // dumpNamedSTLendl ("center ", *bscCentVecPtr);
            // dumpNamedSTLendl ("rms    ", *bscErrVecPtr);
            // dumpNamedSTLendl ("gen    ", *bscGenVecPtr); 
            cout << endl << endl << endl << endl;
         } // if print results
      } // if running covariance PEs

      treePtr->Fill();
   } // for loop
   // write out everything
   outputFilePtr->cd();
   treePtr->Write();
   outputFilePtr->Close();
   cout << "done." << endl;
}

void localInitialize()
{
   using namespace optutl;
   using namespace mrf;
   ns_parser.addOption ("fitData",
                        VariableMapCont::kBool,
                        "Fit data or not",
                        false);
   ns_parser.addOption ("numPEs",
                        VariableMapCont::kInteger,
                        "Number of PEs to throw",
                        0);
   ns_parser.addOption ("output",
                        VariableMapCont::kString,
                        "Name of output root file",
                        "output.root");
   ns_parser.addOption ("idealPE",
                        VariableMapCont::kBool,
                        "Throw ideal PE",
                        false);
   ns_parser.addOption ("savePlots",
                        VariableMapCont::kBool,
                        "Save Plots",
                        false);
   ns_parser.addOption ("printFitResults",
                        VariableMapCont::kBool,
                        "Print fit results for PEs",
                        false);
   ns_parser.addOption ("runFitterTest",
                        VariableMapCont::kBool,
                        "Run testFitter subroutine for debugging",
                        false);
   ns_parser.addOption ("randomSeed",
                        VariableMapCont::kInteger,
                        "Random seed to use for PEs",
                        0);
   ns_parser.addOption ("print",
                        VariableMapCont::kBool,
                        "Print all values after loading config file",
                        false);
   ns_parser.addOption ("covarPEs",
                        VariableMapCont::kInteger,
                        "Number of covariance sets to generate",
                        0);
   ns_parser.addOption ("saveCanvasWithParams",
                        VariableMapCont::kStringVector,
                        "Name and value pairs used to save canvas");
   ns_parser.addOption ("peNames",
                        VariableMapCont::kStringVector,
                        "Names of variables to fluctuate for PEs");
   ns_parser.addOption ("peMeans",
                        VariableMapCont::kDoubleVector,
                        "Means of variables to fluctuate for PEs");
   ns_parser.addOption ("peSigmas",
                        VariableMapCont::kDoubleVector,
                        "Sigmas of variables to fluctuate for PEs");

}
