import FWCore.ParameterSet.Config as cms
import sys
process = cms.Process("ANA")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")


###############################
####### Parameters ############
###############################
print "Parsing these command line arguments: %s" % sys.argv
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')


options.register('usePDFs',
                 0,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Ignore trigger in selection")

options.register('useData',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Use data (1) or MC (0)")

options.register('useLooseElectrons',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Add loose electron collection (1) or not (0)")

options.register('useLooseMuons',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Add loose muon collection (1) or not (0)")

options.register('useMETRes',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Add shifted MET collection (1) or not (0)")

options.register('eleTriggerName',
                 'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "trigger to run for ele paths")

options.register('muTriggerName',
                 'HLT_IsoMu24_eta2p1_v,HLT_Mu30_eta2p1_v,HLT_Mu40_eta2p1_v',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "trigger to run for muon paths")

options.register('runMuons',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Run Muon paths")

options.register('runElectrons',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Run Electron Paths")

options.register('ignoreTrigger',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Ignore trigger bits (1) or not (0)")

options.register('minJets',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Minimum number of jets per event")


options.parseArguments()

print options
 ## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Add MVA triggers
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )

# hopefully pulls in the proper pat trigger
#process.load( "PhysicsTools.PatAlgos.patSequences_cff" )
#from PhysicsTools.PatAlgos.patEventContent_cff import *

inputCutsToIgnore = ["== 1 Tight Lepton","== 1 Lepton","Dilepton Veto","== 1 Tight Lepton, Mu Veto"]
if options.useData == 0 :
    useData = False
    #inputCutsToIgnore = ['Trigger']
    process.GlobalTag.globaltag = 'START53_V7F::All'
    payloads = [
        'Jec12_V3_MC_L1FastJet_AK5PFchs.txt',
        'Jec12_V3_MC_L2Relative_AK5PFchs.txt',
        'Jec12_V3_MC_L3Absolute_AK5PFchs.txt',
        'Jec12_V3_MC_L2L3Residual_AK5PFchs.txt',
        'Jec12_V3_MC_Uncertainty_AK5PFchs.txt',
        ]
else:
    useData = True
    process.GlobalTag.globaltag = 'GR_P_V40_AN1::All'
    payloads = [
        'Jec12_V3_L1FastJet_AK5PFchs.txt',
        'Jec12_V3_L2Relative_AK5PFchs.txt',
        'Jec12_V3_L3Absolute_AK5PFchs.txt',
        'Jec12_V3_L2L3Residual_AK5PFchs.txt',
        'Jec12_V3_Uncertainty_AK5PFchs.txt',
        ]

# triggering's being done in the HLTfilter
#inputCutsToIgnore = ['Trigger']
# run the trigger on the fly
#process.load('PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff')

process.load("Configuration.StandardSequences.MagneticField_cff")

if useData:
    inputFiles = [
            '/store/user/meloam/SingleMu/meloam_feb12_tlbsm53x2_Run2012C_24Aug2012_v1/20130212222033/00000/4463A7D0-9075-E211-95CC-003048F2E8C2.root'
        ]
else :
    inputFiles = [
            '/store/user/gurrola/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/DY4JetsToLL_M-50_SkimPat/405381fe00d9112adafce059de4ce799/skimPat_1000_1_GJN.root',
#'file:///scratch/meloam/SUSHyFT-Scripts/checkouts/pat/src/HighMassAnalysis/Configuration/test/Data_TauTauSkim/skimPat.root']
 #   nothingFiles=[
'/store/user/meloam/ra2tau_pat_8TeV/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014400/0000/skimPat_1.root',
'/store/user/gurrola/DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/DY1JetsToLL_M-50_SkimPat/405381fe00d9112adafce059de4ce799/skimPat_752_1_I91.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_101.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_102.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_103.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_104.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_105.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_106.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_107.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_108.root',
'/store/user/meloam/ra2tau_pat_8TeV/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/ra2tau_skim2015_v5_Summer12_DR53X-PU_S10_START53_V7A-v1/150721_014310/0000/skimPat_109.root',
            #'file:///scratch/meloam/SUSHyFT-Scripts/checkouts/pat/src/HighMassAnalysis/Configuration/test/Data_TauTauSkim/skimPat.root',
            #'/store/user/flanagwh/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/5_3_11_patch6_Vandy/405381fe00d9112adafce059de4ce799/skimPat_100_2_05a.root'
            #'/store/user/gurrola/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/W2JetsToLNu_SkimPat/405381fe00d9112adafce059de4ce799/skimPat_1001_1_nnH.root'
            #'root://se8.accre.vanderbilt.edu:1095//store/user/flanagwh/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/5_3_11_patch6_Vandy/405381fe00d9112adafce059de4ce799/skimPat_688_3_e0l.root',
            #'root://xrootd.unl.edu//store/user/lpctlbsm/meloam/GJets_HT-200To400_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_425_1_CX2.root'
        ]

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring( inputFiles )
)

if len(options.inputFiles) > 0 :
    process.source.fileNames = options.inputFiles


## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

from Analysis.SHyFT.shyftselection_cfi import wplusjetsAnalysis as shyftSelectionInput

if options.usePDFs:
	process.load('Analysis.PdfWeights.pdfWeightProducer_cfi')

process.pfShyftProducerMu = cms.EDFilter('EDSHyFTSelector',
        shyftSelection = shyftSelectionInput.clone(
            muonSrc = cms.InputTag('patMuons'),
            electronSrc = cms.InputTag('patElectrons'),
            metSrc = cms.InputTag('patPfMetT0pcT1Txy'),
            jetSrc = cms.InputTag('selectedPatJets'),
            pvSrc   = cms.InputTag('offlinePrimaryVertices'),
            tauSrc = cms.InputTag('patTaus'),
            ePlusJets = cms.bool( False ),
            muPlusJets = cms.bool( True ),
            matchByHand = cms.bool(False),
            muTrig = cms.string(options.muTriggerName),
            eleTrig = cms.string(options.eleTriggerName),
            jetPtMin = cms.double(10.0),##30
            minJets = cms.int32(options.minJets),
            metMin = cms.double(0.0),
            tauPtMin = cms.double(10.0),
            tauEtaMax = cms.double(3.0),
            muPtMin = cms.double(25.0),##35
            identifier = cms.string('PFMu'),
            cutsToIgnore=cms.vstring( inputCutsToIgnore ),
            useData = cms.bool(useData),
            jetSmear = cms.double(0.0),
            jecPayloads = cms.vstring( payloads ),
            puDataFileName = cms.string("DataPUFile_Full2012.root"),
            puMCFileName = cms.string("S10MC_PUFile.root"),
            ),
            matchByHand = cms.bool(False),
        )
#process.pfShyftProducerMu.userData.userFunctionLabels = cms.vstring('secvtxMass')
process.pfShyftTupleJetsMu = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("pfShyftProducerMu", "jets"),
    pvSrc   = cms.InputTag('offlinePrimaryVertices'),
    lazyParser = cms.untracked.bool(True),
    eventInfo = cms.untracked.bool(False),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("mass"),
            quantity = cms.untracked.string("mass")
            ),
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
            ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
            ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
            ),
        cms.PSet(
            tag = cms.untracked.string("csvhe"),
            quantity = cms.untracked.string("bDiscriminator('combinedSecondaryVertexBJetTags')")
            ),
        cms.PSet(
            tag = cms.untracked.string("secvtxMass"),
            quantity = cms.untracked.string("userFloat('secvtxMass')"),
            ),
        cms.PSet(
            tag = cms.untracked.string("jetArea"),
            quantity = cms.untracked.string("jetArea")
            ),
        cms.PSet(
            tag = cms.untracked.string("tauJetPt"),
            quantity = cms.untracked.string("userFloat('tauJetPt')")
            ),
        cms.PSet(
            tag = cms.untracked.string("tauJetPhi"),
            quantity = cms.untracked.string("userFloat('tauJetPhi')")
            ),
        cms.PSet(
            tag = cms.untracked.string("tauJetEta"),
            quantity = cms.untracked.string("userFloat('tauJetEta')")
            ),
        cms.PSet(
            tag = cms.untracked.string("tauJetMass"),
            quantity = cms.untracked.string("userFloat('tauJetMass')")
            ),
        cms.PSet(
            tag = cms.untracked.string("byLooseCombinedIsolationDeltaBetaCorr3Hits"),
            quantity = cms.untracked.string("userFloat('byLooseCombinedIsolationDeltaBetaCorr3Hits')")
            ),
        cms.PSet(
            tag = cms.untracked.string("byMediumCombinedIsolationDeltaBetaCorr3Hits"),
            quantity = cms.untracked.string("userFloat('byMediumCombinedIsolationDeltaBetaCorr3Hits')")
            ),
        cms.PSet(
            tag = cms.untracked.string("byTightCombinedIsolationDeltaBetaCorr3Hits"),
            quantity = cms.untracked.string("userFloat('byTightCombinedIsolationDeltaBetaCorr3Hits')")
            ),
        cms.PSet(
            tag = cms.untracked.string("byLooseIsolationMVA2"),
            quantity = cms.untracked.string("userFloat('byLooseIsolationMVA2')")
            ),
        cms.PSet(
            tag = cms.untracked.string("byMediumIsolationMVA2"),
            quantity = cms.untracked.string("userFloat('byMediumIsolationMVA2')")
            ),
        cms.PSet(
            tag = cms.untracked.string("byTightIsolationMVA2"),
            quantity = cms.untracked.string("userFloat('byTightIsolationMVA2')")
            ),
        )
    )

if not options.useData:
    process.pfShyftTupleJetsMu.variables.append(
        cms.PSet(
            tag = cms.untracked.string("flavor"),
            quantity = cms.untracked.string("partonFlavour()")
            )
        )

process.pfShyftTupleMuons = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("pfShyftProducerMu", "muons"),
    lazyParser = cms.untracked.bool(True),
    eventInfo = cms.untracked.bool(False),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("tight"),
            quantity = cms.untracked.string("userInt('isTight')")
        ),
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
            ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
            ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
            ),
        cms.PSet(
            tag = cms.untracked.string("pfiso"),
            quantity = cms.untracked.string("userIsolation('pat::PfChargedHadronIso') + " +
                                            "userIsolation('pat::PfNeutralHadronIso') + " +
                                            "userIsolation('pat::PfGammaIso')"
                                            )
            ),
        )
    )

process.pfShyftTupleElectrons = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("pfShyftProducerMu", "electrons"),
    lazyParser = cms.untracked.bool(True),
    eventInfo = cms.untracked.bool(False),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
            ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
            ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
            ),
        cms.PSet(
            tag = cms.untracked.string("pfiso"),
            quantity = cms.untracked.string("userIsolation('pat::PfChargedHadronIso') + " +
                                            "userIsolation('pat::PfNeutralHadronIso') + " +
                                            "userIsolation('pat::PfGammaIso')"
                                            )
            ),
        )
    )


process.pfShyftTupleMETMu = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("pfShyftProducerMu", "MET"),
    lazyParser = cms.untracked.bool(True),
    eventInfo = cms.untracked.bool(False),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
            ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
            )
        )
    )


process.p1 = cms.Path()
if options.runMuons:
    process.p1 = cms.Path(
        process.pfShyftProducerMu*
        process.pfShyftTupleJetsMu*
        process.pfShyftTupleMuons*
        process.pfShyftTupleElectrons*
        process.pfShyftTupleMETMu
        #process.pfShyftTupleTausMu
        )

process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("shyft_edntuple.root"),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p1') ),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*_genpv_*',
                                                                      'keep *_*_pileUp_*',
                                                                      'keep *_pfShyftTuple*_*_*',
                                                                      'keep *_pdfWeightProducer_*_*',
                                                                      'keep *_patTriggerEvent_*_*',
                                                                      'keep edmTriggerResults_TriggerResults__HLT',
                                                                      'keep patTriggerPaths_patTrigger_*_*'
                                                                      )
                               )
process.outpath = cms.EndPath(process.out)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.out.dropMetaData = cms.untracked.string("DROPPED")

if options.useData:
    #suppress the L1 trigger error messages
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
    process.MessageLogger.suppressWarning.append('patTrigger')
    process.MessageLogger.cerr.FwkJob.limit=1
    process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit =
                                                           cms.untracked.int32(0) )

