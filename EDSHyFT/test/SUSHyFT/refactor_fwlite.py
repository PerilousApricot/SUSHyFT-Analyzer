#!/usr/bin/env python

# shyft (py)fwlite step
# Andrew Melo <andrew.m.melo@vanderbilt.edu>

# Make sure ROOT is imported first from our wrapper
from SHyFT.ROOTWrap import ROOT

# All other imports
import copy
from DataFormats.FWLite import Events, Handle
import glob
import itertools
import math
from optparse import OptionParser
import sys
import time

class FWLiteAnalysis:
    def __init__(self, outputFile, maxJets = 9, maxBJets = 9, maxTaus = 9, elePtMin = 35, muPtMin = 35, isData = False):
        self.simpleHists = {}
        self.binnedHists = {}
        self.maxJets = maxJets
        self.maxBJets = maxBJets
        self.maxTaus = maxTaus
        self.isData = isData
        self.elePtMin = elePtMin
        self.muPtMin = muPtMin
        self.jetPtMin = 35
        self.metMin = 20
        self.bTagCut = 0.678
        self.tTagCut = 0.5
        self.ignoreMet = False
        # stupid global variables
        self.output = ROOT.TFile(outputFile + ".root", "recreate")
        print "Openting outputFile %s" % outputFile
        self.output.cd()
        self.binnedHists = []
        if self.isData:
            self.flavorNames = [""]
        else:
            self.flavorNames = ["_b", "_c", "_q"]
        for ijet in range(0,self.maxJets+1):
            self.binnedHists.append([])
            for ibjet in range(0,self.maxBJets+1):
                if ibjet > ijet:
                    continue
                self.binnedHists[ijet].append([])
                for itau in range(0,self.maxTaus+1):
                    if (ibjet + itau) > ijet:
                        continue
                    self.binnedHists[ijet][ibjet].append([])
                    for iflavor,_ in enumerate(self.flavorNames):
                        self.binnedHists[ijet][ibjet][itau].append({})
        self.makeHistograms()

    # http://eli.thegreenplace.net/2009/06/12/safely-using-destructors-in-python
    def close(self):
        if self.output:
            self.output.cd()
            self.output.Write()
            self.output.Close()
            self.output = None

    def addSimpleHist(self, name, title, xSteps, minSteps, maxSteps):
        self.simpleHists[name] = ROOT.TH1F(name, title, xSteps,\
                                            minSteps, maxSteps)
    def fillSimpleHist(self, name, value, weight = 1.0):
        self.simpleHists[name].Fill(value, weight)

    def addBinnedHist(self, name, title, xSteps, minSteps, maxSteps):
        for k,v in iterateDeep(self.binnedHists):
            suffix = self.getSuffix(k)
            self.getBin(k)[name] = ROOT.TH1F(name + suffix, title, xSteps,\
                                            minSteps, maxSteps)
    def getBin(self, binInfo):
        # TODO: Memoize this
        target = self.binnedHists
        for row in binInfo:
            target = target[row]
        return target

    def fillBinnedHist(self, name, binInfo, value, weight = 1.0):
        self.getBin(binInfo)[name].Fill(value, weight)

    def getSuffix(self, binInfo):
        suffix = "_%ij_%ib_%it%s" % (binInfo[0],
                                    binInfo[1],
                                    binInfo[2],
                                    self.flavorNames[binInfo[3]])
        return suffix

    
    def makeHistograms(self):
        p = [
             ["nEvents", "Number of Events;N_{events };Number", 5, 0.5, 4.5],
             ["nEventsNoPU", "Number of Events;N_{events };Number", 5, 0.5, 4.5],
            ]
        for plot in p:
            self.addSimpleHist(*plot)

        names = ['nPV', 'nPVnoPU','lepEta','lepPt', 'centrality','sumEt', 'MET', 'wMT', 'hT', 'std', 'stdt','wrp','nAmbiguous','nJets']
        titles = ['number of Primary Vertices','number of Primary Vertices no PileUp','Lepton #eta', 'Lepton pt', 'Centrality','#sum E_{T}','MET', 'M_{WT}', 'hT', "SingleTopDisc", "SingleTopDiscT",'wJetR', 'Ambiguous Jets', 'nJets']
        bounds = [  [25, -0.5, 24.5],
                    [25,-0.5,24.5],
                    [30,0.,3.0],
                    [100,0.,200],
                    [120,0.,1.2],
                    [100,0.,1000.],
                    [120,0.,300.],
                    [120,0.,300.],
                    [120,0.,1200.],
                    [240,0.,600.],
                    [240,0.,600.],
                    [165,0.,3.3],
                    [5,0.5,4.5],
                    [5,0.5,4.5],
                ]
        for plot in itertools.izip(names,titles,bounds):
            self.addBinnedHist(plot[0],plot[1],*plot[2])

    def getPUWeight(self, event):
        if not self.isData:
            #num of true interations
            nPV = event.getByTitle('puInfo')
            if nPV == None:
                raise RuntimeError, "Couldn't pull PU info from sample"
            genVertices = -1
            for iPV in nPV:
                if iPV.getBunchCrossing() == 0:
                    genVertices = iPV.getPU_NumInteractions()
                    break
            if genVertices == -1:
                raise RuntimeError, "Couldn't pull PU info from sample"
            dataWeight = event.puDataHist.GetBinContent(
                                        event.puDataHist.GetBin(genVertices + 1))
            mcWeight = event.puMCHist.GetBinContent(
                                        event.puMCHist.GetBin(genVertices + 1))
            if mcWeight * event.puDataIntegral == 0:
                PUWeight = 1
            else:
                # I wonder if there's a float-only division operator
                PUWeight = float( float(dataWeight * event.puMCIntegral) /
                                float(mcWeight * event.puDataIntegral) )
        else:
            PUWeight = 1
        return PUWeight

    def handleLepton(self, event, lepStr, minPt):
        lepCount = 0
        lepPts = event.getByTitle(lepStr + 'Pt')
        if not lepPts:
            return 0
        for lepPt in lepPts:
            if lepPt > minPt:
                lepCount += 1
        return lepCount
    
    def leptonCut(self, event):
         # Always require exactly one muon. future should be different leps.
        #eleCount = self.handleLepton(event, "ele", 35)
        eleCount = 0
        muonCount = self.handleLepton(event, "mu", 35)
        return muonCount != 1 or eleCount != 0

    def metCut(self, event):
        if self.ignoreMet:
            return False
        rawMet = event.getByTitle('metPt')
        if not rawMet:
            return True
        return rawMet[0] <= self.metMin

    def shouldCut(self, event):
        # TODO: Calibrate what is the ideal cut-order from the first few events
        # Cut on lep content
        if self.leptonCut(event):
            return True
        # Cut on MET content
        if self.metCut(event):
            return True
        return False

    def countJets(self, event):
        """ returns (njets, nbtags, nttags, ambiguous tags, bIndices) 
        """
        jetPts = event.getByTitle('jetPts')
        jetBDiscs = event.getByTitle('bDisc')
        jetTDiscs = event.getByTitle('tDisc')

        if not jetPts:
            return 0,0,0,0,[]

        nJets = 0
        nBTags = 0
        nTTags = 0
        nAmbiguous = 0
        bTagIndices = []
        for index, jetPt in enumerate(jetPts):
            if jetPt < self.jetPtMin:
                continue
            nJets += 1
            # Here if something is both a B-jet and a Tau, I call it B. Okay?
            if jetBDiscs[index] > self.bTagCut:
                bTagIndices.append(index)
                nBTags += 1
            elif jetTDiscs[index] > self.tTagCut:
                nTTags += 1

            # I'm not convinced this is such a good idea.. count up the number
            # of ambiguous and non-ambiguous jets
            if jetBDiscs[index] > self.bTagCut and jetTDiscs[index] > self.tTagCut:
                nAmbiguous += 1
        assert (nTTags + nBTags) <= nJets, "Don't double-count events"
        return nJets, nBTags, nTTags, nAmbiguous, bTagIndices

    def processEvent(self, event):
        PUWeight = self.getPUWeight(event)
        self.fillSimpleHist("nEvents", 1, PUWeight)
        self.fillSimpleHist("nEventsNoPU", 1)
        if self.shouldCut(event):
            return
        nJets, nBTags, nTTags, nAmbiguous, bTagIndices = self.countJets(event)
        nJets = min(nJets, self.maxJets)
        nBTags = min(nBTags, self.maxBJets)
        bTagIndices = bTagIndices[0:nBTags]
        nTTags = min(nTTags, self.maxTaus)

        # Get kinematics
        metRaw = event.getByTitle("metPt")[0]
        metPhi = event.getByTitle("metPhi")[0]
        jetPts = event.getByTitle("jetPts")
        jetEtas = event.getByTitle("jetEtas")
        jetPhis = event.getByTitle("jetPhis")
        jetMasses = event.getByTitle("jetMasses")
        if not self.isData:
            jetFlavors = event.getByTitle("jetFlavors")
        lepPt = event.getByTitle("lepPts")[0]
        lepEta = event.getByTitle("lepEtas")[0]
        lepPhi = event.getByTitle("lepPhis")[0]
        lepMass = 0.10565837
        allJets = []
        leadingBJetPt = 0
        leadingBJet = None
        leadingLightJetPt = 0
        leadingLightJet = None
        if self.isData:
            eventFlavor = 0
        else:
            eventFlavor = 2
        hT = lepPt + metRaw
        lep_px = lepPt * math.cos( lepPhi )
        lep_py = lepPt * math.sin( lepPhi )
        wPt = lepPt + metRaw
        met_px = metRaw * math.cos( metPhi )
        met_py = metRaw * math.sin( metPhi )
        wPx = lep_px + met_px
        wPy = lep_py + met_py
        wMT = math.sqrt(wPt*wPt-wPx*wPx-wPy*wPy)
        sumEt = 0.
        sumPt = 0.
        sumE = 0.
        deltaR = 5.0
        for ijet, jetPt in enumerate(jetPts):
            ## get the jet 4-vector
            thisJet = ROOT.TLorentzVector()
            thisJet.SetPtEtaPhiM(jetPt,
                                 jetEtas[ijet],
                                 jetPhis[ijet],
                                 jetMasses[ijet])
            allJets.append(thisJet)
            deta = thisJet.Eta() - lepEta
            dphi = thisJet.Phi() - lepPhi
            dphi = normalizedPhi(dphi)
            deltaR = min( math.sqrt(deta*deta + dphi*dphi), deltaR)
            hT    += thisJet.Et()    
            sumEt += thisJet.Et()
            sumPt += thisJet.Pt()
            sumE  += thisJet.E()

            if ijet in bTagIndices:
                if jetPt > leadingBJetPt:
                    leadingBJet = thisJet
                    leadingBJetPt = jetPt
            else:
                if jetPt > leadingLightJetPt:
                    leadingLightJet = thisJet
                    leadingLightJetPt = jetPt
            if not self.isData:
                jetFlavor = jetFlavors[ijet]
                if abs(jetFlavor) == 5:
                    eventFlavor = 0
                elif abs(jetFlavor) == 4:
                    eventFlavor = min(eventFlavor, 1)
                else:
                    eventFlavor = min(eventFlavor, 2)

        lepVector = ROOT.TLorentzVector()
        lepVector.SetPtEtaPhiM(lepPt,
                                lepEta,
                                lepPhi,
                                lepMass)
        metVector = ROOT.TLorentzVector()
        metVector.SetPtEtaPhiM(metRaw,
                                0,
                                metPhi,
                                0)

        vertices = event.getByTitle("vertices")
        if not vertices:
            vertices = []
        if leadingBJet:
            singleTopDiscriminatorT = (leadingBJet + lepVector + metVector).Mt()
            singleTopDiscriminator  = (leadingBJet + lepVector + metVector).M()
            wRLeadingVector = leadingBJet + lepVector + metVector
            if leadingLightJet:
                wRTrailingVector = leadingLightJet
                # TMath::Abs(normalizedPhi(v1.phi()-v2.phi()))
                deltaPhi = wRLeadingVector.Phi() - wRTrailingVector.Phi()
                wRDiscriminator = abs(normalizedPhi(deltaPhi))
        binInfo = (nJets, nBTags, nTTags, eventFlavor)
        self.fillBinnedHist("nJets", binInfo, nJets)
        if nAmbiguous:
            self.fillBinnedHist("nAmbiguous", binInfo, nAmbiguous)
        self.fillBinnedHist("nPV", binInfo, vertices.size(), PUWeight )    
        self.fillBinnedHist("lepEta" , binInfo, lepEta, PUWeight )
        self.fillBinnedHist("lepPt" , binInfo, lepPt, PUWeight )
        if sumE != 0.0:
            self.fillBinnedHist("centrality" , binInfo, sumEt / sumE, PUWeight )
        self.fillBinnedHist("sumEt" , binInfo, sumEt, PUWeight )
        self.fillBinnedHist("MET" , binInfo, metRaw, PUWeight )
        self.fillBinnedHist("wMT" , binInfo, wMT, PUWeight )
        self.fillBinnedHist("hT" , binInfo, hT, PUWeight )
        if leadingBJet:
            self.fillBinnedHist("std" , binInfo, singleTopDiscriminator, PUWeight )
            self.fillBinnedHist("stdt" , binInfo, singleTopDiscriminatorT, PUWeight )
            if leadingLightJet:
                self.fillBinnedHist("wrp" , binInfo, wRDiscriminator, PUWeight )


class SHyFTEventInfo:
    def __init__(self):
        pass

class ROOTEventInfo(SHyFTEventInfo):
    """ thin wrapper around ROOT events """
    def __init__(self, lepStr, postfix, puMCInput, puDataInput):
        self.jetPtHandle         = Handle( "std::vector<float>" )
        self.jetPtLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "pt" )
        self.genJetPtHandle          = Handle( "std::vector<float>" )
        self.genJetPtLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "genJetPt" )
        self.jetEtaHandle         = Handle( "std::vector<float>" )
        self.jetEtaLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "eta" )
        self.jetPhiHandle         = Handle( "std::vector<float>" )
        self.jetPhiLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "phi" )
        self.jetMassHandle         = Handle( "std::vector<float>" )
        self.jetMassLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "mass" )
        self.jetSecvtxMassHandle         = Handle( "std::vector<float>" )
        self.jetSecvtxMassLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "secvtxMass" )
        self.jetSSVHEHandle         = Handle( "std::vector<float>" )
        self.jetSSVHELabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "csvhe" )
        self.jetFlavorHandle         = Handle( "std::vector<float>" )
        self.jetFlavorLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "flavor" )
        self.jetTauPtHandle         = Handle( "std::vector<float>" )
        self.jetTauPtLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "tauJetPt" )
        self.jetTauEtaHandle         = Handle( "std::vector<float>" )
        self.jetTauEtaLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "tauJetEta" )
        self.jetTauPhiHandle         = Handle( "std::vector<float>" )
        self.jetTauPhiLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "tauJetPhi" )
        self.jetTauMassHandle         = Handle( "std::vector<float>" )
        self.jetTauMassLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "tauJetMass" )
        self.tauTagHandle      =Handle( "std::vector<float>" )
        self.tauDB3HitLooseHandle         = Handle( "std::vector<float>" )
        self.tauDB3HitLooseLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "byLooseCombinedIsolationDeltaBetaCorr3Hits" )
        self.tauDB3HitMediumHandle         = Handle( "std::vector<float>" )
        self.tauDB3HitMediumLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "byMediumCombinedIsolationDeltaBetaCorr3Hits" )
        self.tauDB3HitTightHandle         = Handle( "std::vector<float>" )
        self.tauDB3HitTightLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "byTightCombinedIsolationDeltaBetaCorr3Hits" )
        self.tauMVALooseHandle         = Handle( "std::vector<float>" )
        self.tauMVALooseLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "byLooseIsolationMVA2" )
        self.tauMVAMediumHandle         = Handle( "std::vector<float>" )
        self.tauMVAMediumLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "byMediumIsolationMVA2" )
        self.tauMVATightHandle         = Handle( "std::vector<float>" )
        self.tauMVATightLabel    = ( "pfShyftTupleJets" + lepStr +  postfix,   "byTightIsolationMVA2" )
        self.muonPtHandle         = Handle( "std::vector<float>" )
        self.muonPtLabel    = ( "pfShyftTupleMuons"+  postfix,   "pt" )
        self.muonEtaHandle         = Handle( "std::vector<float>" )
        self.muonEtaLabel    = ( "pfShyftTupleMuons"+  postfix,   "eta" )
        self.muonPhiHandle         = Handle( "std::vector<float>" )
        self.muonPhiLabel    = ( "pfShyftTupleMuons"+  postfix,   "phi" )
        self.muonPfisoHandle         = Handle( "std::vector<float>" )
        self.muonPfisoLabel    = ( "pfShyftTupleMuons"+  postfix,   "pfiso" )
        self.electronPtHandle         = Handle( "std::vector<float>" )
        self.electronPtLabel    = ( "pfShyftTupleElectrons"+  postfix,   "pt" )
        self.electronEtaHandle         = Handle( "std::vector<float>" )
        self.electronEtaLabel    = ( "pfShyftTupleElectrons"+  postfix,   "eta" )
        self.electronPhiHandle         = Handle( "std::vector<float>" )
        self.electronPhiLabel    = ( "pfShyftTupleElectrons"+  postfix,   "phi" )
        self.electronPfisoHandle         = Handle( "std::vector<float>" )
        self.electronPfisoLabel    = ( "pfShyftTupleElectrons"+  postfix,   "pfiso" )

        self.metHandle = Handle( "std::vector<float>" )
        self.metLabel = ("pfShyftTupleMET" + lepStr +  postfix,   "pt" )
        self.metPhiHandle = Handle( "std::vector<float>" )
        self.metPhiLabel = ("pfShyftTupleMET" + lepStr +  postfix,   "phi" )

        ### FIXME SHOULD BE SMARTER ###
        self.puMCFile = ROOT.TFile()                                                              
        self.puMCFile = self.puMCFile.Open(puMCInput)
        self.puDataFile = ROOT.TFile()                                                         
        self.puDataFile = self.puDataFile.Open(puDataInput)
        self.puMCHist = self.puMCFile.Get("analyzeHiMassTau/NVertices_0")                             
        self.puDataHist = self.puDataFile.Get("analyzeHiMassTau/NVertices_0") 
        self.puMCIntegral = self.puMCHist.Integral()
        self.puDataIntegral = self.puDataHist.Integral()
        self.puInfoLabel = ("addPileupInfo")
        self.puInfoHandle = Handle("std::vector<PileupSummaryInfo>")
        self.vertH = Handle("std::vector<reco::Vertex>")
        self.vertLabel = ("offlinePrimaryVertices")
        
        self.varLookup = {}
        self.varLookup['muPt'] = (self.muonPtLabel, self.muonPtHandle)
        self.varLookup['elePt'] = (self.electronPtLabel, self.electronPtHandle)
        self.varLookup['puInfo'] = (self.puInfoLabel, self.puInfoHandle)
        self.varLookup['metPt'] = (self.metLabel, self.metHandle)
        self.varLookup['metPhi'] = (self.metPhiLabel, self.metPhiHandle)
        self.varLookup['bDisc'] = (self.jetSSVHELabel, self.jetSSVHEHandle)
        self.varLookup['tDisc'] = (self.tauDB3HitLooseLabel , self.tauDB3HitLooseHandle)
        self.varLookup['jetPts'] = (self.jetPtLabel, self.jetPtHandle)
        self.varLookup['jetEtas'] = (self.jetEtaLabel, self.jetEtaHandle)
        self.varLookup['jetPhis'] = (self.jetPhiLabel, self.jetPhiHandle)
        self.varLookup['jetMasses'] = (self.jetMassLabel, self.jetMassHandle)
        self.varLookup['jetFlavors'] = (self.jetFlavorLabel, self.jetFlavorHandle)
        self.varLookup['puInfo'] = (self.puInfoLabel, self.puInfoHandle) 
        self.varLookup['vertices'] = (self.vertLabel, self.vertH)
        if lepStr == 'Mu':
            self.varLookup['lepPts'] = (self.muonPtLabel, self.muonPtHandle)
            self.varLookup['lepEtas'] = (self.muonEtaLabel, self.muonEtaHandle)
            self.varLookup['lepPhis'] = (self.muonPhiLabel, self.muonPhiHandle)
        elif lepStr == 'Ele':
            self.varLookup['lepPts'] = (self.electronPtLabel, self.electronPtHandle)
            self.varLookup['lepEtas'] = (self.electronEtaLabel, self.electronEtaHandle)
            self.varLookup['lepMasses'] = (self.electronMassLabel, self.electronMassHandle)

    def getByTitle(self, title):
        if title in self.eventCache:
            return self.eventCache[title]
        print "Title %s Label %s Handle %s .. eventID %s -- event %s" % (title, self.varLookup[title][0], id(self.varLookup[title][1]), id(self.event), self.event)
        self.event.getByLabel(*self.varLookup[title])
        if self.varLookup[title][1].isValid():
            product = self.varLookup[title][1].product()
        else:
            product = None
        self.eventCache[title] = product
        return product

    def resetHandles(self):
        self.eventCache = {}
# helper function from c++
def normalizedPhi(phi):
    M_PI = 3.14159265358979323846
    TWO_PI = M_PI * 2
    while ( phi < -M_PI ):
        phi += TWO_PI
    while ( phi >  M_PI ):
        phi -= TWO_PI
    return phi;

def iterateDeep(t, *keys):                                                      
    if isinstance(t, type([])): 
        for k,v in enumerate(t):                                                
            for yieldKey, val in iterateDeep(v,k,*keys):
                yieldKey.insert(0,k)
                yield yieldKey, val                                           
    else:                   
        yield [], t

splitArgs = [list(x[1]) for x in itertools.groupby(sys.argv[1:], lambda x: x=='++') if not x[0]]
def getParser(default=None):
    parser = OptionParser()
    # Input files to use. This is in "glob" format, so you can use wildcards.
    # If you get a "cannot find file" type of error, be sure to use "\*" instead
    # of "*" to make sure you don't confuse the shell. 
    parser.add_option('--files', metavar='F', type='string', action='store',
                      dest='files',
                      help='Input files glob')
    
    parser.add_option('--inputListFile', metavar='F', type='string',
                        action='store', dest='inputListFile',
                        help='Filelist to process')

    # Output name to use. 
    parser.add_option('--outname', metavar='F', type='string', action='store',
                      default='shyft_fwlite' if default else None,
                      dest='outname',
                      help='output name')

    # no MET cut
    parser.add_option('--noMET', metavar='F', action='store_true',
                      default=False if default else None,
                      dest='noMET',
                      help='no MET cut')

    # Using data or not. For MC, truth information is accessed.
    parser.add_option('--useData', metavar='F', action='store_true',
                      default=False if default else None,
                      dest='useData',
                      help='use data')

    # BTag systematics
    parser.add_option('--btagSys', metavar='F', type='string', action='store',
                      default='nominal' if default else None,
                      dest='btagSys',
                      help='BTag Systematic variation. Options are "nominal, up, down, up2, down2"')

    # LFTag systematics
    parser.add_option('--lftagSys', metavar='F', type='string', action='store',
                      default='nominal' if default else None,
                      dest='lftagSys',
                      help='LFTag Systematic variation. Options are "nominal, up, down, up2, down2"')

    #Min Jets
    parser.add_option('--minJets', metavar='F', type='int', action='store',
                      default=1 if default else None,
                      dest='minJets',
                      help='Min number of jets')
    #Jet Pt Cut
    parser.add_option('--jetPt', metavar='F', type='float', action='store',
                      default=35.0 if default else None,
                      dest='jetPt',
                      help='Jet Pt cut')

    # electrons or muons
    parser.add_option('--lepType', metavar='F', type='int', action='store',
                      default=0 if default else None,
                      dest='lepType',
                      help='Lepton type. Options are 0 = muons, 1 = electrons')

    # PU Data distribution
    parser.add_option('--PUDataInput', metavar='P', type='string', action='store',
                      default='DataPUFile_Full2012.root' if default else None,
                      dest='puDataInput',
                      help='Input data PU distribution')

    # PU MC distribution
    parser.add_option('--PUMCInput', metavar='P', type='string', action='store',
                      default='S10MC_PUFile.root' if default else None,
                      dest='puMCInput',
                      help='Input MC PU distribution')

    return parser

analyzerList = []
eventWrapperList = []
# Show the help screen if the user doesn't give any arguments
if len(splitArgs) == 0:
    splitArgs = [["--help"]]

if len(splitArgs) == 1:
    parser = getParser(default = True)
    (options, args) = parser.parse_args(splitArgs[0])
    argList = [[options, args]]
    analyzerList.append(FWLiteAnalysis(outputFile=options.outname,
                                       isData=options.useData))
    eventWrapperList.append(ROOTEventInfo(lepStr='Mu',
                                            postfix='',
                                            puMCInput=options.puMCInput,
                                            puDataInput=options.puDataInput))
else:
    argList = []
    parser = getParser(default = True)
    (mainOptions, mainArgs) = parser.parse_args(splitArgs[0])
    outNames = []
    for oneList in splitArgs[1:]:
        parser = getParser(default = False)
        (subOptions, subArgs) = parser.parse_args(oneList)
        argCopy = mainArgs[:]
        argCopy.extend(subArgs)
        optCopy = copy.copy(mainOptions)
        if getattr(subOptions,'files', None) or getattr(subOptions, 'inputFileList', None):
            raise RuntimeError, "Each subrun must use the same input files"
        # subOpts have priority
        for key in subOptions.__dict__:
            if subOptions.__dict__[key] != None:
                setattr(optCopy, key, subOptions.__dict__[key])
        argList.append([optCopy, argCopy])
        if subOptions.outname in outNames:
            raise RuntimeError, "Two iterations are writing to same output"
        outNames.append(subOptions.outname)
        analyzerList.append(FWLiteAnalysis(outputFile=optCopy.outname, isData=optCopy.useData))
        eventWrapperList.append(ROOTEventInfo(lepStr='mu',
                                            postfix='',
                                            puMCInput=optCopy.puMCInput,
                                            puDataInput=optCopy.puDataInput))

#ROOT.gSystem.Load('libCondFormatsJetMETObjects')

#from CondFormats.JetMETObjects import *
firstOptions = argList[0][0]
if firstOptions.files and firstOptions.inputListFile:
    raise RuntimeError, "You can either specify the files on the command line " +\
                        "OR pass a file list. Not both."

if firstOptions.files:
    print 'Getting files from this glob: ' + firstOptions.files

    # Get the file list. 
    files = glob.glob( firstOptions.files )
elif firstOptions.inputListFile:
    files = []
    fh = open(firstOptions.inputListFile, 'r')
    for line in fh.readlines():
        line = line.rstrip()
        if line:
            files.append(line)
else:
    raise RuntimeError, "No files specified"
print files
events = Events(files)
print "Getting event count"
ntotal = events.size()
print " ... got %s events total" % ntotal
startTime = time.time()
stopTime = time.time()
currentEvent = 0
percentDoneLast = -1
lastTime = time.time()
lastEvent = 0
for event in events:
    percentDone = int(200.0 * currentEvent/ntotal)
    if percentDone != percentDoneLast:
        percentDoneLast = percentDone
        currentTime = time.time()
        if (currentTime - lastTime) != 0:
            eventsPerSecond = (currentEvent - lastEvent) / (currentTime - lastTime)
        else:
            eventsPerSecond = 0
        lastTime = currentTime
        lastEvent = currentEvent
        print '  Processing {0:9.0f}/{1:.0f}: {2:5.1f}% {3:7.0f} ev/sec'.format( 
                    currentEvent, ntotal, percentDone/2.0, eventsPerSecond)
    currentEvent += 1
    for eventWrapper in eventWrapperList:
        eventWrapper.resetHandles()
        eventWrapper.event = None
    for analyze, eventWrapper in zip(analyzerList, eventWrapperList):
        eventWrapper.event = event
        analyze.processEvent(eventWrapper)
    stopTime = time.time()
for analyze in analyzerList:
    analyze.close()
elaspedTime = stopTime - startTime
print "Processed %s events in %.2f seconds %7.0f ev/sec" % (currentEvent, 
                                                            elaspedTime,
                                                            currentEvent/elaspedTime)
