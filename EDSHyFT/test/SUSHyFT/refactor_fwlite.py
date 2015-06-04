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

from Analysis.SHyFTScripts.combinations import EffInfo, EffInfoCombinations
# load JEC uncertainty stuff
jecParStr = ROOT.std.string('Jec12_V3_Uncertainty_AK5PFchs.txt')
jecUnc = ROOT.JetCorrectionUncertainty( jecParStr )


class FWLiteAnalysis:
    def __init__(self, outputFile, 
                        maxJets = 9,
                        maxBJets = 9, 
                        maxTaus = 9,
                        elePtMin = 35,
                        muPtMin = 35,
                        isData = False,
                        noMET = False,
                        metMin = 20,
                        flipTauOrder = False,
                        jetPt = 35,
                        invertTauCut = False,
                        jecScale = 0.0,
                        sfB = 1.0,
                        sfC = 1.0,
                        sfQ = 1.0):
        self.simpleHists = {}
        self.binnedHists = {}
        self.maxJets = maxJets
        self.maxBJets = maxBJets
        self.maxTaus = maxTaus
        self.isData = isData
        self.elePtMin = elePtMin
        self.muPtMin = muPtMin
        self.jetPtMin = int(jetPt)
        self.metMin = int(metMin)
        self.bTagCut = 0.678
        self.tTagCut = 0.5
        self.ignoreMet = noMET
        self.invertTauCut = invertTauCut
        self.flipTauOrder = flipTauOrder
        # stupid global variables
        self.output = ROOT.TFile(outputFile + ".root", "recreate")
        print "Openting outputFile %s" % outputFile
        self.output.cd()
        self.binnedHists = []
        self.totalEvents = 0
        self.eventsPassed = 0
        self.sfB = sfB
        self.sfC = sfC
        self.sfQ = sfQ
        self.jecScale = jecScale
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
        try:
            for row in binInfo:
                target = target[row]
            return target
        except:
            print "Couldn't update bin %s target %s" % (binInfo, target)
            raise

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

        names = ['nPV', 
                 'nPVnoPU',
                 'lepEta',
                 'lepPt',
                 'centrality',
                 'sumEt',
                 'MET',
                 'wMT',
                 'hT',
                 'std',
                 'stdt',
                 'wrp',
                 'nAmbiguous',
                 'nJets',
                 'dPhi2JetMET',
                 'mTLeadMuoMET',
                 '2MuonInvMass',
                 'gt2MuonInvMass',
                 'nMuon',
                 'nElectron']
        titles = ['number of Primary Vertices',
                  'number of Primary Vertices no PileUp',
                  'Lepton #eta',
                  'Lepton pt',
                  'Centrality',
                  '#sum E_{T}',
                  'MET',
                  'M_{WT}',
                  'hT',
                  "SingleTopDisc",
                  "SingleTopDiscT",
                  'wJetR',
                  'Ambiguous Jets',
                  'nJets',
                  'dPhi between 2nd Jet and MET',
                  'mT_{lead muon, MET}',
                  'invariant mass, 2 muons',
                  'invariant mass, #gt muons'
                  'nMuons #gt 10GeV',
                  'nElectrons #gt 10GeV']
        bounds = [  [40, -0.5, 39.5], #npv
                    [40,-0.5,39.5], # npbnopu
                    [30,0.,3.0], #lepeta
                    [100,0.,200], #leppt
                    [120,0.,1.2], # centrality
                    [100,0.,1000.], # sumEt
                    [300,0.,300.], # MET
                    [120,0.,300.], # wmt
                    [120,0.,1200.], #ht
                    [240,0.,600.], #std
                    [240,0.,600.], #stdt
                    [69,-0.05,3.4], #wrp
                    [10,-0.5,10.5], #ambiguous
                    [10,-0.5,10.5], # njet
                    [136,-3.4,3.4], #dphi jet met
                    [200,0,400], # mt
                    [100,0,200], # inv mass
                    [150,0,300], # gt2 muon mass
                    [10,-0.5,10.5],
                    [10,-0.5,10.5],
                ]
        for plot in itertools.izip(names,titles,bounds):
            self.addBinnedHist(plot[0],plot[1],*plot[2])

    def getPUWeight(self, event):
        if not self.isData:
            #num of true interations
            nPV = event.getByTitle('vertices')
            if nPV <= 0:
                raise RuntimeError, "Couldn't pull PU info from sample: %s" % nPV
            nPV = nPV[0]
            dataWeight = event.puDataHist.GetBinContent(
                                        event.puDataHist.GetBin(nPV + 1))
            mcWeight = event.puMCHist.GetBinContent(
                                        event.puMCHist.GetBin(nPV + 1))
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
        eleCount = self.handleLepton(event, "ele", 35)
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
        tauPts = event.getByTitle('tauPts')
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
            if nJets >= self.maxJets:
                break
            nJets += 1
            # Here if something is both a B-jet and a Tau, I call it B. Okay?
            passB = jetBDiscs[index] > self.bTagCut
            if self.invertTauCut:
                passTau = jetTDiscs[index] < self.tTagCut and tauPts[index] > 20
            else:
                passTau = jetTDiscs[index] > self.tTagCut and tauPts[index] > 20

            if passB:
                bTagIndices.append(index)
                nBTags += 1
            elif passTau:
                nTTags += 1

            # I'm not convinced this is such a good idea.. count up the number
            # of ambiguous and non-ambiguous jets
            if passB and passTau:
                nAmbiguous += 1
        assert (nTTags + nBTags) <= nJets, "Don't double-count events"
        return nJets, nBTags, nTTags, nAmbiguous, bTagIndices

    def processEvent(self, event):
        PUWeight = self.getPUWeight(event)
        self.fillSimpleHist("nEvents", 1, PUWeight)
        self.fillSimpleHist("nEventsNoPU", 1, 1.0)
        if self.metCut(event):
            return

        nJets, nBTags, nTTags, nAmbiguous, bTagIndices = self.countJets(event)
        # get necessary kinematics
        if self.isData:
            eventFlavor = 0
        else:
            eventFlavor = 2
            jetFlavors = event.getByTitle('jetFlavors')
        nJets = min(nJets, self.maxJets)
        nBTags = min(nBTags, self.maxBJets)
        bTagIndices = bTagIndices[0:nBTags]
        nTTags = min(nTTags, self.maxTaus)
        binInfo = (nJets, nBTags, nTTags, eventFlavor)
        eleCount2Mu = self.handleLepton(event, "ele", 10)
        muonCount2Mu = self.handleLepton(event, "mu", 10)
        lepPts = event.getByTitle("lepPts")
        lepEtas = event.getByTitle("lepEtas")
        lepPhis = event.getByTitle("lepPhis")
        if muonCount2Mu >= 1:
            lepMass = 0.10565837
            lepVector = ROOT.TLorentzVector()
            lepVector.SetPtEtaPhiM(lepPts[0],
                                    lepEtas[0],
                                    lepPhis[0],
                                    lepMass)

        if (eleCount2Mu == 0 and muonCount2Mu == 2):
            lepVector2 = ROOT.TLorentzVector()
            lepVector2.SetPtEtaPhiM(lepPts[1],
                                    lepEtas[1],
                                    lepPhis[1],
                                    lepMass)

            self.fillBinnedHist("2MuonInvMass", binInfo, (lepVector + lepVector2).M(), PUWeight )    

        if (eleCount2Mu == 0 and muonCount2Mu > 2):
            f = math.factorial
            n = muonCount2Mu
            r = 2
            nCr = f(n) / f(r) / f(n-r)
            for x,y in itertools.combinations(range(muonCount2Mu), 2):
                lepVector1 = ROOT.TLorentzVector()
                lepVector1.SetPtEtaPhiM(lepPts[x],
                                        lepEtas[x],
                                        lepPhis[x],
                                        lepMass)
                lepVector2 = ROOT.TLorentzVector()
                lepVector2.SetPtEtaPhiM(lepPts[y],
                                        lepEtas[y],
                                        lepPhis[y],
                                        lepMass)
                self.fillBinnedHist("gt2MuonInvMass", binInfo,
                                        (lepVector1 + lepVector2).M(),
                                        float(PUWeight)/ float(nCr) )

        # Cut the main kinematic plots
        if self.leptonCut(event):
            return
        # Get kinematics
        metRaw = event.getByTitle("metPt")[0]
        metPhi = event.getByTitle("metPhi")[0]
        jetPts = event.getByTitle("jetPts")
        jetEtas = event.getByTitle("jetEtas")
        jetPhis = event.getByTitle("jetPhis")
        jetMasses = event.getByTitle("jetMasses")
        lepMass = 0.10565837

        allJets = []
        leadingBJetPt = 0
        leadingBJet = None
        leadingLightJetPt = 0
        leadingLightJet = None
        lep_px = lepPts[0] * math.cos( lepPhis[0] )
        lep_py = lepPts[0] * math.sin( lepPhis[0] )
        met_px = metRaw * math.cos( metPhi )
        met_py = metRaw * math.sin( metPhi )
        # dPhi(jet[1], MET)
        if nJets > 1:
            dPhi2JetMET = normalizedPhi(jetPhis[1] - metPhi)
        sumEt = 0.
        sumPt = 0.
        sumE = 0.
        deltaR = 5.0
        isf = 1.0
        hT = 0.0
        effs = []
        for ijet, jetPt in enumerate(jetPts):
            ## get the jet 4-vector
            thisJet = ROOT.TLorentzVector()
            thisJet.SetPtEtaPhiM(jetPt,
                                 jetEtas[ijet],
                                 jetPhis[ijet],
                                 jetMasses[ijet])
            # JEC Systematic
            if not self.isData and self.jecScale != 0.0:
                # Subtract the jet from the met, correct the jet then replace
                # it
                jecUnc.setJetEta( jetEtas[ijet] )
                jecUnc.setJetPt( jetPts[ijet] )
                upOrDown = bool(jecScale > 0.0)
                unc = abs(jecUnc.getUncertainty(upOrDown))
                jetScale = 1 + unc * jecScale
                met_px += thisJet.Px()
                met_py += thisJet.Py()
                thisJet *= jetScale
                met_px -= thisJet.Px()
                met_py -= thisJet.Py()

            # Track btag/lftag systematic
            if not self.isData:
                jetFlavor = jetFlavors[ijet]
                if abs(jetFlavor) == 5:
                    eventFlavor = 0
                    iflavor = 0
                    isf = self.sfB
                elif abs(jetFlavor) == 4:
                    eventFlavor = min(eventFlavor, 1)
                    isf = self.sfC
                    iflavor = 1
                else:
                    eventFlavor = min(eventFlavor, 2)
                    isf = self.sfQ
                    iflavor = 2
                if ijet in bTagIndices:
                    effs.append(EffInfo(ijet, isf, iflavor))
                else:
                    effs.append(EffInfo(ijet, 0.0, iflavor))

            allJets.append(thisJet)
            deta = thisJet.Eta() - lepEtas[0]
            dphi = thisJet.Phi() - lepPhis[0]
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

        metRaw2 = math.sqrt(met_py*met_py + met_px*met_px)
        metPhi2 = math.atan2(met_py,met_px)
        metVector = ROOT.TLorentzVector()
        metVector.SetPtEtaPhiM(metRaw2,
                                0,
                                metPhi2,
                                0)
        wPx = lep_px + met_px
        wPy = lep_py + met_py
        wPt = lepPts[0] + metRaw2
        hT += lepPts[0] + metRaw2
        wMT = math.sqrt(wPt*wPt-wPx*wPx-wPy*wPy)

        if not self.isData:
            effCombos = EffInfoCombinations(effs, verbose=False)
            fillList = []
            for itag in range(0,nJets+1):
                pTag = effCombos.pTag( itag )
                jtag = itag
                if jtag > self.maxBJets:
                    jtag = self.maxBJets
                # use pileup reweighting factor
                pTag*=PUWeight
                newBinInfo = list(binInfo[:])
                newBinInfo[1] = jtag
                # nBjet + nTjet > nJet
                if jtag + binInfo[2] > binInfo[0]:
                    newBinInfo[2] = 0
                fillList.append([newBinInfo, pTag])
        else:
            fillList = [[binInfo, 1]]
        vertices = event.getByTitle("vertices")[0]
        for binInfo, PUWeight in fillList:
            # do the combinatorials
            if leadingBJet:
                singleTopDiscriminatorT = (leadingBJet + lepVector + metVector).Mt()
                singleTopDiscriminator  = (leadingBJet + lepVector + metVector).M()
                wRLeadingVector = leadingBJet + lepVector + metVector
                if leadingLightJet:
                    wRTrailingVector = leadingLightJet
                    # TMath::Abs(normalizedPhi(v1.phi()-v2.phi()))
                    deltaPhi = wRLeadingVector.Phi() - wRTrailingVector.Phi()
                    wRDiscriminator = abs(normalizedPhi(deltaPhi))
            self.fillBinnedHist("nJets", binInfo, nJets)
            if nAmbiguous:
                self.fillBinnedHist("nAmbiguous", binInfo, nAmbiguous)
            self.fillBinnedHist("nPV", binInfo, vertices, PUWeight )    
            self.fillBinnedHist("nPVnoPU", binInfo, vertices, 1 )    
            self.fillBinnedHist("lepEta" , binInfo, lepEtas[0], PUWeight )
            self.fillBinnedHist("lepPt" , binInfo, lepPts[0], PUWeight )
            if sumE != 0.0:
                self.fillBinnedHist("centrality" , binInfo, sumEt / sumE, PUWeight )
            if nJets > 1:
                self.fillBinnedHist("dPhi2JetMET" , binInfo, dPhi2JetMET, PUWeight )
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
        self.vertH = Handle("int")
        self.vertLabel = ("pfShyftProducer" + lepStr +  postfix ,"genpv")
        
        self.varLookup = {}
        self.varLookup['tauPts'] = (self.jetTauPtLabel, self.jetTauPtHandle)
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

    # invertTauCut
    parser.add_option('--invertTauCut', metavar='F', action='store_true',
                      default=False if default else None,
                      dest='invertTauCut',
                      help='Reverses tau discriminator cut')

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

    # JEC systematics
    parser.add_option('--jecSys', metavar='F', type='string', action='store',
                      default='nominal' if default else None,
                      dest='jecSys',
                      help='JEC Systematic variation. Options are "nominal, up, down"')

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

    parser.add_option('--metMin', metavar='P', type='string', action='store',
                      default=20 if default else None,
                      dest='metMin',
                      help='Met cut (GeV)')

    parser.add_option('--flipTauOrder', metavar='P', type='string', action='store',
                      default=False if default else None,
                      dest='flipTauOrder',
                      help='Flip the order of testing for taus and jets')



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
    # Per-jet scale factors:
    sfB = 1.00
    sfC = 1.00
    sfQ = 1.00
    if options.btagSys == 'up':
        sfB = 1.10
        sfC = 1.10
    elif options.btagSys == 'down':
        sfB = 0.90
        sfC = 0.90
    elif options.btagSys == 'up2':
        sfB = 1.20
        sfC = 1.20
    elif options.btagSys == 'down2':
        sfB = 0.80
        sfC = 0.80
        
    if options.lftagSys == 'up':
        sfQ = 1.10
    elif options.lftagSys == 'down':
        sfQ = 0.90
    elif options.lftagSys == 'up2':
        sfQ = 1.20
    elif options.lftagSys == 'down2':
        sfQ = 0.80

    # JEC scales
    jecScale = 0.0
    if options.jecSys == 'up':
        jecScale = 1.0
    elif options.jecSys == 'down':
        jecScale = -1.0


    analyzerList.append(FWLiteAnalysis(outputFile=options.outname,
                                       isData=options.useData,
                                       noMET=options.noMET,
                                       metMin=options.metMin,
                                       flipTauOrder=options.flipTauOrder,
                                       jetPt=options.jetPt,
                                       invertTauCut=options.invertTauCut,
                                       sfB = sfB,
                                       sfC = sfC,
                                       sfQ = sfQ,
                                       jecScale = jecScale))
    eventWrapperList.append(ROOTEventInfo(lepStr='Mu',
                                            postfix='',
                                            puMCInput=options.puMCInput,
                                            puDataInput=options.puDataInput))
else:
    raise Exception, "This doesn't work"
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
        analyzerList.append(FWLiteAnalysis(outputFile=optCopy.outname,
                                            isData=optCopy.useData,
                                            noMET=optCopy.noMET,
                                            metMin=optCopy.metMin,
                                            flipTauOrder=optCopy.flipTauOrder,
                                            jetPt=options.jetPt,
                                            invertTauCut=optCopy.invertTauCut))
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
if events.size():
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
