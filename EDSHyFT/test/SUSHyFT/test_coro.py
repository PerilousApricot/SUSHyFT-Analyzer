#!/usr/bin/env python2.6
import itertools
paramList = [ { 'foo' : True }, { 'foo' : False } ]

# Import everything from ROOT
import sys
import time
# Fuck root for stealing argv just by importing it
oldargv = sys.argv[:]
sys.argv = []
import ROOT
from ROOT import TMath
ROOT.gROOT.Macro("~/rootlogon.C")
sys.argv = oldargv[:]
from DataFormats.FWLite import Events, Handle


def getEvents():
    counter = 0
    for x in Events('/store/user/meloam/MET/melo_v2_METRun2012CPromptReco_/562ac012a02badb8cb519b40e75790e9/shyftDump2_1_1_gRM.root'):
        print "Getting event %s" % counter
        counter = counter + 1
        yield x

def runOnce(params, eventList):
    # put the entire old script in here
    print "Initializing, foo is %s" % params['foo']
    Mpostfix = ''
    muonPtHandle         = Handle( "std::vector<float>" )
    muonPtLabel    = ( "pfShyftTupleMuons"+  Mpostfix,   "pt" )
    muonEtaHandle         = Handle( "std::vector<float>" )
    muonEtaLabel    = ( "pfShyftTupleMuons"+  Mpostfix,   "eta" )
    muonPhiHandle         = Handle( "std::vector<float>" )
    muonPhiLabel    = ( "pfShyftTupleMuons"+  Mpostfix,   "phi" )
    muonPfisoHandle         = Handle( "std::vector<float>" )
    muonPfisoLabel    = ( "pfShyftTupleMuons"+  Mpostfix,   "pfiso" )


    for event in eventList:
        print "Got event, foo is %s" % (params['foo'])
        event.getByLabel (muonPtLabel, muonPtHandle)
        oldtime = time.time()
        if not muonPtHandle.isValid():
            muonPts = None
        else:
            muonPts = muonPtHandle.product()
            print muonPts
        event.getByLabel (muonEtaLabel, muonEtaHandle)
        if not muonEtaHandle.isValid():
            muonPts = None
        else:
            muonPts = muonEtaHandle.product()
            print muonPts
        event.getByLabel(muonPhiLabel, muonPhiHandle)
        if not muonPhiHandle.isValid():
            muonPts = None
        else:
            muonPts = muonPhiHandle.product()
            print muonPts
        event.getByLabel (muonPfisoLabel, muonPfisoHandle)
        if not muonPfisoHandle.isValid():
            muonPts = None
        else:
            muonPts = muonPfisoHandle.product()
            print muonPts
        print "-- time %s" % (time.time() - oldtime)



        yield None
    print "Shutting down, foo is %s" % params['foo']

eventTee = itertools.tee(getEvents())
for x in itertools.izip_longest(runOnce(paramList[0], eventTee[0]), runOnce(paramList[1], eventTee[1])):
    print "Got an outer event %s" % (x,)   
