#!/usr/bin/env python2.6

# Import everything from ROOT
import sys
import os.path
oldargv = sys.argv[:]
sys.argv = []
import ROOT
from ROOT import TMath
ROOT.gROOT.Macro("~/rootlogon.C")
sys.argv = oldargv[:]

from DataFormats.FWLite import Events, Handle

for onePath in sys.argv[1:]:
    junk = ROOT.TFile(onePath)
    print "*** %s" % os.path.basename(onePath)
    print "Total: %s" % junk.Get('nEventsNoPU').Integral()
    for key in sorted(junk.GetListOfKeys()):
        title = key.GetTitle()
        #print key.GetTitle()
