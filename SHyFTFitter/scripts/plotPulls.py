#! /usr/bin/env python

import ROOT
import optparse
import math
import re
import time
import sys
from pprint import pprint

# plots to make
plotsList = [
    # name       hist range       draw string
    ('Pull',     (100, -5, 5),    "(%(meas)s-%(gen)s)/%(err)s", 'gaus'),
    ('AsymPull', (100, -5, 5),  "((%(meas)s-%(gen)s)/%(posErr)s)*(%(meas)s<%(gen)s)-"
     "((%(meas)s-%(gen)s)/%(negErr)s)*(%(meas)s>%(gen)s)", 'gaus'),
    ('ScanPull', (100, -5, 5),  "((%(meas)s-%(gen)s)/%(scanPosErr)s)*(%(meas)s<%(gen)s)-"
     "((%(meas)s-%(gen)s)/%(scanNegErr)s)*(%(meas)s>%(gen)s)", 'gaus'),
    ('Err',     (),              "%(err)s"),
    ('ScanErr', (),              "(%(scanPosErr)s-%(scanNegErr)s)/2"),
    ('Meas',    (),              "%(meas)s-%(gen)s", 'gaus'),
    ('Gen',     (),              "%(gen)s", 'gaus'),
    ('RelUnc',  (),              "%(err)s/abs(%(meas)s)"),
    ]

ggRE = re.compile (r'>>')

def prettyUnit (origEdge, unit, asSmall):
    sign = 1
    if origEdge < 0:
        sign = -1
    edge = abs (origEdge) / unit
    if asSmall:
        return int( math.floor (edge * 0.9) * sign )
    else:      
        return int( math.ceil  (edge * 1.2) * sign )


def calcHistRange (chain, drawString, cuts = ''):
    """Calculates a 'pretty' range for a given draw command"""
    if cuts:
        chain.Draw (drawString, cuts)
    else:
        chain.Draw (drawString, cuts)
    htemp   = ROOT.MakeNullPointer( ROOT.TH1F )
    htemp   = chain.GetHistogram()
    entries = htemp.GetEntries()
    bins    = htemp.GetNbinsX()
    lowerBin = bins
    for bin in range (1, bins):
        percent = htemp.Integral (1, bin) / entries
        if percent > 0.025:
            lowerBin = bin 
            break
    upperBin = 1
    for bin in range (bins, 1, -1):
        percent = htemp.Integral (bin, bins) / entries
        if percent > 0.025:
            upperBin = bin
            break
    lowerEdge = htemp.GetBinLowEdge (lowerBin)
    upperEdge = htemp.GetBinLowEdge (upperBin) + htemp.GetBinWidth (upperBin)
    if upperBin <= 2:
        variable = ggRE.sub ('', drawString)
        cut = "%s < %f" % (variable, upperEdge)
        return calcHistRange (chain, drawString, cut)
    lowAsSmall = False
    uppAsSmall = False
    if lowerEdge > 0:
        lowAsSmall = True
    elif upperEdge < 0:
        uppAsSmall = True
    maxValue = max( abs(lowerEdge), abs(upperEdge))
    unit = 10 ** math.floor( math.log10( maxValue ) )    
    lowUnit = prettyUnit (lowerEdge, unit, lowAsSmall)
    uppUnit = prettyUnit (upperEdge, unit, uppAsSmall)
    diff = lowUnit - uppUnit
    bins = int (math.ceil (100 / diff)) * diff
    return bins, lowUnit * unit, uppUnit * unit


if __name__ == "__main__":
    # Setup options parser
    parser = optparse.OptionParser \
             ("usage: %prog [options] input1.root [input2.root]" \
              "Plots pull distributions, etc.")
    parser.add_option ('--output', dest = 'output', type='string',
                       default='images/',
                       help='output name for plots');
    parser.add_option ('--ext', dest = 'ext', type='string',
                       default='png',
                       help='extention for plots');
    parser.add_option ('--samples', dest = 'samples', type='string',
                       action='append',
                       default=[],
                       help='Only consider these samples');
    options, args = parser.parse_args()
    if not args:
        raise RuntimeError, "Must provide at least one input root file"
    samples = []
    commaRE = re.compile (r'\s*,\s*')    
    for piece in options.samples:
        samples.extend( commaRE.split (piece) )
    ROOT.gROOT.SetStyle ('Plain')
    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptFit (1111)
    files = []
    for arg in args:
        files.append (arg)
    rootfile = ROOT.TFile.Open (files[0])
    namesVec     = ROOT.MakeNullPointer( ROOT.vector("string") )
    bscNamesVec  = ROOT.MakeNullPointer( ROOT.vector("string") )
    scanNamesVec = ROOT.MakeNullPointer( ROOT.vector("string") )
    rootfile.GetObject ("names", namesVec)
    try:
        rootfile.GetObject ("bscNames", bscNamesVec)
    except:
        bscNamesVec = []
    try:
        rootfile.GetObject ("scanNames", scanNamesVec)
    except:
        scanNamesVec = []
    normNames = []
    bscNames  = []
    scanNames = []
    for name in namesVec:
        normNames.append (name)
    for name in bscNamesVec:
        bscNames.append (name)
    for name in scanNamesVec:
        scanNames.append (name)        
    #print 'names', normNames, '\nbsc  ', bscNames,'\nscan ', scanNames
    scanDict = {}
    for scanName in scanNames:
        found = False
        for index, nominalName in enumerate (normNames):
            if scanName == nominalName:
                scanDict[scanName] = index
                found = True
                break
        if not found:
            raise RuntimeError, "Unknown scan name '%s'" % scanName
    #print "scanDict", scanDict
    #sys.exit()
    plotVars = []
    for index, name in enumerate (normNames):
        if samples and name not in samples:
            continue
        scanIndex = scanDict.get (name)
        if scanIndex is None:
            plotVars.append( (name,
                              "measVec[%d]"   % index,
                              "genVec[%d]"    % index,
                              "errVec[%d]"    % index,
                              "posErrVec[%d]" % index,
                              "negErrVec[%d]" % index) )
        else:
            plotVars.append( (name,
                              "measVec[%d]"    % index,
                              "genVec[%d]"     % index,
                              "errVec[%d]"     % index,
                              "posErrVec[%d]"  % index,
                              "negErrVec[%d]"  % index,
                              "scanPosErr[%d]" % scanIndex,
                              "scanNegErr[%d]" % scanIndex) )
    for index, name in enumerate (bscNames):
        plotVars.append( (name,
                          "measBsc[%d]" % index,
                          "genBsc[%d]"  % index,
                          "errBsc[%d]"  % index) )
    #pprint (plotVars)
    c1 = ROOT.TCanvas('c1', 'c1', 600, 400)
    chain = ROOT.TChain ('PEtree')
    for rootfile in files:
        chain.AddFile (rootfile)

    for tup in plotVars:
        ## if tup[0] != 'sing':
        ##     continue
        ## #print tup
        stringDict = { 'name' : tup[0],
                       'meas' : tup[1],
                       'gen'  : tup[2],
                       'err'  : tup[3]}
        if len (tup) > 4:
            stringDict['posErr'] = tup[4]
            stringDict['negErr'] = tup[5]
        if len (tup) > 6:
            stringDict['scanPosErr'] = tup[6]
            stringDict['scanNegErr'] = tup[7]
        for plot in plotsList:
            plotName = stringDict['name'] + plot[0]
            print 'name', plotName
            #drawString = plot[2] % stringDict
            try:
                drawString = plot[2] % stringDict
            except:
                print "skipping '%s'" % plotName
                continue
            if len (plot[1]):
                tup = plot[1]
            else:
                tup = calcHistRange (chain, drawString)
            drawString = plot[2] % stringDict + ">>" + plotName
            print drawString
            hist = ROOT.TH1F (plotName, plotName,
                              tup[0], tup[1], tup[2])
            chain.Draw (drawString)
            if len (plot) > 3:
                hist.Fit( plot[3] )
            c1.Print ('%s%s.%s' % (options.output, plotName,
								   options.ext))
    print "Done."

    del chain
