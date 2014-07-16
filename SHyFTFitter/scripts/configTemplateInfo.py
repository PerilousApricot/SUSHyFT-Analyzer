#! /usr/bin/env python

import ROOT
import optparse, sys, re, pprint, os
from FitterConfig import FitterConfig

# global variables to be filled
histNames  = []
groupNames = []
fileName   = ''
lumi       = 1.
# number of jet and tag bins
totalDict = {}
histList  = []
# REs
commentRE  = re.compile (r'\#.+$')
trailingRE = re.compile (r'\s*$')
sampleRE   = re.compile (r'^\s*\+\s+names\s*=\s*(.+)',        re.IGNORECASE)
groupRE    = re.compile (r'^\s*\+\s+groupNames\s*=\s*(.+)',   re.IGNORECASE)
fileRE     = re.compile (r'^\s*\+\s+templateFile\s*=\s*(.+)', re.IGNORECASE)
lumiRE     = re.compile (r'^\s*\+\s+intLumi\s*=\s*(.+)',      re.IGNORECASE)
commaRE    = re.compile (r'\s*,\s*')
jetRE      = re.compile (r'_(\d+)j')
tagRE      = re.compile (r'_(\d+)t')
htRE       = re.compile (r'_hT', re.IGNORECASE)

colorDict = {
     'Top'  :   2,
     'sing' :  93,
     'Wbb'  :  56,
     'Wcc'  :  62,
     'Wc'   :  65,
     'Wqq'  :  69,
     'EW'   :  89,
     'QCD'  :  33,
     }


if __name__ == "__main__":
    # Setup options parser
    parser = optparse.OptionParser \
             ("usage: %prog [options] templates.root" \
              "Prints out info on templates.")
    parser.add_option ('--lum', dest = 'lum', type='float', default=0.,
                       help='Override integrated luminosity in config file');
    parser.add_option ("--latex", dest='latex',
                       action='store_true',
                       help="Formats output as latex table")
    parser.add_option ("--debug", dest='debug',
                       action='store_true',
                       help="Print out FitterConfig object")
    parser.add_option ('--noData', dest='noData', action='store_true',
                       default=True,
                       help='Do not display data counts')
    parser.add_option ('--Data', dest='noData', action='store_false',
                       help='Display data counts')
    parser.add_option ('--totalMC', dest='totalMC', action='store_true',
                       default=False,
                       help='Display total MC prediction counts')
    parser.add_option ('--file', dest = 'file', type='string', 
                       help='Override root file to use');
    parser.add_option ('--combineGroups', dest = 'combineGroups',
                       action='append', type='string', default=[],
                       help='Groups to combine');
    parser.add_option ('--combineSamples', dest = 'combineSamples',
                       action='append', type='string', default=[],
                       help='Samples to combine');
    parser.add_option ("--groups", dest='groups', action="append",
                       type="string", default=[],
                       help="Which groups to use")
    parser.add_option ("--samples", dest='samples', action="append",
                       type="string", default=[],
                       help="Which samples to use")
    
    ## saveGroup = optparse.OptionGroup (parser, "Save Stacks Options")
    ## saveGroup.add_option ("--saveStacks", dest='saveStacks',
    ##                      action='store_true',
    ##                      help="Saves images of stack of templates")
    ## saveGroup.add_option ("--cms", dest='cms', action='store_true',
    ##                       help="Use CMS titles, etc for plots")
    ## saveGroup.add_option ("--big", dest='big', action='store_true',
    ##                       help="Make big plots")
    ## saveGroup.add_option ("--eps", dest='eps', action='store_true',
    ##                       help='Save .eps files')
    ## parser.add_option_group (saveGroup)
    options, args = parser.parse_args()
    ROOT.gROOT.SetBatch()
    ROOT.gROOT.SetStyle('Plain')    
    if len (args) < 1:
        print "Need to provide configuration file. Aborting."
        sys.exit(1)
    configName = args[0]
    config = FitterConfig (configName, ignoreBinString=True)
    config.noData = options.noData
    config.setValuesFromArgs (args)
    #config.readConfig (configName)
    config.printMCtotal = options.totalMC
    config.latex        = options.latex
    config.setCombineGroups  (options.combineGroups)
    config.setCombineSamples (options.combineSamples)
    samples = []
    for sample in options.samples:
        samples.extend (commaRE.split (sample))
    if samples:
        config.setSamples (samples)
    groups = []
    for group in options.groups:
        groups.extend (commaRE.split (group))
    if groups:
        config.setGroups  (groups)
    if options.file:
         config.fileName = options.file
    if options.lum:
         config.lumi = options.lum
    print "info for %s:" % config.fileName
    config.printInfo()
    if options.debug:
         print "%s" % config
