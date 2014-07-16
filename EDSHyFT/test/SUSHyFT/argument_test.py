#!/usr/bin/env python2.6

from optparse import OptionParser

import sys, pprint, itertools, copy

splitArgs = [list(x[1]) for x in itertools.groupby(sys.argv[1:], lambda x: x=='++') if not x[0]]
pprint.pprint(splitArgs)
def getParser(default=None):
    parser = OptionParser()
    parser.add_option('--files', metavar='F', type='string', action='store',
                      dest='files',
                      help='Input files')

    parser.add_option('--inputListFile',dest='inputListFile',
                      help='Text file listing input files')

    # Output name to use. 
    parser.add_option('--outname', metavar='F', type='string', action='store',
                      default='shyft_fwlite' if default else None,
                      dest='outname',
                      help='output name')

    # Sample name
    parser.add_option('--sampleName', metavar='F', type='string', action='store',
                      default='Top' if default else None,
                      dest='sampleName',
                      help='output name')

    # This will use the loose selections, and negate the isolation
    # criteria
    parser.add_option('--useLoose', metavar='F', action='store_true',
                      default=False if default else None,
                      dest='useLoose',
                      help='use loose leptons (exclusive from tight)')

    # no MET cut
    parser.add_option('--noMET', metavar='F', action='store_true',
                      default=False if default else None,
                      dest='noMET',
                      help='no MET cut')

    # invert MET cut
    parser.add_option('--invMET', metavar='F', action='store_true',
                      default=False if default else None,
                      dest='invMET',
                      help='invert MET cut')

    # Using data or not. For MC, truth information is accessed.
    parser.add_option('--useData', metavar='F', action='store_true',
                      default=False if default else None,
                      dest='useData',
                      help='use data')

    # MET Unclustered energy systematics
    parser.add_option('--metSys', metavar='F', type='string', action='store',
                      default='nominal' if default else None,
                      dest='metSys',
                      help='MET Systematic variation. Options are "nominal, up, down"')
    return parser

if len(splitArgs) == 1:
    parser = getParser(default = True)
    (options, args) = parser.parse_args(splitArgs[0])
    argList = [[options, args]]
else:
    argList = []
    parser = getParser(default = True)
    (mainOptions, mainArgs) = parser.parse_args(splitArgs[0])
    for oneList in splitArgs[1:]:
        print "got one list %s" % oneList
        parser = getParser(default = False)
        (subOptions, subArgs) = parser.parse_args(oneList)
        subOptions = subOptions.__dict__
        argCopy = mainArgs[:]
        argCopy.extend(subArgs)
        optCopy = copy.copy(mainOptions).__dict__
        # subOpts have priority
        for key in subOptions:
            if subOptions[key] != None:
                optCopy[key] = subOptions[key]
        argList.append([optCopy, argCopy])

pprint.pprint(argList)


