#! /usr/bin/env python

import ROOT
import optparse, sys, re, pprint, os, re
from array import array
from FitterConfig import FitterConfig

# What Mark calls everything
histNames = ['EW', 'Wqq', 'Wc', 'Wcc', 'Wbb', 'Top']
# markers and colors
markers = (20, 21, 22, 23, 29)
colors  = ( 2,  4,  8, 1, 6,  88, 15)
legentries = ('1 Jet', '2 Jets', '3 Jets', '4 Jets', '5 Jets')
# number of jet and tag bins
numJets = 5
numTags = 2

if __name__ == "__main__":
    # Setup options parser
    parser = optparse.OptionParser \
             ("usage: %prog [options] templates.root" \
              "Prints out info on templates.")
    parser.add_option ('--writeConfig', dest='writeConfig', type='string',
                       default='',
                       help="Write out config file 'textfile:NAME,mean,width")
    parser.add_option ('--showValues', dest = 'showValues', action='store_true',
                       default=False,
                       help="Show values being fit")
    parser.add_option ('--func', dest = 'func', type='string',
                       default='pol2',
                       help='Function fit name');
    parser.add_option ('--Xaxis', dest = 'Xaxis', type='string',
                       default='shift',
                       help='X Axis title');
    parser.add_option ('--Yaxis', dest = 'Yaxis', type='string',
                       default='Multiplicative factor',
                       help='Y Axis title');
    parser.add_option ('--Xrange', dest = 'Xrange', type='float',
                       default=3.,
                       help='Absolute value of X range');
    parser.add_option ('--Ymax', dest = 'Ymax', type='float',
                       default=2.,
                       help='Maximum y value');
    parser.add_option ('--Ymin', dest = 'Ymin', type='float',
                       default=2.,
                       help='Minimum y value');
    parser.add_option ("--groups", dest='groups', action="append",
                       type="string", default=[],
                       help="Which groups to use")
    parser.add_option ("--samples", dest='samples', action="append",
                       type="string", default=[],
                       help="Which samples to use")
    parser.add_option ('--offset', dest = 'offset', type='int',
                       default=0,
                       help='Offset for colors/markers');
    parser.add_option ("--name", dest='name',
                       type="string", default="btag",
                       help="Name of variation")
    parser.add_option ("--eps", dest='eps',
                       action='store_true', default=False,
                       help="Create eps files")
    parser.add_option ("--pdf", dest='pdf',
                       action='store_true', default=False,
                       help="Create pdf files")
    parser.add_option ("--debug", dest='debug',
                       action='store_true', default=False,
                       help="Print debug information")
    parser.add_option ("--legend", dest='leg',
                       action='store_true', default=False,
                       help="Create a legend on the plots")
    parser.add_option ("--noFit", dest='noFit',
                       action='store_true', default=False,
                       help="Displays no title for plot")
    parser.add_option ("--noTitle", dest='noTitle',
                       action='store_true', default=False,
                       help="Only displays points with no fit")
    parser.add_option ('--combineGroups', dest = 'combineGroups',
                       action='append', type='string', default=[],
                       help='Groups to combine');
    parser.add_option ('--onlyCombinedGroups', dest='onlyCombined',
                       action='store_true',
                       help='Do not fit individual group if used in combined group')
    parser.add_option ('--combineSamples', dest = 'combineSamples',
                       action='append', type='string', default=[],
                       help='Samples to combine');
    parser.add_option ('--text', dest='text', type='string',
                       help='Place text where title would be')
    options, args = parser.parse_args()
    if options.text:
        options.noTitle = True
    ROOT.gROOT.SetStyle('Plain')    
    ROOT.gROOT.SetBatch()
    # Make sure we have at least two files
    if len (args) < 3:
        print "Need to provide a configuration file and at least two " \
              "template rootr file. Aborting."
        sys.exit(1)
    configFile = args.pop(0)
    config = FitterConfig()
    config.noData = True
    config.readConfig (configFile)
    config.useSamples = config.sampleNames
    config.useGroups  = config.groupNames
    config.setCombineGroups  (options.combineGroups)
    config.setCombineSamples (options.combineSamples)
    groups = []
    commaRE = re.compile (r',')
    for group in options.groups:
        groups.extend (commaRE.split (group))
    options.groups = groups
    samples = []
    commaRE = re.compile (r',')
    for sample in options.samples:
        samples.extend (commaRE.split (sample))
    colonRE = re.compile (r':')
    foundZero = False
    for fileval in args:
        pieces = colonRE.split (fileval)
        if len (pieces) != 2:
            raise RuntimeError, "'%s' is not formatted file.root:value." % \
                  fileval
        hfile = ROOT.TFile.Open (pieces[0])
        if not hfile:
            raise RuntimeError, 'Can not open "%s".' % pieces[0]
        config.getNorms( hfile, pieces[1] )
        if pieces[1] == '0':
            foundZero = True
    if not foundZero:
        raise RuntimeError, "Exactly one file must have the value of '0'"
    # now that we've loaded the files and done any combining that we
    # need to do, we can now set the samples and groups.
    if samples:
        config.setSamples (samples)
    if groups:
        config.setGroups  (groups)

    config.renorm (debug=options.debug)
    if options.debug:
        print "%s" % config
        print "renorm"
        pprint.pprint (config.renormDict)

    if options.showValues:
        print '%s' % config
    legend = ROOT.TLegend (0.4, 0.65, 0.6, 0.895)
    configBody = ''
    graphList = []
    fun = ROOT.TF1("fun", options.func)
    fun.FixParameter(0,1)
    for sample, sampleDict in sorted( config.renormDict.iteritems() ):
        if options.debug:
            print "sample", sample
        canvas = ROOT.TCanvas ()
        first = True
        index = -1 + options.offset
        if options.debug:
            pprint.pprint (sampleDict)
        for group, pointsMap in sorted( sampleDict.iteritems() ):
            index += 1
            if options.groups and group not in options.groups:
                if options.debug:
                    print "skipping", group
                continue
            if options.onlyCombined and group not in config.groupNames:
                continue
            if options.debug:
                print "group", group, "index", index
            prettyName = config.prettyGroupName (group)
            xVal = array ('d', [float (x) for x in pointsMap.keys() ] )
            yVal = array ('d', pointsMap.values())
            graph = ROOT.TGraph ( len (xVal), xVal, yVal)
            graph.SetMarkerStyle( markers[index % len (markers) ] )
            graph.SetMarkerColor( colors [index % len (colors)  ]);
            graph.SetLineColor( colors [index % len (colors) ] );
            print "group", group, "color", colors [index % len (colors)  ],\
                  index, index % len (colors)
            graph.SetLineWidth( 2 );
            graph.SetMarkerSize (1.4)
            if options.noTitle:
                graph.SetTitle('')
            else:
                graph.SetTitle ("%s %s variation" % (sample, options.name) )
            legend.AddEntry(graph, prettyName, "lp")
            if not options.noFit:
                 graph.Fit (fun)
                 fit = graph.GetFunction("fun")
                 ## graph.Fit (options.func)
                 ## fit = graph.GetFunction (options.func)
                 fit.SetLineColor( colors [index % len (colors) ] );
                 fit.SetRange (-1 * options.Xrange, options.Xrange)
                 configBody += "-- %-22s : %12s :" % (config.expandGroups (group),
                                                  config.expandSamples (sample))
                 numPar = fit.GetNpar()
                 for parIndex in range (numPar):
                     param = fit.GetParameter (parIndex)
                     if parIndex:
                         configBody += "  %18.15f" % param
                     else:
                         configBody += "  %4.1f"  % param
                 configBody += "\n";
            xAxis = graph.GetXaxis()
            xAxis.SetTitle (options.Xaxis)
            xAxis.Set (1, -3, 3)
            yAxis = graph.GetYaxis()
            yAxis.SetTitle (options.Yaxis)
            graph.SetMinimum (options.Ymin)
            graph.SetMaximum (options.Ymax)
            print options.Ymax
            if first:
                graph.Draw ("AP")
                first = False
            else:
                graph.Draw ("P")
            graphList.append (graph)
        if options.leg:
            #legend.SetY1(1 - legend.GetNRows() * 0.05)
            legend.SetBorderSize(1)
            legend.SetFillColor(10)
            legend.Draw()
        if options.text:
            latex = ROOT.TLatex()
            latex.SetTextSize (0.045)
            latex.DrawLatex (-3 + 0.05 * (3 - -3),
                             1.02 * 4, options.text)
        canvasName = "%s_%s" % (sample, options.name)
        if options.pdf:
            canvasName += ".pdf"
        elif options.eps:
            canvasName += ".eps"
        else:
            canvasName += ".png"
            
        canvas.Print (canvasName)
        legend.Clear()

    #print configBody
        
    if options.writeConfig:
        match = re.match (r'(\S+):(\S+),(\S+),(\S+)', options.writeConfig)
        filename = name = mean = width = ''
        lower, upper, step = -10., 10., 0.001
        if match:
            filename = match.group(1)
            name     = match.group(2)
            mean     = match.group(3)
            width    = match.group(4)
        else:
            match = re.match (r'(\S+):(\S+)', options.writeConfig)
            if match:
                filename = match.group(1)
                name     = match.group(2)
                mean     = 0
                width    = 1
            else:
                print "uh oh"
                sys.exit()
        try:
            filehandle = open (filename, "w")
            filehandle.write ("- %s %s %s %s %s %s\n" % \
                              ( name, mean, width, lower, upper, width) )
            filehandle.write (configBody)
            filehandle.close()
        except:
            print "whoops"
