#! /usr/bin/env python

import optparse
import ROOT
import sys
import re
import os
import signal

if __name__ == "__main__":
    # Setup options parser
    parser = optparse.OptionParser \
             ("usage: %prog [options] file.root [file2.root]")
    parser.add_option ('--mean', dest='mean', action='store_true',
                       help = 'Print mean and RMS of histograms')
    parser.add_option ('--range', dest='range', action='store_true',
                       help = 'Print range info of histograms')
    parser.add_option ('--regex', dest='regex', action='append',
                       type='string', default=[],
                       help='(case insensitive) Regex string to apply to name')
    parser.add_option ('--killAtEnd', dest='killAtEnd', action='store_true',
                       help = 'Tells script to kill self at end to avoid CPU intensive ROOT cleanup when lots of histograms/files are investiaged.')
    parser.add_option ('--types', dest='types', action='store_true',
                       help='Print out name and types only')
    options, args = parser.parse_args()
    regexList = []
    commaRE = re.compile (r'\s*,\*')
    for chunk in options.regex:
        for piece in commaRE.split (chunk):
            regexList.append( re.compile(piece, re.IGNORECASE) )
    if not len (args):
        parser.print_usage()
        sys.exit()
    colonSlashRE = re.compile (r':/(.*)')
    spaces = 24
    hist = None
    for filename in args:
        currentDir = None
        print "\n%s:" % filename
        rootfile = ROOT.TFile.Open (filename)
        dirList = [rootfile]
        outputList = []
        while dirList:
            if currentDir:
                del currentDir
            currentDir = dirList.pop (0)
            path = currentDir.GetPath()
            dirName = ''
            match = colonSlashRE.search (path)
            if match:
                dirName = match.group(1)
            nextKey = ROOT.TIter( currentDir.GetListOfKeys() )
            key = nextKey()
            while key:
                #print key.GetName()
                if hist:
                    del hist
                hist = key.ReadObj()
                if isinstance (hist, ROOT.TDirectory):
                    dirList.append (hist)
                    key = nextKey()
                    continue
                if not options.types and not isinstance (hist, ROOT.TH1):
                    # not a histogram, so don't bother
                    key = nextKey()
                    continue
                name = key.GetName()
                classname = '???'
                try:
                    classname = hist.ClassName()
                except:
                    pass
                    
                if dirName:
                    name = "%s/%s" % (dirName, name)
                if regexList:
                    found = False
                    for regex in regexList:
                        if regex.search(name):
                            found = True
                            break
                    if not found:
                        key = nextKey()
                        continue
                if options.types:
                    payload = classname
                else:
                    payload = "%12.5f" % hist.Integral()
                    if options.mean:
                        payload += "  (%12.5f +- %12.5f) " % \
                                   (hist.GetMean(), hist.GetRMS())
                    if options.range:
                        num = hist.GetNbinsX()
                        payload += "  (%3d: %5s - %5s)" % \
                                   (num,
                                    hist.GetBinLowEdge (1),
                                    hist.GetBinLowEdge (1 + num))
                outputList.append( (name, payload) )
                key = nextKey()
        # while dirList
        for piece in outputList:
            if len(piece[0]) > spaces:
                spaces = len (piece[0])
        if options.types:
            form = '   %%-%ds : %%s' % spaces
        else:
            form = '   %%-%ds   %%s' % spaces
        for piece in outputList:
            print form % piece
        #del rootfile
    #print "Cleaning up.  If this is taking a long time, please hit ctrl-z and then type 'kill %'"
    
    # one shouldn't need to do this, but it just refuses to exit when
    # one processes a lot of files.
    if options.killAtEnd:
        os.kill( os.getpid(), signal.SIGKILL )
 
