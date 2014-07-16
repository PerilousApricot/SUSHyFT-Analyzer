#! /usr/bin/env python

import optparse
import ROOT
import sys
import re
import pprint
import os

commentRE    = re.compile (r'#.+$')
emptyCmdRE   = re.compile (r'^\s*\+\s*(\w+)\s*$')
commandRE    = re.compile (r'^\s*\+\s*(\w+)\s+(\S.*?)\s*$')
colonSlashRE = re.compile (r':/(.*)')
nonSpaceRE   = re.compile (r'\S')
spacesRE     = re.compile (r'\s+')
filenameRE   = re.compile (r'^\s*(\S+)')
doublePipeRE = re.compile (r'\|\|')
rebinRE      = re.compile (r'(\S.+):(\d+),(\d+)')
renormRE     = re.compile (r'(\S.+):(\S+)')
colonRE      = re.compile (r'\s*:\s*')
commaRE      = re.compile (r'\s*,\s*')
carotRE      = re.compile (r'(.+)\^$')
varSetRE     = re.compile (r'(\S+)=(\S+)')
parenRE      = re.compile (r'\((\S+?)\)')

class RebinInfo (object):

    def __init__ (self, maxBin, rebin):
        self.maxBin = maxBin
        self.rebin  = rebin


    def __str__ (self):
        return "max %d rebin %d" % (self.maxBin, self.rebin)

    
    def copyAndModify (self, oldHist, newName):
        if not self.maxBin:
            # just rebin and nothing else
            oldHist.Rebin (self.rebin)
            # this isn't strictly necessary in this piece of code, but
            # if we move it somewhere else, we'll want to make sure
            # the new histogram has the new name.
            oldHist.SetName (newName)
            return oldHist
        # move the old histogram out of the way in case the old
        # histogram and new histogram have the same name
        oldHist.SetName( oldHist.GetName() + '_oldpleaseignore' )
        oldBins = oldHist.GetNbinsX()
        if oldBins < self.maxBin:
            raise RuntimeError, "Number of bins (%d) is less than new edge (%d) for '%s'" \
                  % (oldBins, self.maxBin, oldHist.GetName())
        maximum = oldHist.GetBinLowEdge (self.maxBin) + \
                  oldHist.GetBinWidth   (self.maxBin)
        minimum = oldHist.GetBinLowEdge (1)
        newHist = ROOT.TH1F (newName, newName, self.maxBin, minimum, maximum)
        for binIndex in range (0, self.maxBin): # 0 to maxBin - 1
            newHist.SetBinContent( binIndex, oldHist.GetBinContent (binIndex) )
        total = 0
        # overflow is oldBins + 1 and range is exclusive, so I need to add
        # another '1'
        for binIndex in range (self.maxBin, oldBins + 1 + 1):
            total += oldHist.GetBinContent (binIndex)
        newHist.SetBinContent (self.maxBin, total)
        newHist.Rebin (self.rebin)
        return newHist



class FileObject (object):

    allowedCommands = ['keep', 'sub', 'rebin', 'renorm', 'scale', 'add',
                       'varDefs', 'renormTo']
    gSub      = []
    gRebin    = []
    gAdd      = []
    gAddKeep  = []
    gRenormTo = []
    gScale    = []
    gVarDefs  = {}
    gVarIsSet = set()

    def __init__ (self):
        self.keep     = []
        self.sub      = []
        self.rebin    = []
        self.renorm   = []
        self.scale    = []
        self.add      = []
        self.renormTo = []
        self.varDefs  = []
        self.histDict = {}
        self.nameDict = {}
        self.finished = False

    def __str__ (self):
        return pprint.pformat (self.__dict__)


    def fileExists (self):
        if self.name == 'global':
            return False
        if not os.path.exists (self.name):
            print "Warning: '%s' does not exist" % self.name
            return False
        return True


    def _finish (self):
        # keep
        keep = self.keep
        self.keep = []
        for piece in keep:
            self.keep.append( re.compile (piece, re.IGNORECASE) )
        # sub
        sub = self.sub
        self.sub = []
        for chunk in sub:
            pieces = doublePipeRE.split (chunk)
            if len (pieces) != 2:
                raise RuntimeError, "Sub Error: %s must have exactly two pieces separated by '||'" % chunk
            self.sub.append( ( (re.compile (pieces[0], re.IGNORECASE),
                                pieces[1]) ) )
        # rebin
        rebin = self.rebin
        self.rebin = []
        for chunk in rebin:
            match = rebinRE.search (chunk)
            if not match:
                raise RuntimeError, "Rebin error '%s' does not match a rebin format." % chunk
            self.rebin.append( ( re.compile (match.group(1), re.IGNORECASE),
                                 RebinInfo( int( match.group(2) ),
                                            int( match.group(3) )) ) )
        # renorm
        renorm = self.renorm
        self.renorm = []
        for chunk in renorm:
            match = renormRE.search (chunk)
            if not match:
                raise RuntimeError, "Renorm error '%s' does not match a renorm format." % chunk
            try:
                self.renorm.append( ( re.compile (match.group(1), re.IGNORECASE),
                                      float( match.group(2) ) ) )
            except:
                raise RuntimeError, "Renorm error '%s' is not a valid format" % chunk
        # scale
        scale = self.scale
        self.scale = []
        for chunk in scale:
            match = renormRE.search (chunk)
            if not match:
                raise RuntimeError, "Scale error '%s' does not match a scale format." % chunk
            try:
                self.scale.append( ( re.compile (match.group(1), re.IGNORECASE),
                                      float( match.group(2) ) ) )
            except:
                raise RuntimeError, "Scale error '%s' is not a valid format" % chunk
        # add
        for chunk in self.add:
            two = colonRE.split (chunk)
            if len (two) != 2:
                raise RuntimeError, "Add error '%s' not recognized" % chunk
            target = two[0]
            sources = commaRE.split (two[1])
            if len (sources) < 1:
                raise RuntimeError, "Add error '%s' not enough input pieces" % chunk
            sourceREs = []
            for source in sources:
                match = carotRE.match (source)
                if match:
                    regex = re.compile( match.group(1) )
                    FileObject.gAddKeep.append (regex)
                else:
                    regex = re.compile (source)
                sourceREs.append (regex)
            FileObject.gAdd.append( (target, sourceREs) )
        # renormTo
        for chunk in self.renormTo:
            two = colonRE.split (chunk)
            if len (two) != 2:
                raise RuntimeError, "RenormTo error '%s' not recognized" % chunk
            target = two[0]
            sources = commaRE.split (two[1])
            if len (sources) < 1:
                raise RuntimeError, "RenormTo error '%s' not enough input pieces" % chunk
            sourceREs = []
            for source in sources:
                regex = re.compile (source)
                sourceREs.append (regex)
            FileObject.gRenormTo.append( (target, sourceREs) )
        # varDef
        for chunk in self.varDefs:
            match = varSetRE.search (chunk)
            if not match:
                raise RuntimeError, "varDef error '%s'" % chunk
            name = match.group(1)
            if name not in FileObject.gVarIsSet:
                FileObject.gVarDefs[name] = match.group(2)
        # set global quantities
        if self.name == 'global':
            FileObject.gSub.extend  ( self.sub   )
            FileObject.gRebin.extend( self.rebin )
            FileObject.gScale.extend( self.scale )
        self.finished = True
        


    def copy (self, targetDict, **kwargs):
        if self.loadHistograms (**kwargs):
            self.updateAndWriteHistograms (targetDict, **kwargs)


    def loadHistograms (self, **kwargs):
        verbose = False
        if kwargs.get ('verbose'):
            verbose = True
        if not self.finished:
            self._finish()
        if not self.fileExists():
            #print "skipping\n%s" % self
            # don't bother
            return False
        self.rootfile = ROOT.TFile.Open (self.name)
        if not self.rootfile:
            raise RuntimeError, "Input file '%s' not found." % self.name
        dirList = [self.rootfile]
        outputList = []
        while dirList:
            currentDir = dirList.pop (0)
            path = currentDir.GetPath()
            dirName = ''
            match = colonSlashRE.search (path)
            if match:
                dirName = match.group(1)
            nextKey = ROOT.TIter( currentDir.GetListOfKeys() )
            key = nextKey()
            while key:
                hist = key.ReadObj()
                if isinstance (hist, ROOT.TDirectory):
                    dirList.append (hist)
                    key = nextKey()
                    continue
                if not isinstance (hist, ROOT.TH1):
                    # not a histogram, so don't bother
                    key = nextKey()
                    continue
                name = key.GetName()
                if dirName:
                    name = "%s/%s" % (dirName, name)
                # do we want this one
                for regex in self.keep:
                    if regex.search (name):
                        if verbose:
                            print "   loading %s" % name
                        self.histDict[name] = hist
                        self.nameDict[name] = name
                        break
                key = nextKey()
        return True


    def updateAndWriteHistograms (self, targetDict, **kwargs):
        verbose = False
        if kwargs.get ('verbose'):
            verbose = True
        # first fix any names
        sub = FileObject.gSub[:]
        sub.extend( self.sub[:] )
        for old, new in sorted( self.nameDict.iteritems() ):
            for tup in sub:
                new = tup[0].sub( tup[1], new)
            if old != new:
                if verbose:
                    print"   renaming %s to %s" % (old, new)
                self.nameDict[old] = new
        # rebin if necessary
        rebin = self.rebin[:]
        rebin.extend( FileObject.gRebin )
        for old in sorted( self.histDict.keys() ):
            hist = self.histDict[old]
            new  = self.nameDict[old]
            for tup in rebin:
                if tup[0].search (new):
                    #print "hist %s new %s : %d" % (hist, new, hist.GetNbinsX())
                    self.histDict[old] = tup[1].copyAndModify (hist, new)
                    if verbose:
                        print "   rebin %s as %s" % (new, tup[1])
                    break
        # renorm
        for old, hist in sorted( self.histDict.iteritems() ):
            new = self.nameDict[old]
            for tup in self.renorm:
                if tup[0].search (new):
                    integral = hist.Integral()
                    if not integral:
                        if tup[1]:
                            print "   Skipping renormalization of '%s (%s) because it is empty." % (old, new)
                        continue
                    hist.Scale( tup[1] / integral )
                    if verbose:
                        print "   renorm %s to %f" % (new, tup[1])
        # scale
        scale = self.scale[:]
        scale.extend( FileObject.gScale )
        for old, hist in sorted( self.histDict.iteritems() ):
            new = self.nameDict[old]
            for tup in scale:
                if tup[0].search (new):
                    hist.Scale( tup[1] )
                    if verbose:
                        print "   scale %s to %f" % (new, tup[1])
        
        # rename internally
        for old, hist in sorted (self.histDict.iteritems()):
            new = self.nameDict[old]
            if old != new:
                hist.SetName (new)
            targetDict[new] = hist


    @staticmethod
    def addHistograms (targetDict, noWarning=False):
        # loop over add commands
        for addTup in FileObject.gAdd:
            target = addTup[0]
            sourceList = addTup[1]
            addDict = {}
            # loop over existing histograms
            for name, hist in targetDict.iteritems():
                # loop over source pieces
                for regex in sourceList:
                    newName = regex.sub (target, name)
                    if newName != name:
                        # we have a match
                        addDict.setdefault (newName,[]).append( (name, hist) )
            #pprint.pprint (addDict)
            # add histograms
            for name, combList in addDict.iteritems():
                if len (sourceList) != len (combList):
                    combNames = []
                    for tup in combList:
                        #print "tup", tup
                        combNames.append (tup[0])
                    if not noWarning:
                        print "Warning: %s is not the correct size %d != (%s)" \
                              % (name, len (sourceList), combNames)
                newHist = combList[0][1].Clone (name)
                newHist.Reset()
                for tup in combList:
                    newHist.Add( tup[1] )
                    # should we keep this histogram or delete it?
                    keep = False
                    for regex in FileObject.gAddKeep:
                        if regex.search( tup[0] ):
                            keep = True
                            break
                    if not keep:
                        del targetDict[ tup[0] ]
                targetDict[name] = newHist
                #print "adding '%s'" % name
        #print '\n' * 10


    @staticmethod
    def renormToHistograms (targetDict, noWarning=False):
        # loop over renormTo commands
        for renormToTup in FileObject.gRenormTo:
            target = renormToTup[0]
            sourceList = renormToTup[1]
            renormToDict = {}
            # loop over existing histograms
            for name, hist in targetDict.iteritems():
                # loop over source pieces
                for regex in sourceList:
                    newName = regex.sub (target, name)
                    if newName != name:
                        # we have a match
                        renormToDict.setdefault (newName,[]).append( (name, hist) )
            # pprint.pprint (renormToDict)
            # renormTo histograms
            for name, combList in renormToDict.iteritems():
                if len (sourceList) != len (combList):
                    combNames = []
                    for tup in combList:
                        #print "tup", tup
                        combNames.append (tup[0])
                    if not noWarning:
                        print "Warning: %s is not the correct size %d != (%s)" \
                              % (name, len (sourceList), combNames)
                renormHist = targetDict.get(name)
                if not renormHist:
                    print "Could not find histogram '%s'.  Aborting renormalization" \
                          % name
                    continue
                integral = renormHist.Integral()
                if not integral:
                    if not noWarning:
                        print "'%s' is empty.  Renormalization stopped" \
                              % renormHist.GetName()
                    continue
                total = 0.
                for tup in combList:
                    total += tup[1].Integral()
                renormHist.Scale( total / integral )


    @staticmethod
    def _replaceWithVar (matchObject):
        name = matchObject.group(1)
        if not FileObject.gVarDefs.has_key (name):
            raise RuntimeError, "Replace error: name '%s' is not defined" \
                  % name
        return FileObject.gVarDefs[name]

    

    @staticmethod
    def filename (source):
        return parenRE.sub (FileObject._replaceWithVar, source)
        

if __name__ == "__main__":
    # Setup options parser
    parser = optparse.OptionParser \
             ("usage: %prog [options] config.txt output.root")
    parser.add_option ('--varDef', dest='varDef', action='append',
                       type='string', default=[],
                       help='Comma separated list of variable definitions')
    parser.add_option ('--verbose', dest='verbose',
                       action='store_true',
                       help='Print lots of stuff')
    parser.add_option ('--noRenormWarnings', dest='noRenormWarnings',
                       action='store_true',
                       help='Suppress any renormalization warnings')
    options, fullargs = parser.parse_args()
    args = []
    for arg in fullargs:
        if varSetRE.search (arg):
            options.varDef.append (arg)
        else:
            args.append (arg)
    if len (args) != 2:
        print "args", args
        raise  RuntimeError, "Must provide exactly one configuration file and output file"
    for chunk in options.varDef:
        pieces = commaRE.split (chunk)
        for vardef in pieces:
            match = varSetRE.search(vardef)
            if not match:
                raise RuntimeError, "--varDef error '%s' ('%s')" \
                      % (vardef, chunk)
            name = match.group(1)
            FileObject.gVarDefs[name] = match.group(2)
            FileObject.gVarIsSet.add (name)
    config = open (args[0], 'r')
    ROOT.TH1.AddDirectory (False)
    fileObject = None
    fileList = []
    for line in config:
        orig = line
        #print "::", orig,
        line = commentRE.sub ('', line)
        if not nonSpaceRE.search (line):
            continue
        # is this a command?
        match = commandRE.search (line)        
        if match:
            if not fileObject:
                raise RuntimeError, "Must give filename before giving commands"
            command = match.group(1)
            #print "command '%s'" % command, match.groups()
            arguments = spacesRE.split( match.group(2) )
            if command not in FileObject.allowedCommands:
                raise RuntimeError, "Don't know command '%s' in line '%s'" \
                      % (command, line.strip())
            getattr (fileObject, command).extend( arguments )
            continue
        if emptyCmdRE.search (line):
            print "Warning: '%s' is not a valid command. Skipping" % \
                  orig.strip()
            continue
        # is this a filename?
        match = filenameRE.search (line)
        if match:
            # clean up old fileObject
            if fileObject:
                fileList.append (fileObject)
            fileObject = FileObject()
            #print "name", match.group(1), orig
            fileObject.name = FileObject.filename( match.group(1) )
            continue
    if fileObject:
        fileList.append (fileObject)
    config.close()
    targetDict= {}
    for obj in fileList:
        print obj.name
        obj.copy (targetDict, verbose=options.verbose)
    FileObject.addHistograms (targetDict, options.noRenormWarnings)
    FileObject.renormToHistograms (targetDict, options.noRenormWarnings)
    target = ROOT.TFile.Open (args[1], 'recreate')
    for name, hist in sorted (targetDict.iteritems()):
        hist.Write()
    target.Close()
