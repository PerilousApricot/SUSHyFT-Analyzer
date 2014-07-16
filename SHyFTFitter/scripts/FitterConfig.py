#! /usr/bin/env python

import sys
import re
import pprint
import ROOT
import os


class FitterConfig (object):
     """Class designed to read and parse SHyFT configuration files"""

     # Regexs
     commentRE       = re.compile (r'\#.+$')
     trailingRE      = re.compile (r'\s*$')
     equalRE         = re.compile (r'(\w+)=(\S+)')
     commaRE         = re.compile (r'\s*,\s*')
     colonRE         = re.compile (r'\s*:\s*')
     plusRE          = re.compile (r'\s*\+\s*')
     underRE         = re.compile (r'_')
     jetRE           = re.compile (r'_(\d+)j')
     tagRE           = re.compile (r'_(\d+)t')
     sampleRE        = re.compile (r'^\s*\+\s+names\s*=\s*(.+)',
                                   re.IGNORECASE)
     groupRE         = re.compile (r'^\s*\+\s+groupNames\s*=\s*(.+)',
                                   re.IGNORECASE)
     binsRE          = re.compile (r'^\s*\+\s+binsPerGroup\s*=\s*(.+)',
                                   re.IGNORECASE)
     groupStringsRE  = re.compile (r'^\s*\+\s+groupStrings\s*=\s*(.+)',
                                   re.IGNORECASE)
     fileRE          = re.compile (r'^\s*\+\s+templateFile\s*=\s*(.+)',
                                   re.IGNORECASE)
     lumiRE          = re.compile (r'^\s*\+\s+intLumi\s*=\s*(.+)',
                                   re.IGNORECASE)
     cleanREList     = [ re.compile (r'_hT', re.IGNORECASE),
                         re.compile (r'_met', re.IGNORECASE),
                         re.compile (r'_eta', re.IGNORECASE),
                         re.compile (r'_svm', re.IGNORECASE),
                         re.compile (r'_lepEta', re.IGNORECASE),
                         ]

     # colors
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


     def __init__ (self, configName='', **kwargs):
          """Initializaiton"""
          self.sampleNames     = []
          self.groupNames      = []
          self.groupStrings    = {}
          self.groupBins       = {}
          self.fitterHists     = {}
          self.fileName        = ''
          self.lumi            = 1.
          self.totalDict       = {}
          self.intDict         = {}
          self.histList        = []
          self.useSamples      = []
          self.useGroups       = []
          self.valueDict       = {}
          self.latex           = False
          self.printMCtotal    = False
          self.combineGroups   = []
          self.combineSamples  = []
          self.renormDict      = {}
          self.noData          = False
          self.latex           = False
          self.ignoreBinString = kwargs.get('ignoreBinString')
          if configName:
               self.readConfig (configName)


     def __str__ (self):
          """stringifys self"""
          return pprint.pformat (self.__dict__)


     def readConfig (self, configName):
          """Reads and parses configuration files"""
          try:
               configFile = open (configName, 'r')
          except:
               raise RuntimeError, "Could not open '%s'." % configName
          groupStrings = []
          binStrings = []
          for line in configFile:
               line  = FitterConfig.trailingRE.sub ('', line)
               line  = FitterConfig.commentRE.sub ('', line)
               # histogram names
               match = FitterConfig.sampleRE.search (line)
               if match:
                    self.sampleNames.extend( \
                         FitterConfig.commaRE.split( match.group(1) ) )
                    continue
               # groupNames
               match = FitterConfig.groupRE.search (line)
               if match:
                    self.groupNames.extend( \
                         FitterConfig.commaRE.split( match.group(1) ) )
                    continue
               # binsPerGroup
               match = FitterConfig.binsRE.search (line)
               if match:
                    binStrings.extend( \
                         FitterConfig.commaRE.split( match.group(1) ) )
                    continue
               # groupStrings
               match = FitterConfig.groupStringsRE.search (line)
               if match:
                    groupStrings.extend( \
                         FitterConfig.commaRE.split( match.group(1) ) )
                    continue
               # filename
               match = FitterConfig.fileRE.search (line)
               if match:
                    self.fileName = match.group(1)
                    continue
               # lumi
               match = FitterConfig.lumiRE.search (line)
               if match:
                    self.lumi = float (match.group(1))
                    continue
          configFile.close()

          if len (groupStrings) != len (self.groupNames):
               raise RuntimeError, "Strings and names lengths don't match strings %d != names %d" % (len (groupStrings), len (self.groupNames) )
          for name, pretty in zip (self.groupNames, groupStrings):
               self.groupStrings[name] = FitterConfig.underRE.sub(' ', pretty)
          if len (binStrings) != len (self.groupNames):
               if self.ignoreBinString:
                    return
               raise RuntimeError, "Bins (%d) and names (%d) are different lengtsh." \
                     % (len(binStrings), len(self.groupNames))
          start = 1
          for name, bin in zip (self.groupNames, binStrings):
               try:
                    intBin = int (bin)
                    self.groupBins[name] = (intBin, start, start + intBin - 1)
                    start += intBin
               except:
                    raise RuntimeError, "Bin string '%s' is not an integer" \
                          % bin
          
                   


          # now that we've read everything in, check to make sure
          # various things are kosher:
          for key in self.valueDict.keys():
               if key not in self.sampleNames:
                    raise RuntimeError, "'%s' is not one of the sample names (%s)" % (key, self.sampleNames)


     def setSamples (self, samples):
          '''Set a subset of samples'''
          for sample in samples:
               if sample not in self.sampleNames:
                    raise RuntimeError, "Sample '%s' does not match a given sample" % sample
          self.sampleNames = samples


     def setGroups (self, groups):
          '''Set a subset of groups'''
          for group in groups:
               if group not in self.groupNames:
                    raise RuntimeError, "Group '%s' does not match a given group" % group
          self.groupNames = groups


     def setValuesFromArgs (self, args):
          ''' '''
          for arg in args[1:]:
               match = FitterConfig.equalRE.match (arg)
               if match:
                    self.valueDict[ match.group(1) ] = float( match.group(2) )


     def setCombineGroups (self, listOfStrings):
          '''Expects a list of strings where each string is of the
          format:
          newGroup:oldGroup1+oldGroup2+oldGroup3,newGroup2:oldGroup5+oldGroup6
          '''
          allOldGroups = []
          for chunk in listOfStrings:
               lines = FitterConfig.commaRE.split (chunk)
               for line in lines:
                    twoPieces = FitterConfig.colonRE.split (line)
                    if len(twoPieces) != 2:
                         raise RuntimeError, "setCombineGroups() expects strings to have exactly 1 ':'.  Problem with '%s'" % line
                    newGroup = twoPieces[0]
                    oldGroups = FitterConfig.plusRE.split (twoPieces[1])
                    # let's make sure all of the old groups are actually groups
                    for group in oldGroups:
                         if group not in self.groupNames:
                              raise RuntimeError, "group '%s' not in %s" \
                                    % (group, self.groupNames)
                         if group in allOldGroups:
                              raise RuntimeError, "You can only use each group once. '%s' duplicated" % group
                         allOldGroups.append (group)
                    self.combineGroups.append( (newGroup, oldGroups) )
          


     def setCombineSamples (self, listOfStrings):
          '''Expects a list of strings where each string is of the
          format:
          newSample:oldSample1+oldSample2+oldSample3,newSample2:oldSample5+oldSample6
          '''
          allOldSamples = []
          for chunk in listOfStrings:
               lines = FitterConfig.commaRE.split (chunk)
               for line in lines:
                    twoPieces = FitterConfig.colonRE.split (line)
                    if len(twoPieces) != 2:
                         raise RuntimeError, "setCombineSamples() expects strings to have exactly 1 ':'.  Problem with '%s'" % line
                    newSample = twoPieces[0]
                    oldSamples = FitterConfig.plusRE.split (twoPieces[1])
                    # let's make sure all of the old samples are actually samples
                    for sample in oldSamples:
                         if sample not in self.sampleNames:
                              raise RuntimeError, "sample '%s' not in %s" \
                                    % (sample, self.sampleNames)
                         if sample in allOldSamples:
                              raise RuntimeError, "You can only use each sample once. '%s' duplicated" % sample
                         allOldSamples.append (sample)
                    self.combineSamples.append( (newSample, oldSamples) )


     def loadHistogramsFromFitterOutput (self, filename):
          ROOT.TH1.AddDirectory (False)
          samples = self.useSamples or self.sampleNames
          groups  = self.useGroups  or self.groupNames
          hfile = ROOT.TFile.Open (filename)
          if not hfile:
               raise RuntimeError, "Could not open root file '%s'" % filename
          
          for sample in samples:
               if sample == 'Data':
                    useName = 'newdata'
               else:
                    useName = '%s_updated' % sample
               hist = hfile.Get (useName).Clone(useName + '_')
               if not hist:
                    raise RuntimeError, "could not find '%s'" % useName
               self.fitterHists[sample] = hist
          hfile.Close()
          


     def getNormsFromFitterOutput (self, filename):
          """ """
          self.loadHistogramsFromFitterOutput (filename)
          samples = self.useSamples or self.sampleNames
          groups  = self.useGroups  or self.groupNames
          for sample in samples:
               total = 0.
               hist = self.fitterHists[sample]
               for group in groups:
                    binTup = self.groupBins[group]
                    integral = hist.Integral (binTup[1], binTup[2])
                    self.intDict.setdefault( sample, {}).\
                         setdefault (group, integral)
                    total += integral
               self.totalDict[sample] = total
          self._combineGroupsAndSections()


     def getNorms (self, hfile, value = None):
          """ """
          samples = self.useSamples or self.sampleNames
          groups  = self.useGroups  or self.groupNames
          # print hfile.GetName(), "groups", groups, "samples", samples
          for sample in samples:
               if sample == 'Data' and self.noData:
                    continue
               total = 0.
               for group in groups:
                    histName = sample + group
                    hist = hfile.Get (histName)
                    if not hist:
                         print "Warning: Skipping %s!" % histName
                         continue
                    integral = hist.Integral()
                    integral = hist.Integral()
                    if value is not None:
                         self.intDict.setdefault( sample, {}).\
                              setdefault( group, {}).\
                              setdefault( value, integral)
                    else:
                         total += integral
                         self.intDict.setdefault( sample, {}).\
                              setdefault( group, integral)
               # after looping over all groups, let's take the total
               self.totalDict[sample] = total
          self._combineGroupsAndSections()


     def _combineGroupsAndSections (self):
          ####################
          ## Combine Groups ##
          ####################
          if self.combineGroups:
               # first, combine group counts
               for sample, localDict in self.intDict.iteritems():
                    for combineTup in self.combineGroups:
                         newGroup  = combineTup[0]
                         oldGroups = combineTup[1]
                         oldGroupList  = []
                         for old in oldGroups:
                              oldGroupList.append  ( localDict[old] )
                         if isinstance (oldGroupList[0], dict):
                              localDict[newGroup] = {}
                              FitterConfig._addMapValues (localDict[newGroup],
                                                          oldGroupList)
                         else:
                              localDict[newGroup] = 0.  
                              for groupValue in oldGroupList:
                                   localDict[newGroup] += groupValue
               # now fix group names
               toRemove = []
               for group in self.combineGroups:
                    toRemove.extend( group[1] )
               newGroupNames = []
               for group in self.groupNames:
                    if group not in toRemove:
                         newGroupNames.append (group)
               for group in self.combineGroups:
                    newGroupNames.append( group[0] )
               self.groupNames = newGroupNames
          ######################
          ## Combine Sections ##
          ######################
          if self.combineSamples:
               # copy sample counts
               for sampleTuple in self.combineSamples:
                    newSample  = sampleTuple[0]
                    oldSamples = sampleTuple[1]
                    self.intDict[newSample] = {}
                    oldSampleList = []
                    oldFactorList = []
                    for old in oldSamples:
                         oldSampleList.append( self.intDict[old] )
                         oldFactorList.append( self.valueDict.get (old, 1.) )
                    FitterConfig._addMapValues (self.intDict[newSample],
                                                oldSampleList, oldFactorList)
                    ## # fix totals
                    self.totalDict[newSample] = 0.
                    for sample in oldSamples:
                         self.totalDict[newSample] += \
                             self.totalDict[sample] * \
                             self.valueDict.get (sample, 1.) 
               # now fix sample names
               toRemove = []
               for sample in self.combineSamples:
                    toRemove.extend( sample[1] )
               newSampleNames = []
               for sample in self.sampleNames:
                    if sample not in toRemove:
                         newSampleNames.append (sample)
               for sample in self.combineSamples:
                    newSampleNames.append( sample[0] )
               self.sampleNames = newSampleNames
          ##############
          ## Total MC ##
          ##############
          if self.printMCtotal:
               total = 'Total Pred'
               self.intDict[total] = {}
               sampleDictList = []
               factorList     = []
               self.totalDict[total] = 0.
               for sample in self.sampleNames:
                    if sample == 'Data':
                         continue
                    self.totalDict[total] += self.totalDict[sample] \
                                             * self.valueDict.get (sample, 1.)
                    sampleDictList.append( self.intDict[sample] )
                    factorList.append( self.valueDict.get (sample, 1.) )
               FitterConfig._addMapValues (self.intDict[total],
                                           sampleDictList, factorList)
               if self.sampleNames[0] == 'Data':
                    self.sampleNames.insert (1, total)
               else:
                    self.sampleNames.insert (0, total)



     def renorm (self, renormPoint = '0', **kwargs):
          '''Makes a renormalized copy of intDict'''
          self.renormDict = {}
          for sample in self.sampleNames:
               sampleDict = self.intDict[sample]
               if kwargs.get('debug'):
                    print "renorm sampleDict"
                    pprint.pprint (sampleDict)
               self.renormDict[sample] = {}
               for group, pointsMap in sorted( sampleDict.iteritems() ):
                    renorm = pointsMap[ renormPoint ]                    
                    if not renorm:
                         print "skipping group", group, "sample", sample
                         #raise RuntimeError, "Unshifted value can't be 0."
                         continue
                    for key, value in pointsMap.iteritems():
                         value /= renorm
                         self.renormDict[sample].setdefault (group, {}).\
                              setdefault( key, value )


     def printInfo (self, **options):
          '''Prints out all requested template information'''
          if not options.get ('alreadyLoaded'):
               if not os.path.exists (self.fileName):            
                    print "Can't find '%s'.  Aborting." % self.fileName
                    sys.exit(1)
               hfile = ROOT.TFile.Open (self.fileName)
               self.getNorms (hfile)
          print "Integrated Luminosity: %f" % self.lumi
          numNames = len (self.sampleNames)
          if self.noData:
               numNames -= 1
          if self.latex:
               print "\n\\begin{tabular}{| l |%s}\n\\hline" \
                     % (' r |' * numNames )
                                                   
          # Since we want 1..numTags inclusive, don't forget the '+ 1'
          maxLen = max([len(self.prettyGroupName (name)) \
                        for name in self.groupNames])
          self.prettyFormat = "%%-%ds" % maxLen
          print ' ' * (maxLen),
          if self.latex:
               print '',
          for sample in self.sampleNames:
               if self.noData and 'Data' == sample:
                    continue
               if self.latex:
                    print '& %10s ' % sample,
               else:
                    print '%10s ' % sample,
          if self.latex:
               print "\\\\\n\\hline"
          else:
               print
          for group in self.groupNames:
               self.printLine (group)
          if self.latex:
               print "\\hline"
          print 'Total', ' ' * (maxLen - 6),
          for sample in self.sampleNames:
               factor = self.lumi * self.valueDict.get (sample, 1.)
               if sample == 'Data'and self.noData:
                    continue
               if self.latex:
                    print " & ",
               if sample == 'Data':
                    if self.noData:
                         continue
                    print "%10d " % (self.totalDict [sample]),
               else:
                    print "%10.1f " % (self.totalDict [sample] * factor),
          if self.latex: print "\\\\\n\\hline\n\\end{tabular}\n",
          print
 

     def printLine (self, group, stack = None):
          """Prints line"""
          print self.prettyFormat % self.prettyGroupName (group), 
          first = True
          splitDict = {}
          totalMC = 0
          for sample in self.sampleNames:
               histName = sample + group
               if sample == 'Data' and self.noData:
                    continue
               if self.latex:
                    print " & ",
               # Is this data or MC
               if sample == 'Data':
                    if self.noData:
                         continue
                    # Data
                    factor = 1
                    counts = self.intDict[sample].get (group, -1)
                    if self.noData:
                         continue
                    print "%10d " % counts,
               else:
                    # MC
                    factor = self.lumi * self.valueDict.get (sample, 1.)
                    counts = self.intDict[sample].get (group, -1) * factor
                    print "%10.1f " % counts,
                    ## if stack:
                    ##      hist.SetLineColor = ROOT.kBlack
                    ##      hist.SetLineWidth = 2
                    ##      hist.SetFillColor( colorDict [sample] )
                    ##      hist.Scale (factor)
                    ##      stack.Add (hist)
                    ##      histList.append (hist)
          # for sample
          if self.latex:
               print " \\\\",
          print

     def expandGroups (self, group):
          ''' expands a group if it was combined'''
          for tup in self.combineGroups:
               if group == tup[0]:
                    return " ".join (tup[1])
          return group


     def expandSamples (self, sample):
          ''' expands a sample if it was combined'''
          for tup in self.combineSamples:
               if sample == tup[0]:
                    return " ".join (tup[1])
          return sample


     def prettyGroupName (self, group):
          """Makes group name pretty"""
          if self.groupStrings.has_key (group):
               return self.groupStrings[group]
          prettyGroup = FitterConfig.jetRE.sub (r'\1 jet ', group)
          prettyGroup = FitterConfig.tagRE.sub (r'\1 tag ', prettyGroup)
          for regex in FitterConfig.cleanREList:
               prettyGroup = regex.sub  ('', prettyGroup)
          return prettyGroup


     @staticmethod
     def _addMapValues (newMap, oldMapList, factorList=None):
          '''Assuming the map structrure of the old maps are:

          1) Identical to each other, and
          2) Only contain other maps or addable objects (e.g., floats,
             ints),

          This function will fill the new map entry copying the
          structure where the numbers will be the sums of the old
          maps.'''
          structureMap = oldMapList[0]
          for key, value in structureMap.iteritems():
               # is value a map
               if isinstance (value, dict):
                    # map
                    newMap[key] = {}
                    oldSubList = []
                    for oldMap in oldMapList:
                         oldSubList.append( oldMap[key] )
                    FitterConfig._addMapValues (newMap[key],
                                                oldSubList, factorList)
               else:
                    # addable object
                    newMap[key] = 0
                    if factorList:
                         for oldMap, factor in zip( oldMapList, factorList):
                              newMap[key] += oldMap[key] * factor
                    else:
                         for oldMap in oldMapList:
                              newMap[key] += oldMap[key]


if __name__ == "__main__":
     config = FitterConfig()
     config.readConfig ('config/cms_fitter.config')
     hfile = ROOT.TFile.Open (config.fileName)
     if not hfile:
          print "zut"
          sys.exit()
     config.getNorms (hfile)
     print '%s' % config
