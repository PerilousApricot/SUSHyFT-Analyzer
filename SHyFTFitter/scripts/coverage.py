#! /usr/bin/env python

import ROOT
import optparse
import math
import bisect
import pprint

# formula for scan pull
pull = '((measVec[0]-genVec[0])/scanPosErr[0])*(measVec[0]<genVec[0])-((measVec[0]-genVec[0])/scanNegErr[0])*(measVec[0]>genVec[0])'

# +- 1 sigma
thirtytwo  = 1 - ROOT.TMath.Erf( 1 / math.sqrt(2) )
sixteen    = thirtytwo / 2
eightyfour = 1 - sixteen
epsilon    = 1e-4

class Percentile (object):
    '''Class to hold necessary information to make discovering
    percentile information easy'''

    def __init__ (self, chain):
        self.pull       = pull
        self.chain      = chain
        self.entries    = self.numEntries (999)
        self.sixteen    = sixteen    * self.entries
        self.fifty      = 0.50       * self.entries
        self.eightyfour = eightyfour * self.entries
        self.cache      = [ (-5, self.numEntries (-5)),
                            ( 5, self.numEntries ( 5)) ]

    def __str__ (self):
        return 'entries %d sixteen %f fifty %f eightyfour %f' \
               % (self.entries, self.sixteen, self.fifty,
                  self.eightyfour)


    def numEntries (self, cut):
        '''return numEntries less than cut'''
        return chain.Draw (self.pull,
                           '%s < %f' % (self.pull,
                                        cut),
                           'goff')

    def searchFor (self, counts, start = 0):
        size = len (self.cache)
        prev = self.cache [start]
        found = False
        for index in range (start + 1, size):
            curr = self.cache [index]
            if prev[1] <= counts <= curr[1]:
                found = True
                break
            prev = curr
        if not found:
            raise RuntimeError, "searchFor failed."
        ##########################
        ## Are We Close Enough? ##
        ##########################
        if curr[1] - prev[1] == 1:
            # we're exactly where we want to be
            delta = counts - prev[1]
            return prev[0] + delta * (curr[0] - prev[0])
        xVal = 0.5 * prev[0] + 0.5 * curr[0]
        if curr[0] - prev[0] < epsilon:
            return xVal
        if curr[1] == prev[1]:
            # both of these X values have the same number of counts
            return xVal
        # if we're still here, then we weren't close enough.  Add
        # another point, but do it cleverly so that the list stays
        # sorted
        bisect.insort( self.cache, (xVal, self.numEntries (xVal) ) )
        # Start the search over again.  Since we know the piece we
        # added is larger than the 'prev', just start our search at
        # the 'prev' value instead of at the beginning of the list.
        return self.searchFor (counts, index - 1)

    def return16_50_84 (self):
        return ( self.searchFor (self.sixteen),
                 self.searchFor (self.fifty),
                 self.searchFor (self.eightyfour) )

    def printCache (self):
        pprint.pprint (self.cache)


        
if __name__ == "__main__":
    # Setup options parser
    parser = optparse.OptionParser \
             ("usage: %prog [options] input1.root [input2.root]" \
              "Plots pull distributions, etc.")
    parser.add_option ('--printCache', dest = 'printCache', action='store_true',
                       help='Print integral cache');
    options, args = parser.parse_args()
    if not args:
        raise RuntimeError, "Must provide at least one input root file"
    samples = []
    ROOT.gROOT.SetStyle ('Plain')
    ROOT.gROOT.SetBatch()
    files = []
    for arg in args:
        files.append (arg)
    chain = ROOT.TChain ('PEtree')
    for rootfile in files:
        chain.AddFile (rootfile)
    entries = chain.Draw (pull, '%s > -999' % pull, 'goff')
    percentile = Percentile (chain)
    print '16%% %f 50%% %f 84%% %f' % percentile.return16_50_84()
    if options.printCache:
        percentile.printCache()
