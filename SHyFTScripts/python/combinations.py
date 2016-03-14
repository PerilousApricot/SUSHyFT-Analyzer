import itertools

class EffInfo:
    def __init__(self, index, beff, teff, flavor):
        self.index = index
        self.beff = beff
        self.teff = teff
        self.flavor =flavor
    def __lt__(self, other ):
        return self.index < other.index
    def __gt__(self, other ):
        return self.index > other.index
    def printIt(self):
        print ' [{0:2.0f} : {1:5.2f}, {2:1.0f}] '.format(
            self.index, self.eff, self.flavor ),


class EffInfoCombinations:
    """Helper class to generate probabilities to tag k out of N jets with different probabilities"""
    def __init__(self, k , verbose=False):
        self.k = k
        self.xtot = set( self.k )
        self.ik = 0
        self.efftot = 0.0
        self.verbose = verbose

    def multiplies(self, a ):
        efftot = 1.0
        for ia in a:
            efftot = efftot * ia.beff
        return efftot

    def taumultiplies(self, taus, btagged, nbtagged):
        efftot = 1.0
        for tau in taus:
            if tau in btagged:
                efftot *= tau.teff * ( 1.0 - tau.beff )
            elif tau in nbtagged:
                efftot *= tau.teff
            else:
                raise RuntimeError, "Slippery tau"
        return efftot

    def oneminusmultiplies(self, a ):
        efftot = 1.0
        for ia in a:
            efftot = efftot * (1 - ia.beff)
        return efftot

    def tauoneminusmultiplies(self, taus, btagged, nbtagged):
        efftot = 1.0
        for tau in taus:
            if tau in btagged:
                efftot *= (1 - (tau.teff * ( 1.0 - tau.beff )))
            elif tau in nbtagged:
                efftot *= (1 - tau.teff)
            else:
                raise RuntimeError, "Slippery tau"
        return efftot

    def pTag(self, ijet):
        """Method to call to get the probabiltiy to tag 'ijet' jets"""
        if self.verbose:
            print '-------------- njets = ' + str(ijet) + '----------------'
        a = itertools.combinations(self.k, ijet)
        iefftot = 0.0
        for ia in a:
            xa = set(ia)
            xb = self.xtot.difference(ia)
            eff1 = self.multiplies(xa)
            eff2 = self.oneminusmultiplies(xb)
            iefftot = iefftot + (eff1 * eff2)

            if self.verbose:
                for ixa in xa:
                    ixa.printIt()
                for ixb in xb:
                    ixb.printIt()
                print ', eff1 = {0:6.2f}, eff2 = {1:6.2f}, tot = {2:6.2f}'.format( eff1, eff2, iefftot )
        if self.verbose:
            print '------ Probability to tag {0:3.0f} jets = {1:6.4f}'.format( ijet, iefftot )
        return iefftot

    def pTag2(self, bjet, tjet):
        """Method to call to get the probabiltiy to tag several b/tau jets"""
        if self.verbose:
            print '-------------- njets = ' + str(bjet) + str(tjet) + '----------------'
        bcombo = itertools.combinations(self.k, bjet)
        tcombo = itertools.combinations(self.k, tjet)
        
        efftot = 0.0
        befftot = 0.0
        tefftot = 0.0
        for ib in bcombo:
            btagged = set(ib)
            nbtagged = self.xtot.difference(ib)
            beff1 = self.multiplies(btagged)
            beff2 = self.oneminusmultiplies(nbtagged)
            beff = beff1 * beff2
            befftot += beff

            for it in tcombo:
                ttagged = set(it)
                nttagged = self.xtot.difference(it)
                
                teff1 = self.taumultiplies(ttagged, btagged, nbtagged)
                teff2 = self.tauoneminusmultiplies(nttagged, btagged, nbtagged)
                teff = teff1 * teff2
                tefftot += teff
                
                efftot += beff * teff
                if self.verbose:
                    for ixa in btagged:
                        ixa.printIt()
                    for ixb in nbtagged:
                        ixb.printIt()
                    print ', eff1 = {0:6.2f}, eff2 = {1:6.2f}, tot = {2:6.2f}'.format( beff1, beff2, befftot )
        if self.verbose:
            print '------ Probability to tag {0:3.0f} jets = {1:6.4f}'.format( bjet, befftot )
        return befftot, tefftot, efftot
