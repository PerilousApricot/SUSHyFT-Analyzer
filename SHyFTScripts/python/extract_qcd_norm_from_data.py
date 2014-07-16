#!/usr/bin/python

"""
  extract_qcd_norm_from_data.py - Given a stitched template with un-scaled
  data/QCD and with the rest normalized to their MC predictions, perform a fit
  for each jet/tag bin to extract the QCD normalization

  QCD is pre-tagged so we get the right shape

  Author: Andrew Melo, based on magic_fit.py from Andrew Ivanov
"""
print "butt"
import sys
oldarg = sys.argv[:]
sys.argv = [oldarg[0], '-b']
from ROOT import *
gROOT.Macro("~/rootlogon.C")
sys.argv = oldarg[:]
# ===============
# options
# ===============
from optparse import OptionParser
parser = OptionParser()
# New args added 9/2013
parser.add_option('--stitched-input', type='string', action='append',
        dest='stitched_input',
        default=[],
        help='Input to process')
parser.add_option('--qcd-shape', type='string', action='append',
        dest='qcd_shape',
        default=[],
        help='QCD shapes we want to output')
parser.add_option('--qcd-signal', type='string', action='append',
        dest='qcd_signal',
        default=[],
        help='QCD from signal region (used for normalization)')
parser.add_option('--qcd-norm', type='string', action='append',
        dest='qcd_norm',
        default=[],
        help='QCD to use in mini/MET-fit')




parser.add_option('--verbose',action='store_true',
        default=False,
        dest='verbose',
        help='verbose switch')				  
parser.add_option('--simple',action='store_true',
        default=False,
        dest='simple',
        help='simple_mode')	
parser.add_option('--fit',action='store_true',
        default=False,
        dest='fit',
        help='fitting_mode')
parser.add_option('--Lumi', metavar='D', type='float', action='store',
        default=17000,
        dest='Lumi',
        help='Data Luminosity')
parser.add_option('--var', metavar='T', type='string', action='store',
        default='',
        dest='var',
        help='variable to plot')
parser.add_option('--rebin', metavar='T', type='int', action='store',
        default=1,
        dest='rebin',
        help='rebin x axes to this')
parser.add_option('--minJets', metavar='D', type='int', action='store',
        default=1,
        dest='minJets',
        help='Minimum number of jets for plots')
parser.add_option('--maxJets', metavar='D', type='int', action='store',
        default=5,
        dest='maxJets',
        help='Minimum number of jets for plots')
parser.add_option('--nJets', metavar='D', type='int', action='store',
        default=-1,
        dest='nJets',
        help='Exact number of jets for plots')
parser.add_option('--minTags', metavar='D', type='int', action='store',
        default=0,
        dest='minTags',
        help='Minimum number of tags for plots')
parser.add_option('--maxTags', metavar='D', type='int', action='store',
        default=2,
        dest='maxTags',
        help='Minimum number of tags for plots')
parser.add_option('--nTags', metavar='D', type='int', action='store',
        default=-1,
        dest='nTags',
        help='Exact number of tags for plots')
parser.add_option('--pretagMinTags', type='int', action='store',
        default=0,
        dest='pretagMinTags',
        help='Minimum tags to sweep into pretagged things')
parser.add_option('--pretagMaxTags', type='int', action='store',
        default=2,
        dest='pretagMaxTags',
        help='Maximum tags to sweep into pretagged things')


parser.add_option('--templateDir', metavar='MTD', type='string',
        default='pfShyftAnaNoMET',
        dest='templateDir',
        help='Directory from which to get high statistics templates')

parser.add_option('--qcdDir', metavar='MTD', type='string',
        default='pfShyftAnaQCDWP95NoMET',
        dest='qcdDir',
        help='Directory from which to get qcd MC statistics')

parser.add_option('--subDir', metavar='MTD', type='string',
        default='eleEB',
        dest='subDir',
        help='Directory from which to get EE or EB statistics')

parser.add_option('--outputDir', metavar='MTD', type='string',
        default='plots_772_leg',
        dest='outputDir',
        help='Directory to store output histos')

parser.add_option('--pretagDir', metavar='MTD', type='string',
        default='plots_pre',
        dest='pretagDir',
        help='Directory to store output pretag or some combination of histos')

parser.add_option('--nBin', metavar='D', type='int', action='store',
        default=300,
        dest='nBin',
        help='Number of x-axis bin to display')

parser.add_option('--version', type='string', action='store',
        default='v7',
        dest='version',
        help='root files version')

parser.add_option('--fixBin', metavar='D', type='int', action='store',
        default=1000,
        dest='fixBin',
        help='x-bin range to be fitted')

parser.add_option('--printTable', type='string', action='store',
        default='fitOut.txt',
        dest='printTable',
        help='print event counts to the file')


(options,args) = parser.parse_args()
# ==========end: options =============

# =====================================================
# class that contains all relevant plot quantities		
# =====================================================

class TDistribution:
    def __init__(self, name, var, *filenames, **hists):
        self.name = name
        self.legentry = "legentry"
        self.file = TFile(filenames[0])
        keys = sorted(hists.keys())
                    #print keys
        print "opening %s with key %s" % (filenames[0], keys[0])
        self.hist = None
        for ihist in keys:
            if not self.hist:
                self.hist = self.file.Get(ihist)
                if not self.hist:
                    print "  Histogram (%s) not found, skipping" % ihist
                continue
            histo = self.file.Get(ihist)
            try:	
                self.hist.Add(histo,hists[ihist])
                print "  Added %s to %s" % (ihist, keys[0])
            except (TypeError, AttributeError):
                if options.verbose :
                    print "No histogram ", ihist, "in file", filenames[0], " ... skipping"

                continue					
        for filename in filenames[1:] :
            file = TFile(filename)
            keys = sorted(hists.keys())
            for ihist in keys :
                histo = file.Get(ihist)
                try:
                    self.hist.Add(histo,hists[ihist])
                except TypeError:
                    if options.verbose :
                        print "No histogram ", ihist, "in file", filename, " ... skipping"
                    continue
        if options.rebin > 1 :
            self.hist.Rebin(options.rebin)
        self.integerNorm = self.hist.Integral()
# RooFit quantities:
        self.SF = RooRealVar(name+"SF",name+"SF",1.,0.,10000.)
        self.norm = RooRealVar(name+"Norm",name+"Norm",self.hist.Integral())
        self.N = RooFormulaVar(name+"N",name+"SF*"+name+"Norm",RooArgList(self.norm,self.SF))
        if var != 0 :				
            self.set = RooDataHist(name+"Set",name+"Set",RooArgList(var),self.hist)
            self.pdf = RooHistPdf(name+"Pdf",name+"Pdf",RooArgSet(var), self.set)

# =========== end: class TDistribution =================

if options.verbose :
    print "script to create normalized plots"
dirMain     = "../RootFiles_v5/"
tempDir     = options.templateDir
qcdDir      = options.qcdDir
subDir      = options.subDir
reg         = subDir[-8:]
vx          =options.version

nJet = '{0:1.0f}'.format( options.nJets)
nTag = '{0:1.0f}'.format( options.nTags)
minJ = '{0:1.0f}'.format( options.minJets)
maxJ = '{0:1.0f}'.format( options.maxJets)
minT = '{0:1.0f}'.format( options.minTags)
maxT = '{0:1.0f}'.format( options.maxTags)
lum  = '{0:1.0f}'.format( options.Lumi)

#provide explanatory title for each variable name
##if options.var == "secvtxMass" :
##	xtitle = "Secondary Vertex Mass (GeV)," + nJet+"j_"+nTag+"t"
##elif options.var == "MET" :
##	xtitle = "Missing Transverse Energy (GeV)," + nJet+"j_"+nTag+"t"
##elif options.var == "wMT" :
##        xtitle = "W Transverse Mass (GeV)," + nJet+"j_"+nTag+"t"
##elif options.var == "hT" :
##        xtitle = "hT, #sum (Jet et + MET + lep pt) (GeV)," + nJet+"j_"+nTag+"t"
##elif options.var == "elEta" :
##        xtitle = "electron #eta," + nJet+"j_"+nTag+"t"
##elif options.var == "jetEt" :
##        xtitle = "#sum (jet et) (GeV)" + nJet+"j_"+nTag+"t"        
##else :
##	xtitle =""

if options.var == "secvtxMass" :
    xtitle = "Secondary Vertex Mass (GeV), #geq 3j #geq 1t"
elif options.var == "MET" :
    xtitle = "Missing Transverse Energy (GeV), #geq 3j #geq 1t"
elif options.var == "wMT" :
    xtitle = "W Transverse Mass (GeV), #geq 3j #geq 1t"
elif options.var == "hT" :
    xtitle = "hT, #sum (Jet et + MET + lep pt) (GeV), #geq 3j #geq 1t"
elif options.var == "lepPt" :
    xtitle = "electron pt, #geq 3j #geq 1t"
elif options.var == "lepEta" :
    xtitle = "electron #eta, #geq 3j #geq 1t"
elif options.var == "jetEt" :
    xtitle = "electron jet Et, #geq 3j #geq 1t"
else :
    xtitle =""


if options.nJets >= 0:
    options.minJets = options.nJets
    options.maxJets = options.nJets
    minJ            = nJet
    maxJ            = nJet
if options.nTags >= 0 :
    options.minTags = options.nTags
    options.maxTags = options.nTags
    minT            = nTag
    maxT            = nTag

# add different jets/tags histograms
def add(hists,doPretagged=False):
    retval = {}
    keys = sorted(hists.keys())
    if options.verbose:
        print "Keys: ", keys
        print hists
    myMinTags = options.minTags
    myMaxTags = options.maxTags
    if doPretagged:
        myMinTags = options.pretagMinTags
        myMaxTags = options.pretagMaxTags
    for ikey in keys[:]:
        weight = hists[ikey]
        ikey_tmp = ikey + "_" + options.var
        #print ikey_tmp
        bins = {}
        myMinTag = options.minTags
        for nJets in range(options.minJets, options.maxJets+1):
            for nTags in range(myMinTags, myMaxTags+1):
                temp_key = "_" + str(nJets) + "j_" + str(nTags) + "t"
                bins[temp_key] = weight
        bin_keys = sorted(bins.keys())
        for ibin in bin_keys[:]:
            ikey_new = ikey_tmp + ibin
            hists[ikey_new] = bins[ibin]
        del hists[ikey]	
    if True:#options.verbose :	
        for ikey in hists.iteritems():
            print "In add %s" % [ikey]

def init_var(dist) :
    NBINS = dist.hist.GetNbinsX()
    minX = (dist.hist.GetXaxis()).GetXmin()
    maxX = (dist.hist.GetXaxis()).GetXmax()
    var = RooRealVar(options.var,xtitle,minX,maxX)
    var.setBins(NBINS)
    dist.set = RooDataHist(dist.name+"Set",dist.name+"Set",RooArgList(var),dist.hist)
    return var

#_____________________________data_________________________________
hists = {'Data' : 1}
add(hists)
data = TDistribution("data", 0, *options.stitched_input,**hists)
data.legentry = "Data("+lum+"pb^{-1})"
data.hist.SetMarkerStyle(8)
data.error = "e"
templates = [data]

var = init_var(data)
#___________________________EWK (single-top + WJets)____________________________
hists = {'SingleTop' : 1, 'Wbx' : 1, 'Wcx' : 1, 'Wqq' : 1}
add(hists)
ewk = TDistribution("ewk",var, *options.stitched_input,**hists)
ewk.legentry = "EWK/TOP"
ewk.hist.SetFillColor(4)
templates.append(ewk)

#_________________________________ZJets_________________________________________
hists = {'ZJets' : 1 }
add(hists)
dy = TDistribution("dy", var, *options.stitched_input,**hists)
dy.legentry = "Z+jets"
dy.hist.SetFillColor(2)
templates.append(dy)

#________________________________TTBar_________________________________________
hists = {'Top' : 1 }
add(hists)
top = TDistribution("top", var, *options.stitched_input,**hists)
top.legentry = "t #bar{t}"
top.hist.SetFillColor(206)
templates.append(top)

#____________________________qcd______________________________________
hists = {'QCD' : 1}
hists_tagged = {'QCD' : 1}
add(hists,doPretagged=True)
add(hists_tagged)

# Handle the QCD. First, get the input
if not options.qcd_norm:
    options.qcd_norm = options.stitched_input
if not options.qcd_signal:
    options.qcd_signal = options.stitched_input
if not options.qcd_shape:
    options.qcd_shape = options.stitched_input
qcd = TDistribution("qcd", var, *options.qcd_norm,**hists)
qcd_signal = TDistribution("qcd", var, *options.qcd_signal,**hists_tagged)
qcd_shape  = TDistribution("qcd", var, *options.qcd_shape,**hists)
qcd.legentry = "QCD"
qcd.hist.SetFillColor(6)
templates.append(qcd)


# =====================================================
# ================  FIT =======================
# =====================================================

if options.fit :
    var.setRange("fix",0,options.fixBin)
    sum_pdf = RooAddPdf("SumPdf","SumPdf",RooArgList(qcd.pdf,ewk.pdf,dy.pdf,top.pdf),RooArgList(qcd.N,ewk.N,dy.N,top.N))
    #ewk_constr = RooGaussian("ewk_constr","ewk_constr",ewk.SF,RooFit.RooConst(1.),RooFit.RooConst(0.05))
    zjet_constr = RooGaussian("zjet_constr","zjet_constr",dy.SF,RooFit.RooConst(1.),RooFit.RooConst(0.10))
    #all_constraints = [ zjet_constr.SF ]
    #tot_pdf = RooProdPdf("TotPdf","TotPdf",RooArgList(sum_pdf,ewk_constr))#,zjet_constr))
    tot_pdf = RooProdPdf("TotPdf","TotPdf",RooArgList(sum_pdf,zjet_constr))
    ##sum_pdf = RooAddPdf("SumPdf","SumPdf",RooArgList(qcd.pdf,ewk.pdf,wjets.pdf),RooArgList(qcd.N,ewk.N,wjets.N))
##	ewk_constr = RooGaussian("ewk_constr","ewk_constr",ewk.SF,RooFit.RooConst(1.),RooFit.RooConst(0.05))
##	tot_pdf = RooProdPdf("TotPdf","TotPdf",RooArgList(sum_pdf,ewk_constr))

    #args = [ data.set ]
    #args.extend(all_constraints)
    #args.extend([RooFit.Extended(kTRUE),RooFit.Range("fix"),RooFit.Save()])
    #print "got args %s" % args
    #r = tot_pdf.fitTo(*args)
    r = tot_pdf.fitTo(data.set,RooFit.Constrain(RooArgSet(dy.SF,ewk.SF)),RooFit.Extended(kTRUE),RooFit.Range("fix"),RooFit.Save())
    params = tot_pdf.getVariables()
    print "params are %s" % params
    params.Print("v")


    for idist in templates[1:]:
        print idist.hist.Print()
        print "This is the scale %s" % params.find(idist.name+"SF").getVal()
        idist.hist.Scale(params.find(idist.name+"SF").getVal())

sys.exit(0)
# =====================================================
# ================  PLOT =======================
# =====================================================
for idist in templates:
    print "   ", idist.legentry, "    %5.2f" % idist.hist.Integral()
    NBins = idist.hist.GetNbinsX()
    minB   = idist.hist.GetBinLowEdge(1)
    maxB   = idist.hist.GetBinLowEdge(NBins + 1)
    IMET = int((20 - minB)/(maxB - minB)*float(NBins))
    if idist.hist.GetName() == 'Data_'+options.var+"_"+nJet+"j_"+nTag+"t":
        relErr = 0 
    #elif idist.hist.GetName() != 'Data_'+options.var+"_"+nJet+"j_"+nTag+"t" and options.fit:
        #print  idist.hist.GetName()
    #    relErr = params.find(idist.name+"SF").getError()/params.find(idist.name+"SF").getVal()
    if options.verbose:
        metCalc = minB + (maxB-minB)*float(IMET)/float(NBins);
        print 'MET bin boundary = ', metCalc
hs = THStack("nEvents","nEvents")		
for idist in templates[1:] :
    hs.Add(idist.hist)

# draw
if data.hist.GetMaximum() > hs.GetMaximum() :
    hs.SetMaximum(data.hist.GetMaximum())
hs.Draw("HIST")
data.hist.Draw("esame")

xs = hs.GetXaxis()
xs.SetTitle(xtitle)
#xs.SetRangeUser(0.,options.nBin)
#xs.SetRangeUser(0.,200)
gPad.RedrawAxis()

#legend		
leg = TLegend(0.65,0.8,0.99,0.99)
leg.AddEntry(data.hist,data.legentry,"pl")
for idist in reversed(templates[1:]) :
    opt = ""
    if idist.hist.GetFillColor() :
        opt += "f"
    elif idist.hist.GetLineColor() != 1 :
        opt += "l"
    if idist.legentry != "" :
        leg.AddEntry(idist.hist,idist.legentry,opt)

Ysize = max(4, len(templates))
leg.SetY1(1-0.05*Ysize)
leg.SetBorderSize(1)
leg.SetFillColor(10)
leg.Draw()

c1.SetLogy(1)
if options.fit == 1:
    c1.SaveAs(options.outputDir+"/"+options.var+"_fit_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+"_log.gif")
    c1.SaveAs(options.outputDir+"/"+options.var+"_fit_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+"_log.eps")
    gROOT.ProcessLine(".!epstopdf "+options.outputDir+"/"+options.var+"_fit_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+"_log.eps")    
elif minJ==maxJ and minT==maxT:
    c1.SaveAs(options.outputDir+"/"+options.var+"_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+"_log.gif")
    c1.SaveAs(options.outputDir+"/"+options.var+"_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+"_log.eps")
    gROOT.ProcessLine(".!epstopdf "+options.outputDir+"/"+options.var+"_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+"_log.eps")
else:
    c1.SaveAs(options.pretagDir+"/"+options.var+"_"+subDir+"_"+minJ+"j_"+minT+"t_"+lum+"_log.gif")
    c1.SaveAs(options.pretagDir+"/"+options.var+"_"+subDir+"_"+minJ+"j_"+minT+"t_"+lum+"_log.eps")
    gROOT.ProcessLine(".!epstopdf "+options.pretagDir+"/"+options.var+"_"+subDir+"_"+minJ+"j_"+minT+"t_"+lum+"_log.eps")

c1.SetLogy(0)
if options.fit == 1:
    c1.SaveAs(options.outputDir+"/"+options.var+"_fit_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+".gif")
    c1.SaveAs(options.outputDir+"/"+options.var+"_fit_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+".eps")
    gROOT.ProcessLine(".!epstopdf "+options.outputDir+"/"+options.var+"_fit_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+".eps")
elif minJ==maxJ and minT==maxT:    
    c1.SaveAs(options.outputDir+"/"+options.var+"_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+".gif")
    c1.SaveAs(options.outputDir+"/"+options.var+"_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+".eps")
    gROOT.ProcessLine(".!epstopdf "+options.outputDir+"/"+options.var+"_"+subDir+"_"+nJet+"j_"+nTag+"t_"+lum+".eps")
else:
    c1.SaveAs(options.pretagDir+"/"+options.var+"_"+subDir+"_"+minJ+"j_"+minT+"t_"+lum+".gif")
    c1.SaveAs(options.pretagDir+"/"+options.var+"_"+subDir+"_"+minJ+"j_"+minT+"t_"+lum+".eps")
    gROOT.ProcessLine(".!epstopdf "+options.pretagDir+"/"+options.var+"_"+subDir+"_"+minJ+"j_"+minT+"t_"+lum+".eps")

