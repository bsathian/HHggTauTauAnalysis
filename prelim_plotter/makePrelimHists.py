from __future__ import print_function
import ROOT as r
import numpy as np
import sys,os

ch = r.TChain("Events")
chRuns = r.TChain("Runs")
ch.Add(sys.argv[1]+"/*.root")
chRuns.Add(sys.argv[1]+"/*.root")
xsec = float(sys.argv[3]) * 1000
lumi = 137 #Scale to 3 years!

hists = []
paramList = {"ggMass":["hggMass","(500,0,500)","ggMass > 0"],"tautauMass":["htautauMass","(500,0,500)","tautauMass > 0"],"Photon_pt[0]":["hLeadingPhotonPt","(500,0,500)","passedHPhotons == 1"],"Photon_pt[1]":["hTrailingPhotonPt","(500,0,500)","passedHPhotons == 1"],"Photon_eta[0]":["hLeadingPhotonEta","(1000,-2.5,2.5)","passedHPhotons == 1"],"Photon_eta[1]":["hTrailingPhotonEta","(1000,-2.5,2.5)","passedHPhotons == 1"],"Tau_pt[tauHidx[0]]":["hLeadingTauPt","(500,0,500)","Category_pairs >= 1"], "Tau_pt[tauHidx[1]]":["hTrailingTauPt","(500,0,500)","Category_pairs >=\
    1"],"Tau_eta[tauHidx[0]]":["hLeadingTauEta","(1000,-2.5,2.5)","Category_pairs >= 1"], "Tau_eta[tauHidx[1]]":["hTrailingTauEta","(1000,-2.5,2.5)","Category_pairs >= 1"]}

#Get sum of weights first
sumOfWeights = 0.0
for i in range(chRuns.GetEntries()):
    chRuns.GetEntry(i)
    sumOfWeights += chRuns.genEventSumw

for param,others in paramList.items():
    histName = others[0]
    binning = others[1]
    if len(others) > 2:
        selection = "genWeight*({})".format(others[2])
    else:
        selection = "genWeight"
    scanString = "{}>>{}{}".format(param,histName,binning)
    print(scanString,selection)
    ch.Draw(scanString,selection)
    hists.append(r.gDirectory.Get(histName))
    hists[-1].Scale(xsec/sumOfWeights * lumi)

outputFile = r.TFile(sys.argv[2],"RECREATE")
for hist in hists:
    print(hist)
    hist.SetDirectory(outputFile)

outputFile.Write()



