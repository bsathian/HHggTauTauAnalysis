import plottery.plottery as ply
import ROOT as r
import numpy as np
from collections import defaultdict

output_prefix = "/home/users/bsathian/public_html/HHGGTauTau_plots/prelim_plots/"

def plot_hists(bg_hists,sig_hist,param,legend_labels,extra_output_name = "",options = None):

    if not options:
        options = {
            "yaxis_log":True,
            "xaxis_range":[0,500],
        }

    options["output_name"] = output_prefix+extra_output_name+param+".pdf"
    options["lumi_value"] = 137.2
    options["legend_smart"] = False
    options["legend_percentageinbox"] = False
    options["cms_label"] = "Preliminary"
    options["xaxis_label"] = param



    ply.plot_hist(
            bgs = bg_hists,
            sigs = sig_hist,
            legend_labels = legend_labels,
            sig_labels = ["HH to #gamma#gamma#tau#tau signal"],
            options = options
            )



ggJetsFile = r.TFile("diphoton_80Inf.root")
DYJetsFile = r.TFile("dyjets.root")
GJetsFile = r.TFile("gjets.root")
ZGTo2LGFile = r.TFile("ZGToLLG.root")
signalFile = r.TFile("signal.root")

bg_files = [ggJetsFile,DYJetsFile,GJetsFile,ZGTo2LGFile]
legend_labels = ["#gamma#gamma + Jets","DY","#gamma + Jets","Z#gamma"]
sig_file = [signalFile]

bg_hists = defaultdict(list)
sig_hists = defaultdict(list)

for bg_file in bg_files:
    bg_hists["ggMass"].append(bg_file.Get("hggMass"))
    bg_hists["tautauMass"].append(bg_file.Get("htautauMass"))
    bg_hists["photon_pt_leading"].append(bg_file.Get("hLeadingPhotonPt"))
    bg_hists["photon_pt_trailing"].append(bg_file.Get("hTrailingPhotonPt"))
    bg_hists["tau_pt_leading"].append(bg_file.Get("hLeadingTauPt"))
    bg_hists["tau_pt_trailing"].append(bg_file.Get("hTrailingTauPt"))
    bg_hists["photon_eta_leading"].append(bg_file.Get("hLeadingPhotonEta"))
    bg_hists["photon_eta_trailing"].append(bg_file.Get("hTrailingPhotonEta"))
    bg_hists["tau_eta_leading"].append(bg_file.Get("hLeadingTauEta"))
    bg_hists["tau_eta_trailing"].append(bg_file.Get("hTrailingTauEta"))

sig_hists["ggMass"].append(sig_file[0].Get("hggMass"))
sig_hists["tautauMass"].append(sig_file[0].Get("htautauMass"))
sig_hists["photon_pt_leading"].append(sig_file[0].Get("hLeadingPhotonPt"))
sig_hists["photon_pt_trailing"].append(sig_file[0].Get("hTrailingPhotonPt"))
sig_hists["tau_pt_leading"].append(sig_file[0].Get("hLeadingTauPt"))
sig_hists["tau_pt_trailing"].append(sig_file[0].Get("hTrailingTauPt"))
sig_hists["photon_eta_leading"].append(sig_file[0].Get("hLeadingPhotonEta"))
sig_hists["photon_eta_trailing"].append(sig_file[0].Get("hTrailingPhotonEta"))
sig_hists["tau_eta_leading"].append(sig_file[0].Get("hLeadingTauEta"))
sig_hists["tau_eta_trailing"].append(sig_file[0].Get("hTrailingTauEta"))

rebins = {
    "ggMass" : 10,
    "photon_pt_leading":10,
    "photon_pt_trailing":10,
    "tau_pt_leading":10,
    "tau_pt_trailing":10,
    "photon_eta_leading":10,
    "photon_eta_trailing":10,
    "tau_eta_leading":10,
    "tau_eta_trailing":10,
}

special_options ={
        "ggMass":{"xaxis_range":[0,500], "yaxis_log":True, "yaxis_range":[1,1e6]},
        "tautauMass":{"xaxis_range":[0,500], "yaxis_log":True,"yaxis_range":[1,2e4]},
        "photon_pt_leading":{"xaxis_range":[0,200], "yaxis_log":True,"yaxis_range":[0.1,1e6]},
        "photon_pt_trailing":{"xaxis_range":[0,200], "yaxis_log":True,"yaxis_range":[0.1,1e6]},
        "tau_pt_leading":{"xaxis_range":[0,200], "yaxis_log":True,"yaxis_range":[0.1,1e6]},
        "tau_pt_trailing":{"xaxis_range":[0,200], "yaxis_log":True,"yaxis_range":[0.1,1e6]},
        "photon_eta_leading":{"xaxis_range":[-2.5,2.5], "yaxis_log":True,"yaxis_range":[0.1,1e6]},
        "photon_eta_trailing":{"xaxis_range":[-2.5,2.5], "yaxis_log":True,"yaxis_range":[0.1,1e6]},
        "tau_eta_leading":{"xaxis_range":[-2.5,2.5], "yaxis_log":True,"yaxis_range":[0.1,1e6]},
        "tau_eta_trailing":{"xaxis_range":[-2.5,2.5], "yaxis_log":True,"yaxis_range":[0.1,1e6]},

        }


for param in bg_hists.keys():
    print(param)
    print(bg_hists[param])
    if param in rebins:
        for bg_hist in bg_hists[param]:
            bg_hist.Rebin(rebins[param])
        sig_hists[param].Rebin(rebins[param])
    plot_hists(bg_hists[param],sig_hists[param],param,legend_labels,"",special_options[param])

