#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeFormula.h"
#include "HHGG.h"
#include "tqdm.h"

using namespace std;

int ScanChain(TChain *ch, TString outputFileName, float xsec) {

    TH1F* hggMass = new TH1F("hggMass","mgg",500,0,500);
    TH1F* htautauMass = new TH1F("htautauMass","mtautau",500,0,500);
    TH1F* hLeadingPhotonPt = new TH1F("hLeadingPhotonPt","leadingGammaPt",500,0,500);
    TH1F* hTrailingPhotonPt = new TH1F("hTrailingPhotonPt","trailingGammaPt",500,0,500);
    TH1F* hLeadingPhotonEta = new TH1F("hLeadingPhotonEta","leadingGammaEta",1000,-2.5,2.5);
    TH1F* hTrailingPhotonEta = new TH1F("hTrailingPhotonEta","trailingGammaEta",1000,-2.5,2.5);
    TH1F* hLeadingTauPt = new TH1F("hLeadingTauPt","leadingTauPt",500,0,500);
    TH1F* hTrailingTauPt = new TH1F("hTrailingTauPt","trailingTauPt",500,0,500);
    TH1F* hLeadingTauEta = new TH1F("hLeadingTauEta","leadingTauEta",1000,-2.5,2.5);
    TH1F* hTrailingTauEta = new TH1F("hTrailingTauEta","trailingTauEta",1000,-2.5,2.5);

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();
    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);
    tqdm bar;
    float sumOfWeights = 0;

    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile *file = TFile::Open( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        TString filename(currentFile->GetTitle());
        tree->SetCacheSize(128*1024*1024);
        tree->SetCacheLearnEntries(10);
        hhgg.SetYear(2018);
        hhgg.Init(tree);

        // Can also make custom formulas to eval (but dunno how slow this is)

        for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {
            hhgg.GetEntry(event);
            tree->LoadTree(event);
            nEventsTotal++;
            bar.progress(nEventsTotal, nEventsChain);

            // Analysis code

            // auto jets = hhggtautau::GetVLV("Jet_p4"); if (jets.size() > 0) h1->Fill(jets[0].pt());
            // auto jets = hhggtautau::Jet_p4(); if (jets.size() > 0) h1->Fill(jets[0].pt());
            // auto jetpts = hhggtautau::Jet_pt(); if (jetpts.size() > 0) h1->Fill(jetpts[0]);
            // for (auto& s : {"MET_pt", "CaloMET_pt", "GenMET_pt", "RawMET_pt", "TkMET_pt"}) h1->Fill(hhggtautau::GetF(s));
            // bool highmet = formula.EvalInstance();
            sumOfWeights += hhgg.genWeight();

            if(hhgg.ggMass() > 0)
            {
                hggMass->Fill(hhgg.ggMass(), hhgg.genWeight());
            }
            if(hhgg.tautauMass() > 0)
            {
                htautauMass->Fill(hhgg.tautauMass(), hhgg.genWeight());
            }
            if(hhgg.passedHPhotons() == 1)
            {
                hLeadingPhotonPt->Fill(hhgg.Photon_pt()[hhgg.gHidx()[0]], hhgg.genWeight());
                hLeadingPhotonEta->Fill(hhgg.Photon_eta()[hhgg.gHidx()[0]], hhgg.genWeight());
                hTrailingPhotonPt->Fill(hhgg.Photon_pt()[hhgg.gHidx()[1]], hhgg.genWeight());
                hTrailingPhotonEta->Fill(hhgg.Photon_eta()[hhgg.gHidx()[1]], hhgg.genWeight());
            }
/*            if(hhgg.Category_pairs() >= 1)
            {
                hLeadingTauPt->Fill(hhgg.Tau_pt()[hhgg.tauHidx()[0]]);
                hTrailingTauPt->Fill(hhgg.Tau_pt()[hhgg.tauHidx()[1]]);
                hLeadingTauEta->Fill(hhgg.Tau_eta()[hhgg.tauHidx()[0]]);
                hTrailingTauEta->Fill(hhgg.Tau_eta()[hhgg.tauHidx()[1]]);
            }*/

        } // Event loop
        delete file;
    } // File loop
    bar.finish();

    TFile* f1 = new TFile(outputFileName, "RECREATE");
    hggMass->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hggMass->Write();
    htautauMass->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    htautauMass->Write();
    hLeadingPhotonPt->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hLeadingPhotonPt->Write();
    hTrailingPhotonPt->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hTrailingPhotonPt->Write();
    hLeadingPhotonEta->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hLeadingPhotonEta->Write();
    hTrailingPhotonEta->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hTrailingPhotonEta->Write();
    hLeadingTauPt->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hLeadingTauPt->Write();
    hTrailingTauPt->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hTrailingTauPt->Write();
    hLeadingTauEta->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hLeadingTauEta->Write();
    hTrailingTauEta->Scale(1.0/sumOfWeights * xsec * 1000 * 137.2);
    hTrailingTauEta->Write();
    f1->Write();
    f1->Close();
    return 0;
}
