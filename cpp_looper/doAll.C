{
    gROOT->ProcessLine(".L HHGG.cc+");
    gROOT->ProcessLine(".L ScanChain.C+");
    TChain *ch = new TChain("Events");

    TString dataset_location("/hadoop/cms/store/user/legianni/ProjectMetis/GJets_HT-40To100_17____v1/*.root");
    TString output_name("GJets_40100.root");
    float xsec = 20790;
    //ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/HHggtautau_Era2018_private_prova_vera/*.root");
    ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/GJets_HT-40To100_17____v1/*.root");
    ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/GJets_HT-100To200_17____v1/*.root");
    ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/GJets_HT-200To400_17____v1/*.root");
    ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/GJets_HT-400To600_17____v1/*.root");
    ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/GJets_HT-600ToInf_17____v1/*.root");
    ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/DYJetsToLL_M-50_17____v1/*.root");
    ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/ZGToLLG_01J_17____v1/*.root");
    ch->Add("/hadoop/cms/store/user/legianni/ProjectMetis/DiPhotonJetsBox_MGG-80toInf_17____v1/*.root"); 
    ScanChain(ch, "output.root",1);
}
