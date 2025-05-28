// compareHists.C
// Simple ROOT macro to compare histograms, given their names and the ROOT files they are from.

template <typename T>
T* GetTH(TString fName, TString histName)
{
    TFile *f = new TFile(fName);
    return dynamic_cast<T *>(f->Get(histName));
}

void doComparisonTH1(TString fName1, TString fName2, std::vector<TString> vHistNames)
{
    TCanvas *c;
    for (auto &histName : vHistNames)
    {
        TH1 *h1 = GetTH<TH1>(fName1, histName);
        TH1 *h2 = GetTH<TH1>(fName2, histName);
        c = new TCanvas(Form("c_TH1_%s", histName.Data()), Form("c_TH1_%s", histName.Data()));
        c->cd();
        h1->GetYaxis()->SetTitle("counts / total counts")
        h1->Draw();
        h2->Draw("SAME");
    }
}

void doComparisonTH2(TString fName1, TString fName2, std::vector<TString> vHistNames)
{
    TCanvas *c;
    for (auto &histName : vHistNames)
    {
        TH2 *h1 = GetTH<TH2>(fName1, histName);
        TH2 *h2 = GetTH<TH2>(fName2, histName);
        c = new TCanvas(Form("c_TH2_%s", histName.Data()), Form("c_TH2_%s", histName.Data()));
        c->Divide(2,1);
        c->cd(1);
        h1->Draw();
        c->cd(2);
        h2->Draw();
    }
}

int compareHists()
{
    std::vector<TString> vHistNamesTH1 = {"hMULT", "hNMPIs", "hAvgPtAll", "hAvgPtSoft", "hSphAll", "hSphSoft"};
    std::vector<TString> vHistNamesTH2 = {"hNMPIs_hMult", "hMult_hAvgPt_All", "hMult_hAvgPt_Soft", "hMult_hSph_All", "hMult_hSph_Soft", "hAvgPt_hSph_All", "hAvgPt_hSph_Soft"};
    doComparisonTH1("BplusBminus_minbias_1e5.root", "BplusBminus_hardbbbar_1e5.root", vHistNamesTH1);
    doComparisonTH2("BplusBminus_minbias_1e5.root", "BplusBminus_hardbbbar_1e5.root", vHistNamesTH2);
    return 0;
}