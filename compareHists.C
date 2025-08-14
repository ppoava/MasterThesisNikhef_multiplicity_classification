// compareHists.C
// Simple ROOT macro to compare histograms, given their names and the ROOT files they are from.

template <typename T>
T* GetTH(TString fName, TString histName)
{
    TFile *f = new TFile(fName);
    return dynamic_cast<T *>(f->Get(histName));
}

void doComparisonTH1(TString fName1, TString fLabel1, TString fName2, TString fLabel2, std::vector<TString> vHistNames)
{
    TCanvas *c;
    TLegend *leg;

    TCanvas cBegDummy = TCanvas();
    cBegDummy.SaveAs(Form("Plots/all_plots_%s_%s.pdf(", fLabel1.Data(), fLabel2.Data()));
    for (auto &histName : vHistNames)
    {
        TH1 *h1 = GetTH<TH1>(fName1, histName);
        TH1 *h2 = GetTH<TH1>(fName2, histName);
        h1->Scale(1 / h1->GetEntries());
        h2->Scale(1 / h2->GetEntries());
        c = new TCanvas(Form("c_TH1_%s_%s_%s", histName.Data(), fLabel1.Data(), fLabel2.Data()), Form("c_TH1_%s_%s_%s", histName.Data(), fLabel1.Data(), fLabel2.Data()));
        c->cd();
        h1->GetYaxis()->SetTitle("counts / total counts");
        h1->Draw();
        h2->Draw("SAME");
        h2->SetLineColor(kRed);
        leg = new TLegend;
        leg->AddEntry(h1, fLabel1, "l");
        leg->AddEntry(h2, fLabel2, "l");
        leg->Draw();
        c->SaveAs(Form("Plots/all_plots_%s_%s.pdf", fLabel1.Data(), fLabel2.Data()));
        // c->SaveAs(Form("Plots/cTH1_%s_%s_%s.png", histName.Data(), fLabel1.Data(), fLabel2.Data()));
        // c->SaveAs(Form("Plots/cTH1_%s_%s_%s.pdf", histName.Data(), fLabel1.Data(), fLabel2.Data()));
    }
}

void doComparisonTH2(TString fName1, TString fLabel1, TString fName2, TString fLabel2, std::vector<TString> vHistNames)
{
    TCanvas *c;
    for (auto &histName : vHistNames)
    {
        TH2 *h1 = GetTH<TH2>(fName1, histName);
        TH2 *h2 = GetTH<TH2>(fName2, histName);
        h1->Scale(1 / h1->GetEntries());
        h2->Scale(1 / h2->GetEntries());
        c = new TCanvas(Form("c_TH2_%s_%s_%s", histName.Data(), fLabel1.Data(), fLabel2.Data()), Form("c_TH2_%s_%s_%s", histName.Data(), fLabel1.Data(), fLabel2.Data()));
        c->Divide(2,1);
        c->cd(1);
        h1->SetTitle(Form("%s_%s", fLabel1.Data(), h1->GetTitle()));
        h1->Draw();
        c->cd(2);
        h2->SetTitle(Form("%s_%s", fLabel2.Data(), h2->GetTitle()));
        h2->Draw();
        c->SaveAs(Form("Plots/all_plots_%s_%s.pdf", fLabel1.Data(), fLabel2.Data()));
        // c->SaveAs(Form("Plots/cTH2_%s_%s_%s.png", histName.Data(), fLabel1.Data(), fLabel2.Data()));
        // c->SaveAs(Form("Plots/cTH2_%s_%s_%s.pdf", histName.Data(), fLabel1.Data(), fLabel2.Data()));
    }
    TCanvas cEndDummy = TCanvas();
    cEndDummy.SaveAs(Form("Plots/all_plots_%s_%s.pdf)", fLabel1.Data(), fLabel2.Data()));
}

int compareHists()
{
    std::vector<TString> vHistNamesTH1 = {"hMultAll", "hNMPIs", "hNJets", "hAvgPtAll", "hAvgPtSoft", "hSphAll", "hSphSoft"};
    std::vector<TString> vHistNamesTH2 = {"hNMPIs_hMult", "hMult_hAvgPt_All", "hMult_hAvgPt_Soft", "hMult_hSph_All", "hMult_hSph_Soft", "hMult_hNJets", "hNJets_hSph_All", "hNJets_hSph_Soft", "hNJets_hAvgPt_All", "hNJets_hAvgPt_Soft", "hAvgPt_hSph_All", "hAvgPt_hSph_Soft"};

    doComparisonTH1("BplusBminus_minbias_bb_1e6.root", "minbias-bb", "BplusBminus_hardbbbar_1e5.root", "hardbbbar", vHistNamesTH1);
    doComparisonTH2("BplusBminus_minbias_bb_1e6.root", "minbias-bb", "BplusBminus_hardbbbar_1e5.root", "hardbbbar", vHistNamesTH2);

    doComparisonTH1("BplusBminus_minbias_bb_JUNCTIONS_1e4.root", "minbias-bb-junc", "BplusBminus_hardbbbar_JUNCTIONS_1e4.root", "hardbbbar-junc", vHistNamesTH1);
    doComparisonTH2("BplusBminus_minbias_bb_JUNCTIONS_1e4.root", "minbias-bb-junc", "BplusBminus_hardbbbar_JUNCTIONS_1e4.root", "hardbbbar-junc", vHistNamesTH2);

    doComparisonTH1("DplusDminus_minbias_cc_1e6.root", "minbias-cc", "DplusDminus_hardccbar_1e5.root", "hardccbar", vHistNamesTH1);
    doComparisonTH2("DplusDminus_minbias_cc_1e6.root", "minbias-cc", "DplusDminus_hardccbar_1e5.root", "hardccbar", vHistNamesTH2);

    doComparisonTH1("DplusDminus_minbias_cc_JUNCTIONS_1e4.root", "minbias-cc-junc", "DplusDminus_hardccbar_JUNCTIONS_1e4.root", "hardccbar-junc", vHistNamesTH1);
    doComparisonTH2("DplusDminus_minbias_cc_JUNCTIONS_1e4.root", "minbias-cc-junc", "DplusDminus_hardccbar_JUNCTIONS_1e4.root", "hardccbar-junc", vHistNamesTH2);

    // now the same for hardQCD (including bb and cc bar)
    doComparisonTH1("BplusBminus_minbias_bb_1e6.root", "minbias-bb", "BplusBminus_hardQCD_1e7.root", "hardQCD", vHistNamesTH1);
    doComparisonTH2("BplusBminus_minbias_bb_1e6.root", "minbias-bb", "BplusBminus_hardQCD_1e7.root", "hardQCD", vHistNamesTH2);

    doComparisonTH1("BplusBminus_minbias_bb_JUNCTIONS_1e4.root", "minbias-bb-junc", "BplusBminus_hardQCD_JUNCTIONS_1e4.root", "hardQCD-junc", vHistNamesTH1);
    doComparisonTH2("BplusBminus_minbias_bb_JUNCTIONS_1e4.root", "minbias-bb-junc", "BplusBminus_hardQCD_JUNCTIONS_1e4.root", "hardQCD-junc", vHistNamesTH2);

    doComparisonTH1("DplusDminus_minbias_cc_1e6.root", "minbias-cc", "DplusDminus_hardQCD_1e7.root", "hardQCD", vHistNamesTH1);
    doComparisonTH2("DplusDminus_minbias_cc_1e6.root", "minbias-cc", "DplusDminus_hardQCD_1e7.root", "hardQCD", vHistNamesTH2);

    doComparisonTH1("DplusDminus_minbias_cc_JUNCTIONS_1e4.root", "minbias-cc-junc", "DplusDminus_hardQCD_JUNCTIONS_1e4.root", "hardQCD-junc", vHistNamesTH1);
    doComparisonTH2("DplusDminus_minbias_cc_JUNCTIONS_1e4.root", "minbias-cc-junc", "DplusDminus_hardQCD_JUNCTIONS_1e4.root", "hardQCD-junc", vHistNamesTH2);

    // beauty vs charm
    // okay now?
    doComparisonTH1("BplusBminus_hardbbbar_1e5.root", "bb-only", "DplusDminus_hardccbar_1e5.root", "cc-only", vHistNamesTH1);
    doComparisonTH2("BplusBminus_hardbbbar_1e5.root", "bb-only", "DplusDminus_hardccbar_1e5.root", "cc-only", vHistNamesTH2);

    doComparisonTH1("BplusBminus_hardQCD_1e7.root", "hardQCD B+B-", "DplusDminus_hardQCD_1e7.root", "hardQCD D+D-", vHistNamesTH1);
    doComparisonTH2("BplusBminus_hardQCD_1e7.root", "hardQCD B+B-", "DplusDminus_hardQCD_1e7.root", "hardQCd D+D-", vHistNamesTH2);

    return 0;
}