#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TLine.h"

// Bring in your custom analysis framework headers
#include "../header/Binning.h"
#include "../header/JetHealthHistograms.h"
#include "../header/JetHealthPlotting.h"

const std::vector<double> ptCutsMap = {50.0, 100.0};
const std::vector<double> etaCutsMap = {0.0, 5.1};

void FormatSlidePad(bool isLogX, bool isLogY, bool is2D) {
    gPad->SetTopMargin(0.18); 
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(is2D ? 0.12 : 0.15);
    gPad->SetRightMargin(is2D ? 0.16 : 0.05); 
    if (isLogX) gPad->SetLogx(true);
    if (isLogY) gPad->SetLogy(true);
}

void DrawSlideText(int padNum, bool is2D, const TString& mainTitle, double ptCut, const TString& jetAlgo, const TString& extraLabel = "") {
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    
    double leftX  = is2D ? 0.12 : 0.15;
    double rightX = is2D ? 0.84 : 0.95;
    
    tex.SetTextAlign(22); 
    tex.SetTextSize(0.055); 
    tex.DrawLatex(0.5, 0.87, mainTitle);
    
    if (padNum == 1) {
        tex.SetTextAlign(11); 
        tex.DrawLatex(leftX, 0.94, "#bf{CMS} #it{Internal}");
        tex.SetTextAlign(31); 
        tex.DrawLatex(rightX, 0.94, Form("#bf{p_{T} > %.0f GeV/c}", ptCut));
    } else if (padNum == 2) {
        if (!extraLabel.IsNull()) {
            tex.SetTextAlign(31); 
            tex.DrawLatex(0.5, 0.94, extraLabel); 
        }
    } else if (padNum == 3) {
        tex.SetTextAlign(31); 
        tex.DrawLatex(rightX, 0.94, Form("#bf{%s}", jetAlgo.Data())); 
    }
}

void PlotJetKinematics(TString file1Path = "JetHealth_404350", TString label1="Run 404350", 
                            TString file2Path = "JetHealth_404469_RPmany", TString label2="Run 404469") {

    if (false){
        file1Path = "JetHealth_2026_MC"; label1="2026 MC";
        file2Path = "JetHealth_2026_MC_maskFPIX"; label2="2026 MC FPIX-masked";
    }   
    else if (false){
        file1Path = "JetHealth_2026_MC"; label1="PF Jets";
        file2Path = "JetHealth_2026MC_ak4PF"; label2="Unsubtracted PF Jets";
    }  
    
    //file2Path = "JetHealth_404350_CaloJets_try2"; label2="Calo Jets"; //label2="Run404350 akPu4Calo";
    
    //file2Path = "JetHealth_2025_Data_new"; label2="2025 Data";
      
    file2Path = "JetHealth_404350_PromptRECO"; label2="Run 404350 Prompt RECO";

    gStyle->SetOptStat(0);
    
    TString outDir = file1Path + "_vs_" + file2Path;
    gSystem->mkdir(outDir, true);

    TFile* f1 = TFile::Open(file1Path+".root", "READ");
    TFile* f2 = TFile::Open(file2Path+".root", "READ");
    
    if (!f1 || f1->IsZombie()) { std::cerr << "Cannot open file 1" << std::endl; return; }
    if (!f2 || f2->IsZombie()) { std::cerr << "Cannot open file 2" << std::endl; return; }

    THnSparseF* hnKin1 = (THnSparseF*)f1->Get("hjetkin");
    THnSparseF* hnKin2 = (THnSparseF*)f2->Get("hjetkin");

    BinningStruct bins(50.0);
    struct EtaRange { double lo; double hi; TString name; TString label; };
    std::vector<EtaRange> etaRanges;
    for (size_t i = 0; i < etaCutsMap.size() - 1; ++i) {
        etaRanges.push_back({etaCutsMap[i], etaCutsMap[i+1], Form("eta_%.1f_%.1f", etaCutsMap[i], etaCutsMap[i+1]), 
                             etaCutsMap[i] == 0 ? Form("|#eta| < %.1f", etaCutsMap[i+1]) : Form("%.1f < |#eta| < %.1f", etaCutsMap[i], etaCutsMap[i+1])});
    }
    etaRanges.push_back({0.0, 5.1, "eta_inclusive", "Inclusive #eta"});

    auto StyleTH1Slide = [](TH1D* h, Color_t color, const TString& xTitle) {
        StyleTH1(h, color);
        h->SetLineWidth(1); h->SetMarkerSize(0.01);
        h->SetTitle(""); h->GetYaxis()->SetTitle(""); 
        h->GetXaxis()->SetTitle(xTitle); h->GetXaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleSize(0.05); h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05); h->GetXaxis()->SetTitleOffset(1.1); 
    };

    auto GetKinVectors = [&](int axis, double ptCut, double etaLo, double etaHi, const TString& xTitle) {
        std::vector<TH1D*> v1, v2, vRat;
        for (const auto& hb : bins.hiBins) {
            TH1D *h1, *h2;
            if (etaLo == 0.0) {
                h1 = ProjectTHn1D(hnKin1, axis, {{0, ptCut, 1000.0}, {1, -etaHi, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_1");
                h2 = ProjectTHn1D(hnKin2, axis, {{0, ptCut, 1000.0}, {1, -etaHi, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_2");
            } else {
                TH1D* h1_p = ProjectTHn1D(hnKin1, axis, {{0, ptCut, 1000.0}, {1, etaLo, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_1p");
                TH1D* h1_n = ProjectTHn1D(hnKin1, axis, {{0, ptCut, 1000.0}, {1, -etaHi, -etaLo}, {3, (double)hb.lo, (double)hb.hi}}, "_1n");
                h1_p->Add(h1_n); h1 = h1_p; delete h1_n;
                
                TH1D* h2_p = ProjectTHn1D(hnKin2, axis, {{0, ptCut, 1000.0}, {1, etaLo, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_2p");
                TH1D* h2_n = ProjectTHn1D(hnKin2, axis, {{0, ptCut, 1000.0}, {1, -etaHi, -etaLo}, {3, (double)hb.lo, (double)hb.hi}}, "_2n");
                h2_p->Add(h2_n); h2 = h2_p; delete h2_n;
            }
            NormalizeTH1(h1); NormalizeTH1(h2);
            StyleTH1Slide(h1, hb.color, xTitle); StyleTH1Slide(h2, hb.color, xTitle);
            
            TH1D* hRatioNum = (TH1D*)h1->Clone(); TH1D* hRatioDen = (TH1D*)h2->Clone();
            if (axis == 1 || axis == 2) { hRatioNum->Rebin(4); hRatioDen->Rebin(4); }
            TH1D* hRatio = (TH1D*)hRatioNum->Clone(); hRatio->Divide(hRatioDen);
            
            StyleTH1Slide(hRatio, hb.color, xTitle);
            hRatio->GetYaxis()->CenterTitle(true); hRatio->GetYaxis()->SetTitleSize(0.05); 
            hRatio->GetYaxis()->SetTitleOffset(0.85); hRatio->GetYaxis()->SetNdivisions(505);
            
            v1.push_back(h1); v2.push_back(h2); vRat.push_back(hRatio);
            delete hRatioNum; delete hRatioDen;
        }
        return std::make_tuple(v1, v2, vRat);
    };

    auto DrawAndSaveSlide = [&](const TString& cName, const TString& title, std::vector<TH1D*>& v1, std::vector<TH1D*>& v2, std::vector<TH1D*>& vRat, 
                                double ptCut, const TString& etaLabel, bool isLogY, bool isLogX) {
        TCanvas* cSlide = new TCanvas(cName, title, 3600, 1200); cSlide->Divide(3, 1);
        TLegend* leg = new TLegend(0.40, 0.55, 0.90, 0.80); leg->SetBorderSize(0); leg->SetTextSize(0.045);
        for(size_t i=0; i<bins.hiBins.size(); ++i) leg->AddEntry(v1[i], bins.hiBins[i].title, "lp");
        
        double ymax_main = 0.0, ymax_ratio = -999.0, ymin_ratio = 999.0;
        for(auto h : v1) ymax_main = std::max(ymax_main, h->GetMaximum());
        for(auto h : v2) ymax_main = std::max(ymax_main, h->GetMaximum());
        for(auto h : vRat) {
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                if (h->GetBinContent(b) > 0) {
                    ymax_ratio = std::max(ymax_ratio, h->GetBinContent(b));
                    ymin_ratio = std::min(ymin_ratio, h->GetBinContent(b));
                }
            }
        }
        
        cSlide->cd(1); FormatSlidePad(isLogX, isLogY, false);
        for(size_t i=0; i<v1.size(); ++i) {
            if(i==0) {
                v1[0]->GetYaxis()->SetRangeUser(isLogY ? 1e-6 : 0.0, ymax_main * (isLogY ? 5.0 : 1.35));
                if(isLogX) v1[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                v1[0]->Draw("EP");
            } else v1[i]->Draw("EP SAME");
        }
        DrawSlideText(1, false, label1, ptCut, "", etaLabel);
        
        cSlide->cd(2); FormatSlidePad(isLogX, isLogY, false);
        for(size_t i=0; i<v2.size(); ++i) {
            if(i==0) {
                v2[0]->GetYaxis()->SetRangeUser(isLogY ? 1e-6 : 0.0, ymax_main * (isLogY ? 5.0 : 1.35));
                if(isLogX) v2[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                v2[0]->Draw("EP");
            } else v2[i]->Draw("EP SAME");
        }
        DrawSlideText(2, false, label2, ptCut, "", etaLabel);
        
        cSlide->cd(3); FormatSlidePad(isLogX, false, false);
        for(size_t i=0; i<vRat.size(); ++i) {
            if(i==0) {
                vRat[0]->GetYaxis()->SetRangeUser(ymin_ratio, ymax_ratio);
                if(isLogX) vRat[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                vRat[0]->Draw("EP");
                double lineXm = isLogX ? ptCut : v1[0]->GetXaxis()->GetXmin();
                double lineXM = isLogX ? 500.0 : v1[0]->GetXaxis()->GetXmax();
                TLine* line = new TLine(lineXm, 1, lineXM, 1);
                line->SetLineStyle(2); line->SetLineColor(kBlack); line->Draw();
            } else vRat[i]->Draw("EP SAME");
        }
        DrawSlideText(3, false, Form("Ratio to %s", label2.Data()), ptCut, "", etaLabel);
        
        cSlide->SaveAs(outDir + "/" + cName + ".png"); 
        delete cSlide;
        for(auto h: v1) delete h; for(auto h: v2) delete h; for(auto h: vRat) delete h;
    }; 

    // TASK 1: pT spectra for 1st centrality bin (w/ Integrals)
    {
        TCanvas* cPtCent0 = new TCanvas("cPtCent0", "pT 1st Centrality Bin", 2000, 2000);
        FormatSlidePad(true, true, false);

        TH1D* hPt1 = (TH1D*) hnKin1->Projection(0);
        TH1D* hPt2 = (TH1D*) hnKin2->Projection(0);
        StyleTH1Slide(hPt1, kBlue, "p_{T} (GeV/c)");
        StyleTH1Slide(hPt2, kRed, "p_{T} (GeV/c)");

        double int1 = hPt1->Integral(hPt1->FindBin(50.0), hPt1->GetNbinsX());
        double int2 = hPt2->Integral(hPt2->FindBin(50.0), hPt2->GetNbinsX());

        hPt1->GetYaxis()->SetRangeUser(1, std::max(hPt1->GetMaximum(), hPt2->GetMaximum()) * 200.0);
        hPt1->GetXaxis()->SetRangeUser(1, 500.0);
        
        hPt1->Draw("EP");
        hPt2->Draw("EP SAME");

        TLegend* legPt = new TLegend(0.2, 0.70, 0.8, 0.8);
        legPt->SetBorderSize(0);
        legPt->SetTextSize(0.04);
        TString sub1 = label1; sub1.Remove(0,10);
        TString sub2 = label2; sub2.Remove(0,10);
        legPt->AddEntry(hPt1, Form("%s (Int p_{T}>50: %.2e)", sub1.Data(), int1), "lp");
        legPt->AddEntry(hPt2, Form("%s (Int p_{T}>50: %.2e)", sub2.Data(), int2), "lp");
        legPt->Draw("SAME");
        
        DrawSlideText(3, false, Form("hiBin %.0f-%.0f", bins.hiBins[0].lo, bins.hiBins[0].hi), 50.0, "Run 404350", Form("hiBin %.0f-%.0f", bins.hiBins[0].lo, bins.hiBins[0].hi));
        
        cPtCent0->SaveAs(outDir + "/cPtSpectra_1stCentBin.png");
        delete cPtCent0; delete hPt1; delete hPt2;
    }

    // SLIDES: KINEMATICS
    for (const auto& er : etaRanges) {
        auto [v1, v2, vR] = GetKinVectors(0, bins.ptmin, er.lo, er.hi, "p_{T} (GeV/c)");
        DrawAndSaveSlide(Form("cKinSlide_pt_%s", er.name.Data()), "pT", v1, v2, vR, bins.ptmin, er.label, true, true);
    }
    
    for (double ptC : ptCutsMap) {
        auto [v1, v2, vR] = GetKinVectors(1, ptC, 0.0, 5.1, "#eta");
        DrawAndSaveSlide(Form("cKinSlide_eta_pt%.0f", ptC), "Eta", v1, v2, vR, ptC, "Inclusive #eta", false, false);
    }
    
    {
        auto [v1, v2, vR] = GetKinVectors(2, bins.ptmin, 0.0, 5.1, "#phi (rad)");
        DrawAndSaveSlide("cKinSlide_phi_inclusive", "Phi", v1, v2, vR, bins.ptmin, "Inclusive #eta", false, false);
    }

    std::cout << "Kinematic plots successfully saved to " << outDir << "/" << std::endl;
}