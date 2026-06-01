#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <tuple>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TLine.h"
#include "TProfile.h" 

// Bring in your custom analysis framework headers
#include "../header/Binning.h"
#include "../header/JetHealthHistograms.h"
#include "../header/JetHealthPlotting.h"

// ==============================================================================
// 1. CONSTANTS & CONFIGURATION
// ==============================================================================

// Dynamic Cut Maps
const std::vector<double> ptCutsMap = {50.0, 100.0};
const std::vector<double> etaCutsMap = {0.0, 5.1};

// ==============================================================================
// HELPER FUNCTIONS (Formatting & Organization)
// ==============================================================================

double GetNevents(TFile* f) {
    if (!f || f->IsZombie()) return 1.0;
    TH1F* hvz = (TH1F*)f->Get("hvz");
    if (!hvz) {
        std::cerr << "WARNING: Could not find 'hvz' in " << f->GetName() << ". Normalizing by 1." << std::endl;
        return 1.0;
    }
    return hvz->Integral();
}

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

// ==============================================================================
// MAIN PLOTTING MACRO
// ==============================================================================

void PlotJetHealthComparison(TString file1Path = "JetHealth_404350", TString label1="PF Jets", 
                            TString file2Path = "JetHealth_404350_ak4PF", TString label2="Unsubtracted PF Jets") {

    if (false){
        file1Path = "JetHealth_2026_MC"; label1="2026 MC";
        file2Path = "JetHealth_2026_MC_maskFPIX"; label2="2026 MC FPIX-masked";
    }   
    else if (true){
        file1Path = "JetHealth_2026_MC"; label1="PF Jets";
        file2Path = "JetHealth_2026MC_ak4PF"; label2="Unsubtracted PF Jets";
    }  
    else if (false){
        file2Path = "JetHealth_404350_CaloJets_try2"; label2="Calo Jets"; //label2="Run404350 akPu4Calo";
    }    
    else if (false){
        file2Path = "JetHealth_2025_Data_new"; label2="2025 Data";
    }                   
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRedBlue);
    
    TString outDir = file1Path + "_vs_" + file2Path;
    gSystem->mkdir(outDir, true);

    TFile* f1 = TFile::Open(file1Path+".root", "READ");
    TFile* f2 = TFile::Open(file2Path+".root", "READ");
    
    if (!f1 || f1->IsZombie()) { std::cerr << "Cannot open file 1: " << file1Path << std::endl; return; }
    if (!f2 || f2->IsZombie()) { std::cerr << "Cannot open file 2: " << file2Path << std::endl; return; }
    
    double nEvt1 = 1.; //GetNevents(f1);
    double nEvt2 = 1.; //GetNevents(f2);

    THnSparseF* hnKin1 = (THnSparseF*)f1->Get("hjetkin");
    THnSparseF* hnKin2 = (THnSparseF*)f2->Get("hjetkin");
    THnSparseF* hnPF1  = (THnSparseF*)f1->Get("hjetpf");
    THnSparseF* hnPF2  = (THnSparseF*)f2->Get("hjetpf");

    BinningStruct bins(50.0);

    // Build Eta Ranges directly from etaCutsMap
    struct EtaRange { double lo; double hi; TString name; TString label; };
    std::vector<EtaRange> etaRanges;
    for (size_t i = 0; i < etaCutsMap.size() - 1; ++i) {
        etaRanges.push_back({etaCutsMap[i], etaCutsMap[i+1], Form("eta_%.1f_%.1f", etaCutsMap[i], etaCutsMap[i+1]), 
                             etaCutsMap[i] == 0 ? Form("|#eta| < %.1f", etaCutsMap[i+1]) : Form("%.1f < |#eta| < %.1f", etaCutsMap[i], etaCutsMap[i+1])});
    }
    etaRanges.push_back({0.0, 5.1, "eta_inclusive", "Inclusive #eta"});

    // ==============================================================================
    // CANVAS 2: 2D ETA-PHI MAPS
    // ==============================================================================
    for (double ptC : ptCutsMap) {
        for (const auto& hb : bins.hiBins) { 
            TString cName = Form("cMap_pt%.0f_hb%.0f_%.0f", ptC, hb.lo, hb.hi);
            TCanvas* cMap = new TCanvas(cName, "Eta-Phi Map", 1600, 800); cMap->Divide(2, 1);
            
            cMap->cd(1); FormatSlidePad(false, false, true);
            TH2D* h1 = ProjectTHn2D(hnKin1, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "map1");
            //h1->Scale(1.0 / nEvt1); 
            h1->SetTitle(";#eta;#phi (rad)");
            h1->GetXaxis()->SetRangeUser(-2.4999,2.4999);
            h1->GetXaxis()->SetTitleOffset(0.85); h1->GetYaxis()->SetTitleOffset(0.85); h1->Draw("colz");
            DrawSlideText(1, true, label1, ptC, "Run 404350");
            
            cMap->cd(2); FormatSlidePad(false, false, true);
            TH2D* h2 = ProjectTHn2D(hnKin2, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "map2");
            //h2->Scale(1.0 / nEvt2); 
            h2->SetTitle(";#eta;#phi (rad)");
            h2->GetXaxis()->SetRangeUser(-2.4999,2.4999);
            h2->GetXaxis()->SetTitleOffset(0.85); h2->GetYaxis()->SetTitleOffset(0.85); h2->Draw("colz");
            DrawSlideText(2, true, label2, ptC, "Run 404350", Form("#bf{MC 2026 hiBin %.0f-%.0f}", hb.lo, hb.hi));
            
            cMap->SaveAs(outDir + "/" + cName + ".png");
            delete h1; delete h2; delete cMap;
        }
    }

    // ------------------------------------------------------------------------------
    // REUSABLE LAMBDAS FOR 1D PROJECTIONS & DRAWING
    // ------------------------------------------------------------------------------
    
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

    auto GetPFVectors = [&](int pf, double etaLo, double etaHi, const TString& xTitle) {
        std::vector<TH1D*> v1, v2, vRat;
        for (const auto& hb : bins.hiBins) {
            TH1D* h1 = ProjectTHn1D(hnPF1, 0, {{1, (double)pf, (double)pf + 1.0}, {2, etaLo, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_1");
            TH1D* h2 = ProjectTHn1D(hnPF2, 0, {{1, (double)pf, (double)pf + 1.0}, {2, etaLo, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_2");
            
            NormalizeTH1(h1); NormalizeTH1(h2);
            StyleTH1Slide(h1, hb.color, xTitle); StyleTH1Slide(h2, hb.color, xTitle);
            
            TH1D* hRatio = (TH1D*)h1->Clone(); hRatio->Divide(h2);
            StyleTH1Slide(hRatio, hb.color, xTitle);
            hRatio->GetYaxis()->CenterTitle(true); hRatio->GetYaxis()->SetTitleSize(0.05); 
            hRatio->GetYaxis()->SetTitleOffset(0.85); hRatio->GetYaxis()->SetNdivisions(505);
            
            v1.push_back(h1); v2.push_back(h2); vRat.push_back(hRatio);
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
        double maxDev = std::max(std::abs(ymax_ratio - 1.0), std::abs(1.0 - ymin_ratio));
        double zoomLim = std::min(maxDev * 1.2, 1.0); 
        for(size_t i=0; i<vRat.size(); ++i) {
            if(i==0) {
                vRat[0]->GetYaxis()->SetRangeUser(ymin_ratio, ymax_ratio);
                if(isLogX) vRat[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                vRat[0]->Draw("EP");
                double lineXm = isLogX ? ptCut : v1[0]->GetXaxis()->GetXmin();
                double lineXM = isLogX ? 500.0 : v1[0]->GetXaxis()->GetXmax();
                if (cName.Contains("cPFVsEta")) {lineXm = -3; lineXM = 3;}
                TLine* line = new TLine(lineXm, 1, lineXM, 1);
                line->SetLineStyle(2); line->SetLineColor(kBlack); line->Draw();
            } else vRat[i]->Draw("EP SAME");
        }
        DrawSlideText(3, false, Form("Ratio to %s", label2.Data()), ptCut, "", etaLabel);
        
        // --- Draw specific vertical lines for ALL PF vs Eta plots ---
        // UPDATED: Triggers for any canvas name containing "cPFVsEta"
        if (cName.Contains("cPFVsEta")) {
            std::vector<double> targetEtas = {-2.3, -1.6, 1.5, 2.1};
            
            // Grab the Y-limits you already calculated for the pads
            double yminMain = isLogY ? 1e-6 : 0.0;
            double ymaxMain = ymax_main * (isLogY ? 5.0 : 1.35);
            double yminRat = ymin_ratio;
            double ymaxRat = ymax_ratio;

            for (double x : targetEtas) {
                TLine* line = new TLine();
                line->SetLineStyle(3); // 3 = dotted
                line->SetLineColor(kBlack);
                line->SetLineWidth(2);
                
                cSlide->cd(1); line->DrawLine(x, yminMain, x, ymaxMain);
                cSlide->cd(2); line->DrawLine(x, yminMain, x, ymaxMain);
                cSlide->cd(3); line->DrawLine(x, yminRat, x, ymaxRat);
            }
        }
        // ------------------------------------------------------------

        // NOW you can save and clean up memory
        cSlide->SaveAs(outDir + "/" + cName + ".png"); 
        delete cSlide;
        for(auto h: v1) delete h; 
        for(auto h: v2) delete h; 
        for(auto h: vRat) delete h;
    }; // <--- THIS IS THE ONLY CLOSING BRACE FOR THE LAMBDA

    // ==============================================================================
    // TASK 1: pT spectra for 1st centrality bin (w/ Integrals)
    // ==============================================================================
    {
        TCanvas* cPtCent0 = new TCanvas("cPtCent0", "pT 1st Centrality Bin", 2000, 2000);
        FormatSlidePad(true, true, false);

        // Project pT (axis 0) for the 1st hiBin, inclusive in eta. 
        // Using custom ProjectTHn1D to respect the specific centralities correctly rather than THnSparse Projection.
        TH1D* hPt1 = (TH1D*) hnKin1->Projection(0);
        TH1D* hPt2 = (TH1D*) hnKin2->Projection(0);
        StyleTH1Slide(hPt1, kBlue, "p_{T} (GeV/c)");
        StyleTH1Slide(hPt2, kRed, "p_{T} (GeV/c)");

        // Calculate integrals for pT >= 50
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

    // ==============================================================================
    // SLIDES: KINEMATICS
    // ==============================================================================

    // 1. pT: For each eta cut (using inclusive pt 50.0 cut base)
    for (const auto& er : etaRanges) {
        auto [v1, v2, vR] = GetKinVectors(0, bins.ptmin, er.lo, er.hi, "p_{T} (GeV/c)");
        DrawAndSaveSlide(Form("cKinSlide_pt_%s", er.name.Data()), "pT", v1, v2, vR, bins.ptmin, er.label, true, true);
    }
    
    // 2. Eta: For each pT cut (using inclusive eta base)
    for (double ptC : ptCutsMap) {
        auto [v1, v2, vR] = GetKinVectors(1, ptC, 0.0, 5.1, "#eta");
        DrawAndSaveSlide(Form("cKinSlide_eta_pt%.0f", ptC), "Eta", v1, v2, vR, ptC, "Inclusive #eta", false, false);
    }
    
    // 3. Phi: Only one inclusive plot
    {
        auto [v1, v2, vR] = GetKinVectors(2, bins.ptmin, 0.0, 5.1, "#phi (rad)");
        DrawAndSaveSlide("cKinSlide_phi_inclusive", "Phi", v1, v2, vR, bins.ptmin, "Inclusive #eta", false, false);
    }

    // ==============================================================================
    // SLIDES: PF FRACTIONS 
    // ==============================================================================
    for (const auto& er : etaRanges) {
        for (int pf = 0; pf < 5; ++pf) {
            auto [v1, v2, vR] = GetPFVectors(pf, er.lo, er.hi, PFTypeTitles.at(pf));
            DrawAndSaveSlide(Form("cPFSlide_%s_%s", PFTypeNames.at(pf), er.name.Data()), PFTypeTitles.at(pf), 
                             v1, v2, vR, bins.ptmin, er.label, PFTypeLogY(pf), false);
        }
    }
    
    // ==============================================================================
    // TASK 2: All PF Fractions vs Eta (3x1 Canvas)
    // ==============================================================================
    for (int pf = 0; pf < 5; ++pf) {
        std::vector<TH1D*> v1_eta, v2_eta, vRat_eta;
        
        for (const auto& hb : bins.hiBins) {
            // Project X=axis 2 (eta), Y=axis 0 (PF fraction) dynamically for PF Type 'pf'
            TH2D* h2D_1 = ProjectTHn2D(hnPF1, 2, 0, {{1, (double)pf, (double)pf + 1.0}, {3, (double)hb.lo, (double)hb.hi}}, Form("pf_2d_1_%d_type%d", (int)hb.lo, pf));
            TH2D* h2D_2 = ProjectTHn2D(hnPF2, 2, 0, {{1, (double)pf, (double)pf + 1.0}, {3, (double)hb.lo, (double)hb.hi}}, Form("pf_2d_2_%d_type%d", (int)hb.lo, pf));
            
            // Extract the profile (mean of fraction in each eta bin)
            TProfile* p1 = h2D_1->ProfileX(Form("prof1_%d_type%d", (int)hb.lo, pf));
            TProfile* p2 = h2D_2->ProfileX(Form("prof2_%d_type%d", (int)hb.lo, pf));
            
            // Cast down to TH1D for uniform styling via standard drawing struct
            TH1D* h1 = p1->ProjectionX(Form("hpf_eta_1_%d_type%d", (int)hb.lo, pf));
            TH1D* h2 = p2->ProjectionX(Form("hpf_eta_2_%d_type%d", (int)hb.lo, pf));

            if (pf == 0) printf("hiLow %.1f hiHi %.1f entries %.2f mean %.2f \n",(double)hb.lo, (double)hb.hi, h1->GetEntries(), h1->GetMean());
            
            h1->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h2->GetXaxis()->SetRangeUser(-3.0, 3.0);
            
            StyleTH1Slide(h1, hb.color, "#eta");
            StyleTH1Slide(h2, hb.color, "#eta");
            
            TH1D* hRatio = (TH1D*)h1->Clone(Form("hrat_%d_type%d", (int)hb.lo, pf));
            hRatio->Divide(h2);
            
            StyleTH1Slide(hRatio, hb.color, "#eta");
            hRatio->GetYaxis()->CenterTitle(true); 
            hRatio->GetYaxis()->SetTitleSize(0.05); 
            hRatio->GetYaxis()->SetTitleOffset(0.85); 
            hRatio->GetYaxis()->SetNdivisions(505);

            hRatio->GetXaxis()->SetRangeUser(-3.0, 3.0);
            
            v1_eta.push_back(h1); v2_eta.push_back(h2); vRat_eta.push_back(hRatio);
            
            delete h2D_1; delete h2D_2; delete p1; delete p2;
        }
        
        // Pass arrays to standard slide formatter, dynamically constructing names and titles via mapping
        TString canvasName = TString("cPFVsEta_") + PFTypeNames.at(pf);
        TString canvasTitle = TString(PFTypeTitles.at(pf)) + " vs #eta";
        DrawAndSaveSlide(canvasName, canvasTitle, v1_eta, v2_eta, vRat_eta, bins.ptmin, PFTypeTitles.at(pf), false, false);
    }

    std::cout << "All plots successfully saved to " << outDir << "/" << std::endl;
}