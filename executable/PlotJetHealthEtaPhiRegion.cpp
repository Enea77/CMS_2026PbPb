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
#include "TBox.h"

// Bring in your custom analysis framework headers
#include "../header/Binning.h"
#include "../header/JetHealthHistograms.h"
#include "../header/JetHealthPlotting.h"

// ==============================================================================
// 1. CONSTANTS & CONFIGURATION
// ==============================================================================

// Dynamic Cut Maps
const std::vector<double> ptCutsMap = {50.0};

// Custom Eta-Phi Target Regions
struct Region { 
    double etaLo; double etaHi; 
    double phiLo; double phiHi; 
    TString name; TString label; 
};

const std::vector<Region> targetRegions = {
    {-2.3, -1.6, -2.15, -1.75, "Region1", "-2.3 < #eta < -1.6, -2.15 < #phi < -1.75"},
    { 1.5,  2.1, -3.0, -2.5, "Region2",   "1.5 < #eta < 2.1, -3.0 < #phi < -2.5"}
};

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
            tex.SetTextAlign(11); 
            tex.DrawLatex(leftX, 0.94, extraLabel); 
        }
        if (is2D){
            tex.SetTextAlign(31); 
            tex.DrawLatex(rightX, 0.94, Form("#bf{%s}", jetAlgo.Data())); 
        }
    } else if (padNum == 3) {
        tex.SetTextAlign(31); 
        tex.DrawLatex(rightX, 0.94, Form("#bf{%s}", jetAlgo.Data())); 
    }
}

// ==============================================================================
// MAIN PLOTTING MACRO
// ==============================================================================

void PlotJetHealthEtaPhiRegion(TString file1Path = "JetHealth_404350", TString label1="Run404350", 
                               TString file2Path = "JetHealth_2026_MC", TString label2="2026 MC") {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRedBlue);
    
    TString outDir = file1Path + "_vs_" + file2Path + "_Regions";
    gSystem->mkdir(outDir, true);
    TString jetAlgo = "akCs4PF";

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

    // ------------------------------------------------------------------------------
    // GENERATE AND SAVE COMMON LEGEND SEPARATELY
    // ------------------------------------------------------------------------------
    TCanvas* cLegend = new TCanvas("cLegend", "Legend", 1600, 1200);
    TLegend* commonLeg = new TLegend(0.1, 0.1, 0.9, 0.9);
    commonLeg->SetBorderSize(0); 
    commonLeg->SetTextSize(0.06); // Larger text for standalone image
    
    // Create dummy histograms just to populate the legend markers correctly
    std::vector<TH1D*> dummyHists;
    for(size_t i=0; i<bins.hiBins.size(); ++i) {
        TH1D* hDummy = new TH1D(Form("dummy_%zu", i), "", 1, 0, 1);
        hDummy->SetMarkerColor(bins.hiBins[i].color);
        hDummy->SetLineColor(bins.hiBins[i].color);
        hDummy->SetFillColor(bins.hiBins[i].color);
        hDummy->SetFillStyle(1001);
        hDummy->SetMarkerStyle(20);
        dummyHists.push_back(hDummy);
        commonLeg->AddEntry(hDummy, bins.hiBins[i].title, "f");
    }
    
    commonLeg->Draw();
    cLegend->SaveAs(outDir + "/Common_Legend.png");
    delete commonLeg;
    delete cLegend;
    for(auto h : dummyHists) delete h;


    // ------------------------------------------------------------------------------
    // REUSABLE LAMBDAS FOR 1D PROJECTIONS & DRAWING
    // ------------------------------------------------------------------------------
    
    auto StyleTH1Slide = [](TH1D* h, Color_t color, const TString& xTitle) {
        StyleTH1(h, color);
        h->SetFillColor(color);
        h->SetLineWidth(4); h->SetMarkerSize(0.01);
        h->SetTitle(""); h->GetYaxis()->SetTitle(""); 
        h->GetXaxis()->SetTitle(xTitle); h->GetXaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleSize(0.05); h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05); h->GetXaxis()->SetTitleOffset(0.85); 
    };

    auto GetKinVectors = [&](int axis, double ptCut, const Region& reg, const TString& xTitle) {
        std::vector<TH1D*> v1, v2, vRat;
        for (const auto& hb : bins.hiBins) {
            std::vector<SparseRange> cuts = {
                {0, ptCut, 1000.0}, 
                {1, reg.etaLo, reg.etaHi}, 
                {2, reg.phiLo, reg.phiHi}, 
                {3, (double)hb.lo, (double)hb.hi}
            };
            
            TH1D* h1 = ProjectTHn1D(hnKin1, axis, cuts, "_1");
            TH1D* h2 = ProjectTHn1D(hnKin2, axis, cuts, "_2");
            
            NormalizeTH1(h1); NormalizeTH1(h2);
            StyleTH1Slide(h1, hb.color, xTitle); StyleTH1Slide(h2, hb.color, xTitle);
            
            TH1D* hRatioNum = (TH1D*)h1->Clone(); TH1D* hRatioDen = (TH1D*)h2->Clone();
            TH1D* hRatio = (TH1D*)hRatioNum->Clone(); hRatio->Divide(hRatioDen);
            
            StyleTH1Slide(hRatio, hb.color, xTitle);
            hRatio->GetYaxis()->CenterTitle(true); hRatio->GetYaxis()->SetTitleSize(0.05); 
            hRatio->GetYaxis()->SetTitleOffset(0.85); hRatio->GetYaxis()->SetNdivisions(505);
            
            v1.push_back(h1); v2.push_back(h2); vRat.push_back(hRatio);
            delete hRatioNum; delete hRatioDen;
        }
        return std::make_tuple(v1, v2, vRat);
    };

    auto GetPFVectors = [&](int pf, double ptCut, const Region& reg, const TString& xTitle) {
        std::vector<TH1D*> v1, v2, vRat;
        for (const auto& hb : bins.hiBins) {
            std::vector<SparseRange> cuts = {
                {1, (double)pf, (double)pf + 1.0}, 
                {2, reg.etaLo, reg.etaHi}, 
                {3, reg.phiLo, reg.phiHi}, 
                {4, ptCut, 1000.0},                
                {5, (double)hb.lo, (double)hb.hi}  
            };
            
            TH1D* h1 = ProjectTHn1D(hnPF1, 0, cuts, "_1");
            TH1D* h2 = ProjectTHn1D(hnPF2, 0, cuts, "_2");
            
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
                                double ptCut, const TString& regionLabel, bool isLogY, bool isLogX) {
        TCanvas* cSlide = new TCanvas(cName, title, 3600, 1200); cSlide->Divide(3, 1);
        
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
        
        // Tighter framing: 1.05 for linear, 5.0 for log to leave just a sliver of space
        double yFrameMult = isLogY ? 5.0 : 1.05;

        cSlide->cd(1); FormatSlidePad(isLogX, isLogY, false);
        for(size_t i=0; i<v1.size(); ++i) {
            if(i==0) {
                v1[0]->GetYaxis()->SetRangeUser(isLogY ? 1e-6 : 0.0, ymax_main * yFrameMult);
                if(isLogX) v1[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                v1[0]->Draw("EP");
            } else v1[i]->Draw("EP SAME");
        }
        DrawSlideText(1, false, label1, ptCut, jetAlgo, regionLabel);
        
        cSlide->cd(2); FormatSlidePad(isLogX, isLogY, false);
        for(size_t i=0; i<v2.size(); ++i) {
            if(i==0) {
                v2[0]->GetYaxis()->SetRangeUser(isLogY ? 1e-6 : 0.0, ymax_main * yFrameMult);
                if(isLogX) v2[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                v2[0]->Draw("EP");
            } else v2[i]->Draw("EP SAME");
        }
        DrawSlideText(2, false, label2, ptCut, jetAlgo, regionLabel);
        
        cSlide->cd(3); FormatSlidePad(isLogX, false, false);
        double maxDev = std::max(std::abs(ymax_ratio - 1.0), std::abs(1.0 - ymin_ratio));
        double zoomLim = std::min(maxDev * 1.2, 1.0); 
        TLine* line = nullptr; // Track memory for deletion
        
        for(size_t i=0; i<vRat.size(); ++i) {
            if(i==0) {
                vRat[0]->GetYaxis()->SetRangeUser(std::max(0.0, 1.0 - zoomLim), std::min(2.0, 1.0 + zoomLim));
                if(isLogX) vRat[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                vRat[0]->Draw("EP");
                double lineXm = isLogX ? ptCut : vRat[0]->GetXaxis()->GetXmin();
                double lineXM = isLogX ? 500.0 : vRat[0]->GetXaxis()->GetXmax();
                line = new TLine(lineXm, 1, lineXM, 1);
                line->SetLineStyle(2); line->SetLineColor(kBlack); line->Draw();
            } else vRat[i]->Draw("EP SAME");
        }
        DrawSlideText(3, false, Form("Ratio to %s", label2.Data()), ptCut, jetAlgo, regionLabel);
        
        cSlide->SaveAs(outDir + "/" + cName + ".png"); 
        
        // Cleanup memory
        delete line; 
        delete cSlide;
        for(auto h: v1) delete h; for(auto h: v2) delete h; for(auto h: vRat) delete h;
    };

    // ==============================================================================
    // SLIDES: KINEMATICS (Per Region)
    // ==============================================================================

    for (const auto& reg : targetRegions) {
        auto [v1_pt, v2_pt, vR_pt] = GetKinVectors(0, bins.ptmin, reg, "p_{T} (GeV/c)");
        DrawAndSaveSlide(Form("cKinSlide_pt_%s", reg.name.Data()), "pT", v1_pt, v2_pt, vR_pt, bins.ptmin, reg.label, true, true);
        
        auto [v1_eta, v2_eta, vR_eta] = GetKinVectors(1, bins.ptmin, reg, "#eta");
        DrawAndSaveSlide(Form("cKinSlide_eta_%s", reg.name.Data()), "Eta", v1_eta, v2_eta, vR_eta, bins.ptmin, reg.label, false, false);

        auto [v1_phi, v2_phi, vR_phi] = GetKinVectors(2, bins.ptmin, reg, "#phi (rad)");
        DrawAndSaveSlide(Form("cKinSlide_phi_%s", reg.name.Data()), "Phi", v1_phi, v2_phi, vR_phi, bins.ptmin, reg.label, false, false);
    }

    // ==============================================================================
    // SLIDES: PF FRACTIONS (Per Region)
    // ==============================================================================
    
    for (const auto& reg : targetRegions) {
        for (int pf = 0; pf < 5; ++pf) {
            auto [v1, v2, vR] = GetPFVectors(pf, bins.ptmin, reg, PFTypeTitles.at(pf));
            DrawAndSaveSlide(Form("cPFSlide_%s_%s", PFTypeNames.at(pf), reg.name.Data()), PFTypeTitles.at(pf), 
                             v1, v2, vR, bins.ptmin, reg.label, PFTypeLogY(pf), false);
        }
    }

    // ==============================================================================
    // CANVAS 2: 2D ETA-PHI MAPS (Inclusive kinematics + Highlight Boxes)
    // ==============================================================================
    for (double ptC : ptCutsMap) {
        for (const auto& hb : bins.hiBins) { 
            TString cName = Form("cMap_pt%.0f_hb%.0f_%.0f", ptC, hb.lo, hb.hi);
            TCanvas* cMap = new TCanvas(cName, "Eta-Phi Map", 3200, 1600); cMap->Divide(2, 1);
            
            std::vector<TBox*> boxesToCleanup; // Memory tracker for dynamically allocated TBoxes
            
            // File 1 Map
            cMap->cd(1); FormatSlidePad(false, false, true);
            TH2D* h1 = ProjectTHn2D(hnKin1, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "map1");
            h1->Scale(1.0 / nEvt1); h1->SetTitle(";#eta;#phi (rad)");
            h1->GetXaxis()->SetTitleOffset(0.85); h1->GetYaxis()->SetTitleOffset(0.85); h1->Draw("colz");
            
            for (const auto& reg : targetRegions) {
                TBox* box = new TBox(reg.etaLo, reg.phiLo, reg.etaHi, reg.phiHi);
                box->SetFillStyle(0);
                box->SetLineColor(kRed);
                box->SetLineWidth(3);
                box->Draw("SAME");
                boxesToCleanup.push_back(box);
            }
            DrawSlideText(1, true, label1, ptC, jetAlgo);
            
            // File 2 Map
            cMap->cd(2); FormatSlidePad(false, false, true);
            TH2D* h2 = ProjectTHn2D(hnKin2, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "map2");
            h2->Scale(1.0 / nEvt2); h2->SetTitle(";#eta;#phi (rad)");
            h2->GetXaxis()->SetTitleOffset(0.85); h2->GetYaxis()->SetTitleOffset(0.85); h2->Draw("colz");
            
            for (const auto& reg : targetRegions) {
                TBox* box = new TBox(reg.etaLo, reg.phiLo, reg.etaHi, reg.phiHi);
                box->SetFillStyle(0);
                box->SetLineColor(kRed);
                box->SetLineWidth(3);
                box->Draw("SAME");
                boxesToCleanup.push_back(box);
            }
            DrawSlideText(2, true, label2, ptC, jetAlgo, Form("#bf{hiBin %.0f-%.0f}", hb.lo, hb.hi));
            
            cMap->SaveAs(outDir + "/" + cName + ".png");
            
            // Explicit cleanup
            for (auto b : boxesToCleanup) delete b; 
            delete h1; delete h2; delete cMap;
        }
    }

    std::cout << "All plots successfully saved to " << outDir << "/" << std::endl;
}