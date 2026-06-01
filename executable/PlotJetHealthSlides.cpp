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

// Bring in your custom analysis framework headers
#include "../header/Binning.h"
#include "../header/JetHealthHistograms.h"
#include "../header/JetHealthPlotting.h"

// ==============================================================================
// 1. CONSTANTS & CONFIGURATION
// ==============================================================================

const TString file2026Data = "JetHealth_404350_CaloJets.root"; 
const TString file2025Data = "JetHealth_2025_Data.root";
const TString file2026MC   = "JetHealth_2026_MC.root";

const TString runNumber = "404359";
const TString jetAlgo   = "akPu4Calo";

// Which ratio to plot? (true = 2026 MC, false = 2025 Data)
const bool useMCForRatio = false; // Defaults to 2025 Data

// Dynamic Cut Maps
const std::vector<double> ptCutsMap = {50.0, 100.0, 200.0};
const std::vector<double> etaCutsMap = {0.0, 1.6, 2.4, 5.1};

const TString outDir = "JetHealthSlides_" + runNumber + "_" + jetAlgo;

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
            tex.SetTextAlign(22); 
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

void PlotJetHealthSlides() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRedBlue);
    gSystem->mkdir(outDir, true);

    TFile* f26Data = TFile::Open(file2026Data, "READ");
    TFile* f25Data = TFile::Open(file2025Data, "READ");
    TFile* f26MC   = TFile::Open(file2026MC, "READ");
    
    if (!f26Data || f26Data->IsZombie()) { std::cerr << "Cannot open 2026 Data file!" << std::endl; return; }
    
    double nEvt26Data = GetNevents(f26Data);
    double nEvtComp   = 1.0;
    TFile* fComp      = nullptr;
    TString compLabel = "";

    if (useMCForRatio) {
        nEvtComp = GetNevents(f26MC); fComp = f26MC; compLabel = "2026 MC";
    } else {
        nEvtComp = GetNevents(f25Data); fComp = f25Data; compLabel = "2025 Data";
    }

    THnSparseF* hnKin26Data = (THnSparseF*)f26Data->Get("hjetkin");
    THnSparseF* hnKin25Data = (THnSparseF*)f25Data->Get("hjetkin");
    THnSparseF* hnKin26MC = (THnSparseF*)f26MC->Get("hjetkin");
    THnSparseF* hnKinComp   = (THnSparseF*)fComp->Get("hjetkin");
    THnSparseF* hnPF26Data  = (THnSparseF*)f26Data->Get("hjetpf");
    THnSparseF* hnPFComp    = (THnSparseF*)fComp->Get("hjetpf");

    BinningStruct bins(50.0);

    // Build Eta Ranges directly from etaCutsMap
    struct EtaRange { double lo; double hi; TString name; TString label; };
    std::vector<EtaRange> etaRanges;
    for (size_t i = 0; i < etaCutsMap.size() - 1; ++i) {
        etaRanges.push_back({etaCutsMap[i], etaCutsMap[i+1], Form("eta_%.1f_%.1f", etaCutsMap[i], etaCutsMap[i+1]), 
                             etaCutsMap[i] == 0 ? Form("|#eta| < %.1f", etaCutsMap[i+1]) : Form("%.1f < |#eta| < %.1f", etaCutsMap[i], etaCutsMap[i+1])});
    }
    etaRanges.push_back({0.0, 5.1, "eta_inclusive", "Inclusive #eta"});

    // ------------------------------------------------------------------------------
    // REUSABLE LAMBDAS FOR 1D PROJECTIONS & DRAWING
    // ------------------------------------------------------------------------------
    
    auto StyleTH1Slide = [](TH1D* h, Color_t color, const TString& xTitle) {
        StyleTH1(h, color);
        h->SetLineWidth(1); h->SetMarkerSize(0.01);
        h->SetTitle(""); h->GetYaxis()->SetTitle(""); 
        h->GetXaxis()->SetTitle(xTitle); h->GetXaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleSize(0.05); h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05); h->GetXaxis()->SetTitleOffset(0.85); 
    };

    auto GetKinVectors = [&](int axis, double ptCut, double etaLo, double etaHi, const TString& xTitle) {
        std::vector<TH1D*> v26, vComp, vRat;
        for (const auto& hb : bins.hiBins) {
            TH1D *h26, *hComp;
            if (etaLo == 0.0) {
                h26   = ProjectTHn1D(hnKin26Data, axis, {{0, ptCut, 1000.0}, {1, -etaHi, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_26");
                hComp = ProjectTHn1D(hnKinComp,   axis, {{0, ptCut, 1000.0}, {1, -etaHi, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_comp");
            } else {
                TH1D* h26_p = ProjectTHn1D(hnKin26Data, axis, {{0, ptCut, 1000.0}, {1, etaLo, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_26p");
                TH1D* h26_n = ProjectTHn1D(hnKin26Data, axis, {{0, ptCut, 1000.0}, {1, -etaHi, -etaLo}, {3, (double)hb.lo, (double)hb.hi}}, "_26n");
                h26_p->Add(h26_n); h26 = h26_p; delete h26_n;
                
                TH1D* hCp = ProjectTHn1D(hnKinComp, axis, {{0, ptCut, 1000.0}, {1, etaLo, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_Cp");
                TH1D* hCn = ProjectTHn1D(hnKinComp, axis, {{0, ptCut, 1000.0}, {1, -etaHi, -etaLo}, {3, (double)hb.lo, (double)hb.hi}}, "_Cn");
                hCp->Add(hCn); hComp = hCp; delete hCn;
            }
            NormalizeTH1(h26); NormalizeTH1(hComp);
            StyleTH1Slide(h26, hb.color, xTitle); StyleTH1Slide(hComp, hb.color, xTitle);
            
            TH1D* hRatioNum = (TH1D*)h26->Clone(); TH1D* hRatioDen = (TH1D*)hComp->Clone();
            if (axis == 1 || axis == 2) { hRatioNum->Rebin(4); hRatioDen->Rebin(4); }
            TH1D* hRatio = (TH1D*)hRatioNum->Clone(); hRatio->Divide(hRatioDen);
            
            StyleTH1Slide(hRatio, hb.color, xTitle);
            //hRatio->GetYaxis()->SetTitle(Form("Ratio to %s", compLabel.Data()));
            hRatio->GetYaxis()->CenterTitle(true); hRatio->GetYaxis()->SetTitleSize(0.05); 
            hRatio->GetYaxis()->SetTitleOffset(0.85); hRatio->GetYaxis()->SetNdivisions(505);
            
            v26.push_back(h26); vComp.push_back(hComp); vRat.push_back(hRatio);
            delete hRatioNum; delete hRatioDen;
        }
        return std::make_tuple(v26, vComp, vRat);
    };

    auto GetPFVectors = [&](int pf, double etaLo, double etaHi, const TString& xTitle) {
        std::vector<TH1D*> v26, vComp, vRat;
        for (const auto& hb : bins.hiBins) {
            TH1D* h26   = ProjectTHn1D(hnPF26Data, 0, {{1, (double)pf, (double)pf + 1.0}, {2, etaLo, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_26");
            TH1D* hComp = ProjectTHn1D(hnPFComp,   0, {{1, (double)pf, (double)pf + 1.0}, {2, etaLo, etaHi}, {3, (double)hb.lo, (double)hb.hi}}, "_comp");
            
            NormalizeTH1(h26); NormalizeTH1(hComp);
            StyleTH1Slide(h26, hb.color, xTitle); StyleTH1Slide(hComp, hb.color, xTitle);
            
            TH1D* hRatio = (TH1D*)h26->Clone(); hRatio->Divide(hComp);
            StyleTH1Slide(hRatio, hb.color, xTitle);
            //hRatio->GetYaxis()->SetTitle(Form("Ratio to %s", compLabel.Data()));
            hRatio->GetYaxis()->CenterTitle(true); hRatio->GetYaxis()->SetTitleSize(0.05); 
            hRatio->GetYaxis()->SetTitleOffset(0.85); hRatio->GetYaxis()->SetNdivisions(505);
            
            v26.push_back(h26); vComp.push_back(hComp); vRat.push_back(hRatio);
        }
        return std::make_tuple(v26, vComp, vRat);
    };

    auto DrawAndSaveSlide = [&](const TString& cName, const TString& title, std::vector<TH1D*>& v26, std::vector<TH1D*>& vComp, std::vector<TH1D*>& vRat, 
                                double ptCut, const TString& etaLabel, bool isLogY, bool isLogX) {
        TCanvas* cSlide = new TCanvas(cName, title, 3600, 1200); cSlide->Divide(3, 1);
        TLegend* leg = new TLegend(0.40, 0.55, 0.90, 0.80); leg->SetBorderSize(0); leg->SetTextSize(0.045);
        for(size_t i=0; i<bins.hiBins.size(); ++i) leg->AddEntry(v26[i], bins.hiBins[i].title, "lp");
        
        double ymax_main = 0.0, ymax_ratio = -999.0, ymin_ratio = 999.0;
        for(auto h : v26) ymax_main = std::max(ymax_main, h->GetMaximum());
        for(auto h : vComp) ymax_main = std::max(ymax_main, h->GetMaximum());
        for(auto h : vRat) {
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                if (h->GetBinContent(b) > 0) {
                    ymax_ratio = std::max(ymax_ratio, h->GetBinContent(b));
                    ymin_ratio = std::min(ymin_ratio, h->GetBinContent(b));
                }
            }
        }
        
        cSlide->cd(1); FormatSlidePad(isLogX, isLogY, false);
        for(size_t i=0; i<v26.size(); ++i) {
            if(i==0) {
                v26[0]->GetYaxis()->SetRangeUser(isLogY ? 1e-6 : 0.0, ymax_main * (isLogY ? 5.0 : 1.35));
                if(isLogX) v26[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                v26[0]->Draw("EP");
            } else v26[i]->Draw("EP SAME");
        }
        leg->Draw("SAME"); 
        DrawSlideText(1, false, Form("2026 Data (Run %s)", runNumber.Data()), ptCut, jetAlgo, etaLabel);
        
        cSlide->cd(2); FormatSlidePad(isLogX, isLogY, false);
        for(size_t i=0; i<vComp.size(); ++i) {
            if(i==0) {
                vComp[0]->GetYaxis()->SetRangeUser(isLogY ? 1e-6 : 0.0, ymax_main * (isLogY ? 5.0 : 1.35));
                if(isLogX) vComp[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                vComp[0]->Draw("EP");
            } else vComp[i]->Draw("EP SAME");
        }
        DrawSlideText(2, false, compLabel, ptCut, jetAlgo, etaLabel);
        
        cSlide->cd(3); FormatSlidePad(isLogX, false, false);
        double maxDev = std::max(std::abs(ymax_ratio - 1.0), std::abs(1.0 - ymin_ratio));
        double zoomLim = std::min(maxDev * 1.2, 1.0); 
        for(size_t i=0; i<vRat.size(); ++i) {
            if(i==0) {
                vRat[0]->GetYaxis()->SetRangeUser(std::max(0.0, 1.0 - zoomLim), std::min(2.0, 1.0 + zoomLim));
                if(isLogX) vRat[0]->GetXaxis()->SetRangeUser(ptCut, 500);
                vRat[0]->Draw("EP");
                double lineXm = isLogX ? ptCut : vRat[0]->GetXaxis()->GetXmin();
                double lineXM = isLogX ? 500.0 : vRat[0]->GetXaxis()->GetXmax();
                TLine* line = new TLine(lineXm, 1, lineXM, 1);
                line->SetLineStyle(2); line->SetLineColor(kBlack); line->Draw();
            } else vRat[i]->Draw("EP SAME");
        }
        DrawSlideText(3, false, Form("Ratio to %s", compLabel.Data()), ptCut, jetAlgo, etaLabel);
        
        cSlide->SaveAs(outDir + "/" + cName + ".png"); delete cSlide;
        for(auto h: v26) delete h; for(auto h: vComp) delete h; for(auto h: vRat) delete h;
    };

    // ==============================================================================
    // SLIDES: KINEMATICS
    // ==============================================================================

    // 1. pT: For each eta cut (using inclusive pt 50.0 cut base)
    for (const auto& er : etaRanges) {
        auto [v26, vC, vR] = GetKinVectors(0, bins.ptmin, er.lo, er.hi, "p_{T} (GeV/c)");
        DrawAndSaveSlide(Form("cKinSlide_pt_%s", er.name.Data()), "pT", v26, vC, vR, bins.ptmin, er.label, true, true);
    }
    
    // 2. Eta: For each pT cut (using inclusive eta base)
    for (double ptC : ptCutsMap) {
        auto [v26, vC, vR] = GetKinVectors(1, ptC, 0.0, 5.1, "#eta");
        DrawAndSaveSlide(Form("cKinSlide_eta_pt%.0f", ptC), "Eta", v26, vC, vR, ptC, "Inclusive #eta", false, false);
    }
    
    // 3. Phi: Only one inclusive plot
    {
        auto [v26, vC, vR] = GetKinVectors(2, bins.ptmin, 0.0, 5.1, "#phi (rad)");
        DrawAndSaveSlide("cKinSlide_phi_inclusive", "Phi", v26, vC, vR, bins.ptmin, "Inclusive #eta", false, false);
    }

    // ==============================================================================
    // SLIDES: PF FRACTIONS 
    // ==============================================================================
    for (const auto& er : etaRanges) {
        for (int pf = 0; pf < 5; ++pf) {
            auto [v26, vC, vR] = GetPFVectors(pf, er.lo, er.hi, PFTypeTitles.at(pf));
            DrawAndSaveSlide(Form("cPFSlide_%s_%s", PFTypeNames.at(pf), er.name.Data()), PFTypeTitles.at(pf), 
                             v26, vC, vR, bins.ptmin, er.label, PFTypeLogY(pf), false);
        }
    }

    // ==============================================================================
    // CANVAS 2: 2D ETA-PHI MAPS
    // ==============================================================================
    if(f25Data && !f25Data->IsZombie()) {
        for (double ptC : ptCutsMap) {
            for (const auto& hb : bins.hiBins) { 
                TString cName = Form("cMap_pt%.0f_hb%.0f_%.0f", ptC, hb.lo, hb.hi);
                TCanvas* cMap = new TCanvas(cName, "Eta-Phi Map", 2400, 1000); cMap->Divide(3, 1);
                
                cMap->cd(1); FormatSlidePad(false, false, true);
                TH2D* h26 = ProjectTHn2D(hnKin26Data, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "map26");
                h26->Scale(1.0 / nEvt26Data); h26->SetTitle(";#eta;#phi (rad)");
                h26->GetXaxis()->SetTitleOffset(0.85); h26->GetYaxis()->SetTitleOffset(0.85); h26->Draw("colz");
                DrawSlideText(1, true, Form("2026 Data (Run %s)", runNumber.Data()), ptC, jetAlgo);
                
                cMap->cd(2); FormatSlidePad(false, false, true);
                TH2D* h25 = ProjectTHn2D(hnKinComp, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "map25");
                h25->Scale(1.0 / nEvtComp); h25->SetTitle(";#eta;#phi (rad)");
                h25->GetXaxis()->SetTitleOffset(0.85); h25->GetYaxis()->SetTitleOffset(0.85); h25->Draw("colz");
                DrawSlideText(2, true, compLabel, ptC, jetAlgo, Form("#bf{hiBin %.0f-%.0f}", hb.lo, hb.hi));
                
                cMap->cd(3); FormatSlidePad(false, false, true);
                TH2D* hMC = ProjectTHn2D(hnKin26MC, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "mapMC");
                hMC->Scale(1.0 / nEvtComp); hMC->SetTitle(";#eta;#phi (rad)");
                hMC->GetXaxis()->SetTitleOffset(0.85); hMC->GetYaxis()->SetTitleOffset(0.85); hMC->Draw("colz");
                DrawSlideText(3, true, "2026 MC", ptC, jetAlgo);
                
                cMap->SaveAs(outDir + "/" + cName + ".png");
                delete h26; delete h25; delete hMC; delete cMap;
            }
        }
    }

    std::cout << "All plots successfully saved to " << outDir << "/" << std::endl;
}

/*
==============================================================================
MACRO DOCUMENTATION (SLIDES FORMAT)
==============================================================================
What this code does:
This macro generates highly-uniform 3x1 slide-ready canvases directly from THnSparse.
It scales layouts dynamically, automatically generating plots for all combinations of 
ptCutsMap and etaCutsMap provided.

Normalizations:
- 1D Plots (Kinematics & PF): Shape Normalized (`1.0/Integral()`) via `NormalizeTH1()`.
- 2D Maps (Eta-Phi): Yield Normalized by 1.0/N_events (calculated from hvz vertex).

Constraints & Notes:
- Kinematic symmetric cuts: When an eta boundary > 0 is selected (e.g. 1.6 < |eta| < 2.4), 
  the macro seamlessly splices [-2.4, -1.6] and [1.6, 2.4] from `hjetkin` into a single histogram.
- PF Fractions pT restriction: `hjetpf` only possesses 4 axes (Frac, Type, absEta, hiBin).
  Therefore, pT cut filtering cannot be dynamically applied in post. The plot headers
  will default to referencing the inclusive pT filter executed inside the analyzer (50.0).
==============================================================================
*/