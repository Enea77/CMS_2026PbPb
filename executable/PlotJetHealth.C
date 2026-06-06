#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

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

const TString file2026Data = "output_JetHealth_404350.root"; 
const TString file2025Data = "JetHealth_2025_Data.root";
const TString file2026MC   = "JetHealth_2026_MC.root";

const TString runNumber = "404350";
const TString jetAlgo   = "akCs4PF";

// Which ratio to plot? (true = 2026 MC, false = 2025 Data)
const bool useMCForRatio = false;

// pt cuts for the 2D Eta-Phi maps
const std::vector<double> ptCutsMap = {50.0, 100.0, 200.0};

// Output directory
const TString outDir = "JetHealthPlots_"+ runNumber;

// ==============================================================================
// HELPER FUNCTIONS
// ==============================================================================

// Extract number of events for the un-normalized 2D maps
double GetNevents(TFile* f) {
    if (!f || f->IsZombie()) return 1.0;
    TH1F* hvz = (TH1F*)f->Get("hvz");
    if (!hvz) {
        std::cerr << "WARNING: Could not find 'hvz' in " << f->GetName() << ". Normalizing by 1." << std::endl;
        return 1.0;
    }
    return hvz->Integral();
}

// ==============================================================================
// MAIN PLOTTING MACRO
// ==============================================================================

void PlotJetHealth() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRedBlue);
    gSystem->mkdir(outDir, true);

    // 1. Open files
    TFile* f26Data = TFile::Open(file2026Data, "READ");
    TFile* f25Data = TFile::Open(file2025Data, "READ");
    TFile* f26MC   = TFile::Open(file2026MC, "READ");
    
    if (!f26Data || f26Data->IsZombie()) { std::cerr << "Cannot open 2026 Data file!" << std::endl; return; }
    
    double nEvt26Data = GetNevents(f26Data);
    double nEvtComp   = 1.0;
    TFile* fComp      = nullptr;
    TString compLabel = "";

    if (useMCForRatio) {
        if (!f26MC || f26MC->IsZombie()) { std::cerr << "Cannot open 2026 MC file!" << std::endl; return; }
        nEvtComp = GetNevents(f26MC);
        fComp = f26MC;
        compLabel = "2026 MC";
    } else {
        if (!f25Data || f25Data->IsZombie()) { std::cerr << "Cannot open 2025 Data file!" << std::endl; return; }
        nEvtComp = GetNevents(f25Data);
        fComp = f25Data;
        compLabel = "2025 Data";
    }

    THnSparseF* hnKin26Data = (THnSparseF*)f26Data->Get("hjetkin");
    THnSparseF* hnKinComp   = (THnSparseF*)fComp->Get("hjetkin");
    THnSparseF* hnPF26Data  = (THnSparseF*)f26Data->Get("hjetpf");
    THnSparseF* hnPFComp    = (THnSparseF*)fComp->Get("hjetpf");

    // Init Binning config from header
    BinningStruct bins(50.0);

    // ==============================================================================
    // CANVAS 1: KINEMATICS (3x2)
    // ==============================================================================
    TCanvas* cKin = new TCanvas("cKin", "Jet Kinematics", 2400, 1600);
    cKin->Divide(3, 2);
    
    const char* kinTitles[3] = {"p_{T} (GeV/c)", "#eta", "#phi (rad)"};
    
    TLatex texKin;
    texKin.SetNDC();
    texKin.SetTextFont(42);
    texKin.SetTextSize(0.05);
    
    for (int axis = 0; axis < 3; axis++) {
        // --- TOP ROW: MAIN DISTRIBUTIONS ---
        cKin->cd(axis + 1);
        gPad->SetTopMargin(0.12);
        gPad->SetBottomMargin(0.12);
        gPad->SetLeftMargin(0.12);
        
        if (axis == 0) {
            gPad->SetLogy(true);
            gPad->SetLogx(true);
        }
        
        // Legend shifted slightly left to avoid plotting area cutoff
        TLegend* leg = new TLegend(0.30, 0.55, 0.80, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.045);
        
        std::vector<TH1D*> hMainVec;
        std::vector<TH1D*> hRatioVec;
        
        double ymax_main  = 0.0;
        double ymax_ratio = -999.0;
        double ymin_ratio = 999.0;
        
        // Loop over the hiBins defined in BinningStruct
        for (size_t i = 0; i < bins.hiBins.size(); ++i) {
            const auto& hb = bins.hiBins[i];
            
            TString suf26   = Form("_26_ax%d_hb%zu", axis, i);
            TString sufComp = Form("_comp_ax%d_hb%zu", axis, i);
            
            // Using your native ProjectTHn1D
            TH1D* h26   = ProjectTHn1D(hnKin26Data, axis, {{3, (double)hb.lo, (double)hb.hi}}, suf26);
            TH1D* hComp = ProjectTHn1D(hnKinComp,   axis, {{3, (double)hb.lo, (double)hb.hi}}, sufComp);
            
            // Using your native NormalizeTH1 and StyleTH1
            NormalizeTH1(h26);
            NormalizeTH1(hComp);
            StyleTH1(h26, hb.color);
            h26->SetLineWidth(1);
            h26->SetMarkerSize(0.01);
            
            // Formatting
            h26->SetTitle(""); 
            h26->GetYaxis()->SetTitle(""); 
            h26->GetXaxis()->SetTitle(kinTitles[axis]);
            h26->GetXaxis()->CenterTitle(true);
            h26->GetXaxis()->SetTitleSize(0.05);
            h26->GetXaxis()->SetLabelSize(0.05);
            h26->GetYaxis()->SetLabelSize(0.05);
            
            ymax_main = std::max(ymax_main, h26->GetMaximum());
            hMainVec.push_back(h26);
            leg->AddEntry(h26, hb.title, "lp");
            
            // --- Ratio Calculation with Rebinning ---
            TH1D* hRatioNum = (TH1D*)h26->Clone(Form("ratioNum_ax%d_hb%zu", axis, i));
            TH1D* hRatioDen = (TH1D*)hComp->Clone(Form("ratioDen_ax%d_hb%zu", axis, i));
            
            // Rebin by 4 for eta and phi
            if (axis == 1 || axis == 2) {
                hRatioNum->Rebin(4);
                hRatioDen->Rebin(4);
            }
            
            TH1D* hRatio = (TH1D*)hRatioNum->Clone(Form("ratio_ax%d_hb%zu", axis, i));
            hRatio->Divide(hRatioDen);
            
            StyleTH1(hRatio, hb.color); // Apply native styling to ratio too
            hRatio->SetLineWidth(1);
            hRatio->SetMarkerSize(0.01);
            
            hRatio->SetTitle("");
            hRatio->GetYaxis()->SetTitle(Form("Ratio to %s", compLabel.Data()));
            hRatio->GetYaxis()->CenterTitle(true);
            hRatio->GetYaxis()->SetTitleSize(0.05);
            hRatio->GetYaxis()->SetLabelSize(0.05);
            hRatio->GetYaxis()->SetTitleOffset(1.2);
            hRatio->GetYaxis()->SetNdivisions(505);
            hRatio->GetXaxis()->SetTitleSize(0.05);
            hRatio->GetXaxis()->SetLabelSize(0.05);
            
            for (int b = 1; b <= hRatio->GetNbinsX(); ++b) {
                double val = hRatio->GetBinContent(b);
                if (val > 0) {
                    ymax_ratio = std::max(ymax_ratio, val);
                    ymin_ratio = std::min(ymin_ratio, val);
                }
            }
            hRatioVec.push_back(hRatio);
            
            delete hRatioDen; 
            delete hRatioNum;
        }
        
        cKin->cd(axis + 1);
        for (size_t i = 0; i < hMainVec.size(); ++i) {
            if (i == 0) {
                if (axis == 0) {
                    hMainVec[0]->GetYaxis()->SetRangeUser(1e-6, ymax_main * 5.0);
                    hMainVec[0]->GetXaxis()->SetRangeUser(bins.ptmin, 500);
                } else {
                    hMainVec[0]->GetYaxis()->SetRangeUser(0.0, ymax_main * 1.35);
                }
                hMainVec[i]->Draw("EP");
            } else {
                hMainVec[i]->Draw("EP SAME");
            }
        }
        // Only draw legend on the rightmost pad
        if (axis == 2) leg->Draw("SAME");
        
        // Global Labels
        if (axis == 0) {
            texKin.SetTextAlign(11);
            texKin.DrawLatex(0.12, 0.92, "#bf{CMS} #it{Internal}");
            texKin.DrawLatex(0.55, 0.92, Form("#bf{p_{T} > %.0f GeV/c}", bins.ptmin));
        } else if (axis == 1) {
            texKin.SetTextAlign(22);
            texKin.DrawLatex(0.5, 0.92, Form("Run %s", runNumber.Data()));
        } else if (axis == 2) {
            texKin.SetTextAlign(31);
            texKin.DrawLatex(0.88, 0.92, Form("#bf{%s}", jetAlgo.Data()));
        }
        
        // --- BOTTOM ROW: RATIO ---
        cKin->cd(axis + 4);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.12);
        if (axis == 0) gPad->SetLogx(true); 
        
        // Dynamic, centered zoom logic bounded by [0, 2]
        double maxDeviation = std::max(std::abs(ymax_ratio - 1.0), std::abs(1.0 - ymin_ratio));
        double zoomLimit = std::min(maxDeviation * 1.2, 1.0); 
        
        double ratioRangeMin = std::max(0.0, 1.0 - zoomLimit);
        double ratioRangeMax = std::min(2.0, 1.0 + zoomLimit);
        
        for (size_t i = 0; i < hRatioVec.size(); ++i) {
            if (i == 0) {
                hRatioVec[0]->GetYaxis()->SetRangeUser(ratioRangeMin, ratioRangeMax);
                if (axis == 0) hRatioVec[0]->GetXaxis()->SetRangeUser(bins.ptmin, 500);
                hRatioVec[0]->Draw("EP");
                
                // Match the TLine bounds strictly to the axis plotting bounds
                double lineXmin = (axis == 0) ? bins.ptmin : hRatioVec[0]->GetXaxis()->GetXmin();
                double lineXmax = (axis == 0) ? 500.0 : hRatioVec[0]->GetXaxis()->GetXmax();
                TLine* line = new TLine(lineXmin, 1, lineXmax, 1);
                
                line->SetLineStyle(2);
                line->SetLineColor(kBlack);
                line->Draw();
            } else {
                hRatioVec[i]->Draw("EP SAME");
            }
        }
    }
    cKin->SaveAs(outDir + "/Canvas1_Kinematics.png");

    // ==============================================================================
    // CANVAS 2: ETA-PHI MAPS (Multiple 3x1)
    // ==============================================================================
    if(f25Data && !f25Data->IsZombie()) {
        THnSparseF* hnKin25Data = (THnSparseF*)f25Data->Get("hjetkin");
        THnSparseF* hnKin26MC   = (THnSparseF*)f26MC->Get("hjetkin");
        double nEvt25 = GetNevents(f25Data);
        double nEvtMC = GetNevents(f26MC);
        
        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(42);
        
        auto formatPadAndHist = [](TH2D* h) {
            gPad->SetTopMargin(0.18);
            gPad->SetRightMargin(0.16);
            gPad->SetBottomMargin(0.12);
            gPad->SetLeftMargin(0.12);
            
            h->GetXaxis()->SetTitleOffset(0.9);
            h->GetYaxis()->SetTitleOffset(0.9);
            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);
        };

        for (double ptC : ptCutsMap) {
            for (const auto& hb : bins.hiBins) { // Looping directly over BinningStruct's hiBins
                TString cName = Form("cMap_pt%.0f_hb%.0f_%.0f", ptC, hb.lo, hb.hi);
                TCanvas* cMap = new TCanvas(cName, "Eta-Phi Map", 2400, 1000);
                cMap->Divide(3, 1);
                
                TString hiBinLabel = Form("#bf{hiBin %.0f-%.0f}", hb.lo, hb.hi);
                TString ptLabel    = Form("#bf{p_{T} > %.0f GeV/c}", ptC);
                TString algoLabel  = Form("#bf{%s}", jetAlgo.Data());
                
                // --- Pad 1: 2026 Data ---
                cMap->cd(1);
                TH2D* h26 = ProjectTHn2D(hnKin26Data, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "map26");
                h26->Scale(1.0 / nEvt26Data);
                h26->SetTitle(";#eta;#phi (rad)");
                formatPadAndHist(h26);
                h26->Draw("colz");
                
                tex.SetTextAlign(22); tex.SetTextSize(0.055); 
                tex.DrawLatex(0.5, 0.87, Form("2026 Data (Run %s)", runNumber.Data()));
                tex.SetTextAlign(11); tex.SetTextSize(0.045);
                tex.DrawLatex(0.12, 0.94, "#bf{CMS} #it{Preliminary}");
                tex.SetTextAlign(31); 
                tex.DrawLatex(0.84, 0.94, ptLabel); 
                
                // --- Pad 2: 2025 Data ---
                cMap->cd(2);
                TH2D* h25 = ProjectTHn2D(hnKin25Data, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "map25");
                h25->Scale(1.0 / nEvt25);
                h25->SetTitle(";#eta;#phi (rad)");
                formatPadAndHist(h25);
                h25->Draw("colz");
                
                tex.SetTextAlign(22); tex.SetTextSize(0.055); 
                tex.DrawLatex(0.5, 0.87, "2025 Data");
                tex.SetTextAlign(22); tex.SetTextSize(0.045);
                tex.DrawLatex(0.5, 0.94, hiBinLabel); 
                
                // --- Pad 3: 2026 MC ---
                cMap->cd(3);
                TH2D* hMC = ProjectTHn2D(hnKin26MC, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, "mapMC");
                hMC->Scale(1.0 / nEvtMC);
                hMC->SetTitle(";#eta;#phi (rad)");
                formatPadAndHist(hMC);
                hMC->Draw("colz");
                
                tex.SetTextAlign(22); tex.SetTextSize(0.055); 
                tex.DrawLatex(0.5, 0.87, "2026 MC");
                tex.SetTextAlign(31); tex.SetTextSize(0.045); 
                tex.DrawLatex(0.84, 0.94, algoLabel); 
                
                cMap->SaveAs(outDir + "/" + cName + ".png");
                
                delete h26; delete h25; delete hMC; delete cMap;
            }
        }
    } else {
        std::cerr << "Skipping Canvas 2: Missing either 2025 Data or 2026 MC file." << std::endl;
    }

    // ==============================================================================
    // CANVAS 3: PF FRACTIONS (5x2)
    // ==============================================================================
    TCanvas* cPF = new TCanvas("cPF", "PF Fractions", 4000, 1600);
    cPF->Divide(5, 2);
    
    for (int pf = 0; pf < 5; ++pf) {
        cPF->cd(pf + 1);
        gPad->SetTopMargin(0.12);
        gPad->SetBottomMargin(0.12);
        gPad->SetLeftMargin(0.12);
        
        // Uses boolean logic from JetHealthPlotting.h
        if (PFTypeLogY(pf)) gPad->SetLogy();
        
        // Legend shifted slightly left to avoid plotting area cutoff
        TLegend* leg = new TLegend(0.30, 0.55, 0.80, 0.85);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.045);
        
        std::vector<TH1D*> hMainVec;
        std::vector<TH1D*> hRatioVec;
        
        double ymax_main  = 0.0;
        double ymax_ratio = -999.0;
        double ymin_ratio = 999.0;
        
        for (size_t i = 0; i < bins.hiBins.size(); ++i) {
            const auto& hb = bins.hiBins[i];
            
            TString suf26   = Form("_pf26_%d_hb%zu", pf, i);
            TString sufComp = Form("_pfComp_%d_hb%zu", pf, i);
            
            // Projecting PF Fractions
            TH1D* h26   = ProjectTHn1D(hnPF26Data, 0, {{1, (double)pf, (double)pf + 1.0}, {3, (double)hb.lo, (double)hb.hi}}, suf26);
            TH1D* hComp = ProjectTHn1D(hnPFComp,   0, {{1, (double)pf, (double)pf + 1.0}, {3, (double)hb.lo, (double)hb.hi}}, sufComp);
            
            NormalizeTH1(h26);
            NormalizeTH1(hComp);
            StyleTH1(h26, hb.color);
            h26->SetLineWidth(1);
            h26->SetMarkerSize(0.01);
            
            h26->SetTitle(""); 
            h26->GetYaxis()->SetTitle(""); 
            h26->GetXaxis()->SetTitle(PFTypeTitles.at(pf)); // Sourced from JetHealthHistograms.h
            h26->GetXaxis()->CenterTitle(true);
            h26->GetXaxis()->SetTitleSize(0.05);
            h26->GetXaxis()->SetLabelSize(0.05);
            h26->GetYaxis()->SetLabelSize(0.05);
            
            ymax_main = std::max(ymax_main, h26->GetMaximum());
            hMainVec.push_back(h26);
            leg->AddEntry(h26, hb.title, "lp");
            
            // Ratio
            TH1D* hRatio = (TH1D*)h26->Clone(Form("ratiopf_%d_hb%zu", pf, i));
            hRatio->Divide(hComp);
            StyleTH1(hRatio, hb.color);
            hRatio->SetLineWidth(1);
            hRatio->SetMarkerSize(0.01);
            
            hRatio->SetTitle("");
            hRatio->GetYaxis()->SetTitle(Form("Ratio to %s", compLabel.Data()));
            hRatio->GetYaxis()->CenterTitle(true);
            hRatio->GetYaxis()->SetTitleSize(0.05);
            hRatio->GetYaxis()->SetLabelSize(0.05);
            hRatio->GetYaxis()->SetTitleOffset(1.2);
            hRatio->GetYaxis()->SetNdivisions(505);
            hRatio->GetXaxis()->SetTitleSize(0.05);
            hRatio->GetXaxis()->SetLabelSize(0.05);
            
            for (int b = 1; b <= hRatio->GetNbinsX(); ++b) {
                double val = hRatio->GetBinContent(b);
                if (val > 0) {
                    ymax_ratio = std::max(ymax_ratio, val);
                    ymin_ratio = std::min(ymin_ratio, val);
                }
            }
            hRatioVec.push_back(hRatio);
            delete hComp;
        }
        
        cPF->cd(pf + 1);
        for (size_t i = 0; i < hMainVec.size(); ++i) {
            if (i == 0) {
                hMainVec[0]->GetYaxis()->SetRangeUser(PFTypeLogY(pf) ? 1e-6 : 0.0, ymax_main * (PFTypeLogY(pf) ? 5.0 : 1.35));
                hMainVec[i]->Draw("EP");
            } else {
                hMainVec[i]->Draw("EP SAME");
            }
        }
        // Only draw legend on the rightmost pad
        if (pf == 4) leg->Draw("SAME");
        
        // Spread global text across the 5 top pads
        if (pf == 0) {
            texKin.SetTextAlign(11);
            texKin.DrawLatex(0.12, 0.92, "#bf{CMS} #it{Internal}");
        } else if (pf == 2) {
            texKin.SetTextAlign(22);
            texKin.DrawLatex(0.5, 0.92, Form("Run %s", runNumber.Data()));
        } else if (pf == 4) {
            texKin.SetTextAlign(31);
            texKin.DrawLatex(0.88, 0.92, Form("#bf{%s}", jetAlgo.Data()));
        }
        
        // --- BOTTOM ROW: PF RATIO ---
        cPF->cd(pf + 6);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.12);
        
        double maxDeviation = std::max(std::abs(ymax_ratio - 1.0), std::abs(1.0 - ymin_ratio));
        double zoomLimit = std::min(maxDeviation * 1.2, 1.0); 
        double ratioRangeMin = std::max(0.0, 1.0 - zoomLimit);
        double ratioRangeMax = std::min(2.0, 1.0 + zoomLimit);
        
        for (size_t i = 0; i < hRatioVec.size(); ++i) {
            if (i == 0) {
                hRatioVec[0]->GetYaxis()->SetRangeUser(ratioRangeMin, ratioRangeMax);
                hRatioVec[0]->Draw("EP");
                
                TLine* line = new TLine(hRatioVec[0]->GetXaxis()->GetXmin(), 1, 
                                        hRatioVec[0]->GetXaxis()->GetXmax(), 1);
                line->SetLineStyle(2);
                line->SetLineColor(kBlack);
                line->Draw();
            } else {
                hRatioVec[i]->Draw("EP SAME");
            }
        }
    }
    cPF->SaveAs(outDir + "/Canvas3_PFFractions.png");

    std::cout << "All plots successfully saved to " << outDir << "/" << std::endl;
}

/*
==============================================================================
MACRO DOCUMENTATION
==============================================================================
What this code does:
This macro reads output ROOT files generated by a jet analysis framework containing
THnSparse histograms for jet kinematics and Particle Flow (PF) fractions. It 
projects these multi-dimensional histograms into 1D and 2D distributions, applies
necessary normalizations, and computes Data vs. MC (or Data vs. Data) ratios. 
It generates three main sets of canvases:
1. 1D Jet Kinematics (pT, eta, phi) split by centrality (hiBin).
2. 2D Eta-Phi heatmaps of jet yields for various pT and hiBin thresholds.
3. 1D PF Fractions (CHF, NHF, CEF, NEF, MUF) split by centrality.

Cuts Applied & How to Change Them:
- The pT cut is globally defined by `ptCutsMap` for the 2D maps, and the minimum
  pT bounds are governed by `bins.ptmin` (default 50.0 GeV/c) instantiated from 
  the BinningStruct header.
- Centrality (hiBin) cuts are iteratively applied via the BinningStruct `hiBins`.
  To change these centrality cuts or their associative colors, modify the 
  `Binning.h` header directly.
- To toggle the ratio reference between 2026 MC and 2025 Data, simply flip the 
  `useMCForRatio` boolean at the top of the file.

Normalizations Used:
- Canvas 1 & 3 (1D plots): Area normalized. The main distributions are scaled by 
  `1.0 / Integral()` via `NormalizeTH1()` to directly compare shapes between 
  Data and MC (or Data vs. Data).
- Canvas 2 (2D maps): Event normalized. Distributions are scaled by `1.0 / N_events`
  where N_events is extracted from the integral of the corresponding `hvz` vertex 
  histogram, providing a true per-event yield representation.
==============================================================================
*/