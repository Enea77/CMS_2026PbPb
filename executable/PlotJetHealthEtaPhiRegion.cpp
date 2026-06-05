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
            tex.DrawLatex(0.85, 0.94, extraLabel); 
        }
    } else if (padNum == 3) {
        tex.SetTextAlign(31); 
        tex.DrawLatex(rightX, 0.94, Form("#bf{%s}", jetAlgo.Data())); 
    }
}

// ==============================================================================
// MAIN PLOTTING MACRO
// ==============================================================================

void PlotJetHealthEtaPhiRegion(TString file1Path = "JetHealth_404350", TString label1="FPIX out (Run 404350)", 
                            TString file2Path = "JetHealth_404469_RPmany", TString label2="FPIX in (Run 404469)") {

    TString extratext = "akCs4PF";

    if (false){
        file1Path = "JetHealth_2026_MC"; label1="2026 MC";
        file2Path = "JetHealth_2026_MC_maskFPIX"; label2="2026 MC FPIX-masked";
    }   
    else if (true){
        file1Path = "JetHealth_404469_RPmany"; label1="PF Jets";
        file2Path = "JetHealth_404469_RPmany_ak4PF"; label2="Unsubtracted PF Jets";
        extratext = "FPIX in (Run 404469)";
    }  
    
    //file2Path = "JetHealth_404350_CaloJets_try2"; label2="Calo Jets"; //label2="Run404350 akPu4Calo";
    
    //file1Path = "JetHealth_2025_Data_new"; label1="2025 Data";
      
    //file2Path = "JetHealth_404350_PromptRECO"; label2="Run 404350 Prompt RECO";

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRedBlue);
    
    TString outDir = file1Path + "_vs_" + file2Path;
    gSystem->mkdir(outDir, true);

    TFile* f1 = TFile::Open(file1Path+".root", "READ");
    TFile* f2 = TFile::Open(file2Path+".root", "READ");
    
    if (!f1 || f1->IsZombie()) { std::cerr << "Cannot open file 1: " << file1Path << std::endl; return; }
    if (!f2 || f2->IsZombie()) { std::cerr << "Cannot open file 2: " << file2Path << std::endl; return; }
    
    double nEvt1 = GetNevents(f1);
    double nEvt2 = GetNevents(f2);

    THnSparseF* hnKin1 = (THnSparseF*)f1->Get("hjetkin");
    THnSparseF* hnKin2 = (THnSparseF*)f2->Get("hjetkin");
    THnSparseF* hnPF1  = (THnSparseF*)f1->Get("hjetpf");
    THnSparseF* hnPF2  = (THnSparseF*)f2->Get("hjetpf");

    BinningStruct bins(50.0);

    // Define the eta regions to loop over
    std::vector<double> etaCuts = {1.5, 2.0, 3.0, 5.1};

    for (double etaCut : etaCuts) {
        for (double ptC : ptCutsMap) {
            for (const auto& hb : bins.hiBins) { 
                // Include eta in the canvas name to prevent files from overwriting each other
                TString cName = Form("cMap_eta%.1f_pt%.0f_hb%.0f_%.0f", etaCut, ptC, hb.lo, hb.hi);
                TCanvas* cMap = new TCanvas(cName, "Eta-Phi Map", 1600, 800); cMap->Divide(2, 1);
                
                // 1. Project both histograms first (using unique internal names for ROOT memory safety)
                TH2D* h1 = ProjectTHn2D(hnKin1, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, Form("map1_eta%.1f", etaCut));
                //h1->Scale(1.0 / nEvt1); 
                h1->Scale(1.0/h1->GetEntries());

                TH2D* h2 = ProjectTHn2D(hnKin2, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, Form("map2_eta%.1f", etaCut));
                //h2->Scale(1.0 / nEvt2); 
                h2->Scale(1.0/h2->GetEntries());

                h1->GetXaxis()->SetRangeUser(-etaCut + 0.0001, etaCut - 0.0001);
                h2->GetXaxis()->SetRangeUser(-etaCut + 0.0001, etaCut - 0.0001);
                               
                // =================================================================
                // OPTION: Match Z-axis scales (Comment out this block to decouple them)
                double maxZ = std::max(h1->GetMaximum(), h2->GetMaximum());
                double minZ = std::min(h1->GetMinimum(), h2->GetMinimum());
                h1->SetMaximum(maxZ); h1->SetMinimum(minZ);
                h2->SetMaximum(maxZ); h2->SetMinimum(minZ);
                // =================================================================
                
                // 2. Format and Draw Pad 1
                cMap->cd(1); FormatSlidePad(false, false, true);
                h1->SetTitle(";#eta;#phi (rad)");
                // Dynamically scales the range (subtracting 0.0001 keeps your original 1.5 -> 1.4999 logic)
                h1->GetXaxis()->SetTitleOffset(0.85); h1->GetYaxis()->SetTitleOffset(0.85); h1->Draw("colz");
                DrawSlideText(1, true, label1, ptC, "Run 404350");
                
                // 3. Format and Draw Pad 2
                cMap->cd(2); FormatSlidePad(false, false, true);
                h2->SetTitle(";#eta;#phi (rad)");
                h2->GetXaxis()->SetTitleOffset(0.85); h2->GetYaxis()->SetTitleOffset(0.85); h2->Draw("colz");
                DrawSlideText(2, true, label2, ptC, "Run 404350", Form("#bf{%s hiBin %.0f-%.0f}",extratext.Data(), hb.lo, hb.hi));

                // 4. Save and Clean up
                cMap->SaveAs(outDir + "/" + cName + ".png");
                delete h1; delete h2; delete cMap;
            }
        }
    }
}