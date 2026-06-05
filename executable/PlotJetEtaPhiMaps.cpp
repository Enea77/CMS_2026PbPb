#include <iostream>
#include <vector>
#include <utility>

#include "TFile.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"

// Bring in your custom analysis framework headers
#include "../header/Binning.h"
#include "../header/JetHealthHistograms.h"
#include "../header/JetHealthPlotting.h"

const std::vector<double> ptCutsMap = {50.0, 100.0};

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

void PlotJetEtaPhiMaps(TString file1Path = "JetHealth_404350", TString label1="PF Jets", 
                       TString file2Path = "JetHealth_404350_ak4PF", TString label2="Unsubtracted PF Jets") {

    if (false){
        file1Path = "JetHealth_2026_MC"; label1="PF Jets";
        file2Path = "JetHealth_2026MC_ak4PF"; label2="Unsubtracted PF Jets";
    }  

    if (true){
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

    THnSparseF* hnKin1 = (THnSparseF*)f1->Get("hjetkin");
    THnSparseF* hnKin2 = (THnSparseF*)f2->Get("hjetkin");
    
    BinningStruct bins(50.0);

    // SAFETY: Clear axis ranges in case previous macros left the THnSparse axes restricted
    for (int i = 0; i < hnKin1->GetNdimensions(); ++i) {
        hnKin1->GetAxis(i)->SetRange(0, -1);
        hnKin2->GetAxis(i)->SetRange(0, -1);
    }

    for (double ptC : ptCutsMap) {
        for (const auto& hb : bins.hiBins) { 
            TString cName = Form("cMap_pt%.0f_hb%.0f_%.0f", ptC, hb.lo, hb.hi);
            TCanvas* cMap = new TCanvas(cName, "Eta-Phi Map", 1600, 800); cMap->Divide(2, 1);
            
            cMap->cd(1); FormatSlidePad(false, false, true);
            
            // 1. Create a truly unique name for this specific iteration
            TString n1 = Form("map1_pt%.0f_hb%.0f_%.0f", ptC, hb.lo, hb.hi);
            
            // 2. Pass that unique name into your new robust ProjectTHn2D function
            // Note: Indices (1, 2) for Eta/Phi and (0, 3) for pT/hiBin remain correct for hjetkin!
            TH2D* h1 = ProjectTHn2D(hnKin1, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, n1);
            
            h1->SetTitle(";#eta;#phi (rad)");
            h1->GetXaxis()->SetRangeUser(-1.4999,1.4999);
            h1->GetXaxis()->SetTitleOffset(0.85); h1->GetYaxis()->SetTitleOffset(0.85); h1->Draw("colz");
            DrawSlideText(1, true, label1, ptC, "Run 404350");
            
            cMap->cd(2); FormatSlidePad(false, false, true);
            
            // 1. Unique name for the second file's map
            TString n2 = Form("map2_pt%.0f_hb%.0f_%.0f", ptC, hb.lo, hb.hi);
            
            // 2. Project
            TH2D* h2 = ProjectTHn2D(hnKin2, 1, 2, {{0, ptC, 1000.0}, {3, (double)hb.lo, (double)hb.hi}}, n2);
            
            h2->SetTitle(";#eta;#phi (rad)");
            h2->GetXaxis()->SetRangeUser(-1.4999,1.4999);
            h2->GetXaxis()->SetTitleOffset(0.85); h2->GetYaxis()->SetTitleOffset(0.85); h2->Draw("colz");
            DrawSlideText(2, true, label2, ptC, "Run 404350", Form("#bf{MC 2026 hiBin %.0f-%.0f}", hb.lo, hb.hi));
            
            cMap->SaveAs(outDir + "/" + cName + ".png");
            
            // Cleanup memory
            delete h1; delete h2; delete cMap;
        }
    }
    std::cout << "Eta-Phi maps successfully saved to " << outDir << "/" << std::endl;
}