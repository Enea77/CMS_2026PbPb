#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
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
            tex.DrawLatex(0.6, 0.94, extraLabel); 
        }
    } else if (padNum == 3) {
        tex.SetTextAlign(31); 
        tex.DrawLatex(rightX, 0.94, Form("#bf{%s}", jetAlgo.Data())); 
    }
}

void PlotJetFractions(TString file1Path = "JetHealth_404350", TString label1="Run 404350", 
                      TString file2Path = "JetHealth_404350_ak4PF", TString label2="Unsubtracted PF Jets") {

    if (false){
        file1Path = "JetHealth_2026_MC"; label1="PF Jets";
        file2Path = "JetHealth_2026MC_ak4PF"; label2="Unsubtracted PF Jets";
    }  

    if (true){
        file2Path = "JetHealth_2025_Data_new"; label2="2025 Data";
    }  

    gStyle->SetOptStat(0);
    
    TString outDir = file1Path + "_vs_" + file2Path;
    gSystem->mkdir(outDir, true);

    TFile* f1 = TFile::Open(file1Path+".root", "READ");
    TFile* f2 = TFile::Open(file2Path+".root", "READ");
    
    if (!f1 || f1->IsZombie()) { std::cerr << "Cannot open file 1" << std::endl; return; }
    if (!f2 || f2->IsZombie()) { std::cerr << "Cannot open file 2" << std::endl; return; }

    THnSparseF* hnPF1  = (THnSparseF*)f1->Get("hjetpf");
    THnSparseF* hnPF2  = (THnSparseF*)f2->Get("hjetpf");

    for(int i=0; i<hnPF1->GetNdimensions(); ++i) std::cout << "Axis " << i << ": " << hnPF1->GetAxis(i)->GetTitle() << std::endl;
    
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

    auto GetPFVectors = [&](int pf, double etaLo, double etaHi, const TString& xTitle) {
        std::vector<TH1D*> v1, v2, vRat;
        for (const auto& hb : bins.hiBins) {
            // Give fully unique names to avoid directory collisions
            TString n1 = Form("h1_pf%d_eta%.1f_%.1f_hb%.0f_%.0f", pf, etaLo, etaHi, hb.lo, hb.hi);
            TString n2 = Form("h2_pf%d_eta%.1f_%.1f_hb%.0f_%.0f", pf, etaLo, etaHi, hb.lo, hb.hi);
            
            TH1D* h1 = ProjectTHn1D(hnPF1, 0, {{1, (double)pf, (double)pf + 1.0}, {2, etaLo, etaHi}, {5, (double)hb.lo, (double)hb.hi}}, n1);
            TH1D* h2 = ProjectTHn1D(hnPF2, 0, {{1, (double)pf, (double)pf + 1.0}, {2, etaLo, etaHi}, {5, (double)hb.lo, (double)hb.hi}}, n2);

            NormalizeTH1(h1); NormalizeTH1(h2);
            StyleTH1Slide(h1, hb.color, xTitle); StyleTH1Slide(h2, hb.color, xTitle);
            
            TString nRat = Form("hrat_pf%d_eta%.1f_%.1f_hb%.0f_%.0f", pf, etaLo, etaHi, hb.lo, hb.hi);
            TH1D* hRatio = (TH1D*)h1->Clone(nRat); hRatio->Divide(h2);
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
        
        if (cName.Contains("cPFVsEta")) {
            std::vector<double> targetEtas = {-2.3, -1.6, 1.5, 2.1};
            double yminMain = isLogY ? 1e-6 : 0.0;
            double ymaxMain = ymax_main * (isLogY ? 5.0 : 1.35);

            for (double x : targetEtas) {
                TLine* line = new TLine();
                line->SetLineStyle(3); 
                line->SetLineColor(kBlack);
                line->SetLineWidth(2);
                
                cSlide->cd(1); line->DrawLine(x, yminMain, x, ymaxMain);
                cSlide->cd(2); line->DrawLine(x, yminMain, x, ymaxMain);
                cSlide->cd(3); line->DrawLine(x, ymin_ratio, x, ymax_ratio);
            }
        }
        
        cSlide->SaveAs(outDir + "/" + cName + ".png"); 
        delete cSlide;
        for(auto h: v1) delete h; for(auto h: v2) delete h; for(auto h: vRat) delete h;
    }; 

    // SLIDES: PF FRACTIONS 
    /*
    for (const auto& er : etaRanges) {
        for (int pf = 0; pf < 5; ++pf) {
            auto [v1, v2, vR] = GetPFVectors(pf, er.lo, er.hi, PFTypeTitles.at(pf));
            DrawAndSaveSlide(Form("cPFSlide_%s_%s", PFTypeNames.at(pf), er.name.Data()), PFTypeTitles.at(pf), 
                             v1, v2, vR, bins.ptmin, er.label, PFTypeLogY(pf), false);
        }
    }
        */
    
    // SAFETY: Clear axis ranges in case previous loops left the THnSparse axes restricted 
    // (This prevents caching bugs inside ROOT's projection engine)
    for (int i = 0; i < hnPF1->GetNdimensions(); ++i) {
        hnPF1->GetAxis(i)->SetRange(0, -1);
        hnPF2->GetAxis(i)->SetRange(0, -1);
    }

    // TASK 2: All PF Fractions vs Eta (3x1 Canvas)
    for (int pf = 0; pf < 5; ++pf) {
        std::vector<TH1D*> v1_eta, v2_eta, vRat_eta;
        
        for (const auto& hb : bins.hiBins) {
            
            // Generate COMPLETELY unique names using both hb.lo AND hb.hi so ROOT 
            // doesn't clobber the pointers from the 0-20 and 0-200 bins
            TString n2D_1 = Form("pf_2d_1_%.0f_%.0f_type%d", hb.lo, hb.hi, pf);
            TString n2D_2 = Form("pf_2d_2_%.0f_%.0f_type%d", hb.lo, hb.hi, pf);

            // Added {4, bins.ptmin, 1000.0} to apply the pT > 50 cut
            TH2D* h2D_1 = ProjectTHn2D(hnPF1, 2, 0, {{1, (double)pf, (double)pf + 1.0}, {4, bins.ptmin, 1000.0}, {5, (double)hb.lo, (double)hb.hi}}, n2D_1);
            TH2D* h2D_2 = ProjectTHn2D(hnPF2, 2, 0, {{1, (double)pf, (double)pf + 1.0}, {4, bins.ptmin, 1000.0}, {5, (double)hb.lo, (double)hb.hi}}, n2D_2);

            TString nProf_1 = Form("prof1_%.0f_%.0f_type%d", hb.lo, hb.hi, pf);
            TString nProf_2 = Form("prof2_%.0f_%.0f_type%d", hb.lo, hb.hi, pf);

            TProfile* p1 = h2D_1->ProfileX(nProf_1);
            TProfile* p2 = h2D_2->ProfileX(nProf_2);
            
            TString n1 = Form("hpf_eta_1_%.0f_%.0f_type%d", hb.lo, hb.hi, pf);
            TString n2 = Form("hpf_eta_2_%.0f_%.0f_type%d", hb.lo, hb.hi, pf);

            TH1D* h1 = p1->ProjectionX(n1);
            TH1D* h2 = p2->ProjectionX(n2);

            if (pf == 0) printf("hiLow %.1f hiHi %.1f entries %.2f mean %.2f \n",(double)hb.lo, (double)hb.hi, h1->GetEntries(), h1->GetMean());
            
            h1->GetXaxis()->SetRangeUser(-3.0, 3.0);
            h2->GetXaxis()->SetRangeUser(-3.0, 3.0);
            
            StyleTH1Slide(h1, hb.color, "#eta");
            StyleTH1Slide(h2, hb.color, "#eta");
            
            TString nRat = Form("hrat_%.0f_%.0f_type%d", hb.lo, hb.hi, pf);
            TH1D* hRatio = (TH1D*)h1->Clone(nRat);
            hRatio->Divide(h2);
            
            StyleTH1Slide(hRatio, hb.color, "#eta");
            hRatio->GetYaxis()->CenterTitle(true); hRatio->GetYaxis()->SetTitleSize(0.05); 
            hRatio->GetYaxis()->SetTitleOffset(0.85); hRatio->GetYaxis()->SetNdivisions(505);
            hRatio->GetXaxis()->SetRangeUser(-3.0, 3.0);
            
            v1_eta.push_back(h1); v2_eta.push_back(h2); vRat_eta.push_back(hRatio);
            
            delete h2D_1; delete h2D_2; delete p1; delete p2;
        }
        
        TString canvasName = TString("cPFVsEta_") + PFTypeNames.at(pf);
        TString canvasTitle = TString(PFTypeTitles.at(pf)) + " vs #eta";
        DrawAndSaveSlide(canvasName, canvasTitle, v1_eta, v2_eta, vRat_eta, bins.ptmin, PFTypeTitles.at(pf), false, false);
    }

    std::cout << "PF Fraction plots successfully saved to " << outDir << "/" << std::endl;
}