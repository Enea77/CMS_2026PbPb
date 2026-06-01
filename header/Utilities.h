#ifndef UTILITIES_H
#define UTILITIES_H

// Utility functions for histogram creation, styling, and plotting

#include "TH1F.h"
#include "TH1I.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TString.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <initializer_list>
#include <vector>

#include "Binning.h"

// Histograms //

template<typename T>
void WriteAll(T* h){if(h)h->Write();}
 
template<typename T>
void WriteAll(std::vector<T>& v){for(auto& entry : v) WriteAll(entry);}

// THnSparse

template<typename T>
T* MakeTHnSparse(const TString& name, const TString& title, std::initializer_list<AxisBins> axes){
    const std::vector<AxisBins> axVec(axes);
    const Int_t ndim = axVec.size();
    std::vector<Int_t> nbins;
    std::vector<Double_t> xmin, xmax;
    for (const auto& ax : axVec) {
        nbins.push_back(ax.nBins);
        xmin.push_back(ax.lo);
        xmax.push_back(ax.hi);
    }
    T* h = new T(name, title, ndim, nbins.data(), xmin.data(), xmax.data());
    for(Int_t i=0; i<ndim; i++){h->GetAxis(i)->SetTitle(axVec[i].title);}
    return h;
}

struct SparseRange{
    Int_t axis;
    Double_t lo;
    Double_t hi;
};

inline void SetAxisRange(THnSparse* h, Int_t axis, Double_t lo, Double_t hi){h->GetAxis(axis)->SetRangeUser(lo, hi);}
inline void ResetAxisRange(THnSparse* h, Int_t axis){h->GetAxis(axis)->SetRange(0, -1);}

// 1D Projection
inline TH1D* ProjectTHn1D(THnSparse* h, Int_t axis, const std::vector<SparseRange>& ranges, const TString& newName) {
    for (const auto& r : ranges) { h->GetAxis(r.axis)->SetRangeUser(r.lo, r.hi - 0.0001); }
    TH1D* tmp = h->Projection(axis, "E"); 
    TH1D* out = (TH1D*)tmp->Clone(newName);
    out->SetTitle("");
    delete tmp; 
    for (const auto& r : ranges) { h->GetAxis(r.axis)->SetRange(0, -1); }
    return out;
}

// 2D Projection
inline TH2D* ProjectTHn2D(THnSparse* h, Int_t xAxis, Int_t yAxis, const std::vector<SparseRange>& ranges, const TString& newName) {
    for (const auto& r : ranges) { h->GetAxis(r.axis)->SetRangeUser(r.lo, r.hi - 0.0001); }
    TH2D* tmp = h->Projection(yAxis, xAxis, "E"); 
    TH2D* out = (TH2D*)tmp->Clone(newName);
    out->SetTitle("");
    delete tmp; 
    for (const auto& r : ranges) { h->GetAxis(r.axis)->SetRange(0, -1); }
    return out;
}

// TH1

template<typename T>
T* MakeTH1(const TString& hname, const AxisBins& bins){
    T* h = new T(hname, hname, bins.nBins, bins.lo, bins.hi);
    h->SetTitle(bins.title);
    return h;
}

inline void NormalizeTH1(TH1* h){
    const double integral = h->Integral();
    if(integral > 0.0){h->Scale(1.0 / integral);}
}

// Plotting //

struct PlotConfig{
    TString runNumber = "404338";
    TString globalTag = "";
    TString jetAlgo   = "";
    
    float xmin  = 20.0;
    float xmax  = 200.0;
    float ymax = 1.1;
    float ymin = 0.0;

    int canvasSize = 2400;
    float imageScaling = 1.0;

    float legfrac = 0.2;
};

// Canvas

inline TCanvas* MakeSinglePadCanvas(const TString& name, const PlotConfig& cfg, bool grid = false){
    TCanvas* c = new TCanvas(name, "", cfg.canvasSize, cfg.canvasSize);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetRealAspectRatio();
    if(grid) {
        c->SetGridx();
        c->SetGridy();
    }
    return c;
}

inline TCanvas* MakeColzCanvas(const TString& name, const PlotConfig& cfg){
    TCanvas* c = new TCanvas(name, "", cfg.canvasSize, cfg.canvasSize);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.12);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.15);
    c->SetRealAspectRatio();
    return c;
}

struct SplitCanvas{
    TCanvas* c = nullptr;
    TPad* plotpad = nullptr;
    TPad* legpad = nullptr;
};

inline SplitCanvas MakeSplitPadCanvas(const TString& name, const PlotConfig& cfg){
    SplitCanvas sc;
    sc.c = new TCanvas(name, "", cfg.canvasSize, cfg.canvasSize);
    sc.c->cd();

    sc.legpad = new TPad(name + "_legpad", "", 0.0, 1.0 - cfg.legfrac, 1.0, 1.0);
    sc.legpad->SetTopMargin(0.05);
    sc.legpad->SetBottomMargin(0.0);
    sc.legpad->SetLeftMargin(0.0);
    sc.legpad->SetRightMargin(0.0);
    sc.legpad->SetFillStyle(0);
    sc.legpad->Draw();

    sc.plotpad = new TPad(name + "_plotpad", "", 0.0, 0.0, 1.0, 1.0 - cfg.legfrac);
    sc.plotpad->SetTopMargin(0.02);
    sc.plotpad->SetBottomMargin(0.12);
    sc.plotpad->SetLeftMargin(0.12);
    sc.plotpad->SetRightMargin(0.05);
    sc.plotpad->SetGridx();
    sc.plotpad->SetGridy();
    sc.plotpad->Draw();

    return sc;
}

// Legend

inline TLegend* MakeLegend(float xmin = 0.575, float ymin = 0.75, float xmax = 0.8, float ymax = 0.95){
    TLegend* l = new TLegend(xmin, ymin, xmax, ymax);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(0.030);
    return l;
}

inline void AddInfoEntries(TLegend* l, const PlotConfig& cfg){
    if(!cfg.runNumber.IsNull()){l->AddEntry((TObject*)nullptr, cfg.runNumber,  "");}
    if(!cfg.globalTag.IsNull()){l->AddEntry((TObject*)nullptr, cfg.globalTag,  "");}
    if(!cfg.jetAlgo.IsNull()){l->AddEntry((TObject*)nullptr, cfg.jetAlgo,    "");}
}

inline TString MakePlotDir(const TString& prefix = "plots"){
    TDatime now;
    TString timestamp = Form("%04d-%02d-%02d_%02d-%02d-%02d", now.GetYear(), now.GetMonth(), now.GetDay(), now.GetHour(), now.GetMinute(), now.GetSecond());
    TString dir = prefix + "_" + timestamp;
    gSystem->mkdir(dir, true);
    return dir;
}

// Style

inline void SetPlotStyle(const PlotConfig& cfg){
    gStyle->SetImageScaling(cfg.imageScaling);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.04, "XY");
    gStyle->SetLabelSize(0.04, "XY");
    gStyle->SetGridColor(kGray+2);
    gStyle->SetGridWidth(1);
    gStyle->SetGridStyle(3);
    gStyle->SetPalette(kBird);
}

inline void StyleTH1(TH1* h, Color_t color){
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(4);
}

// TLatex

inline void DrawLabel(const TString& text, float x, float y, float textSize = 0.035, Int_t align = 11){
    if(text.IsNull()){return;}
    TLatex* tex = new TLatex(x, y, text);
    tex->SetNDC();
    tex->SetTextAlign(align);
    tex->SetTextSize(textSize);
    tex->Draw();
}

inline void DrawCMSLabel(float x = 0.12, float y = 0.965, float textSize = 0.035) {
    DrawLabel("CMS #bf{#it{Internal}}", x, y, textSize);
}

inline void DrawJetAlgoLabel(const TString& jetAlgo, float x = 0.8, float y = 0.04, float textSize = 0.035){
    DrawLabel(Form("#bf{%s}", jetAlgo.Data()), x, y, textSize);
}


#endif