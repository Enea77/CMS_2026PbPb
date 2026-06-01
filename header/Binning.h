#ifndef BINNING_H
#define BINNING_H

#include "TString.h"
#include "Rtypes.h"

#include "MyColors.h"

#include <vector>

struct AxisBins{
    Int_t nBins;
    Float_t lo;
    Float_t hi;
    TString title = "";
};

struct RangeBin{
    Float_t lo;
    Float_t hi;
    TString title;
    TString shortName;
    Color_t color = kBlack;
};

struct BinningStruct{
    Float_t ptmin;

    // histogram axes
    AxisBins vz = {40, -20.0, 20.0, "v_{z} (cm)"};
    AxisBins pt = {50, 50.0, 1000.0, "p_{T} (GeV/c)"};
    AxisBins eta = {102, -5.1, 5.1, "#eta"};
    AxisBins abseta = {51, 0, 5.1, "|#eta|"}; 
    AxisBins phi = {64, -3.2, 3.2, "#phi (rad)"};
    AxisBins nref = {30, 0, 30, "nref"};
    AxisBins trig = {2, 0, 2, "trigger decision"};
    AxisBins pfFrac = {20, 0.0, 1.0, "Particle FLow energy fraction"};
    AxisBins pfMult = {50, 0, 50, "Particle Flow candidate multiplicity"};
    AxisBins hiBin = {200, 0, 200, "hiBin"};

    // centrality
    std::vector<RangeBin> hiBins = {
        {0, 200, "Inclusive", "_hiBin_0_to_200", KlimtPink},
        {0, 20, "0% < Centrality < 10%", "_hiBin_0_to_20", KlimtRed},
        {20, 60, "10% < Centrality < 30%", "_hiBin_20_to_60", KlimtYellow},
        {60, 100, "30% < Centrality < 50%", "_hiBin_60_to_100", KlimtGreen},
        {100, 160, "50% < Centrality < 80%", "_hiBin_100_to_160", KlimtBlue},
        {160, 200, "80% < Centrality < 100%", "_hiBin_160_to_200", KlimtPurple}
    };

    // eta
    std::vector<RangeBin> etaBins = {
        {0.0, 1.2, "|#eta| < 1.2", "_eta_0p0_to_1p2"},
        {0.0, 2.4, "|#eta| < 2.4", "_eta_0p0_to_2p4"},
        {0.0, 5.1, "|#eta| < 5.1", "_eta_0p0_to_5p1"},
        {2.7, 5.1, "2.7 < |#eta| < 5.1", "_eta_2p7_to_5p1"}
    };

    BinningStruct(Float_t ptcut = 50.0) : ptmin(ptcut), pt({50, ptcut, 1000.0, "p_{T} (GeV/c)"}){}
};

#endif
