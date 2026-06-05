#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "THnSparse.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"

// Bring in your analysis framework headers to handle translation mappings
#include "../header/Binning.h"
#include "../header/JetSelection_PbPb.h"
#include "../header/JetStruct.h"
#include "../header/EventStructs_PbPb.h" 

// Match the array limits defined in your subsystem configuration
static constexpr Int_t local_maxnref = 999;

void MakePFSparse(const std::string& inputListPath = "files_404350.txt", const std::string& outputFilePath="PFcandidates_Run404350.root", bool isMC = false) {
    
    // 1. Initialize Synchronized Chains across all 5 framework modules
    TChain* chainPF   = new TChain("particleFlowAnalyser/pftree");
    TChain* chainJet  = new TChain("akCs4PFJetAnalyzer/t");
    TChain* chainEvt  = new TChain("hiEvtAnalyzer/HiTree");
    TChain* chainSkim = new TChain("skimanalysis/HltTree");
    TChain* chainHLT  = new TChain("hltanalysis/HltTree");
    
    // 2. Read the file list text file
    std::ifstream fileList(inputListPath);
    std::string rootFileName;
    if (!fileList.is_open()) {
        std::cerr << "Error: Cannot open input file list: " << inputListPath << std::endl;
        return;
    }
    
    std::cout << "Chaining target files into synchronized frameworks..." << std::endl;
    int fileCount = 0;
    while (std::getline(fileList, rootFileName)) {
        if (!rootFileName.empty()) {
            chainPF->Add(rootFileName.c_str());
            chainJet->Add(rootFileName.c_str());
            chainEvt->Add(rootFileName.c_str());
            chainSkim->Add(rootFileName.c_str());
            chainHLT->Add(rootFileName.c_str());
            fileCount++;
        }
    }
    std::cout << "Successfully tied " << fileCount << " files together." << std::endl;

    // 3. Set up standard PF candidate branch pointers
    std::vector<float>* pfPt  = nullptr;
    std::vector<float>* pfEta = nullptr;
    std::vector<float>* pfPhi = nullptr;
    std::vector<int>* pfId  = nullptr;
    chainPF->SetBranchAddress("pfPt",  &pfPt);
    chainPF->SetBranchAddress("pfEta", &pfEta);
    chainPF->SetBranchAddress("pfPhi", &pfPhi);
    chainPF->SetBranchAddress("pfId",  &pfId);
    
    // 4. Instantiate framework classes for Event, Skim Filters, and Jet Arrays
    EventStruct evt;
    FiltersStruct fltr;
    JetStruct<local_maxnref> jt;
    JetSelect js; // Loads default JECDatabase veto map automatically

    // 5. Map internal tracking fields to true branch strings using framework maps
    for (const auto& bMap : evt.BranchMap(isMC)) {
        chainEvt->SetBranchAddress(bMap.first, bMap.second);
    }
    for (const auto& bMap : fltr.BranchMap()) {
        chainSkim->SetBranchAddress(bMap.first, bMap.second);
    }
    for (const auto& bMap : jt.BranchMap(isMC)) {
        chainJet->SetBranchAddress(bMap.first, bMap.second);
    }

    // 6. Bind Minimum Bias hardware trigger decision
    Int_t L1minBias = 0;
    chainHLT->SetBranchStatus("*", 0);
    chainHLT->SetBranchStatus("L1_MinimumBiasHF1_AND_BptxAND", 1);
    chainHLT->SetBranchAddress("L1_MinimumBiasHF1_AND_BptxAND", &L1minBias);

    // 7. Define 10D THnSparse axes and phase-space bin boundaries
    const int ndims = 10; 
    int bins[ndims]     = { 100,  100,  64,             8,     100,   100,   64,            50,   50,  200   };
    double xmin[ndims]  = { 0.0, -5.0, -TMath::Pi(),   -0.5,   0.0,  -5.0,  -TMath::Pi(),   0.0,  0.0,   0.0  };
    double xmax[ndims]  = { 50.0, 5.0,  TMath::Pi(),    7.5,   500.0, 5.0,   TMath::Pi(),   0.4,  1.0, 200.0  };
    
    THnSparseF* hnPF = new THnSparseF("hnPF", 
        "10D Centrality-Correlated PF Sparse;p_{T}^{PF};#eta^{PF};#phi^{PF};ID;p_{T}^{Jet};#eta^{Jet};#phi^{Jet};#DeltaR;z;hiBin", 
        ndims, bins, xmin, xmax);
        
    hnPF->Sumw2();

    // 8. Synchronized Execution Loop
    Long64_t nEntries = chainPF->GetEntries();
    std::cout << "Starting entry processing loop over " << nEntries << " events..." << std::endl;
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        chainPF->GetEntry(i);
        chainJet->GetEntry(i);
        chainEvt->GetEntry(i);
        chainSkim->GetEntry(i);
        chainHLT->GetEntry(i);
        
        if (i % 100000 == 0) {
            std::cout << "Processed " << i << " / " << nEntries << " entries..." << std::endl;
        }

        // Apply primary event tracking vertex limits and skim filtration flags
        if ((fltr.ppvF == 0) || (fltr.pclustF == 0) || (fltr.pphfF == 0) || (TMath::Abs(evt.vz) > 15.0)) {
            continue;
        }
        
        // Skip events failing the Minimum Bias trigger condition or containing 0 valid jets
        if ((L1minBias == 0) || (jt.reco.nref == 0)) {
            continue;
        }

        if (!pfPt || !pfEta || !pfPhi || !pfId) continue;
        
        size_t nPF = pfPt->size();
        
        // Loop over individual candidate configurations
        for (size_t j = 0; j < nPF; ++j) {
            double pEta = (*pfEta)[j];
            double pPhi = (*pfPhi)[j];
            double pPt  = (*pfPt)[j];
            
            int matchedJetIdx = -1;
            double minDeltaR = 0.4; 
            
            // Geometrical matching to authorized jet cores
            for (int k = 0; k < jt.reco.nref; ++k) {
                if (k >= local_maxnref) break; 

                // Skip background fluctuations or jets failing loose baseline parameters
                if (jt.reco.pt[k] < 10.0) continue; 
                
                // Validate via your custom tight Jet Identification criteria and Veto Maps
                if (!js.JetSelection(jt.reco.eta[k], jt.reco.phi[k], jt.reco.pf.CEF[k], jt.reco.pf.NEF[k], jt.reco.pf.MUF[k])) {
                    continue;
                }

                double deta = pEta - jt.reco.eta[k];
                double dphi = pPhi - jt.reco.phi[k];
                
                while (dphi > TMath::Pi())  dphi -= 2.0 * TMath::Pi();
                while (dphi <= -TMath::Pi()) dphi += 2.0 * TMath::Pi();
                
                double deltaR = TMath::Sqrt(deta*deta + dphi*dphi);
                if (deltaR < minDeltaR) {
                    minDeltaR = deltaR;
                    matchedJetIdx = k;
                }
            }
            
            // Populate the sparse matrix array upon successful matching verification
            if (matchedJetIdx != -1) {
                double jetPtVal = jt.reco.pt[matchedJetIdx];
                double zVal = (jetPtVal > 0) ? (pPt / jetPtVal) : 0.0;
                if (zVal > 1.0) zVal = 1.0; 
                
                double fillArray[ndims] = { 
                    (double)pPt, 
                    (double)pEta, 
                    (double)pPhi, 
                    (double)(*pfId)[j],
                    (double)jetPtVal,
                    (double)jt.reco.eta[matchedJetIdx],
                    (double)jt.reco.phi[matchedJetIdx],
                    (double)minDeltaR,
                    (double)zVal,
                    (double)evt.hiBin
                };
                hnPF->Fill(fillArray);
            }
        }
    }

    // 9. Consolidate and deploy out files
    TFile* outFile = TFile::Open(outputFilePath.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Output path creation target failed." << std::endl;
        return;
    }
    
    outFile->cd();
    hnPF->Write();
    outFile->Close();
    
    std::cout << "Consolidated 10D Sparse written to file: " << outputFilePath << std::endl;
    
    delete chainPF;
    delete chainJet;
    delete chainEvt;
    delete chainSkim;
    delete chainHLT;
}