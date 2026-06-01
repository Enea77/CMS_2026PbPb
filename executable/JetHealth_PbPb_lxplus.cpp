#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TDirectory.h"

#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>

#include "../header/Binning.h"
#include "../header/BranchMapping.h"
#include "../header/EventStructs_PbPb.h"
#include "../header/JetSelection_PbPb.h"
#include "../header/JetStruct.h"
#include "../header/JetHealthHistograms.h"
#include "../header/JetHealthPlotting.h"
#include "../header/JetTriggers_2026PbPb_MC.h"
// #include "../header/JetTriggers_2025PbPb.h"

// minimum pt of jets to include
static constexpr Float_t ptcut = 10.0;

// minimum pt of jets to include
const TString& run_number = "";

// more than number of jets in any event being processed
static constexpr Int_t maxnref = 999;

void run(const TString& input_filelist, const TString& output, bool isMC);

int main(int argc, char* argv[]){
    if(argc < 4){
        std::cerr << "Usage: ./JetHLT <filelist.txt> <output.root> <isMC>" << std::endl;
        return 1;
    }
    std::string isMCArg = argv[3];
    if(isMCArg != "true" && isMCArg != "false" && isMCArg != "1" && isMCArg != "0"){
        std::cerr << "ERROR: isMC must be true or false" << std::endl;
        return 1;
    }
    bool isMC = (isMCArg == "true" || isMCArg == "1");
    run(argv[1], argv[2], isMC);
    return 0;
}

void JetHealth_PbPb_lxplus(const TString& input_filelist="/afs/cern.ch/user/n/nbarnett/public/HLT/5_22_2026/filelist_2025PbPb_HP2.txt", const TString& output="JetHealth_2025_Data.root", bool isMC=false){run(input_filelist, output, isMC);}
 
void run(const TString& input_filelist, const TString& output, bool isMC){

    // turning off stat box for histograms
    gStyle->SetOptStat(0);

    // ttree names
    const Int_t nTTrees = 4;
    TString sTTrees[nTTrees] = {"hiEvtAnalyzer/HiTree", "skimanalysis/HltTree", "akCs4PFJetAnalyzer/t", "hltanalysis/HltTree"}; //akPu4CaloJetAnalyzer/caloJetTree   ak4PFJetAnalyzer/t

    // for applying Jet Selections from headers
    JetSelect js;

    // event objects
    EventStruct evt;
    FiltersStruct fltr;

    // jet object
    JetStruct<maxnref> jt;

    // binning
    BinningStruct bins(ptcut);

    // initializing histograms
    JetHealthStruct<maxnref> hists(bins);

    // minbias trigger
    Int_t L1minBias;

    // getting list of root files to process
    std::ifstream myfile(input_filelist);
    std::string filename;

    // keeping track of how many files have been processed
    int filenumber = 0;

    // looping over files
    while(getline(myfile, filename)){
        filenumber+=1;

        // reading and staging input file
        TString input = filename;
        TFile *fi = TFile::Open(input,"read");
        if(!fi || fi->IsZombie()){throw std::runtime_error("ERROR: Could not open input file " + std::string(input.Data()));}
        fi->cd();

        // showing the file being processed in the terminal
        std::cout<<"processing file "<< filenumber <<": "<<input<<std::endl;

        // getting ttrees, printing iff the ttree isn't in the input file
        TTree *ttrees[nTTrees];
        for(int t=0; t<nTTrees; t++){
            ttrees[t] = (TTree*)fi->Get(sTTrees[t]);
            if(!ttrees[t]){throw std::runtime_error(Form("ERROR: Could not find TTree '%s' in file",sTTrees[t].Data()));}
        }

        // assigning variables to branches
        SetBranches(ttrees[0],evt.BranchMap(isMC));
        SetBranches(ttrees[1],fltr.BranchMap());
        SetBranches(ttrees[2],jt.BranchMap(isMC));

        // getting minbias trigger decision
        ttrees[3]->SetBranchStatus("*",0);
        ttrees[3]->SetBranchStatus("L1_MinimumBiasHF1_AND_BptxAND",1);
        ttrees[3]->SetBranchAddress("L1_MinimumBiasHF1_AND_BptxAND", &L1minBias);

        // looping over events
        Long64_t nentries = ttrees[0]->GetEntries();
        for(Long64_t i=0; i<nentries; i++){
            
            // getting event and event filter info from each event
            ttrees[0]->GetEntry(i);
            ttrees[1]->GetEntry(i);

            // // applying JSON //
            // if(!dcs.isGood(run, lumi)){continue;}

            // filling event histograms
            hists.vz_unpassed->Fill(evt.vz, evt.w);
            hists.pclustF->Fill(fltr.pclustF, evt.w);
            hists.ppvF->Fill(fltr.ppvF, evt.w);
            hists.pphfF->Fill(fltr.pphfF, evt.w);

            // checking if the event passes event filters and |vz| < 15
            if((fltr.ppvF==0)||(fltr.pclustF==0)||(fltr.pphfF==0)||(TMath::Abs(evt.vz)>15)){continue;}

            // getting jet and trigger ttree info
            ttrees[2]->GetEntry(i);
            ttrees[3]->GetEntry(i);

            /// Skipping events without jets or that don't fire teh minbias l1 trigger
            if((jt.reco.nref==0)||(L1minBias==0)){continue;}

            // filling event histograms
            hists.vz->Fill(evt.vz, evt.w);
            hists.hiBin->Fill(evt.hiBin, evt.w);
            hists.nref->Fill(jt.reco.nref, evt.w);

            // looping through all jets in each event
            for(unsigned int j=0; j<jt.reco.nref; j++){

                // pt cut on jets
                if(jt.reco.pt[j]<ptcut){continue;}
                
                // using header file for jet identification
                if(!js.JetSelection(jt.reco.eta[j], jt.reco.phi[j], jt.reco.pf.CEF[j], jt.reco.pf.NEF[j],jt.reco.pf.MUF[j])){continue;}
                
                // filling kinematic and PF fraction histograms //
                hists.FillKin(jt.reco, j, evt.hiBin, evt.w);
                hists.FillPF(jt.reco, j, evt.hiBin, evt.w);
            }
        }
        fi->Close();
    }
    
    JetHealthPlotConfig cfg;
    cfg.jetAlgo = "akCs4PF";
    cfg.runNumber = run_number;
    // making output file and storing histograms
    TFile *fo = new TFile(output,"recreate");
    //SaveJetHealthPlots(hists, bins, cfg, fo);
    hists.Write(fo);
    fo->Close();
}
