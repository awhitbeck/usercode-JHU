#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>

using namespace std;

// = = = = = = = = = = 
// global variables
// = = = = = = = = = = 

bool debug_=false;


// = = = = = = = = = = 
// main function 
// = = = = = = = = = = 

void addDifferences(TString inputFile="TrackSplitting_Split_Alignment.root",TString outFileName = "")
{
  TFile* fin = TFile::Open(inputFile);

  if (outFileName == "")
  {
    outFileName = inputFile;
    outFileName.ReplaceAll(".root","_wDiffs.root");
  }

  fin->cd("cosmicValidation");
  TTree* ch=(TTree*)fin->Get("splitterTree"); 

  if (ch==0x0){
    cout << "ERROR: no tree found" << endl;
    return; 
  }

  TFile *newfile = TFile::Open(outFileName,"recreate");

  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  

  double dxy1,dxy2,dz1,dz2;
  double theta1,theta2,eta1,eta2,phi1,phi2;
  double pt1,pt2,qoverpt1,qoverpt2;

  ch->SetBranchAddress("dxy1_spl",&dxy1);
  ch->SetBranchAddress("dxy2_spl",&dxy2);

  ch->SetBranchAddress("dz1_spl",&dz1);
  ch->SetBranchAddress("dz2_spl",&dz2);

  ch->SetBranchAddress("theta1_spl",&theta1);
  ch->SetBranchAddress("theta2_spl",&theta2);

  ch->SetBranchAddress("eta1_spl",&eta1);
  ch->SetBranchAddress("eta2_spl",&eta2);

  ch->SetBranchAddress("phi1_spl",&phi1);
  ch->SetBranchAddress("phi2_spl",&phi2);

  ch->SetBranchAddress("pt1_spl",&pt1);
  ch->SetBranchAddress("pt2_spl",&pt2);

  ch->SetBranchAddress("qoverpt1_spl",&qoverpt1);
  ch->SetBranchAddress("qoverpt2_spl",&qoverpt2);

  double delta_dxy, delta_dz, delta_theta, delta_eta, delta_phi, delta_pt, delta_qoverpt;

  evt_tree->Branch("Delta_dxy",&delta_dxy,"Delta_dxy/D");
  evt_tree->Branch("Delta_dz",&delta_dz,"Delta_dz/D");
  evt_tree->Branch("Delta_theta",&delta_theta,"Delta_theta/D");
  evt_tree->Branch("Delta_eta",&delta_eta,"Delta_eta/D");
  evt_tree->Branch("Delta_phi",&delta_phi,"Delta_phi/D");
  evt_tree->Branch("Delta_pt",&delta_pt,"Delta_pt/D");
  evt_tree->Branch("Delta_qoverpt",&delta_qoverpt,"Delta_qoverpt/D");

  for(int i=0; i<ch->GetEntries(); i++){

    ch->GetEntry(i);

    if(debug_){
      cout << "pt1: " << pt1 << " pt2: " << pt2 << endl;
      cout << "delta pt: " << pt1-pt2 << endl;
      cout << "dz1: " << dz1 << " dz2: " << dz2 << endl;
      cout << "delta dz: " << dz1-dz2 << endl;
      cout << "dxy1: " << dxy1 << " dxy2: " << dxy2 << endl;
      cout << "delta dxy: " << dxy1-dxy2 << endl;
    }
    delta_dxy= dxy1-dxy2;
    delta_dz = dz1-dz2;
    delta_theta = theta1-theta2;
    delta_eta= eta1-eta2;
    delta_phi= phi1-phi2;
    delta_pt = pt1-pt2;
    delta_qoverpt = qoverpt1-qoverpt2;

    evt_tree->Fill();

  }

  fin->Close();
  evt_tree->Write("splitterTree");
  newfile->Close();

}
