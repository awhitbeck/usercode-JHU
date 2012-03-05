#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TString.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TString.h"
#include <sstream>
#include <string>
#include <vector>
#include "../src/AngularPdfFactory.cc"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

/*
//be sure to compile/load everything:
 gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L ../PDFs/RooXZsZs_5D.cxx+
.L ../src/AngularPdfFactory.cc+
.L MELA.C+
*/

using namespace RooFit ;

 TFile *f = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");

vector<double> my8DTemplate(bool normalized,double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1){
  //read from a file the 3D and 2D template
  TH1F *h_mzz= (TH1F*)(f->Get("h_mzz"));
  TH3F *h_mzzm1m2= (TH3F*)(f->Get("h_mzzm1m2"));
  TH2F *h_mzzcosthetastar= (TH2F*)(f->Get("h_mzzcosthetastar"));
  TH2F *h_mzzcostheta1= (TH2F*)(f->Get("h_mzzcostheta1"));
  TH2F *h_mzzcostheta2= (TH2F*)(f->Get("h_mzzcostheta2"));
  TH2F *h_mzzphi1= (TH2F*)(f->Get("h_mzzphi1"));
  TH2F *h_mzzphi= (TH2F*)(f->Get("h_mzzphi"));

  //multiply the P values
  double n = h_mzz->GetBinContent(h_mzz->FindBin(mZZ));
  double Pmzzm1m2 = h_mzzm1m2->GetBinContent(h_mzzm1m2->FindBin(mZZ,m1,m2));

  // - - - - - - - - - - - - - - - whitbeck
  // if bin has no events: add 1
  // safety feature to prevent LD = 1 as a
  // result of low statistics

  if(Pmzzm1m2==0){
    Pmzzm1m2++;
    }
  // - - - - - - - - - - - - - - - 

  double Pmzzcosthetastar = h_mzzcosthetastar->GetBinContent(h_mzzcosthetastar->FindBin(mZZ,costhetastar));
  double Pmzzcostheta2 = h_mzzcostheta2->GetBinContent(h_mzzcostheta2->FindBin(mZZ,costheta2));
  double Pmzzcostheta1 = h_mzzcostheta1->GetBinContent(h_mzzcostheta1->FindBin(mZZ,costheta1));
  double Pmzzphi1 = h_mzzphi1->GetBinContent(h_mzzphi1->FindBin(mZZ,phi1));
  double Pmzzphi = h_mzzphi->GetBinContent(h_mzzphi->FindBin(mZZ,phi));

  //normalization
  double binwidth_mzzm1m2 = h_mzzm1m2->GetYaxis()->GetBinWidth(1) * h_mzzm1m2->GetZaxis()->GetBinWidth(1);
  double binwidth_mzzcosthetastar = h_mzzcosthetastar->GetYaxis()->GetBinWidth(1);
  double binwidth_mzzcostheta1 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
  double binwidth_mzzcostheta2 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
  double binwidth_mzzphi1 = h_mzzphi1->GetYaxis()->GetBinWidth(1);
  double binwidth_mzzphi = h_mzzphi->GetYaxis()->GetBinWidth(1);

  double Pmzzm1m2_norm = Pmzzm1m2/(n*binwidth_mzzm1m2); 
  double Pmzzcosthetastar_norm = Pmzzcosthetastar/(n*binwidth_mzzcosthetastar);
  double Pmzzcostheta1_norm = Pmzzcostheta1/(n*binwidth_mzzcostheta1);
  double Pmzzcostheta2_norm = Pmzzcostheta2/(n*binwidth_mzzcostheta2);
  double Pmzzphi1_norm = Pmzzphi1/(n*binwidth_mzzphi1);
  double Pmzzphi_norm = Pmzzphi/(n*binwidth_mzzphi);

  vector <double> P;
  P.push_back(Pmzzm1m2);
  P.push_back(Pmzzcosthetastar);
  P.push_back(Pmzzcostheta1);
  P.push_back(Pmzzcostheta2);
  P.push_back(Pmzzphi);
  P.push_back(Pmzzphi1);

  vector <double> P_norm;
  P_norm.push_back(Pmzzm1m2_norm);
  P_norm.push_back(Pmzzcosthetastar_norm);
  P_norm.push_back(Pmzzcostheta1_norm);
  P_norm.push_back(Pmzzcostheta2_norm);
  P_norm.push_back(Pmzzphi_norm);
  P_norm.push_back(Pmzzphi1_norm);

  
  if(normalized)
    return P_norm;
  else
    return P;
}


pair<double,double> likelihoodDiscriminant (double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1,double scaleFactor=5.0){

  RooRealVar* z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0,180);
  RooRealVar* z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0,120); 
  RooRealVar* costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1,1);  
  RooRealVar* costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1,1);
  RooRealVar* phi_rrv= new RooRealVar("phi","#Phi",-3.1415,3.1415);
  RooRealVar* mzz_rrv= new RooRealVar("mzz","mZZ",110,180);
  AngularPdfFactory *SMHiggs = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  SMHiggs->makeSMHiggs();
  SMHiggs->makeParamsConst(true);
  
  z1mass_rrv->setVal(m1);  
  z2mass_rrv->setVal(m2);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  mzz_rrv->setVal(mZZ);
  double Psig = SMHiggs->getVal(mZZ);

  z1mass_rrv->setVal(m1);  
 
  vector <double> P=my8DTemplate(1, mZZ,  m1,  m2,  costhetastar,  costheta1,  costheta2,  phi,  phi1);

  double Pbackg = P[0]*P[1]*P[2]*P[3]*P[4]*P[5]*scaleFactor; 

  // - - - - - - - - - - - - - - - - - - - - - Whitbeck 
  // check whether P[i] is zero and print warning
  // message if so

  char* varName[6]={"m1/m2","costhetastar","costheta1","coshteta2","phi","phi1"};
  for(int iVar=0; iVar<6; iVar++){

    if(P[iVar]==0 && (m1+m2)<mZZ && m2>20 && mZZ>110 && mZZ<180)
	cout << " uh oh... Probability of " << varName[iVar] << " is zero." << endl;
  }
  // - - - - - - - - - - - - - - - - - - - - - 
 
  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete mzz_rrv; 

  delete SMHiggs;

  return make_pair(Psig,Pbackg);
}

//=======================================================================

void addDtoTree(char* inputFile){

  char inputFileName[100];
  char outputFileName[150];
  sprintf(inputFileName,"%s.root",inputFile);
  sprintf(outputFileName,"%s_withDiscriminants.root",inputFile);

  TFile* sigFile = new TFile(inputFileName);
  TTree* sigTree;
    if(sigFile)
        sigTree = (TTree*) sigFile->Get("angles");
    if(!sigTree)
        return;

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = new TTree("newTree","angles"); 

  double m1,m2,mzz,h1,h2,hs,phi,phi1,D;
  sigTree->SetBranchAddress("z1mass",&m1);
  sigTree->SetBranchAddress("z2mass",&m2);
  sigTree->SetBranchAddress("zzmass",&mzz);
  sigTree->SetBranchAddress("costheta1",&h1); 
  sigTree->SetBranchAddress("costheta2",&h2);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("phi",&phi);  
  sigTree->SetBranchAddress("phistar1",&phi1);

  newTree->Branch("z1mass",&m1,"z1mass/D");
  newTree->Branch("z2mass",&m2,"z2mass/D");
  newTree->Branch("zzmass",&mzz,"zzmass/D");
  newTree->Branch("costheta1",&h1,"costheta1/D"); 
  newTree->Branch("costheta2",&h2,"costheta2/D");
  newTree->Branch("costhetastar",&hs,"costhetastar/D");
  newTree->Branch("phi",&phi,"phi/D");  
  newTree->Branch("phistar1",&phi1,"phistar1/D");
  newTree->Branch("melaLD",&D,"melaLD/D");  
  
  for(int iEvt=0; iEvt<sigTree->GetEntries(); iEvt++){

    if(iEvt%5000==0) cout << "event: " << iEvt << endl;

    sigTree->GetEntry(iEvt);

    if(mzz>110. && mzz<180. && m2>20) 
      {

      //MELA LD
      pair<double,double> P = likelihoodDiscriminant(mzz, m1, m2, hs, h1, h2, phi, phi1);
      D=P.first/(P.first+P.second);
      newTree->Fill();

    }
   }

  newFile->cd();
  newTree->Write("angles"); 
  newFile->Close();

}

