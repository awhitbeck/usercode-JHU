#include "RooRealVar.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include <sstream>
#include <vector>
#include "../src/AngularPdfFactory.cc"
#include "../PDFs/RooqqZZ_JHU.h"

/* - - - - - - - - - - - - - - - - - - - - - - - 
================================================

be sure to compile/load everything:

root -l -n 
gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L ../PDFs/RooXZsZs_5D.cxx+
.L ../src/AngularPdfFactory.cc+
.L ../PDFs/RooqqZZ_JHU.cxx+
.L MELA.C+
 - -  - - - - - - - - - - - - - - -  - - -  -  -
===============================================*/

using namespace RooFit ;

int mZZbins=350;
int lowMzz=100;
int highMzz=800;
int lowM2=12;

TFile *tempf = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");

template <typename U>
void checkZorder(U& z1mass, U& z2mass,
                 U& costhetastar, U& costheta1,
                 U& costheta2, U& phi,
                 U& phistar1){

  U tempZ1mass=z1mass;
  U tempZ2mass=z2mass;
  U tempH1=costheta1;
  U tempH2=costheta2;
  U tempHs=costhetastar;
  U tempPhi1=phistar1;
  U tempPhi=phi;

  if(z2mass>z1mass){

    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    phistar1=-tempPhi1-tempPhi;
    if(phistar1>3.1415)
      phistar1=phistar1-2*3.1415;
    if(phistar1<-3.1415)
      phistar1=phistar1+2*3.1415;

  }else
    return;

}

vector<double> my8DTemplate(bool normalized,double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1){

  //read from a file the 3D and 2D template
  TH1F *h_mzz= (TH1F*)(tempf->Get("h_mzz"));
  TH3F *h_mzzm1m2= (TH3F*)(tempf->Get("h_mzzm1m2"));
  TH2F *h_mzzcosthetastar= (TH2F*)(tempf->Get("h_mzzcosthetastar"));
  TH2F *h_mzzcostheta1= (TH2F*)(tempf->Get("h_mzzcostheta1"));
  TH2F *h_mzzcostheta2= (TH2F*)(tempf->Get("h_mzzcostheta2"));
  TH2F *h_mzzphi1= (TH2F*)(tempf->Get("h_mzzphi1"));
  TH2F *h_mzzphi= (TH2F*)(tempf->Get("h_mzzphi"));

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

//=======================================================================

pair<double,double> likelihoodDiscriminant (double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1,double scaleFactor=5.0){

  RooRealVar* z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0,180);
  RooRealVar* z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0,120); 
  RooRealVar* costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1,1);  
  RooRealVar* costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1,1);
  RooRealVar* phi_rrv= new RooRealVar("phi","#Phi",-3.1415,3.1415);
  RooRealVar* costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1,1); 
  RooRealVar* phi1_rrv= new RooRealVar("phi1","#Phi^{*}_{1}",-3.1415,3.1415);
  RooRealVar* mzz_rrv= new RooRealVar("mzz","mZZ",80,1000);

  AngularPdfFactory *SMHiggs = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  SMHiggs->makeSMHiggs();
  SMHiggs->makeParamsConst(true);
  RooqqZZ_JHU* SMZZ = new RooqqZZ_JHU("SMZZ","SMZZ",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv);

  checkZorder<double>(m1,m2,costhetastar,costheta1,costheta2,phi,phi1);

  z1mass_rrv->setVal(m1);  
  z2mass_rrv->setVal(m2);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  costhetastar_rrv->setVal(costhetastar);
  phi1_rrv->setVal(phi1);
  mzz_rrv->setVal(mZZ);

  vector <double> P=my8DTemplate(1, mZZ,  m1,  m2,  costhetastar,  costheta1,  costheta2,  phi,  phi1);

  double Pbackg=0.;
  double Psig=0.;

  if(mZZ>80 && mZZ<180){
    Pbackg = P[0]*P[1]*P[2]*P[3]*P[4]*P[5]*5.0;
    Psig = SMHiggs->getVal(mZZ);
  }if(mZZ>180&&mZZ<=2*91.188){
    z1mass_rrv->setVal(mZZ/2. - 1e-9);  // turns out signal norm will by zero if m1+m2==mzz
    z2mass_rrv->setVal(mZZ/2. - 1e-9);  // adding tiny amount to avoid this
    Pbackg = SMZZ->getVal()/(SMZZ->createIntegral(RooArgSet(*costhetastar_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*phi1_rrv))->getVal())*10.0;
    Psig = SMHiggs->PDF->getVal()/(SMHiggs->PDF->createIntegral(RooArgSet(*costheta1_rrv,*costheta2_rrv,*phi_rrv))->getVal());
  }if(mZZ>2*91.188){
    z1mass_rrv->setVal(91.188);
    z2mass_rrv->setVal(91.188);
    Pbackg = SMZZ->getVal()/(SMZZ->createIntegral(RooArgSet(*costhetastar_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*phi1_rrv))->getVal())*10.0;
    Psig = SMHiggs->PDF->getVal()/(SMHiggs->PDF->createIntegral(RooArgSet(*costheta1_rrv,*costheta2_rrv,*phi_rrv))->getVal());
  }

  // - - - - - - - - - - - - - - - - - - - - - Whitbeck 
  // check whether P[i] is zero and print warning
  // message if so

  char* varName[6]={"m1/m2","costhetastar","costheta1","coshteta2","phi","phi1"};
  for(int iVar=0; iVar<6; iVar++){

    if(P[iVar]==0 && (m1+m2)<mZZ && m2>4 && mZZ>80 && mZZ<180)
	cout << " uh oh... Probability of " << varName[iVar] << " is zero." << endl;
  }
  // - - - - - - - - - - - - - - - - - - - - - 
 
  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete costhetastar_rrv;
  delete phi1_rrv;
  delete mzz_rrv; 
  delete SMZZ;
  delete SMHiggs;

  return make_pair(Psig,Pbackg);
}

//=======================================================================

void addDtoTree(char* inputFile){

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

  char inputFileName[100];
  char outputFileName[150];
  sprintf(inputFileName,"%s.root",inputFile);
  sprintf(outputFileName,"%s_withDiscriminants.root",inputFile);

  TFile* sigFile = new TFile(inputFileName);
  TTree* sigTree=0;
    if(sigFile)
        sigTree = (TTree*) sigFile->Get("SelectedTree");
    if(!sigTree){
      cout<<"ERROR could not find the tree!"<<endl;
      return;
    }

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = new TTree("newTree","SelectedTree"); 

  float m1,m2,mzz,h1,h2,hs,phi,phi1,D;
  float ZZPt, w;
  
  sigTree->SetBranchAddress("Z1Mass",&m1);
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  sigTree->SetBranchAddress("costheta1",&h1); 
  sigTree->SetBranchAddress("costheta2",&h2);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("phi",&phi);  
  sigTree->SetBranchAddress("phistar1",&phi1);
  sigTree->SetBranchAddress("ZZPt",&ZZPt);
  sigTree->SetBranchAddress("MC_weight",&w);

  newTree->Branch("z1mass",&m1,"z1mass/F");
  newTree->Branch("z2mass",&m2,"z2mass/F");
  newTree->Branch("zzmass",&mzz,"zzmass/F");
  newTree->Branch("costheta1",&h1,"costheta1/F"); 
  newTree->Branch("costheta2",&h2,"costheta2/F");
  newTree->Branch("costhetastar",&hs,"costhetastar/F");
  newTree->Branch("phi",&phi,"phi/F");  
  newTree->Branch("phistar1",&phi1,"phistar1/F");
  newTree->Branch("ZZPt",&ZZPt,"ZZPt/F");
  newTree->Branch("melaLD",&D,"melaLD/F");
  newTree->Branch("MC_weight",&w,"MC_weight/F");

  for(int iEvt=0; iEvt<sigTree->GetEntries(); iEvt++){

    if(iEvt%5000==0) 
      cout << "event: " << iEvt << endl;

    sigTree->GetEntry(iEvt);

    checkZorder<float>(m1,m2,hs,h1,h2,phi,phi1);

    if(mzz>100. && mzz<1000. && m2>12. ) 
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

//=======================================================================

vector<TH1F*> LDDistributionBackground(char* fileName="../datafiles/7TeV/training/EWKZZ4l_Powheg_mZ4GeV_total_v25_wResolution_withDiscriminants.root"){

  TChain* chain = new TChain("SelectedTree");

  chain->Add(fileName);

  float mZZ, m2, m1, costhetastar, costheta1, costheta2, phi, phi1;
  float MC_weight=1;
  float mela=-99;
  chain->SetBranchAddress("ZZMass",&mZZ);
  chain->SetBranchAddress("Z2Mass",&m2);
  chain->SetBranchAddress("Z1Mass",&m1);
  chain->SetBranchAddress("costhetastar",&costhetastar);
  chain->SetBranchAddress("phi",&phi);
  chain->SetBranchAddress("costheta1",&costheta1);
  chain->SetBranchAddress("costheta2",&costheta2);
  chain->SetBranchAddress("phistar1",&phi1);
  chain->SetBranchAddress("MC_weight_noxsec",&MC_weight);  
  chain->SetBranchAddress("ZZLD",&mela);

  //TFile *f = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");

  TH1F *h_LDbackground= new TH1F("LD_background","LD_background",31,0,1.01);
  vector<TH1F*> vh_LDbackground;
  for (int i=1; i<mZZbins+1; i++){
    std::string s;
    std::stringstream out;
    out << i;
    s = out.str();
    TString name = "h_LDbackground_"+s;
    vh_LDbackground.push_back((new TH1F(name,name,30,0,1)));
    vh_LDbackground[i-1]->Sumw2();
  }

  for (Int_t i=0; i<chain->GetEntries();i++) {

    mela=-99;

    if(i%100000 ==0)
      cout<<"event "<<i<<endl;
      
    chain->GetEvent(i); 

    if( !(mZZ>lowMzz && mZZ<highMzz) ) // && m2>lowM2) )
      continue;
          
    if(mela<0){

      checkZorder<float>(m1,m2,costhetastar,costheta1,costheta2,phi,phi1);
      
      pair<double,double> P =  likelihoodDiscriminant(mZZ, m1, m2, costhetastar, costheta1, costheta2, phi, phi1);
      
      mela=P.first/(P.first+P.second);
    }
      
    h_LDbackground->Fill(mela,MC_weight);
    (vh_LDbackground[(int)((mZZ-100.)/2.)])->Fill(mela,MC_weight);

  }

  TCanvas *LD = new TCanvas("LD_background","LD_background",400,400);
  h_LDbackground->Draw();
  LD->Print("LD_background.eps");

  vh_LDbackground.push_back(h_LDbackground);
 
  return vh_LDbackground;
}


//=======================================================================
vector<TH1F*> LDDistributionSignal(char* fileName="../datafiles/7TeV/testBuildModel/SMHiggs_115_JHU_v0_wResolution_withDiscriminants.root"){

  TChain* chain = new TChain("SelectedTree");
  chain->Add(fileName);
 
  float mZZ, m2, m1, costhetastar, costheta1, costheta2, phi, phi1;
  float MC_weight=1;
  float mela=-99;
  chain->SetBranchAddress("ZZMass",&mZZ);
  chain->SetBranchAddress("Z2Mass",&m2);
  chain->SetBranchAddress("Z1Mass",&m1);
  chain->SetBranchAddress("costhetastar",&costhetastar);
  chain->SetBranchAddress("phi",&phi);
  chain->SetBranchAddress("costheta1",&costheta1);
  chain->SetBranchAddress("costheta2",&costheta2);
  chain->SetBranchAddress("phistar1",&phi1);
  chain->SetBranchAddress("MC_weight_noxsec",&MC_weight);
  chain->SetBranchAddress("ZZLD",&mela);

  //TFile *f = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");

  TH1F *h_LDsignal= new TH1F("LD_signal","LD_signal",31,0,1.01);
  vector<TH1F*> vh_LDsignal;
  for (int i=1; i<mZZbins+1; i++){
    std::string s;
     std::stringstream out;
     out << i;
     s = out.str();
     TString name = "  %h_LDsignal_"+s;
     vh_LDsignal.push_back((new TH1F(name,name,30,0,1)));
     vh_LDsignal[i-1]->Sumw2();
  }

  for (Int_t i=0; i<chain->GetEntries();i++) {

    mela=-99;

    if(i%100000 ==0)
      cout<<"event "<<i<<endl;
    
    chain->GetEvent(i); 

    if( !(mZZ>lowMzz && mZZ<highMzz ) ) //&& m2>lowM2) )
      continue;
      
    if(mela<0){

      checkZorder<float>(m1,m2,costhetastar,costheta1,costheta2,phi,phi1);
      
      pair<double,double> P =  likelihoodDiscriminant(mZZ, m1, m2, costhetastar, costheta1, costheta2, phi, phi1);

      mela=P.first/(P.first+P.second);

    }
    
    h_LDsignal->Fill(mela,MC_weight);
    (vh_LDsignal[(int)((mZZ-100.)/2.)])->Fill(mela,MC_weight);
     
  }
    
  cout << "Drawing" << endl;

  TCanvas *LD = new TCanvas("LD_signal","LD_signal",400,400);
  h_LDsignal->Draw();
  LD->Print("LD_signal.eps");

  vh_LDsignal.push_back(h_LDsignal);

  cout << "done" << endl;

  return vh_LDsignal;
}

//=======================================================================

TH2F* smoothTemplate(TH2F* oldTemp, TH2F* numEvents){

  int effectiveArea=3;

  TH2F* newTemp = new TH2F(*oldTemp);

  double average=0;
  int nBins=0;

  for(int i=1; i<=350; i++){
    
    //if(i<40)
      //cout << "bin: " << i ;

    for(int j=1; j<=30; j++){
      
      //if( i<40 && numEvents->GetBinContent(i,j)!=0 &&  numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)>.5) 
      //cout << " " <<  numEvents->GetBinError(i,j)/ numEvents->GetBinContent(i,j) ;

      if( numEvents->GetBinContent(i,j)!=0 && numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)<.5 )  continue;

      if( (numEvents->GetBinContent(i,j)==0 || numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)>.5)&& i>40 && i<100 && j>15){
	effectiveArea=2;
	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  for(int b=-effectiveArea; b<=effectiveArea; b++){
	    if( i+a<41 || i+a>350 || j+b<1 || j+b>30 ) continue;
	    average+=oldTemp->GetBinContent(i+a,j+b);
	    nBins++;
	  }
	}
      }else if( (numEvents->GetBinContent(i,j)==0 || numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)>.5)&& i>100 ){
	effectiveArea=3;

	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  for(int b=-effectiveArea; b<=effectiveArea; b++){
	    if( i+a<41 || i+a>350 || j+b<1 || j+b>30 ) continue;
	    average+=oldTemp->GetBinContent(i+a,j+b);
	    nBins++;
	  }
	}
      }else if( (numEvents->GetBinContent(i,j)==0 || numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)>.5)&& i<=40 ){
	effectiveArea=1;
	
	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  for(int b=-effectiveArea; b<=effectiveArea; b++){
	    if( i+a<1 || i+a>40 || j+b<1 || j+b>30 ) continue;
	    average+=oldTemp->GetBinContent(i+a,j+b);
	    nBins++;
	  }
	}
      }else continue;

      newTemp->SetBinContent(i,j,average/nBins);
      average=0;
      nBins=0;
    }// end loop over D bins
    //if(i<40) cout << endl;
  }// end loop over mZZ bins

  double norm=0;

  for(int i=1; i<=350; i++){
    for(int j=1; j<=30; j++){
      norm+=newTemp->GetBinContent(i,j);
    }

    for(int j=1; j<=30; j++){
      newTemp->SetBinContent(i,j,newTemp->GetBinContent(i,j)/norm);
    }

    norm=0;

  }

  return newTemp;

}

//=======================================================================

pair<TH2F*,TH2F*> reweightForInterference(TH2F* temp){

  TH2F* tempUp = new TH2F(*temp);
  TH2F* tempDn = new TH2F(*temp);
  
  pair<TH2F*,TH2F*> histoPair(0,0);

  // ---------------------
  // functions for scaling
  // ---------------------
  
  double oldTempValue=0;
  double newTempValue=0;
  int point=-1;

  const int numPoints=9;

  double low[numPoints]   ={100.,        122.,        128.,        135.,        145.,        155.,        165.,         175.,        185.0}; 
  double high[numPoints]  ={122.,        128.,        135.,        145.,        155.,        165.,        175.,         185.,        1000.0}; 
  double slope[numPoints] ={3.41629e-01, 3.02312e-01, 2.41973e-01, 1.16892e-01, 4.91500e-02, 3.08193e-02, -6.77412e-02, 1.13210e-02, 0.0 }; // R = yIntr+slope*D
  double yIntr[numPoints] ={8.39142e-01, 8.55536e-01, 8.87408e-01, 9.39391e-01, 9.75200e-01, 9.90251e-01,  1.03383e+00, 9.39441e-01, 1.0 }; //

  for(int i=1; i<=temp->GetNbinsX(); i++){
    point = -1;

      // choose correct scale factor
      for(int p=0; p<numPoints; p++){
	if( (i*2.+101.)>=low[p] && (i*2.+101.)<high[p] ){
	  point = p;
	}
      }
      if(point == -1){
	cout << "ERROR: could not find correct scale factor"<< endl;
	return histoPair;
      }

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      oldTempValue = temp->GetBinContent(i,j);
      newTempValue = oldTempValue*(slope[point]*(double)j/30.+yIntr[point]);
      tempUp->SetBinContent(i,j,newTempValue);
      newTempValue = oldTempValue*(-slope[point]*(double)j/30.+2.-yIntr[point]);
      tempDn->SetBinContent(i,j,newTempValue);

    }// end loop over Y bins

    // -------------- normalize mZZ slice ----------------

    double norm_up=(tempUp->ProjectionY("temp",i,i))->Integral();
    double norm_dn=(tempDn->ProjectionY("temp",i,i))->Integral();

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      tempUp->SetBinContent(i,j,tempUp->GetBinContent(i,j)/norm_up);
      tempDn->SetBinContent(i,j,tempDn->GetBinContent(i,j)/norm_dn);

    }

    // ---------------------------------------------------

  }// end loop over X bins

  histoPair.first  = tempUp;
  histoPair.second = tempDn;

  return histoPair;

}
//=======================================================================

pair<TH2F*,TH2F*> reweightForCRunc(TH2F* temp){

  cout << "reweightForCRunc" << endl;

  TH2F* tempUp = new TH2F(*temp);
  TH2F* tempDn = new TH2F(*temp);

  pair<TH2F*,TH2F*> histoPair(0,0);

  // ---------------------
  // functions for scaling
  // ---------------------
  
  double oldTempValue=0;
  double newTempValue=0;
  int point=-1;

  const int numPoints=8;

  double low[numPoints]   ={100.,        120.,        140.,         160.,     180.,     220.,     260.,     300. }; 
  double high[numPoints]  ={120.,        140.,        160.,         180.,     220.,     260.,     300.,     1000.};
  double slope[numPoints] ={4.71836e-01, 1.17671e-01, -3.81680e-01, -1.20481, -1.21944, -2.06928, -1.35337, 0.0  };
  double yIntr[numPoints] ={6.83860e-01, 9.38454e-01, 1.12690,      1.24502,  1.72764,  2.11050,  1.52771,  1.0  }; 

  for(int i=1; i<=temp->GetNbinsX(); i++){
      point = -1;

      // choose correct scale factor
      for(int p=0; p<numPoints; p++){
	if( (i*2.+101.)>=low[p] && (i*2.+101.)<high[p] ){
	  point = p;
	}
      }
      if(point == -1){
	cout << "ERROR: could not find correct scale factor"<< endl;
	return histoPair;
      }

    for(int j=1; j<=temp->GetNbinsY(); j++){

      oldTempValue = temp->GetBinContent(i,j);
      newTempValue = oldTempValue*(slope[point]*(double)j/30.+yIntr[point]);
      tempUp->SetBinContent(i,j,newTempValue);
      newTempValue = oldTempValue*(-slope[point]*(double)j/30.+2.-yIntr[point]);
      tempDn->SetBinContent(i,j,newTempValue);

    }// end loop over Y bins

    // -------------- normalize mZZ slice ----------------

    double norm_up=(tempUp->ProjectionY("temp",i,i))->Integral();
    double norm_dn=(tempDn->ProjectionY("temp",i,i))->Integral();


    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      tempUp->SetBinContent(i,j,tempUp->GetBinContent(i,j)/norm_up);
      tempDn->SetBinContent(i,j,tempDn->GetBinContent(i,j)/norm_dn);

    }

    // ---------------------------------------------------

  }// end loop over X bins

  histoPair.first  = tempUp;
  histoPair.second = tempDn;

  return histoPair;

}

//=======================================================================

void storeLDDistribution(bool signal,char* fileName, char* tag,bool smooth=true){

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    
  vector<TH1F*> vh_LD;

  TFile *file;
  char temp[50];

  if(signal){
    vh_LD = LDDistributionSignal(fileName);
    sprintf(temp,"../datafiles/Dsignal_%s.root",tag);
    file= new TFile(temp,"recreate");
  }else{
    vh_LD = LDDistributionBackground(fileName);
    sprintf(temp,"../datafiles/Dbackground_%s.root",tag);
    file= new TFile(temp,"recreate");
  }

  cout<<"LD computed. Vector size "<<vh_LD.size()<<endl;

  TH2F* h_numEvents = new TH2F("numEvents","numEvents",mZZbins,lowMzz,highMzz,vh_LD[0]->GetNbinsX(),vh_LD[0]->GetXaxis()->GetXmin(),vh_LD[0]->GetXaxis()->GetXmax());

  for (int i=1; i<mZZbins+1; i++){
    for(int j=1; j<=vh_LD[0]->GetNbinsX(); j++){
      //cout << vh_LD[i-1]->GetBinContent(j) << endl;
      h_numEvents->SetBinContent(i,j,vh_LD[i-1]->GetBinContent(j));
      h_numEvents->SetBinError(i,j,vh_LD[i-1]->GetBinError(j));
    }
  }

  for (int i=1; i<(mZZbins+2); i++){ //the last one is integrated over mzz
    if(vh_LD[i-1]->Integral()>0)
      vh_LD[i-1]->Scale(1./vh_LD[i-1]->Integral());
  }

  TH2F* h_mzzD = new TH2F("h_mzzD","h_mzzD",mZZbins,lowMzz,highMzz,vh_LD[0]->GetNbinsX(),vh_LD[0]->GetXaxis()->GetXmin(),vh_LD[0]->GetXaxis()->GetXmax());

  for (int i=1; i<mZZbins+1; i++){
    for(int j=1; j<=vh_LD[0]->GetNbinsX(); j++){
      //cout << vh_LD[i-1]->GetBinContent(j) << endl;
      h_mzzD->SetBinContent(i,j,vh_LD[i-1]->GetBinContent(j));
    }
  }

  TH2F* oldTemp = new TH2F(*h_mzzD);

  if(smooth)
    h_mzzD = smoothTemplate(h_mzzD,h_numEvents);

  pair<TH2F*,TH2F*> histoPair;
  if(!signal)
    histoPair = reweightForCRunc(h_mzzD);
  else
    histoPair = reweightForInterference(h_mzzD);

  file->cd();
  h_mzzD->Write();
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  file->Close();

}

//=======================================================================
void genMELApdf(bool isSig=true){

  TFile* inputTempFile;
  if(isSig)
    inputTempFile = new TFile("../datafiles/Dsignal.root");
  else
    inputTempFile = new TFile("../datafiles/Dbackground.root");

  TH2F* mzzDTemplate;
  if(inputTempFile)
    mzzDTemplate = (TH2F*) inputTempFile->Get("h_mzzD");

  if(!inputTempFile)
    return;

  double integral[mZZbins];
  double checkInt;
  string lines="";
  stringstream convert;

  //=====================================                                                                                                    
  // apply norm and write to file                                                                                                            
  //=====================================                                                                                                    

  for(int i=0; i<mZZbins; i++){
    checkInt=0;
    for(int j=0; j<100; j++){

      lines+=" PmzzD[";
      convert << i;
      lines+=convert.str()+"][";
      convert.str("");
      convert << j;
      lines+=convert.str()+"]=";
      convert.str("");
      if(integral[i]>0){
        convert << (double)mzzDTemplate->GetBinContent(i+1,j+1);
        checkInt+=(double)mzzDTemplate->GetBinContent(i+1,j+1)/integral[i];
      }
      else
        convert << 0.00;
      lines+=convert.str()+";  ";
      convert.str("");

      if(j%10==0)
        lines+="\n";

    }

  }
  ifstream templateFile;
  if(isSig)
    templateFile.open("../PDFs/RooMELAModelSig.tpl");
  else
    templateFile.open("../PDFs/RooMELAModelBkg.tpl");
  string templateLine="";
  ofstream newFile;
  if(isSig)
    newFile.open("../PDFs/RooMELAModelSig.cc",ios::out);
  else
    newFile.open("../PDFs/RooMELAModelBkg.cc",ios::out);

  while(templateFile.good() && !templateFile.eof()){

    getline(templateFile,templateLine);
    size_t found=templateLine.find("<INSERT ARRAY>");

    if(found!=string::npos){
      newFile << lines  << endl;
    }else{
      newFile << templateLine << endl;
    }

  }

  templateFile.close();
  newFile.close();

}

//=======================================================================

void calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4Lep11, TLorentzVector thep4Lep12, TLorentzVector thep4Z2, TLorentzVector thep4Lep21, TLorentzVector thep4Lep22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2){
	
	
	//std::cout << "In calculate angles..." << std::endl;
	
	double norm;
	
	TVector3 boostX = -(thep4H.BoostVector());
	TLorentzVector thep4Z1inXFrame( thep4Z1 );
	TLorentzVector thep4Z2inXFrame( thep4Z2 );	
	thep4Z1inXFrame.Boost( boostX );
	thep4Z2inXFrame.Boost( boostX );
	TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
	TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
	
	// calculate phi1, phi2, costhetastar
	phi1 = theZ1X_p3.Phi();
	phi2 = theZ2X_p3.Phi();
	
	///////////////////////////////////////////////
	// check for z1/z2 convention, redefine all 4 vectors with convention
	///////////////////////////////////////////////	
	TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;

	/* old convention of choosing Z1 ------------------------------
	p4H = thep4H;
	if ((phi1 < 0)&&(phi1 >= -TMath::Pi())){
		p4Z1 = thep4Z2; p4M11 = thep4Lep21; p4M12 = thep4Lep22;
		p4Z2 = thep4Z1; p4M21 = thep4Lep11; p4M22 = thep4Lep12;		
		costhetastar = theZ2X_p3.CosTheta();
	}
	else{
		p4Z1 = thep4Z1; p4M11 = thep4Lep11; p4M12 = thep4Lep12;
		p4Z2 = thep4Z2; p4M21 = thep4Lep21; p4M22 = thep4Lep22;
		costhetastar = theZ1X_p3.CosTheta();
	} ---------------------------------------------- */

	p4Z1 = thep4Z1; p4M11 = thep4Lep11; p4M12 = thep4Lep12;
	p4Z2 = thep4Z2; p4M21 = thep4Lep21; p4M22 = thep4Lep22;
	costhetastar = theZ1X_p3.CosTheta();
	
	//std::cout << "phi1: " << phi1 << ", phi2: " << phi2 << std::endl;
	
	// now helicity angles................................
	// ...................................................
	TVector3 boostZ1 = -(p4Z1.BoostVector());
	TLorentzVector p4Z2Z1(p4Z2);
	p4Z2Z1.Boost(boostZ1);
	//find the decay axis
	/////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
	TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
	norm = 1/(unitx_1.Mag());
	unitx_1*=norm;
	//boost daughters of z2
	TLorentzVector p4M21Z1(p4M21);
	TLorentzVector p4M22Z1(p4M22);
	p4M21Z1.Boost(boostZ1);
	p4M22Z1.Boost(boostZ1);
	//create z and y axes
	/////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
	TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
	TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
	TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
	norm = 1/(unitz_1.Mag());
	unitz_1 *= norm;
	TVector3 unity_1 = unitz_1.Cross(unitx_1);
	
	//caculate theta1
	TLorentzVector p4M11Z1(p4M11);
	p4M11Z1.Boost(boostZ1);
	TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
	TVector3 unitM11 = p3M11.Unit();
	double x_m11 = unitM11.Dot(unitx_1); double y_m11 = unitM11.Dot(unity_1); double z_m11 = unitM11.Dot(unitz_1);
	TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
	costheta1 = M11_Z1frame.CosTheta();
	//std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
	//////-----------------------old way of calculating phi---------------/////////
	phi = M11_Z1frame.Phi();
	
	//set axes for other system
	TVector3 boostZ2 = -(p4Z2.BoostVector());
	TLorentzVector p4Z1Z2(p4Z1);
	p4Z1Z2.Boost(boostZ2);
	TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
	norm = 1/(unitx_2.Mag());
	unitx_2*=norm;
	//boost daughters of z2
	TLorentzVector p4M11Z2(p4M11);
	TLorentzVector p4M12Z2(p4M12);
	p4M11Z2.Boost(boostZ2);
	p4M12Z2.Boost(boostZ2);
	TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
	TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
	TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
	norm = 1/(unitz_2.Mag());
	unitz_2*=norm;
	TVector3 unity_2 = unitz_2.Cross(unitx_2);
	//calcuate theta2
	TLorentzVector p4M21Z2(p4M21);
	p4M21Z2.Boost(boostZ2);
	TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
	TVector3 unitM21 = p3M21.Unit();
	double x_m21 = unitM21.Dot(unitx_2); double y_m21 = unitM21.Dot(unity_2); double z_m21 = unitM21.Dot(unitz_2);
	TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
	costheta2 = M21_Z2frame.CosTheta();
	
	// calculate phi
	//calculating phi_n
	TLorentzVector n_p4Z1inXFrame( p4Z1 );
	TLorentzVector n_p4M11inXFrame( p4M11 );
	n_p4Z1inXFrame.Boost( boostX );
	n_p4M11inXFrame.Boost( boostX );        
	TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
	TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
	TVector3 n_unitz_1( n_p4Z1inXFrame_unit );
	//// y-axis is defined by neg lepton cross z-axis
	//// the subtle part is here...
	//////////TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross( n_unitz_1 );
	TVector3 n_unity_1 = n_unitz_1.Cross( n_p4M11inXFrame_unit );
	TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );
	
	TLorentzVector n_p4M21inXFrame( p4M21 );
	n_p4M21inXFrame.Boost( boostX );
	TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
	//rotate into other plane
	TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );
	
	///////-----------------new way of calculating phi-----------------///////
	//double phi_n =  n_p4M21inXFrame_unitprime.Phi();
	/// and then calculate phistar1
	TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
	TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );
	// negative sign is for arrow convention in paper
	phistar1 = (n_p4PartoninXFrame_unitprime.Phi());
	
	// and the calculate phistar2
	TLorentzVector n_p4Z2inXFrame( p4Z2 );
	n_p4Z2inXFrame.Boost( boostX );
	TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
	///////TLorentzVector n_p4M21inXFrame( p4M21 );
	//////n_p4M21inXFrame.Boost( boostX );        
	////TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();  
	TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
	//// y-axis is defined by neg lepton cross z-axis
	//// the subtle part is here...
	//////TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross( n_unitz_2 );
	TVector3 n_unity_2 = n_unitz_2.Cross( n_p4M21inXFrame_unit );
	TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );
	TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
	phistar2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());
	
	double phistar12_0 = phistar1 + phistar2;
	if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
	else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
	else phistar12 = phistar12_0;
	
}

