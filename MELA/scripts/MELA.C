 #include "RooRealVar.h"
#include "RooWorkspace.h"
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
#include "../src/pseudoMELA.cc"
#include "../PDFs/RooqqZZ_JHU.h"
#include "../PDFs/RooTsallis.h"
#include "../PDFs/RooTsallisExp.h"
#include "../PDFs/RooRapidityBkg.h"
#include "../PDFs/RooRapiditySig.h"

/* - - - - - - - - - - - - - - - - - - - - - - - 
================================================

be sure to compile/load everything:

root -l -n 
gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L ../PDFs/RooXZsZs_5D.cxx+
.L ../src/AngularPdfFactory.cc+
.L ../PDFs/RooqqZZ_JHU.cxx+
.L ../PDFs/RooTsallis.cxx+              <-- ***RC***
.L ../PDFs/RooTsallisExp.cxx+              <-- ***RC***
.L ../PDFs/RooRapidityBkg.cxx+	       <-- ***CM&CY***
.L ../PDFs/RooRapiditySig.cxx+	       <-- ***CM&CY***
.L MELA.C+
 - -  - - - - - - - - - - - - - - -  - - -  -  -
===============================================*/

using namespace RooFit ;

int lowMzz=100.;
int highMzz=800.;
int lowM2=12.;

TFile *tempf = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");

void setTemplate(char* file){
  if (tempf !=0)
    delete tempf;

  tempf= new TFile(file,"READ");
  tempf->Print("v");
}

void initAllFunctions(RooWorkspace *ws, int LHCsqrts) {

  RooRealVar* z1mass_rrv = new RooRealVar("z1mass_rrv","m_{Z1}",0,180);
  RooRealVar* z2mass_rrv = new RooRealVar("z2mass_rrv","m_{Z2}",0,120); 
  RooRealVar* costheta1_rrv = new RooRealVar("costheta1_rrv","cos#theta_{1}",-1,1);  
  RooRealVar* costheta2_rrv = new RooRealVar("costheta2_rrv","cos#theta_{2}",-1,1);
  RooRealVar* phi_rrv= new RooRealVar("phi_rrv","#Phi",-3.1415,3.1415);
  RooRealVar* costhetastar_rrv = new RooRealVar("costhetastar_rrv","cos#theta^{*}",-1,1); 
  RooRealVar* phi1_rrv= new RooRealVar("phi1_rrv","#Phi^{*}_{1}",-3.1415,3.1415);
  RooRealVar* mzz_rrv= new RooRealVar("mzz_rrv","mZZ",80,1000);
  RooRealVar* y_rrv= new RooRealVar("y_rrv","y_{ZZ}",-4.0,4.0);
  RooRealVar* pt_rrv= new RooRealVar("pt_rrv","p_{T,ZZ}",0.0,1000.);
  RooRealVar* sqrtS_rrv = new RooRealVar("sqrtS_rrv","#sqrt{s}",1000.,14000.);
  ws->import(*z1mass_rrv);
  ws->import(*z2mass_rrv);
  ws->import(*costheta1_rrv);
  ws->import(*costheta2_rrv);
  ws->import(*phi_rrv);
  ws->import(*costhetastar_rrv);
  ws->import(*phi1_rrv);
  ws->import(*mzz_rrv);
  ws->import(*y_rrv);
  ws->import(*pt_rrv);
  ws->import(*sqrtS_rrv);
 
  static const int NptparamsS = 17;
  static const int NptparamsB = 11;

  char* rrvnamesB[NptparamsB] = {"m","n0","n1","n2","ndue","bb0","bb1","bb2","T0","T1","T2"};
  RooRealVar *ptparamsB[NptparamsB];
  RooArgSet* allparamsB = new RooArgSet();
  for (int i = 0; i < NptparamsB; i++) {
    ptparamsB[i] = new RooRealVar(rrvnamesB[i],rrvnamesB[i],-10000.,10000.);
    allparamsB->add(*ptparamsB[i]);
  }

  char* rrvnamesS[NptparamsS] = {"ms","ns0","ns1","ns2","ndues","bbs0","bbs1","bbs2","Ts0","Ts1","Ts2","bbdues0","bbdues1","bbdues2","fexps0","fexps1","fexps2"};
  RooRealVar *ptparamsS[NptparamsS];
  RooArgSet* allparamsS = new RooArgSet();
  for (int i = 0; i < NptparamsS; i++) {
    ptparamsS[i] = new RooRealVar(rrvnamesS[i],rrvnamesS[i],-10000.,10000.);
    allparamsS->add(*ptparamsS[i]);
  }
 
  RooqqZZ_JHU* SMZZ = new RooqqZZ_JHU("SMZZ","SMZZ",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv);
  ws->import(*SMZZ);
 
  char fileName[200];
  sprintf(fileName,"../datafiles/allParamsSig_%dTeV.txt",LHCsqrts);

  RooTsallisExp* sigPt = new RooTsallisExp("sigPt","sigPt",*pt_rrv,*mzz_rrv,
					   *ptparamsS[0],*ptparamsS[1],*ptparamsS[2],
					   *ptparamsS[3],*ptparamsS[4],*ptparamsS[5],
					   *ptparamsS[6],*ptparamsS[7],*ptparamsS[8],
					   *ptparamsS[9],*ptparamsS[10],*ptparamsS[11],
					   *ptparamsS[12],*ptparamsS[13],*ptparamsS[14],
					   *ptparamsS[15],*ptparamsS[16]);

  allparamsS->readFromFile(fileName,0);
  ws->import(*sigPt);

  sprintf(fileName,"../datafiles/allParamsBkg_%dTeV.txt",LHCsqrts);
  RooTsallis* bkgPt = new RooTsallis("bkgPt","bkgPt",*pt_rrv,*mzz_rrv,
				     *ptparamsB[0],*ptparamsB[1],*ptparamsB[2],
				     *ptparamsB[3],*ptparamsB[4],*ptparamsB[5],
				     *ptparamsB[6],*ptparamsB[7],*ptparamsB[8],
				     *ptparamsB[9],*ptparamsB[10]);
  allparamsB->readFromFile(fileName,0);
  ws->import(*bkgPt);

  RooRapiditySig* sigY = new RooRapiditySig("sigY", "sigY", *y_rrv, *mzz_rrv, *sqrtS_rrv);
  ws->import(*sigY);

  RooRapidityBkg* bkgY = new RooRapidityBkg("bkgY", "bkgY", *y_rrv, *mzz_rrv, *sqrtS_rrv);
  ws->import(*bkgY);

}


template <typename U>
void checkZorder(U& z1mass, U& z2mass,
                 U& costhetastar, U& costheta1,
                 U& costheta2, U& phi,
                 U& phistar1, U& pt, U& ipsilon){

  U tempZ1mass=z1mass;
  U tempZ2mass=z2mass;
  U tempH1=costheta1;
  U tempH2=costheta2;
  U tempHs=costhetastar;
  U tempPhi1=phistar1;
  U tempPhi=phi;
  U tempPt=pt;
  U tempY=ipsilon;

  if(z2mass>z1mass){

    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    pt=tempPt;
    ipsilon=tempY;
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

pair<double,double> likelihoodDiscriminant (RooWorkspace *ws, 
					    double mZZ, double m1, double m2, 
					    double costhetastar, double costheta1, double costheta2, double phi, double phi1, 
					    int LHCsqrts, 
					    bool withPt = false, double pt = 0.0, 
					    bool withY = false, double y = 0.0,
					    double scaleFactor = 5.0){
  
  ws->var("sqrtS_rrv")->setVal(LHCsqrts*1000.);
  ws->var("sqrtS_rrv")->setConstant(kTRUE);

  checkZorder<double>(m1,m2,costhetastar,costheta1,costheta2,phi,phi1,pt,y);

  ws->var("z1mass_rrv")->setVal(m1);  
  ws->var("z2mass_rrv")->setVal(m2);
  ws->var("costheta1_rrv")->setVal(costheta1);
  ws->var("costheta2_rrv")->setVal(costheta2);
  ws->var("phi_rrv")->setVal(phi);
  ws->var("costhetastar_rrv")->setVal(costhetastar);
  ws->var("phi1_rrv")->setVal(phi1);
  ws->var("mzz_rrv")->setVal(mZZ);
  ws->var("pt_rrv")->setVal(pt);
  ws->var("y_rrv")->setVal(y);

  vector <double> P=my8DTemplate(1, mZZ,  m1,  m2,  costhetastar,  costheta1,  costheta2,  phi,  phi1);

  AngularPdfFactory *SMHiggs = new AngularPdfFactory(ws->var("z1mass_rrv"),ws->var("z2mass_rrv"),ws->var("costheta1_rrv"),ws->var("costheta2_rrv"),ws->var("phi_rrv"),ws->var("mzz_rrv"));
  SMHiggs->makeSMHiggs();
  SMHiggs->makeParamsConst(true);

  double Pbackg=0.;
  double Psig=0.;

  if(mZZ>80 && mZZ<180){
    Pbackg = P[0]*P[1]*P[2]*P[3]*P[4]*P[5]*5.0;
    Psig = SMHiggs->getVal(mZZ);
  }if(mZZ>180&&mZZ<=2*91.188){
    ws->var("z1mass_rrv")->setVal(mZZ/2. - 1e-9);  // turns out signal norm will by zero if m1+m2==mzz
    ws->var("z2mass_rrv")->setVal(mZZ/2. - 1e-9);  // adding tiny amount to avoid this
    Pbackg = ws->pdf("SMZZ")->getVal()/( ws->pdf("SMZZ")->createIntegral(RooArgSet(*ws->var("costhetastar_rrv"),*ws->var("costheta1_rrv"),*ws->var("costheta2_rrv"),*ws->var("phi_rrv"),*ws->var("phi1_rrv")))->getVal())*10.0;
    Psig = SMHiggs->PDF->getVal()/(SMHiggs->PDF->createIntegral(RooArgSet(*ws->var("costheta1_rrv"),*ws->var("costheta2_rrv"),*ws->var("phi_rrv")))->getVal());
  }if(mZZ>2*91.188){
    ws->var("z1mass_rrv")->setVal(91.188);
    ws->var("z2mass_rrv")->setVal(91.188);
    Pbackg = ws->pdf("SMZZ")->getVal()/( ws->pdf("SMZZ")->createIntegral(RooArgSet(*ws->var("costhetastar_rrv"),*ws->var("costheta1_rrv"),*ws->var("costheta2_rrv"),*ws->var("phi_rrv"),*ws->var("phi1_rrv")))->getVal())*10.0;
    Psig = SMHiggs->PDF->getVal()/(SMHiggs->PDF->createIntegral(RooArgSet(*ws->var("costheta1_rrv"),*ws->var("costheta2_rrv"),*ws->var("phi_rrv")))->getVal());
  }

  if (withPt) {
    Pbackg *= ws->pdf("bkgPt")->getVal()/(ws->pdf("bkgPt")->createIntegral(RooArgSet(*ws->var("pt_rrv")))->getVal());
    Psig *= ws->pdf("sigPt")->getVal()/(ws->pdf("sigPt")->createIntegral(RooArgSet(*ws->var("pt_rrv")))->getVal());
  }
  if(withY) {
    Pbackg *= ws->pdf("bkgY")->getVal()/(ws->pdf("bkgY")->createIntegral(RooArgSet(*ws->var("y_rrv")))->getVal());
    Psig *= ws->pdf("sigY")->getVal()/(ws->pdf("sigY")->createIntegral(RooArgSet(*ws->var("y_rrv")))->getVal());
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
 
  delete SMHiggs;

  return make_pair(Psig,Pbackg);
}

//=======================================================================

void addDtoTree(char* inputFile, int LHCsqrts = 7, 
		float minMzz = 100., float maxMzz = 800., 
		bool containsPt = false, bool containsY = false){

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

  if (LHCsqrts < 7 || LHCsqrts > 8) {
    cout << "LHC energy must be 7 or 8 TeV!" << endl;
    return;
  }

  pseudoMELA pseudoMELAcalculator;

  char inputFileName[100];
  char outputFileName[150];
  sprintf(inputFileName,"%s.root",inputFile);

  if(!containsPt && !containsY)
    {
      sprintf(outputFileName,"%s_%d-%d_withDiscriminants.root",inputFile,int(minMzz),int(maxMzz));
    }
  else if(containsPt && !containsY)
    {
      sprintf(outputFileName,"%s_%d-%d_withDiscriminantsPt.root",inputFile,int(minMzz),int(maxMzz));
    }
  else if(!containsPt && containsY)
    {
      sprintf(outputFileName,"%s_%d-%d_withDiscriminantsY.root",inputFile,int(minMzz),int(maxMzz));
    }
  else if(containsPt && containsY)
    {
      sprintf(outputFileName,"%s_%d-%d_withDiscriminantsPtY.root",inputFile,int(minMzz),int(maxMzz));
    }

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

  float m1,m2,mzz,h1,h2,hs,phi,phi1,D,psD;
  float ZZPt, ZZY, w;
  
  sigTree->SetBranchAddress("Z1Mass",&m1);
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  sigTree->SetBranchAddress("helcosthetaZ1",&h1); 
  sigTree->SetBranchAddress("helcosthetaZ2",&h2);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("helphi",&phi);  
  sigTree->SetBranchAddress("phistarZ1",&phi1);

  if (containsPt) sigTree->SetBranchAddress("ZZPt",&ZZPt);
  if (containsY) sigTree->SetBranchAddress("ZZRapidity", &ZZY);
  sigTree->SetBranchAddress("MC_weight",&w);

  newTree->Branch("z1mass",&m1,"z1mass/F");
  newTree->Branch("z2mass",&m2,"z2mass/F");
  newTree->Branch("zzmass",&mzz,"zzmass/F");
  newTree->Branch("costheta1",&h1,"costheta1/F"); 
  newTree->Branch("costheta2",&h2,"costheta2/F");
  newTree->Branch("costhetastar",&hs,"costhetastar/F");
  newTree->Branch("phi",&phi,"phi/F");  
  newTree->Branch("phistar1",&phi1,"phistar1/F");
  if (containsPt) {
    newTree->Branch("ZZPt",&ZZPt,"ZZPt/F");
    if(!containsY)
      newTree->Branch("melaLDWithPt",&D,"melaLDWithPt/F");
    else newTree->Branch("melaLDWithPtY", &D,"melaLDWithPtY/F");
  }
  if (containsY) {
    newTree->Branch("ZZY", &ZZY, "ZZY/F");
    if(!containsPt) newTree->Branch("melaLDWithY", &D, "melaLDWithY/F");
  }
  if (!containsPt && !containsY) newTree->Branch("melaLD",&D,"melaLD/F");
  newTree->Branch("pseudomelaLD",&psD,"pseudomelaLD/F");  
  newTree->Branch("MC_weight",&w,"MC_weight/F");

  RooWorkspace *ws = new RooWorkspace("ws");
  initAllFunctions(ws,LHCsqrts);

  for(int iEvt=0; iEvt<sigTree->GetEntries(); iEvt++){

    if(iEvt%5000==0) 
      cout << "event: " << iEvt << endl;

    sigTree->GetEntry(iEvt);
    
    checkZorder<float>(m1,m2,hs,h1,h2,phi,phi1,ZZPt,ZZY);
    
    if(mzz>minMzz && mzz<maxMzz && m2>lowM2 ) 
      {
	
	//MELA LD
	pair<double,double> P;
	if(!containsPt && !containsY)
	  P = likelihoodDiscriminant(ws, mzz, m1, m2, hs, h1, h2, phi, phi1, LHCsqrts);
	else if (containsPt && !containsY) 
	  P = likelihoodDiscriminant(ws, mzz, m1, m2, hs, h1, h2, phi, phi1, LHCsqrts, true, ZZPt, false);
	else if (!containsPt && containsY)
	  P = likelihoodDiscriminant(ws, mzz, m1, m2, hs, h1, h2, phi, phi1, LHCsqrts, false, 0., true, ZZY);
	else if(containsPt && containsY)
	  P = likelihoodDiscriminant(ws, mzz, m1, m2, hs, h1, h2, phi, phi1, LHCsqrts, true, ZZPt, true, ZZY);
	
	D=P.first/(P.first+P.second);

        psD = pseudoMELAcalculator.eval(mzz, m1, m2, hs, h1, h2, phi, phi1);

	newTree->Fill();
	
      }
  }

  newFile->cd();
  newTree->Write("angles"); 
  newFile->Close();
  
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
