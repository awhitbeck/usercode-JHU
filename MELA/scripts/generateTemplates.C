#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include <sstream>
#include <vector>

/* - - - - - - - - - - - - - - - - - - - - - - - 
================================================

be sure to compile/load everything:

root -l -n 
gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L generateTEmplates.C+
- -  - - - - - - - - - - - - - - -  - - -  -  -
===============================================*/

const int mZZbins=350;
int lowMzz=100;
int highMzz=800;
int lowM2=12;

char* discriminantName="ZZLD";

TH1F* mzzBinning = new TH1F("mzzBinning","mzzBinning",350,100,800);

TFile *tempf = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");

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

  double low[numPoints]   ={100.,        120.,        140.,         160.,     180.  };
  double high[numPoints]  ={120.,        140.,        160.,         180.,     1000. };

  /* ================ systematics for pseudoMELA ==========================
  double slope[numPoints] ={-3.32705e-01, -1.90814e-01, -9.77189e-01, -3.81680e-01, 0.0 };
  double yIntr[numPoints] ={ 9.05727e-01, 9.95995e-01,  1.40367e+00,  1.12690,      1.0 }; 
  ==================================================================*/

  // ================ systematics for MELA ==========================
  double slope[numPoints] ={4.71836e-01, 1.17671e-01, -3.81680e-01, -1.20481, -1.21944, -2.06928, -1.35337, 0.0 };
  double yIntr[numPoints] ={6.83860e-01, 9.38454e-01, 1.12690,      1.24502,  1.72764,  2.11050,  1.52771,  1.0 }; 
  //==================================================================


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
TH2F* reweightForInterference(TH2F* temp){

  // for interference reweighting of MELA
  TF1* reweightFunc = new TF1("reweightFunc","gaus",100,1000);

  reweightFunc->SetParameter(0,0.354258);
  reweightFunc->SetParameter(1,114.909);
  reweightFunc->SetParameter(2,17.1512);

  /* ===================================================
  // for interference reweighting of pseudo-MELA
  TF1* reweightFunc = new TF1("reweightFunc","([0]+[1]*(x-110) )*0.5*(1 + TMath::Erf([2]*(x -[3]) ))*0.5*(1 + TMath::Erf([4]*([5]-x) ))  ",100,200);

  reweightFunc->SetParameter(0,-5.66409e-01);
  reweightFunc->SetParameter(1, 1.22591e-02);
  reweightFunc->SetParameter(2, 1.64942e+00);
  reweightFunc->SetParameter(3, 1.10080e+02);
  reweightFunc->SetParameter(4, 2.10905e+00);
  reweightFunc->SetParameter(5, 1.78529e+02);
  ==================================================== */

  TH2F* newTemp = new TH2F(*temp);
  
  // ---------------------
  // functions for scaling
  // ---------------------
  
  double oldTempValue=0;
  double newTempValue=0;

  double slope;

  for(int i=1; i<=temp->GetNbinsX(); i++){

    // choose correct scale factor

    // for reweighting MELA
    if(i<8){
      slope=.354;
    }else{
      slope=reweightFunc->Eval((double)((i-1)*2+101));
    }

    /* ==============================================
    // for reweighting pseudo-MELA
    slope = reweightFunc->Eval((double)((i-1)*2+101));
    ============================================== */

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      oldTempValue = temp->GetBinContent(i,j);
      newTempValue = oldTempValue*(1+slope*((double)j/30.-.5));
      newTemp->SetBinContent(i,j,newTempValue);

    }// end loop over Y bins

    // -------------- normalize mZZ slice ----------------

    double norm=(newTemp->ProjectionY("temp",i,i))->Integral();

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      newTemp->SetBinContent(i,j,newTemp->GetBinContent(i,j)/norm);

    }

    // ---------------------------------------------------

  }// end loop over X bins

  return newTemp;

}

//=======================================================================

TH2F* fillTemplate(char* channel="4mu", int sampleIndex=0,bool isLowMass=true){

  string sample[4]={"H*","ZZTo*","ggZZ*","H*Pse"};
  string sampleName[4]={"signal","qqZZ","ggZZ","signal_PS"};


  TChain* bkgMC = new TChain("SelectedTree");
  char temp[100];
  if(isLowMass){

    if(sampleIndex==3){
      sprintf(temp,"CJLSTtree_Jun25_2012/JHUsignal/HZZ%sTree_%s.root",channel,sample[sampleIndex].c_str());
      bkgMC->Add(temp);
    }else{
      sprintf(temp,"CJLSTtrees_June27_2012/7plus8TeV_FSR/HZZ%sTree_%s.root",channel,sample[sampleIndex].c_str());
      bkgMC->Add(temp);
    }
  }else{
    if(sampleIndex==3){
      sprintf(temp,"CJLSTtree_Jun25_2012/JHUsignal/HZZ*Tree_%s.root",sample[sampleIndex].c_str());
      bkgMC->Add(temp);
    }else{
      sprintf(temp,"CJLSTtrees_June27_2012/7plus8TeV_FSR/HZZ*Tree_%s.root",sample[sampleIndex].c_str());
      bkgMC->Add(temp);
    }
  }
  
  
  float mzz,mela,w;
  
  bkgMC->SetBranchAddress("ZZMass",&mzz);
  bkgMC->SetBranchAddress("ZZLD",&mela);
  bkgMC->SetBranchAddress("MC_weight_noxsec",&w);
  
  TH2F* bkgHist;
  if(!isLowMass)
    bkgHist = new TH2F("bkgHisto","bkgHisto",310,180,800,30,0,1);
  else
    bkgHist = new TH2F("bkgHisto","bkgHisto",40,100,180,30,0,1);

  bkgHist->Sumw2();

  // fill histogram

  for(int i=0; i<bkgMC->GetEntries(); i++){

    bkgMC->GetEntry(i);

    if(w<.0015 /*&& mela>.5 */ ){
      
      bkgHist->Fill(mzz,mela,w);

    }

  }

  int nXbins=(isLowMass)?40:310;
  int nYbins=30;
    
  // normalize slices

  double norm;
  TH1F* tempProj;
  
  for(int i=1; i<=nXbins; i++){
    
    tempProj = (TH1F*) bkgHist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();

    for(int j=1; j<=nYbins; j++){
      bkgHist->SetBinContent(i,j, bkgHist->GetBinContent(i,j)/norm   );
    }

  }
  
  // average 

  TH2F* notSmooth = new TH2F(*bkgHist);

  if(!isLowMass){
    
    int binMzz;
    int effectiveArea=1;
    double average=0,binsUsed=0;

    for(int i=1; i<=nXbins; i++){
      for(int j=1; j<=nYbins; j++){
	
	binMzz=(i-1)*2+181;

	if( binMzz<300 ) continue;
	if( binMzz>=300 && binMzz<350 ) effectiveArea=1;
	if( binMzz>=350 && binMzz<500 ) effectiveArea=3;
	if( binMzz>=500 && binMzz<600 ) effectiveArea=5;
	if( binMzz>=600 ) effectiveArea=7;
	
	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  if(a+i<1 || a+i>nXbins || j>nYbins || j<1) continue;
	  average+= notSmooth->GetBinContent(a+i,j);
	  binsUsed++;
	}
	
	bkgHist->SetBinContent(i,j,average/binsUsed);

	average=0;
	binsUsed=0;
	
      } // end loop over D
    } // end loop over mZZ
  } // end of horizontal averaging
  
  // smooth
  
  bkgHist->Smooth();
  if(!isLowMass)
    bkgHist->Smooth();
  
  for(int i=1; i<=nXbins; i++){
    for(int j=1; j<=nYbins; j++){
      if(bkgHist->GetBinContent(i,j)==0)
	bkgHist->SetBinContent(i,j,.00001);
    }// for(int j=1; j<=nYbins; j++){
  }// for(int i=1; i<=nXbins; i++){

  return bkgHist;
  
}

//=======================================================================

TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp){
  
  TH2F* h_mzzD = new TH2F("h_mzzD","h_mzzD",350,100,800,30,0,1);
  
  // copy lowmass into h_mzzD
  for(int i=1; i<=40; i++){
    for(int j=1; j<=30; j++){
      h_mzzD->SetBinContent(i,j, lowTemp->GetBinContent(i,j)  );
    }// end loop over D
  }// end loop over mZZ

  // copy high mass into h_mzzD
  for(int i=1; i<=310; i++){
    for(int j=1; j<=30; j++){
      h_mzzD->SetBinContent(i+40,j, highTemp->GetBinContent(i,j)  );
    }// end loop over D
  }// end loop over mZZ

  return h_mzzD;

}

//=======================================================================

void makeTemplate(char* channel="4mu"){

  char temp[150];

  sprintf(temp,"../datafiles/Dsignal_%s.root",channel);
  TFile* fsig = new TFile(temp,"RECREATE");
  sprintf(temp,"../datafiles/Dsignal_ALT_%s.root",channel);
  TFile* fpssig = new TFile(temp,"RECREATE");
  sprintf(temp,"../datafiles/Dbackground_qqZZ_%s.root",channel);
  TFile* fqqZZ = new TFile(temp,"RECREATE");
  sprintf(temp,"../datafiles/Dbackground_ggZZ_%s.root",channel);
  TFile* fggZZ = new TFile(temp,"RECREATE");
  TH2F* oldTemp;

  pair<TH2F*,TH2F*> histoPair;

  TH2F* low,*high,*h_mzzD;
  
  // ========================================
  // SM Higgs template

  low = fillTemplate(channel,0,true);
  high = fillTemplate(channel,0,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "correcting for interference and adding syst" << endl;
  if(strcmp(channel,"2e2mu"))
    h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference

  cout << "h_mzzD: " << h_mzzD << endl;

  // --------------------------------------------------

  fsig->cd();

  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  h_mzzD->Write("h_mzzD_up");
  h_mzzD->Write("h_mzzD_dn");
  fsig->Close();

  /* ========================================
  // pseudo scalar template

  low = fillTemplate(channel,3,true);
  high = fillTemplate(channel,3,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "correcting for interference and adding syst" << endl;
  if(strcmp(channel,"2e2mu"))
    h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference

  cout << "h_mzzD: " << h_mzzD << endl;

  // --------------------------------------------------

  fpssig->cd();
  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  fpssig->Close();

  ======================================= */
  // qqZZ template

  low = fillTemplate(channel,1,true);
  high = fillTemplate(channel,1,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "apply systematics for zjets control region" << endl;
  
  histoPair = reweightForCRunc(h_mzzD);

  // --------------------------------------------------

  fqqZZ->cd();
  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  fqqZZ->Close();

  // ==========================
  // ggZZ templates
  
  low = fillTemplate(channel,2,true);
  high = fillTemplate(channel,2,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "apply systematics for zjets control region" << endl;
  
  histoPair = reweightForCRunc(h_mzzD);

  // --------------------------------------------------

  fggZZ->cd();
  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  fggZZ->Close();

}

//=======================================================================

void storeLDDistribution(){

  makeTemplate("4mu");
  makeTemplate("4e");
  makeTemplate("2e2mu");

}
