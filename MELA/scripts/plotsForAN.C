#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include <iostream>
#include <vector>

using namespace RooFit;
using namespace std;

void checkAlternativeTemplates(bool signal=true, char* channel="4mu", int i=1, int j=40){

  char temp[100];
  if(signal)
    sprintf(temp,"../datafiles/Dsignal_%s.root",channel);
  else
    sprintf(temp,"../datafiles/Dbackground_qqZZ_%s.root",channel);
 
  TFile* f = new TFile(temp);

  TH2F* h_mzzD = (TH2F*) f->Get("h_mzzD");
  TH2F* h_mzzD_up = (TH2F*) f->Get("h_mzzD_up");
  TH2F* h_mzzD_dn = (TH2F*) f->Get("h_mzzD_dn");

  TH1F* nom = (TH1F*) h_mzzD->ProjectionY("nom",i,j);
  nom->SetLineWidth(2);
  TH1F* up  = (TH1F*) h_mzzD_up->ProjectionY("up",i,j);
  up->SetLineWidth(2);
  up->SetLineColor(2);
  TH1F* dn  = (TH1F*) h_mzzD_dn->ProjectionY("dn",i,j);
  dn->SetLineWidth(2);
  dn->SetLineColor(4);

  TCanvas* can = new TCanvas("can","can",400,400);

  if(up->GetMaximum()>dn->GetMaximum()){
    up->Draw();
    dn->Draw("SAME");
    nom->Draw("SAME");
  }else{
    up->Draw();
    dn->Draw("SAME");
    nom->Draw("SAME");
  }

  if(signal){
    sprintf(temp,"altTempXcheck_signal_%s_%i-%i.png",channel,(int)(i*2+98),(int)(j*2+101));
    can->SaveAs(temp);
    sprintf(temp,"altTempXcheck_signal_%s_%i-%i.eps",channel,(int)(i*2+98),(int)(j*2+101));
    can->SaveAs(temp);
  }else{
    sprintf(temp,"altTempXcheck_background_%s_%i-%i.png",channel,(int)(i*2+98),(int)(j*2+101));
    can->SaveAs(temp);
    sprintf(temp,"altTempXcheck_background_%s_%i-%i.eps",channel,(int)(i*2+98),(int)(j*2+101));
    can->SaveAs(temp);
  }

}

void checkAllAltTemp(){

  for(int i=1; i<=40; i+=5){
    checkAlternativeTemplates(true,"4mu",i,i+4);  
    checkAlternativeTemplates(true,"4e",i,i+4);  
    checkAlternativeTemplates(true,"2e2mu",i,i+4);  
    checkAlternativeTemplates(false,"4mu",i,i+4);  
    checkAlternativeTemplates(false,"4e",i,i+4);  
    checkAlternativeTemplates(false,"2e2mu",i,i+4);  
  }

}

void anglesAndMasses(int index=0){

  double lowM[3] = {100,180,300};
  double highM[3]= {135,220,500};
  int mH[3]   = {120,200,400};

  char fileName[150];
  sprintf(fileName,"../datafiles/ZZ*AnalysisTree_H%i_withDiscriminants.root",mH[index]);
    

  TChain* sigTree =new TChain("angles");
  sigTree->Add(fileName);
  TChain* bkgTree =new TChain("angles");
  bkgTree->Add("../datafiles/ZZ*AnalysisTree_ZZTo*_withDiscriminants.root");

  if(!sigTree || sigTree->GetEntries()<=0 || !bkgTree || bkgTree->GetEntries()<=0) 
    return;

  TH1F* hsig_m1 = new TH1F("hsig_m1",";m_{Z1};",35,50,120);
  TH1F* hbkg_m1 = new TH1F("hbkg_m1",";m_{Z1};",35,50,120);
  hsig_m1->Sumw2();
  hbkg_m1->Sumw2();
  TH1F* hsig_m2 = new TH1F("hsig_m2",";m_{Z2};",55,10,120);
  TH1F* hbkg_m2 = new TH1F("hbkg_m2",";m_{Z2};",55,10,120);
  hsig_m2->Sumw2();
  hbkg_m2->Sumw2();
  TH1F* hsig_h2 = new TH1F("hsig_h2",";cos#theta_{2};",20,-1,1);
  TH1F* hbkg_h2 = new TH1F("hbkg_h2",";cos#theta_{2};",20,-1,1);
  hsig_h2->Sumw2();
  hbkg_h2->Sumw2();
  TH1F* hsig_h1 = new TH1F("hsig_h1",";cos#theta_{1};",20,-1,1);
  TH1F* hbkg_h1 = new TH1F("hbkg_h1",";cos#theta_{1};",20,-1,1);
  hsig_h1->Sumw2();
  hbkg_h1->Sumw2();
  TH1F* hsig_hs = new TH1F("hsig_hs",";cos#theta^{*};",20,-1,1);
  TH1F* hbkg_hs = new TH1F("hbkg_hs",";cos#theta^{*};",20,-1,1);
  hsig_hs->Sumw2();
  hbkg_hs->Sumw2();
  TH1F* hsig_phi = new TH1F("hsig_phi",";#Phi;",20,-3.1415,3.1415);
  TH1F* hbkg_phi = new TH1F("hbkg_phi",";#Phi;",20,-3.1415,3.1415);
  hsig_phi->Sumw2();
  hbkg_phi->Sumw2();
  TH1F* hsig_phi1 = new TH1F("hsig_phi1",";#Phi_{1};",20,-3.1415,3.1415);
  TH1F* hbkg_phi1 = new TH1F("hbkg_phi1",";#Phi_{1};",20,-3.1415,3.1415);
  hsig_phi1->Sumw2();
  hbkg_phi1->Sumw2();

  float mzz,h1,h2,hs,phi,phi1,m1,m2,w;

  sigTree->SetBranchAddress("zzmass",&mzz);
  sigTree->SetBranchAddress("z1mass",&m1);
  sigTree->SetBranchAddress("z2mass",&m2);
  sigTree->SetBranchAddress("costheta1",&h1);
  sigTree->SetBranchAddress("costheta2",&h2);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("phi",&phi);
  sigTree->SetBranchAddress("phistar1",&phi1);
  sigTree->SetBranchAddress("MC_weight",&w);

  bkgTree->SetBranchAddress("zzmass",&mzz);
  bkgTree->SetBranchAddress("z1mass",&m1);
  bkgTree->SetBranchAddress("z2mass",&m2);
  bkgTree->SetBranchAddress("costheta1",&h1);
  bkgTree->SetBranchAddress("costheta2",&h2);
  bkgTree->SetBranchAddress("costhetastar",&hs);
  bkgTree->SetBranchAddress("phi",&phi);
  bkgTree->SetBranchAddress("phistar1",&phi1);
  bkgTree->SetBranchAddress("MC_weight",&w);
  
  for(int i=0; i<sigTree->GetEntries(); i++){

    sigTree->GetEntry(i);
    
    //cout << mzz << " " << m2 << " " << m1 << " " << h1 << " " << h2 << " " << hs << endl;

    if(mzz>lowM[index]&&mzz<highM[index]){

      hsig_m1->Fill(m1,w);
      hsig_m2->Fill(m2,w);
      hsig_h1->Fill(h1,w);
      hsig_h2->Fill(h2,w);
      hsig_hs->Fill(hs,w);
      hsig_phi->Fill(phi,w);
      hsig_phi1->Fill(phi1,w);

    }

  }

  vector<TH1F*> histoSig;
  histoSig.push_back(hsig_m1);
  histoSig.push_back(hsig_m2);
  histoSig.push_back(hsig_h1);
  histoSig.push_back(hsig_h2);
  histoSig.push_back(hsig_hs);
  histoSig.push_back(hsig_phi);
  histoSig.push_back(hsig_phi1);
  
  for(int i=0; i<bkgTree->GetEntries(); i++){
    
    bkgTree->GetEntry(i);

      //cout << mzz << " " << m2 << " " << m1 << " " << h1 << " " << h2 << " " << hs << endl;
      
      if(mzz>lowM[index]&&mzz<highM[index]){
	
	hbkg_m1->Fill(m1,w);
	hbkg_m2->Fill(m2,w);
	hbkg_h1->Fill(h1,w);
	hbkg_h2->Fill(h2,w);
	hbkg_hs->Fill(hs,w);
	hbkg_phi->Fill(phi,w);
	hbkg_phi1->Fill(phi1,w);
	
      }

    }

  vector<TH1F*> histoBkg;
  histoBkg.push_back(hbkg_m1);
  histoBkg.push_back(hbkg_m2);
  histoBkg.push_back(hbkg_h1);
  histoBkg.push_back(hbkg_h2);
  histoBkg.push_back(hbkg_hs);
  histoBkg.push_back(hbkg_phi);
  histoBkg.push_back(hbkg_phi1);

  for(int i=0; i<histoBkg.size(); i++){

    histoSig[i]->Scale(1/histoSig[i]->Integral());
    histoSig[i]->GetYaxis()->SetRangeUser(0,1.3*(histoSig[i]->GetMaximum()+histoSig[i]->GetBinError(histoSig[i]->GetMaximumBin())));
    histoSig[i]->SetLineColor(2);
    histoSig[i]->SetLineWidth(2);
    
    histoBkg[i]->Scale(1/histoBkg[i]->Integral());
    histoBkg[i]->GetYaxis()->SetRangeUser(0,1.3*(histoSig[i]->GetMaximum()+histoSig[i]->GetBinError(histoSig[i]->GetMaximumBin())));
    histoBkg[i]->SetLineColor(4);
    histoBkg[i]->SetLineWidth(2);
    histoBkg[i]->SetLineStyle(2);

  }

  TCanvas *can_m1 = new TCanvas("can_m1","can_m1",400,400);
  TCanvas *can_m2 = new TCanvas("can_m2","can_m2",400,400);
  TCanvas *can_h1 = new TCanvas("can_h1","can_h1",400,400);
  TCanvas *can_h2 = new TCanvas("can_h2","can_h2",400,400);
  TCanvas *can_hs = new TCanvas("can_hs","can_hs",400,400);
  TCanvas *can_phi = new TCanvas("can_phi","can_phi",400,400);
  TCanvas *can_phi1 = new TCanvas("can_phi1","can_phi1",400,400);
  
  vector<TCanvas*> canvases;
  canvases.push_back(can_m1);
  canvases.push_back(can_m2);
  canvases.push_back(can_h1);
  canvases.push_back(can_h2);
  canvases.push_back(can_hs);
  canvases.push_back(can_phi);
  canvases.push_back(can_phi1);
  
  vector<string> names;
  names.push_back("m1");
  names.push_back("m2");
  names.push_back("h1");
  names.push_back("h2");
  names.push_back("hs");
  names.push_back("phi");
  names.push_back("phi1");

  char temp[50];

  for(int i=0; i<canvases.size(); i++){

    canvases[i]->cd();
    
    if(histoBkg[i]->GetMaximum()>histoSig[i]->GetMaximum()){
      histoBkg[i]->Draw("HISTE");
      histoSig[i]->Draw("SAMEHISTE");
    }else{
      histoSig[i]->Draw("HISTE");
      histoBkg[i]->Draw("SAMEHISTE");
    }
    
    sprintf(temp,"compareSignalvsBackground_H%i_%s.eps",mH[index],names[i].c_str());
    canvases[i]->SaveAs(temp);

    sprintf(temp,"compareSignalvsBackground_H%i_%s.png",mH[index],names[i].c_str());
    canvases[i]->SaveAs(temp);
    
  }

}

void MELAtemplate(char* channel="4mu",bool lowMass=true){

  char fileName[150];
  sprintf(fileName,"../datafiles/Dsignal_%s.root",channel);
  TFile* sigFile = new TFile(fileName);
  sprintf(fileName,"../datafiles/Dbackground_%s.root",channel);
  TFile* bkgFile = new TFile(fileName);

  TH2F* sigTemplate = (TH2F*) sigFile->Get("h_mzzD");
  TH2F* bkgTemplate = (TH2F*) bkgFile->Get("h_mzzD");
  
  if(lowMass){
    sigTemplate->GetXaxis()->SetRangeUser(100,180);
    bkgTemplate->GetXaxis()->SetRangeUser(100,180);
  }else{
    sigTemplate->GetXaxis()->SetRangeUser(180,800);
    bkgTemplate->GetXaxis()->SetRangeUser(180,800);
  }
 
  sigTemplate->GetXaxis()->SetTitle("m_{4l}");
  sigTemplate->GetYaxis()->SetTitle("D");
  bkgTemplate->GetXaxis()->SetTitle("m_{4l}");
  bkgTemplate->GetYaxis()->SetTitle("D");
  

  TCanvas* canSig = new TCanvas("canSig","canSig",400,400);
  sigTemplate->Draw("COL");
  TCanvas* canBkg = new TCanvas("canBkg","canBkg",400,400);
  bkgTemplate->Draw("COL");

  if(lowMass)
    sprintf(fileName,"MELAtemplateSmooth_signal_%s_lowMass.eps",channel);
  else
    sprintf(fileName,"MELAtemplateSmooth_signal_%s_highMass.eps",channel);

  canSig->SaveAs(fileName);

  if(lowMass)
    sprintf(fileName,"MELAtemplateSmooth_background_%s_lowMass.eps",channel);
  else
    sprintf(fileName,"MELAtemplateSmooth_background_%s_highMass.eps",channel);

  canBkg->SaveAs(fileName);

}


// 2D datacards:
// /scratch0/hep/whitbeck/4lHelicity/Combination/2012higgsReview2D/workspaces/

void cocktailSyst2(int mass, char* channel, char* tempFileName, double lowM=100, double highM=1000){

  gSystem->Load("libHiggsAnalysisCombinedLimit.so");

  // get workspace

  char fileName[100];  
  sprintf(fileName,"/scratch0/hep/whitbeck/4lHelicity/Combination/2012higgsReview2D/workspaces/hzz4l_%sS.%i.0.input.root",channel,mass);
  
  cout << fileName << endl;

  TFile* wspFile = new TFile(fileName);  
  RooWorkspace* w = (RooWorkspace*) wspFile->Get("w");

  if(!w) 
    return;

  w->var("CMS_zz4l_mass")->SetName("zzmass");

  cout << "changed names" << w->var("zzmass") << endl;

  RooRealVar MC_weight("MC_weight","MC_weight",0,1000);

  RooArgSet obs(*(w->var("zzmass")),*(w->var("melaLD")));
  RooArgSet obsW(*(w->var("zzmass")),*(w->var("melaLD")),MC_weight);
  
  cout << " initialized argset " << endl;

  // load MC 
  if(strcmp(channel,"2e2mu")==0)
    channel="2mu2e";

  sprintf(fileName,"../datafiles/ZZ%sAnalysisTree_H%i_withDiscriminants.root",channel,mass);
  cout << fileName << endl;
  TChain* treeMC = new TChain("angles");
 
  treeMC->Add(fileName);
  
   if(!treeMC || treeMC->GetEntries()<=0 )
    return ;

   cout << treeMC->GetEntries() << endl;

   char cutString[100]="";
   sprintf(cutString,"zzmass>%i&&zzmass<%i",(int)lowM,(int)highM);

   RooDataSet* dataMC = new RooDataSet("dataMC","dataMC",treeMC,obsW,cutString,"MC_weight");  // need to first load workspace and change names
   cout << "dataset size: " << dataMC->sumEntries() << endl;

   // make 2D model
   
   TFile* tempFile = new TFile(tempFileName);
   TH2F* temp = (TH2F*) tempFile->Get("h_mzzD");
   TH2F* oldTemp = (TH2F*) tempFile->Get("oldTemp");
   
   cout << "loaded templates " << temp << " " << oldTemp << endl;
   
   RooDataHist* sigTempDataHist = new RooDataHist("sigTempDataHist","sigTempDataHist",obs,temp);
   RooDataHist* oldSigTempDataHist = new RooDataHist("oldSigTempDataHist","oldSigTempDataHist",obs,oldTemp);
  
   cout << "initialized RooDataHist" << endl;
   
   RooHistPdf* sigTemplatePdf_ggH = new RooHistPdf("sigTemplatePdf_ggH","sigTemplatePdf_ggH",obs,*sigTempDataHist);
   RooHistPdf* oldSigTemplatePdf_ggH = new RooHistPdf("oldSigTemplatePdf_ggH","oldSigTemplatePdf_ggH",obs,*sigTempDataHist);
   
   cout << "initialized RooHistPdf " << sigTemplatePdf_ggH << " " << oldSigTemplatePdf_ggH << endl;
   
   RooProdPdf *sig2d_ggH;
   RooProdPdf *oldSig2d_ggH;
   if(mass<170){
     sig2d_ggH =  new RooProdPdf("sig2d_ggH","sig2d_ggH",*(w->pdf("signalCB_ggH")),RooFit::Conditional(*sigTemplatePdf_ggH,RooArgSet(*w->var("melaLD"))));
     oldSig2d_ggH =  new RooProdPdf("oldSig2d_ggH","oldSig2d_ggH",*(w->pdf("signalCB_ggH")),RooFit::Conditional(*sigTemplatePdf_ggH,RooArgSet(*w->var("melaLD"))));
   }
   else{
     sig2d_ggH =  new RooProdPdf("sig2d_ggH","sig2d_ggH",*(w->pdf("sig_ggH")),RooFit::Conditional(*sigTemplatePdf_ggH,RooArgSet(*w->var("melaLD"))));
     oldSig2d_ggH =  new RooProdPdf("oldSig2d_ggH","oldSig2d_ggH",*(w->pdf("sig_ggH")),RooFit::Conditional(*sigTemplatePdf_ggH,RooArgSet(*w->var("melaLD"))));
   }
   
   cout << "initialized 2D PDFs" << endl;
   
   RooDataSet *sigPDF = sig2d_ggH->generate(obs,100000);
   sigPDF->reduce(cutString);
   RooDataSet *oldSigPDF = oldSig2d_ggH->generate(obs,100000);
   sigPDF->reduce(cutString);
   
   // plot MELA distribution
   
   RooPlot* melaPlot = w->var("melaLD")->frame(30);
   dataMC->plotOn(melaPlot,Rescale(sigPDF->sumEntries()/dataMC->sumEntries()));
   sigPDF->plotOn(melaPlot,MarkerStyle(2),MarkerColor(4));
   oldSigPDF->plotOn(melaPlot,MarkerStyle(3),MarkerColor(2));
   //sig2d_ggH->plotOn(melaPlot);
   
   TCanvas* can = new TCanvas("can","can",400,400);
   melaPlot->Draw();
   sprintf(fileName,"2DpdfProjxCheck_mH%i_%s_%i-%i.png",mass,channel,(int)lowM,(int)highM);
   can->SaveAs(fileName);
   sprintf(fileName,"2DpdfProjxCheck_mH%i_%s_%i-%i.eps",mass,channel,(int)lowM,(int)highM);
   can->SaveAs(fileName);
}

void runAllcocktailSyst(char* channel="4mu"){

  int mH[15]={120,130,140,150,160,170,180,190,200,250,300,350,400,500,600};
  
  char temp[100];
  
  sprintf(temp,"../datafiles/Dsignal_%s.root",channel);
  
  for(int i=0; i<15; i++){
      cocktailSyst2(mH[i],channel,temp);    
  }

}

void cocktailSyst(int index,char* channel){


  int mH[15]={120,130,140,150,160,170,180,190,200,250,300,350,400,500,600};
  double wH[15]={1.,1.,1.,1.,1.,1.,1.,1.,1.43,4.04,8.43,15.2,29.2,68.,123.};

  char temp[150];
  sprintf(temp,"../datafiles/ZZ%sAnalysisTree_H%i_withDiscriminants.root",channel,mH[index]);

  TChain* tree = new TChain("angles");
  tree->Add(temp);

  if(!tree || tree->GetEntries()<=0 )
    return;

  sprintf(temp,"MC_weight*(zzmass>%i&&zzmass<%i)",int(mH[index]-wH[index]*2),int(mH[index]+wH[index]*2));
  
  cout << "drawing tree: " << temp << endl;
  tree->Draw("melaLD>>hMC(30,0,1)",temp);
  TH1F* histoMC = (TH1F*) gDirectory->Get("hMC");
  histoMC->SetLineColor(4);
  histoMC->SetLineStyle(2);
  histoMC->SetLineWidth(2);

  sprintf(temp,"../datafiles/Dsignal_%s.root",channel);
  cout << "loading template: " << temp << endl;
  TFile* tempFile = new TFile(temp);
  TH2F* Htemp = (TH2F*) tempFile->Get("h_mzzD");

  cout << "lower index: " << int(mH[index]/2.-wH[index]-50+1) << " higher index: " << int(mH[index]/2.+wH[index]-50+1)<< endl;

  TH1F* hTemp = (TH1F*) Htemp->ProjectionY("hTemp",int(mH[index]/2.-wH[index]-50+1),int(mH[index]/2.+wH[index]-50+1));
  hTemp->SetLineColor(2);
  hTemp->SetLineWidth(2);
  hTemp->GetXaxis()->SetTitle("D");
  hTemp->GetYaxis()->SetRangeUser(0,hTemp->GetMaximum()>histoMC->GetMaximum()?1.3*hTemp->GetMaximum():histoMC->GetMaximum()*1.3);

  TCanvas* can = new TCanvas("can","can",400,400);

  hTemp->DrawNormalized();
  histoMC->DrawNormalized("SAME");

  TLegend* leg = new TLegend(.2,.7,.5,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  
  leg->AddEntry(hTemp,"cocktail","l");
  sprintf(temp,"H%i",mH[index]);
  leg->AddEntry(histoMC,temp,"l");
  leg->Draw();

  sprintf(temp,"cocktailSyst_mH%i_%s.eps",mH[index],channel);
  can->SaveAs(temp);
  sprintf(temp,"cocktailSyst_mH%i_%s.png",mH[index],channel);
  can->SaveAs(temp);
}

void smoothTemplate(char* fileName="../datafiles/Dbackground_4mu.root",  int effectiveArea=2){

  TFile* file = new TFile(fileName);

  TH2F* oldTemp = (TH2F*) file->Get("h_mzzD");
  oldTemp->GetXaxis()->SetRangeUser(180,800);
  TH2F* newTemp = new TH2F(*oldTemp);

  TCanvas* oldCan = new TCanvas("oldCan","oldCan",400,400);
  oldTemp->Draw("COLZ");
  
  TCanvas* newCan = new TCanvas("newCan","newCan",400,400);
    
  double average=0;
  int nBins=0;

  for(int i=41; i<=400; i++){
    for(int j=1; j<=30; j++){

      for(int a=-effectiveArea; a<=effectiveArea; a++){
	for(int b=-effectiveArea; b<=effectiveArea; b++){
	  if( i+a<41 || i+a>400 || j+b<1 || j+b>30 ) continue;
	  average+=oldTemp->GetBinContent(i+a,j+b);
	  nBins++;
	}
      }

      newTemp->SetBinContent(i,j,average/nBins);
      average=0;
      nBins=0;

    }
  }

  double norm=0;

  for(int i=41; i<=400; i++){
    for(int j=1; j<=30; j++){
      norm+=newTemp->GetBinContent(i,j);
    }

    for(int j=1; j<=30; j++){
      newTemp->SetBinContent(i,j,newTemp->GetBinContent(i,j)/norm);
    }

    norm=0;

  }

  newTemp->Draw("COLZ");

}

void crossCheckSmoothSlices(char* channel="4mu", bool isSig=true, int start=340, int end=350){

  char fileName[100];
  if(isSig)
    sprintf(fileName,"../datafiles/Dsignal_%s.root",channel);
  else
    sprintf(fileName,"../datafiles/Dbackground_%s.root",channel);

  TFile* file = new TFile(fileName);

  TH2F* oldTemp = (TH2F*) file->Get("oldTemp");
  TH2F* newTemp = (TH2F*) file->Get("h_mzzD");
  
  TH1F* oldHisto = (TH1F*) oldTemp->ProjectionY("old",start,end);
  oldHisto->SetLineColor(4);
  oldHisto->SetLineWidth(2);
  TH1F* newHisto = (TH1F*) newTemp->ProjectionY("new",start,end);
  newHisto->SetLineColor(2);
  newHisto->SetLineWidth(2);
  newHisto->SetLineStyle(2);

  oldHisto->Draw();
  newHisto->Draw("SAME");

}
