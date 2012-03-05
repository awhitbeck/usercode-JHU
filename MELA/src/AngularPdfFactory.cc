#include "RooRealVar.h"
#include "../PDFs/RooXZsZs_5D.cxx"
#include <cmath>

class AngularPdfFactory{

public:

  RooRealVar* mZ;     
  RooRealVar* mX;     
  RooRealVar* gamZ;   
    
  RooRealVar* a1Val;  
  RooRealVar* phi1Val;
  RooRealVar* a2Val;  
  RooRealVar* phi2Val;
  RooRealVar* a3Val;  
  RooRealVar* phi3Val;
    
  RooRealVar* R1Val;  
  RooRealVar* R2Val;  
  
  RooXZsZs_5D *PDF;

  int modelIndex;  //0 - SM Higgs, 1 - PS Higgs

  AngularPdfFactory(RooRealVar* m1,RooRealVar* m2,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* mZZ){

    // Parameters
    mZ     = new RooRealVar("mZ","mZ",91.188);
    gamZ   = new RooRealVar("gamZ","gamZ",2.5);
           
    a1Val  = new RooRealVar("a1Val","a1Val",1);
    phi1Val= new RooRealVar("phi1Val","phi1Val",0);
    a2Val  = new RooRealVar("a2Val","a2Val",0);
    phi2Val= new RooRealVar("phi2Val","phi2Val",0);
    a3Val  = new RooRealVar("a3Val","a3Val",0);
    phi3Val= new RooRealVar("phi3Val","phi3Val",0);
           
    R1Val  = new RooRealVar("R1Val","R1Val",0.15);
    R2Val  = new RooRealVar("R2Val","R2Val",0.15);

    PDF = new RooXZsZs_5D("PDF","PDF",*m1,*m2,*h1,*h2,*Phi,*a1Val,*phi1Val,*a2Val,*phi2Val,*a3Val,*phi3Val,*mZ,*gamZ,*mZZ,*R1Val,*R2Val);

  };

  ~AngularPdfFactory(){

    delete mZ;
    delete gamZ;
    delete a1Val;
    delete phi1Val;
    delete a2Val;
    delete phi2Val;
    delete a3Val;
    delete phi3Val;
    delete R1Val;
    delete R2Val;
    
    delete PDF;

  };

  void makeSMHiggs(){
    a1Val->setVal(1.0);
    phi1Val->setVal(0.0);
    a2Val->setVal(0.0);
    phi2Val->setVal(0.0);
    a3Val->setVal(0.0);
    phi3Val->setVal(0.0);
    modelIndex=0;
  };

    void makePSHiggs(){
    a1Val->setVal(0.0);
    phi1Val->setVal(0.0);
    a2Val->setVal(0.0);
    phi2Val->setVal(0.0);
    a3Val->setVal(1.0);
    phi3Val->setVal(0.0);
    modelIndex=1;
  };

  void makeParamsConst(bool yesNo=true){
    if(yesNo){
      a1Val->setConstant(kTRUE);
      phi1Val->setConstant(kTRUE);
      a2Val->setConstant(kTRUE);
      phi2Val->setConstant(kTRUE);
      a3Val->setConstant(kTRUE);
      phi3Val->setConstant(kTRUE);
      gamZ->setConstant(kTRUE);
      mZ->setConstant(kTRUE);
      R1Val->setConstant(kTRUE);
      R2Val->setConstant(kTRUE);
    }else{
      a1Val->setConstant(kFALSE);
      phi1Val->setConstant(kFALSE);
      a2Val->setConstant(kFALSE);
      phi2Val->setConstant(kFALSE);
      a3Val->setConstant(kFALSE);
      phi3Val->setConstant(kFALSE);
      gamZ->setConstant(kFALSE);
      mZ->setConstant(kFALSE);
      R1Val->setConstant(kFALSE);
      R2Val->setConstant(kFALSE);
    }
  };

  double getVal(double mZZ){

    double Norm[70];
    
    if(modelIndex==0){  // SMHiggs

      Norm[0] =4.314 ;
      Norm[1] =5.10233 ;
      Norm[2] =6.02799 ;
      Norm[3] =7.11033 ;
      Norm[4] =8.37067 ;
      Norm[5] =9.83238 ;
      Norm[6] =11.5211 ;
      Norm[7] =13.4647 ;
      Norm[8] =15.6938 ;
      Norm[9] =18.2416 ;
      Norm[10] =21.1442 ;
      Norm[11] =24.4408 ;
      Norm[12] =28.1742 ;
      Norm[13] =32.3904 ;
      Norm[14] =37.1394 ;
      Norm[15] =42.4756 ;
      Norm[16] =48.4575 ;
      Norm[17] =55.1486 ;
      Norm[18] =62.6178 ;
      Norm[19] =70.9392 ;
      Norm[20] =80.1935 ;
      Norm[21] =90.4677 ;
      Norm[22] =101.856 ;
      Norm[23] =114.461 ;
      Norm[24] =128.393 ;
      Norm[25] =143.771 ;
      Norm[26] =160.728 ;
      Norm[27] =179.402 ;
      Norm[28] =199.95 ;
      Norm[29] =222.537 ;
      Norm[30] =247.346 ;
      Norm[31] =274.576 ;
      Norm[32] =304.443 ;
      Norm[33] =337.186 ;
      Norm[34] =373.064 ;
      Norm[35] =412.361 ;
      Norm[36] =455.393 ;
      Norm[37] =502.504 ;
      Norm[38] =554.076 ;
      Norm[39] =610.53 ;
      Norm[40] =672.335 ;
      Norm[41] =740.01 ;
      Norm[42] =814.137 ;
      Norm[43] =895.364 ;
      Norm[44] =984.421 ;
      Norm[45] =1082.13 ;
      Norm[46] =1189.42 ;
      Norm[47] =1307.34 ;
      Norm[48] =1437.1 ;
      Norm[49] =1580.07 ;
      Norm[50] =1737.85 ;
      Norm[51] =1912.28 ;
      Norm[52] =2105.5 ;
      Norm[53] =2320.04 ;
      Norm[54] =2558.88 ;
      Norm[55] =2825.59 ;
      Norm[56] =3124.46 ;
      Norm[57] =3460.69 ;
      Norm[58] =3840.71 ;
      Norm[59] =4272.51 ;
      Norm[60] =4766.18 ;
      Norm[61] =5334.7 ;
      Norm[62] =5995. ;
      Norm[63] =6769.71 ;
      Norm[64] =7689.76 ;
      Norm[65] =8798.6 ;
      Norm[66] =10159.2 ;
      Norm[67] =11866.3 ;
      Norm[68] =14067.6 ;
      Norm[69] =17002.5 ;

    }else if(modelIndex==1){  //PS Higgs

      Norm[0] =0.0441173 ;
      Norm[1] =0.0513195 ;
      Norm[2] =0.0598708 ;
      Norm[3] =0.070012 ;
      Norm[4] =0.0820181 ;
      Norm[5] =0.0962008 ;
      Norm[6] =0.112912 ;
      Norm[7] =0.132548 ;
      Norm[8] =0.155552 ;
      Norm[9] =0.182418 ;
      Norm[10] =0.213695 ;
      Norm[11] =0.249993 ;
      Norm[12] =0.291985 ;
      Norm[13] =0.340414 ;
      Norm[14] =0.396098 ;
      Norm[15] =0.459934 ;
      Norm[16] =0.532906 ;
      Norm[17] =0.61609 ;
      Norm[18] =0.710665 ;
      Norm[19] =0.817915 ;
      Norm[20] =0.939243 ;
      Norm[21] =1.07618 ;
      Norm[22] =1.23038 ;
      Norm[23] =1.40366 ;
      Norm[24] =1.598 ;
      Norm[25] =1.81553 ;
      Norm[26] =2.05859 ;
      Norm[27] =2.32971 ;
      Norm[28] =2.63166 ;
      Norm[29] =2.96743 ;
      Norm[30] =3.34029 ;
      Norm[31] =3.75381 ;
      Norm[32] =4.21186 ;
      Norm[33] =4.71867 ;
      Norm[34] =5.27887 ;
      Norm[35] =5.8975 ;
      Norm[36] =6.58008 ;
      Norm[37] =7.33266 ;
      Norm[38] =8.16187 ;
      Norm[39] =9.075 ;      
      Norm[40] =10.0801 ;
      Norm[41] =11.1859 ;
      Norm[42] =12.4023 ;
      Norm[43] =13.7401 ;
      Norm[44] =15.2112 ;
      Norm[45] =16.829 ;
      Norm[46] =18.6084 ;
      Norm[47] =20.566 ;
      Norm[48] =22.7206 ;
      Norm[49] =25.0933 ;
      Norm[50] =27.7078 ;
      Norm[51] =30.5913 ;
      Norm[52] =33.7745 ;
      Norm[53] =37.293 ;
      Norm[54] =41.1876 ;
      Norm[55] =45.5056 ;
      Norm[56] =50.3023 ;
      Norm[57] =55.6428 ;
      Norm[58] =61.6043 ;
      Norm[59] =68.2793 ;
      Norm[60] =75.7794 ;
      Norm[61] =84.2417 ;
      Norm[62] =93.836 ;
      Norm[63] =104.776 ;
      Norm[64] =117.338 ;
      Norm[65] =131.88 ;
      Norm[66] =148.887 ;
      Norm[67] =169.028 ;
      Norm[68] =193.257 ;
      Norm[69] =222.994 ;

    }else{

      cout << "ERROR: modelIndex not set... Did you make sure to pick a specific model??" << endl;
      return PDF->getVal()/1.0;

    }

    if((int)floor(mZZ-110)>69){

      cout << "Normalization is not available for this value of mZZ: " << (int)floor(mZZ-110) << endl;
      return PDF->getVal()/Norm[69];

    }if((int)floor(mZZ-110)<0){

      cout << "Normalization is not available for this value of mZZ: " << (int)floor(mZZ-110) << endl;
      return PDF->getVal()/Norm[0];

    }

    return PDF->getVal()/Norm[(int)floor(mZZ-110)];

  };

  double getValIntegrOutAngles(RooRealVar* m1,RooRealVar* m2,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* mZZ){
    RooAbsPdf* PDFIntegratedOut =PDF->createProjection(RooArgSet(*h1,*h2,*Phi));
    double norm = PDFIntegratedOut->getNorm(RooArgSet(*m1, *m2, *mZZ));
    cout<<"norm "<<norm<<endl;
    double val = PDFIntegratedOut->getVal()/norm;
    cout<<"val "<<val<<endl;
   return val;
  }

};
