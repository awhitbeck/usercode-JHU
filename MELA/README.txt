Content:
--------
PDFs:
RooXZsZs_5D.cxx  RooXZsZs_5D.h 
-> 5D PDF for signal

src:
AngularPdfFactory.cc  
-> Utility class to initialize properly the 5D signal PDF

datafiles:
my8DTemplateNotNorm.root
-> 8D template PDF for qq->ZZ background (m2>20 GeV, mZZ 110-180 GeV)

scripts:
MELA.C  
-> script to evaluate MELA on a given sample (root tree)

Instructions:
------------
root 
 gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L ../PDFs/RooXZsZs_5D.cxx+
.L ../src/AngularPdfFactory.cc+
.L MELA.C+

addDtoTree("nameOfYourFile")

Where nameOfYourFiles is the name of the file which contain
your tree. A new file will be created containing
a new tree where the value of LD has been added.

Notice the tree format should be:
   sigTree->SetBranchAddress("z1mass",&m1);
  sigTree->SetBranchAddress("z2mass",&m2);
  sigTree->SetBranchAddress("zzmass",&mzz);
  sigTree->SetBranchAddress("costheta1",&h1); 
  sigTree->SetBranchAddress("costheta2",&h2);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("phi",&phi);  
  sigTree->SetBranchAddress("phistar1",&phi1);

Or you may change the macro addDtoTree to adapt to the
format of your tree