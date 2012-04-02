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
.L ../PDFs/RooqqZZ_JHU.cxx+
.L MELA.C+

addDtoTree("nameOfYourFile")

Where nameOfYourFiles.root is the name of the file which contain
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

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MELA.C contains also a macro ("calculateAngles") which computes the 5 angles starting from the 4-momenta of the Higgs, the 2 Zeds and the 4 Leptons

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MELA.C can also be used to generate a 2D pdf: mZZ vs D.
The mZZ projection is hard coded in RooMELAModel*.tpl. 
This can be changed to your favorite mZZ projection.  The
D portion is configured such that the projection onto 
mZZ is exactly the function you input.  

To generate signal and background PDFs:

root
 gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L ../PDFs/RooXZsZs_5D.cxx+
.L ../src/AngularPdfFactory.cc+
.L ../PDFs/RooqqZZ_JHU.cxx+
.L MELA.C+
//store D vs mZZ template for signal
storeLDDistribution(true,"mySignalFile.root")  
//store D vs mZZ template for background
storeLDDistribution(false,"myBackgroundFile.root")  
genMELApdf(true)
genMELApdf(false)

/* Note the input file names are passed to 
TChain::Add(), so wild cards can be used.  The
tree name is assumed to be angles but this can
be changed in LDDistributionSignal and 
LDDistributionBackground */

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~