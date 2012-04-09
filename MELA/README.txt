==================
version notes:  ||
==================

V00-00-01 - 

first stable version.  Low mass LD only.  

-------------------

V00-00-02 - 

corrected bug in calculateAngles() and added protection against 
possibility of m2>m1 (pdfs implicitly assume m1>m2 below threshold)

LD included for highmass (>180), PDFs still only generated for 
100<mZZ<180.

-------------------

V00-00-03 - 

8D template PDF for qq->ZZ background has been extended to 
80<mZZ<185

PDF can now be generated for arbitrary values of mZZ by setting 
global variables mZZbins, lowMzz, and highMzz.  Low m2 cut can 
be set with lowM2.  

==================

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
-> 8D template PDF for qq->ZZ background (m2>4 GeV, mZZ 80-185 GeV)

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

note: by defualt this function only write events which pass the 
following (loose) cuts:

80<mZZ<1000
mz2>4

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MELA.C contains also a macro ("calculateAngles") which computes the 5 angles starting from the 4-momenta of the Higgs, the 2 Zeds and the 4 Leptons

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generating 2D PDF - 

MELA.C can also be used to generate a 2D pdf: mZZ vs D.
The mZZ projection is hard coded in RooMELAModel*.tpl. 
This can be changed to your favorite mZZ projection.  The
D portion is configured such that the projection onto 
mZZ is exactly the function you input. Ranges to be used 
for mZZ and mZ2 can be set via the following global variables:

mZZbins
lowMzz
highMzz
lowM2

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