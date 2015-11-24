#include <iostream>
#include "HHKinFit2/HHKinFit2/HHKinFit.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraintEHardM.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraintLikelihood.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectEConstBeta.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectMET.h"
#include "HHKinFit2/HHKinFit2/HHFitObject.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectComposite.h"
#include "HHKinFit2/HHKinFit2/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace HHKinFit2;

int main(int argc, char* argv[])
{
  double mh1 = 111.0;
  double mh2 = 222.0;

  //prepare tau objects
  HHFitObjectE* tau1 = new HHFitObjectEConstM(HHLorentzVector(71,0,0,80));
  HHFitObjectE* tau2 = new HHFitObjectEConstM(HHLorentzVector(0,58,0,60));
  tau1->setFitLimitsE(tau1->getInitial4Vector(),mh1,tau2->getInitial4Vector());
  tau2->setFitLimitsE(tau2->getInitial4Vector(),mh1,tau1->getInitial4Vector());

  //prepare bjet objects
  HHFitObjectE* b1 = new HHFitObjectEConstBeta(HHLorentzVector(71,0,0,80));
  HHFitObjectE* b2 = new HHFitObjectEConstBeta(HHLorentzVector(0,58,0,60));
  b1->setFitLimitsE(0,400);
  b2->setFitLimitsE(0,400);
  b1->setCovMatrix(100);
  b2->setCovMatrix(100);

  //prepare MET object
  HHFitObjectMET* met = new HHFitObjectMET(TVector2(10,20));
  met->setCovMatrix(100,100,0);

  //prepare composite object: Higgs
  HHFitObject* higgs  = new HHFitObjectComposite(tau1,tau2,b1,b2,met);
  HHFitObject* higgs1  = new HHFitObjectComposite(tau1,tau2);
  HHFitObject* higgs2  = new HHFitObjectComposite(b1,b2);

  //prepare constraints
  HHFitConstraint* c_invmh1 = new HHFitConstraintEHardM(tau1, tau2, mh1);
  HHFitConstraint* c_invmh2 = new HHFitConstraintEHardM(b1, b2, mh2);
  HHFitConstraint* c_b1 = new HHFitConstraint4Vector(b1, false, false, false, true);
  HHFitConstraint* c_b2 = new HHFitConstraint4Vector(b2, false, false, false, true);
  HHFitConstraint* c_balance = new HHFitConstraint4Vector(higgs, true, true, false, false);

  //fit
  HHKinFit* singlefit = new HHKinFit();

  singlefit->addFitObjectE(tau1);
  singlefit->addFitObjectE(b1);

  singlefit->addConstraint(c_invmh1);
  singlefit->addConstraint(c_invmh2);
  singlefit->addConstraint(c_b1);
  singlefit->addConstraint(c_b2);
  singlefit->addConstraint(c_balance);

  singlefit->fit();

//  higgs->print();
  higgs1->print();
  higgs2->print();

  std::cout << "final chi2: " << singlefit->getChi2() << std::endl;

//  TGraph gr(singlefit->getChi2Function(100));
//  gr.Draw();
//
//  TFile f("out.root","UPDATE");
//  gr.Write();
//  f.Close();
//
//
//  try{
//    std::cout << singlefit->getChi2() << std::endl;
//  }
//  catch(const HHCovarianceMatrixException& e){
//    std::cout << e.what() << std::endl;
//    std::cout << "will fix it manually" << std::endl;
//    met->setCovMatrix(100,100,0);
//    std::cout << fit->getChi2() << std::endl;
//  }
//  catch(...){
//    std::cout << "caught an unexpected exception" << std::endl;
//  }

  return(0);
}
