#include "HHKinFit2/HHKinFit2/HHTauTauEventGenerator.h"

#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "HHKinFit2/HHKinFit2/HHKinFit.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraintEHardM.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraintLikelihood.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectMET.h"
#include "HHKinFit2/HHKinFit2/HHFitObject.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectComposite.h"
#include "HHKinFit2/HHKinFit2/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"

#include "HHKinFit2/HHKinFit2/HHKinFitSingleHMaster.h"

#include <time.h>

using HHKinFit2::HHTauTauEventGenerator;
using HHKinFit2::HHEnergyRangeException;
using HHKinFit2::HHFitObjectE;
using HHKinFit2::HHFitObjectEConstM;
using HHKinFit2::HHFitObjectMET;
using HHKinFit2::HHFitObjectComposite;
using HHKinFit2::HHFitObject;
using HHKinFit2::HHFitConstraint;
using HHKinFit2::HHFitConstraintEHardM;
using HHKinFit2::HHFitConstraint4Vector;
using HHKinFit2::HHFitConstraintLikelihood;
using HHKinFit2::HHLorentzVector;

#include <iostream>
#include <chrono>

int main(int argc, char* argv[])
{
  auto start = std::chrono::high_resolution_clock::now();

  double timer0=0;
  double timer1=0;
  double timer2=0;

  TF1* PDF1 = new TF1("PDF1","2*x",0,2);
  TF1* PDF2 = new TF1("PDF2","2*x",0,2);

  TMatrixD covarmatrix(2,2);
  covarmatrix[0][0]=130;
  covarmatrix[0][1]=0;
  covarmatrix[1][0]=0;
  covarmatrix[1][1]=130;
  HHTauTauEventGenerator testgenerator(PDF1,PDF2,covarmatrix);

  int events = 1000000;
  TH1D h_deltachi2("h_dchi2","deltachi2",200,-1,1);
  TH1D h_deltaE("h_dE","deltaE",200,-0.1,0.1);
  TH1D h_chi2_1("h_chi2_1","chi2 HHKinFit1",200,0,20);
  TH1D h_chi2_2("h_chi2_2","chi2 HHKinFit2",200,0,20);
  TH1D h_prob_1("h_prob_1","P(chi2) HHKinFit1",50,0,1);
  TH1D h_prob_2("h_prob_2","P(chi2) HHKinFit2",50,0,1);


  for(unsigned int i=0; i<events; i++){
//    std::cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    try{
      testgenerator.generateEvent();
    }
    catch(const HHEnergyRangeException& e){
      //std::cout << e.what() << std::endl;
      i--;
      continue;
    }
    timer0+= std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count();


    //KinFit:
    double mass = testgenerator.getMhiggs();
    HHLorentzVector gentauvis1 = testgenerator.getTau1Vis();
    HHLorentzVector gentauvis2 = testgenerator.getTau2Vis();

    TVector2 ptmiss (testgenerator.getMETwithsigma()[0],testgenerator.getMETwithsigma()[1]);
    TLorentzVector ptmissLV  = TLorentzVector(ptmiss.Px(),ptmiss.Py(),0,sqrt(ptmiss.Px()*ptmiss.Px()+ptmiss.Py()*ptmiss.Py()));

    // std::cout << "input energies:" << std::endl;
    // std::cout << gentauvis1.E() << std::endl;
    // std::cout << gentauvis2.E() << std::endl;
    // std::cout << "........................" << std::endl;

    //prepare tau objects
    start = std::chrono::high_resolution_clock::now();
    HHFitObjectE* tau1 = new HHFitObjectEConstM(gentauvis1);//  visible Tau1-component from HHTauTauEventGenerator
    HHFitObjectE* tau2 = new HHFitObjectEConstM(gentauvis2);//  visible Tau2-component from HHTauTauEventGenerator
    tau1->setLowerFitLimitE(tau1->getInitial4Vector());
    tau1->setUpperFitLimitE(mass,tau2->getInitial4Vector());
    tau2->setLowerFitLimitE(tau2->getInitial4Vector());
    tau2->setUpperFitLimitE(mass,tau1->getInitial4Vector());

    //prepare MET object
    HHFitObjectMET* met = new HHFitObjectMET(ptmiss);//Use Met components from HHTauTauEventGenerator
    met->setCovMatrix(covarmatrix[0][0],covarmatrix[1][1],covarmatrix[1][0]); // set Covarmatrix with Matrix inserted in HHTauTauEventGenerator
    
    //prepare composite object: Higgs
    HHFitObject* higgs  = new HHFitObjectComposite(tau1, tau2, met);
    
    //prepare constraints
    HHFitConstraint* invm = new HHFitConstraintEHardM(tau1, tau2, mass);
    HHFitConstraint* balance = new HHFitConstraint4Vector(higgs, true, true, false, false);
    
    //fit
    HHKinFit2::HHKinFit singlefit;
    singlefit.addFitObjectE(tau1);
    singlefit.addConstraint(invm);
    singlefit.addConstraint(balance);

    try{
      singlefit.fit();
    }
    catch(HHKinFit2::HHEnergyRangeException const& e){ // can happen due to derivative calculation...
      std::cout << e.what() << std::endl;
      i--;
      continue;
    }
    timer1+= std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count();

    double chi2 = singlefit.getChi2();
    double e = tau1->getFit4Vector().E();
    delete(tau1);
    delete(tau2);
    delete(met);
    delete(higgs);
    delete(invm);
    delete(balance);

    if (!((singlefit.getConvergence()==1)||(singlefit.getConvergence()==2))) {
      i--;
      continue;
    }

    //std::cout << "........................" << std::endl;

    //-----------------------------------------------------------
    //intance of fitter master class
    start = std::chrono::high_resolution_clock::now();
    HHKinFitSingleHMaster oldkinfit = HHKinFitSingleHMaster(&gentauvis1,&gentauvis2);
    oldkinfit.setAdvancedBalance(&ptmissLV,covarmatrix);
    oldkinfit.addMhHypothesis(mass);
    oldkinfit.doFullFit();
    timer2+= std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count();

    double oldchi2 = oldkinfit.getBestChi2FullFit();
    double olde = oldkinfit.m_tau1_fitted.E();
    //-----------------------------------------------------------
    // std::cout << "........................" << std::endl;
    // std::cout << "new chi2:  " << chi2 << std::endl;
    // std::cout << "old chi2:  " << oldchi2 << std::endl;
    // std::cout << "new Etau1: " << tau1->getFit4Vector().E() << std::endl;
    // std::cout << "old Etau1: " << oldkinfit.m_tau1_fitted.E() << std::endl;
    // std::cout << "new Etau2: " << tau2->getFit4Vector().E() << std::endl;
    // std::cout << "old Etau2: " << oldkinfit.m_tau2_fitted.E() << std::endl;

    h_deltachi2.Fill((chi2-oldchi2)/oldchi2);
    h_deltaE.Fill((e-olde)/olde);
    h_chi2_1.Fill(oldchi2);
    h_chi2_2.Fill(chi2);
    h_prob_1.Fill(TMath::Prob(oldchi2,1));
    h_prob_2.Fill(TMath::Prob(chi2,1));
    if (i%10000==0) std::cout << i << std::endl;
  }

  TFile f("out.root","RECREATE");
  h_deltachi2.Write();
  h_deltaE.Write();
  h_chi2_1.Write();
  h_chi2_2.Write();
  h_prob_1.Write();
  h_prob_2.Write();

  f.Close();

  std::cout << "event generation: " << timer0/events << "us/event." << std::endl;
  std::cout << "new kinfit: " << timer1/events << "us/event." << std::endl;
  std::cout << "old kinfit: " << timer2/events << "us/event." << std::endl;

  return(0);
}
