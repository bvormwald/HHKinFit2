#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

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
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHLimitSettingException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"

#include "HHKinFit2/HHKinFit2/HHKinFitMaster.h"

#include <time.h>
#include <chrono>

using HHKinFit2::HHEnergyRangeException;
using HHKinFit2::HHLimitSettingException;
using HHKinFit2::HHFitObjectE;
using HHKinFit2::HHFitObjectEConstM;
using HHKinFit2::HHFitObjectEConstBeta;
using HHKinFit2::HHFitObjectMET;
using HHKinFit2::HHFitObjectComposite;
using HHKinFit2::HHFitObject;
using HHKinFit2::HHFitConstraint;
using HHKinFit2::HHFitConstraintEHardM;
using HHKinFit2::HHFitConstraint4Vector;
using HHKinFit2::HHFitConstraintLikelihood;
using HHKinFit2::HHLorentzVector;

double GetBjetResolution(double eta, double et);

int main(int argc, char* argv[])
{
  double timer1=0;
  double timer2=0;
  double events=0;
  auto start = std::chrono::high_resolution_clock::now();

  TFile* f = new TFile("/afs/desy.de/user/v/vormwald/workspace/HHKinFit2/examples/ntuple.root","READ");
  if (f->IsZombie()) return(0);
  TTree* t = (TTree*) f->Get("EleTauKinFit");
  
  TH1F* h_mH_2    = new TH1F("h_mH_2",   "HHKinFit2: distribution of mH_{fit};mH_{fit} [GeV];entries/5GeV", 40,200,400);
  TH1F* h_chi2_2  = new TH1F("h_chi2_2", "HHKinFit2: distribution of #chi^{2}_{fit};#chi^{2}_{fit};entries/0.25", 100,0,25);
  TH1F* h_prob_2  = new TH1F("h_prob_2", "HHKinFit2: distribution of P(#chi^{2}_{fit});P;entries/0.02", 50,0,1);

  TH1F* h_mH_1    = new TH1F("h_mH_1",   "HHKinFit1: distribution of mH_{fit};mH_{fit} [GeV];entries/5GeV", 40,200,400);
  TH1F* h_chi2_1  = new TH1F("h_chi2_1", "HHKinFit1: distribution of #chi^{2}_{fit};#chi^{2}_{fit};entries/0.25", 100,0,25);
  TH1F* h_prob_1  = new TH1F("h_prob_1", "HHKinFit1: distribution of P(#chi^{2}_{fit});P;entries/0.02", 50,0,1);
  
  TH1F* h_dmH    = new TH1F("h_dmH",   "relative deviation mH;mH_({HHKinFit1}-mH_{HHKinFit2})/mH_{HHKinFit1};entries", 200,-0.1,0.1);
  TH1F* h_dchi2  = new TH1F("h_dchi2", "relative deviation #chi^{2};#chi^{2}_({HHKinFit1}-#chi^{2}_{HHKinFit2})/#chi^{2}_{HHKinFit1};entries",200,-1,1);
  TH1F* h_dprob  = new TH1F("h_dprib", "relative deviation P(#chi^{2});P(#chi^{2}_({HHKinFit1})-P(#chi^{2}_{HHKinFit2}))/P(#chi^{2}_{HHKinFit1});entries",200,-1,1);


  //Leaf types  
   Int_t           njets;
   Double_t        b1_dR;
   Double_t        b1_csv;
   Double_t        b1_px;
   Double_t        b1_py;
   Double_t        b1_pz;
   Double_t        b1_E;
   Double_t        b2_dR;
   Double_t        b2_csv;
   Double_t        b2_px;
   Double_t        b2_py;
   Double_t        b2_pz;
   Double_t        b2_E;
   Double_t        tauvis1_dR;
   Double_t        tauvis1_px;
   Double_t        tauvis1_py;
   Double_t        tauvis1_pz;
   Double_t        tauvis1_E;
   Double_t        tauvis2_dR;
   Double_t        tauvis2_px;
   Double_t        tauvis2_py;
   Double_t        tauvis2_pz;
   Double_t        tauvis2_E;
   Double_t        met_pt;
   Double_t        met_px;
   Double_t        met_py;
   Double_t        met_cov00;
   Double_t        met_cov10;
   Double_t        met_cov01;
   Double_t        met_cov11;


  // List of branches
   TBranch        *b_njets;
   TBranch        *b_b1_dR;
   TBranch        *b_b1_csv;
   TBranch        *b_b1_px;
   TBranch        *b_b1_py;
   TBranch        *b_b1_pz;
   TBranch        *b_b1_E;
   TBranch        *b_b2_dR;
   TBranch        *b_b2_csv;
   TBranch        *b_b2_px;
   TBranch        *b_b2_py;
   TBranch        *b_b2_pz;
   TBranch        *b_b2_E;
   TBranch        *b_tauvis1_dR;
   TBranch        *b_tauvis1_px;
   TBranch        *b_tauvis1_py;
   TBranch        *b_tauvis1_pz;
   TBranch        *b_tauvis1_E;
   TBranch        *b_tauvis2_dR;
   TBranch        *b_tauvis2_px;
   TBranch        *b_tauvis2_py;
   TBranch        *b_tauvis2_pz;
   TBranch        *b_tauvis2_E;
   TBranch        *b_met_pt;
   TBranch        *b_met_px;
   TBranch        *b_met_py;
   TBranch        *b_met_cov00;
   TBranch        *b_met_cov10;
   TBranch        *b_met_cov01;
   TBranch        *b_met_cov11;
     
   t->SetBranchAddress("Njets", &njets, &b_njets);
   t->SetBranchAddress("b1_dR", &b1_dR, &b_b1_dR);
   t->SetBranchAddress("b1_csv", &b1_csv, &b_b1_csv);
   t->SetBranchAddress("b1_px", &b1_px, &b_b1_px);
   t->SetBranchAddress("b1_py", &b1_py, &b_b1_py);
   t->SetBranchAddress("b1_pz", &b1_pz, &b_b1_pz);
   t->SetBranchAddress("b1_E", &b1_E, &b_b1_E);
   t->SetBranchAddress("b2_dR", &b2_dR, &b_b2_dR);
   t->SetBranchAddress("b2_csv", &b2_csv, &b_b2_csv);
   t->SetBranchAddress("b2_px", &b2_px, &b_b2_px);
   t->SetBranchAddress("b2_py", &b2_py, &b_b2_py);
   t->SetBranchAddress("b2_pz", &b2_pz, &b_b2_pz);
   t->SetBranchAddress("b2_E", &b2_E, &b_b2_E);
   t->SetBranchAddress("tauvis1_dR", &tauvis1_dR, &b_tauvis1_dR);
   t->SetBranchAddress("tauvis1_px", &tauvis1_px, &b_tauvis1_px);
   t->SetBranchAddress("tauvis1_py", &tauvis1_py, &b_tauvis1_py);
   t->SetBranchAddress("tauvis1_pz", &tauvis1_pz, &b_tauvis1_pz);
   t->SetBranchAddress("tauvis1_E", &tauvis1_E, &b_tauvis1_E);
   t->SetBranchAddress("tauvis2_dR", &tauvis2_dR, &b_tauvis2_dR);
   t->SetBranchAddress("tauvis2_px", &tauvis2_px, &b_tauvis2_px);
   t->SetBranchAddress("tauvis2_py", &tauvis2_py, &b_tauvis2_py);
   t->SetBranchAddress("tauvis2_pz", &tauvis2_pz, &b_tauvis2_pz);
   t->SetBranchAddress("tauvis2_E", &tauvis2_E, &b_tauvis2_E);
   t->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   t->SetBranchAddress("met_px", &met_px, &b_met_px);
   t->SetBranchAddress("met_py", &met_py, &b_met_py);
   t->SetBranchAddress("met_cov00", &met_cov00, &b_met_cov00);
   t->SetBranchAddress("met_cov10", &met_cov10, &b_met_cov10);
   t->SetBranchAddress("met_cov01", &met_cov01, &b_met_cov01);
   t->SetBranchAddress("met_cov11", &met_cov11, &b_met_cov11);
  
   Int_t max_nevent = t->GetEntriesFast();
//   Int_t max_nevent = 10;

  // event loop
  for (int i = 0; i < max_nevent; i++) {
    t->GetEntry(i);
    
    //use only btaged jets and check for correct gen match
    if (b1_csv<0.681 || b2_csv<0.681 || njets < 2 || b1_dR > 0.2 || b2_dR > 0.2 || tauvis1_dR > 0.1 || tauvis2_dR > 0.1)   continue;

    double mh1 = 125.0;
    double mh2 = 125.0;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //prepare tau objects
    start = std::chrono::high_resolution_clock::now();

    HHLorentzVector tauin1(tauvis1_px,tauvis1_py,tauvis1_pz,tauvis1_E);
    HHLorentzVector tauin2(tauvis2_px,tauvis2_py,tauvis2_pz,tauvis2_E);
    tauin1.SetMkeepE(1.77682);
    tauin2.SetMkeepE(1.77682);
    HHFitObjectE* tau1 = new HHFitObjectEConstM(tauin1);
    HHFitObjectE* tau2 = new HHFitObjectEConstM(tauin2);

    try{
      tau1->setFitLimitsE(tau1->getInitial4Vector(),mh1,tau2->getInitial4Vector());
      tau2->setFitLimitsE(tau2->getInitial4Vector(),mh1,tau1->getInitial4Vector());
    }
    catch(HHLimitSettingException const& e){
      std::cout << "Exception while setting tau limits:" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "Tau energies are not compatible with invariant mass constraint." << std::endl;
      continue;
    }

    //prepare bjet objects
    HHLorentzVector bin1(b1_px,b1_py,b1_pz,b1_E);
    HHLorentzVector bin2(b2_px,b2_py,b2_pz,b2_E);
    double b1sig = GetBjetResolution(bin1.Eta(),bin1.Et());
    double b2sig = GetBjetResolution(bin2.Eta(),bin2.Et());
//    std::cout << "sigma b1 (inv cov): " << b1sig << " " << 1./pow(b1sig,2) << std::endl;
//    std::cout << "sigma b2 (inv cov): " << b2sig << " " << 1./pow(b2sig,2) << std::endl;
    HHLorentzVector bin1min = bin1;
    if (b1_E-5.0*b1sig<=0) bin1min.SetEkeepBeta(10);
    else bin1min.SetEkeepBeta(b1_E-5.0*b1sig);
    HHLorentzVector bin2min = bin2;
    if (b2_E-5.0*b2sig<=0) bin2min.SetEkeepBeta(10);
    else bin2min.SetEkeepBeta(b2_E-5.0*b2sig);
    HHFitObjectE* b1 = new HHFitObjectEConstBeta(bin1);
    HHFitObjectE* b2 = new HHFitObjectEConstBeta(bin2);
    try{
      b1->setFitLimitsE(bin1min, mh2, bin2min);
      b2->setFitLimitsE(bin2min, mh2, bin1min);
    }
    catch(HHLimitSettingException const& e){
      std::cout << "Exception while setting b-jet limits" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "Bjet energies are not compatible within 5 sigma with invariant mass constraint." << std::endl;
      continue;
    }
    b1->setCovMatrix(pow(b1sig,2));
    b2->setCovMatrix(pow(b2sig,2));

    //prepare MET object
    HHFitObjectMET* met = new HHFitObjectMET(TVector2(met_px,met_py));
    met->setCovMatrix(met_cov00,met_cov11,met_cov10);

    //prepare composite object: Higgs
    HHFitObject* higgs  = new HHFitObjectComposite(tau1,tau2,b1,b2,met);
    HHFitObject* higgs1  = new HHFitObjectComposite(tau1,tau2);
    HHFitObject* higgs2  = new HHFitObjectComposite(b1,b2);

    //prepare constraints
    HHFitConstraint* c_invmh1 = new HHFitConstraintEHardM(tau1, tau2, mh1);
    HHFitConstraint* c_invmh2 = new HHFitConstraintEHardM(b1, b2, mh2);
    HHFitConstraint* c_b1 = new HHFitConstraint4Vector(b1, false, false, false, true);
    HHFitConstraint* c_b2 = new HHFitConstraint4Vector(b2, false, false, false, true);
    HHFitConstraint* c_balance;
    try{
      c_balance = new HHFitConstraint4Vector(higgs, true, true, false, false);
    }
    catch(HHKinFit2::HHCovarianceMatrixException const& e){
      std::cout << e.what() << std::endl;
      std::cout << "replace covariance matrix by diag(100,100)" << std::endl;
      TMatrixD replace(4,4);
      replace(0,0)=100;
      replace(1,1)=100;
      replace.Print();
      higgs->setCovMatrix(replace);
      higgs->print();
      c_balance = new HHFitConstraint4Vector(higgs, true, true, false, false);
    }

    //fit
    HHKinFit2::HHKinFit* singlefit = new HHKinFit2::HHKinFit();

    tau1->setInitStart((tau1->getUpperFitLimitE()+tau1->getLowerFitLimitE())/2);
    tau1->setInitPrecision(0.1);
    tau1->setInitStepWidth(0.1*(tau1->getUpperFitLimitE()-tau1->getLowerFitLimitE()));
    tau1->setInitDirection(1.0);

    b1->setInitStart(b1->getInitial4Vector().E());
    b1->setInitPrecision(0.002*b1->getInitial4Vector().E());
    b1->setInitStepWidth(0.5*b1sig);
    b1->setInitDirection(1.0);

    singlefit->addFitObjectE(tau1);
    singlefit->addFitObjectE(b1);

    singlefit->addConstraint(c_invmh1);
    singlefit->addConstraint(c_invmh2);
    singlefit->addConstraint(c_b1);
    singlefit->addConstraint(c_b2);
    singlefit->addConstraint(c_balance);

    singlefit->fit();

    timer1+= std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count();

    //obtain results from different hypotheses
    Double_t chi2_best = singlefit->getChi2();
    Double_t mh_best = higgs->getFit4Vector().M();
    //singlefit->printChi2();


    /////////////////////////////////////////////////////////////////////////////////////////////////
    start = std::chrono::high_resolution_clock::now();

    TLorentzVector tauin1OLD(tauvis1_px,tauvis1_py,tauvis1_pz,tauvis1_E);
    TLorentzVector tauin2OLD(tauvis2_px,tauvis2_py,tauvis2_pz,tauvis2_E);
    TLorentzVector bin1OLD(b1_px,b1_py,b1_pz,b1_E);
    TLorentzVector bin2OLD(b2_px,b2_py,b2_pz,b2_E);
    TLorentzVector ptmiss  = TLorentzVector(met_px,met_py,0,met_pt);
    TMatrixD metcov(2,2);
    metcov(0,0)=met_cov00;
    metcov(1,0)=met_cov10;
    metcov(0,1)=met_cov01;
    metcov(1,1)=met_cov11;

    HHKinFitMaster kinFits = HHKinFitMaster(&bin1OLD,&bin2OLD,&tauin1OLD,&tauin2OLD);
    kinFits.setAdvancedBalance(&ptmiss,metcov);
    //kinFits.setSimpleBalance(ptmiss.Pt(),10); //alternative which uses only the absolute value of ptmiss in the fit
    kinFits.addMh1Hypothesis(mh1);
    kinFits.addMh2Hypothesis(mh2);
    kinFits.doFullFit();

    timer2+= std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count();

    //obtain results from different hypotheses
    Double_t chi2_best_old = kinFits.getBestChi2FullFit();
    Double_t mh_best_old = kinFits.getBestMHFullFit();

    h_mH_2->Fill(mh_best);
    h_chi2_2->Fill(chi2_best);
    double prob = TMath::Prob(chi2_best, 2);
    h_prob_2->Fill(prob);

    h_mH_1->Fill(mh_best_old);
    h_chi2_1->Fill(chi2_best_old);
    double prob_old = TMath::Prob(chi2_best_old, 2);
    h_prob_1->Fill(prob_old);

    h_dmH->Fill((mh_best_old-mh_best)/mh_best_old);
    h_dchi2->Fill((chi2_best_old-chi2_best)/chi2_best_old);
    h_dprob->Fill((prob_old-prob)/prob_old);
    events++;
  }

  TFile fout("result.root","RECREATE");
  h_mH_1->Write();
  h_chi2_1->Write();
  h_prob_1->Write();
  h_mH_2->Write();
  h_chi2_2->Write();
  h_prob_2->Write();
  h_dmH->Write();
  h_dchi2->Write();
  h_dprob->Write();
  
  std::cout << "new kinfit: " << timer1/events << "us/event." << std::endl;
  std::cout << "old kinfit: " << timer2/events << "us/event." << std::endl;

  return (0);
}


double
GetBjetResolution(double eta, double et){
  double det=0;
  double de=10;

  if(0.000<=abs(eta) && abs(eta)<0.087){
  det = et * (sqrt(0.0686*0.0686 + (1.03/sqrt(et))*(1.03/sqrt(et)) + (1.68/et)*(1.68/et)));
  de = 1.0/sin(2 * atan(exp(-(0.000+0.087)/2))) * det;
  }

  if(0.087<=abs(eta) && abs(eta)<0.174){
  det = et * (sqrt(0.0737*0.0737 + (1.01/sqrt(et))*(1.01/sqrt(et)) + (1.74/et)*(1.74/et)));
  de = 1.0/sin(2 * atan(exp(-(0.087+0.174)/2))) * det;
  }

  if(0.174<=abs(eta) && abs(eta)<0.261){
  det = et * (sqrt(0.0657*0.0657 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (5.16e-06/et)*(5.16e-06/et)));
  de = 1.0/sin(2 * atan(exp(-(0.174+0.261)/2))) * det;
  }

  if(0.261<=abs(eta) && abs(eta)<0.348){
  det = et * (sqrt(0.062*0.062 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (0.000134/et)*(0.000134/et)));
  de = 1.0/sin(2 * atan(exp(-(0.261+0.348)/2))) * det;
  }

  if(0.348<=abs(eta) && abs(eta)<0.435){
  det = et * (sqrt(0.0605*0.0605 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (1.84e-07/et)*(1.84e-07/et)));
  de = 1.0/sin(2 * atan(exp(-(0.348+0.435)/2))) * det;
  }

  if(0.435<=abs(eta) && abs(eta)<0.522){
  det = et * (sqrt(0.059*0.059 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (9.06e-09/et)*(9.06e-09/et)));
  de = 1.0/sin(2 * atan(exp(-(0.435+0.522)/2))) * det;
  }

  if(0.522<=abs(eta) && abs(eta)<0.609){
  det = et * (sqrt(0.0577*0.0577 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (5.46e-06/et)*(5.46e-06/et)));
  de = 1.0/sin(2 * atan(exp(-(0.522+0.609)/2))) * det;
  }

  if(0.609<=abs(eta) && abs(eta)<0.696){
  det = et * (sqrt(0.0525*0.0525 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (4.05e-05/et)*(4.05e-05/et)));
  de = 1.0/sin(2 * atan(exp(-(0.609+0.696)/2))) * det;
  }

  if(0.696<=abs(eta) && abs(eta)<0.783){
  det = et * (sqrt(0.0582*0.0582 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.17e-05/et)*(1.17e-05/et)));
  de = 1.0/sin(2 * atan(exp(-(0.696+0.783)/2))) * det;
  }

  if(0.783<=abs(eta) && abs(eta)<0.870){
  det = et * (sqrt(0.0649*0.0649 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (7.85e-06/et)*(7.85e-06/et)));
  de = 1.0/sin(2 * atan(exp(-(0.783+0.870)/2))) * det;
  }

  if(0.870<=abs(eta) && abs(eta)<0.957){
  det = et * (sqrt(0.0654*0.0654 + (1.1/sqrt(et))*(1.1/sqrt(et)) + (1.09e-07/et)*(1.09e-07/et)));
  de = 1.0/sin(2 * atan(exp(-(0.870+0.957)/2))) * det;
  }

  if(0.957<=abs(eta) && abs(eta)<1.044){
  det = et * (sqrt(0.0669*0.0669 + (1.11/sqrt(et))*(1.11/sqrt(et)) + (1.87e-06/et)*(1.87e-06/et)));
  de = 1.0/sin(2 * atan(exp(-(0.957+1.044)/2))) * det;
  }

  if(1.044<=abs(eta) && abs(eta)<1.131){
  det = et * (sqrt(0.0643*0.0643 + (1.15/sqrt(et))*(1.15/sqrt(et)) + (2.76e-05/et)*(2.76e-05/et)));
  de = 1.0/sin(2 * atan(exp(-(1.044+1.131)/2))) * det;
  }

  if(1.131<=abs(eta) && abs(eta)<1.218){
  det = et * (sqrt(0.0645*0.0645 + (1.16/sqrt(et))*(1.16/sqrt(et)) + (1.04e-06/et)*(1.04e-06/et)));
  de = 1.0/sin(2 * atan(exp(-(1.131+1.218)/2))) * det;
  }

  if(1.218<=abs(eta) && abs(eta)<1.305){
  det = et * (sqrt(0.0637*0.0637 + (1.19/sqrt(et))*(1.19/sqrt(et)) + (1.08e-07/et)*(1.08e-07/et)));
  de = 1.0/sin(2 * atan(exp(-(1.218+1.305)/2))) * det;
  }

  if(1.305<=abs(eta) && abs(eta)<1.392){
  det = et * (sqrt(0.0695*0.0695 + (1.21/sqrt(et))*(1.21/sqrt(et)) + (5.75e-06/et)*(5.75e-06/et)));
  de = 1.0/sin(2 * atan(exp(-(1.305+1.392)/2))) * det;
  }

  if(1.392<=abs(eta) && abs(eta)<1.479){
  det = et * (sqrt(0.0748*0.0748 + (1.2/sqrt(et))*(1.2/sqrt(et)) + (5.15e-08/et)*(5.15e-08/et)));
  de = 1.0/sin(2 * atan(exp(-(1.392+1.479)/2))) * det;
  }

  if(1.479<=abs(eta) && abs(eta)<1.566){
  det = et * (sqrt(0.0624*0.0624 + (1.23/sqrt(et))*(1.23/sqrt(et)) + (2.28e-05/et)*(2.28e-05/et)));
  de = 1.0/sin(2 * atan(exp(-(1.479+1.566)/2))) * det;
  }

  if(1.566<=abs(eta) && abs(eta)<1.653){
  det = et * (sqrt(0.0283*0.0283 + (1.25/sqrt(et))*(1.25/sqrt(et)) + (4.79e-07/et)*(4.79e-07/et)));
  de = 1.0/sin(2 * atan(exp(-(1.566+1.653)/2))) * det;
  }

  if(1.653<=abs(eta) && abs(eta)<1.740){
  det = et * (sqrt(0.0316*0.0316 + (1.21/sqrt(et))*(1.21/sqrt(et)) + (5e-05/et)*(5e-05/et)));
  de = 1.0/sin(2 * atan(exp(-(1.653+1.740)/2))) * det;
  }

  if(1.740<=abs(eta) && abs(eta)<1.830){
  det = et * (sqrt(2.29e-07*2.29e-07 + (1.2/sqrt(et))*(1.2/sqrt(et)) + (1.71e-05/et)*(1.71e-05/et)));
  de = 1.0/sin(2 * atan(exp(-(1.740+1.830)/2))) * det;
  }

  if(1.830<=abs(eta) && abs(eta)<1.930){
  det = et * (sqrt(5.18e-09*5.18e-09 + (1.14/sqrt(et))*(1.14/sqrt(et)) + (1.7/et)*(1.7/et)));
  de = 1.0/sin(2 * atan(exp(-(1.830+1.930)/2))) * det;
  }

  if(1.930<=abs(eta) && abs(eta)<2.043){
  det = et * (sqrt(2.17e-07*2.17e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (2.08/et)*(2.08/et)));
  de = 1.0/sin(2 * atan(exp(-(1.930+2.043)/2))) * det;
  }

  if(2.043<=abs(eta) && abs(eta)<2.172){
  det = et * (sqrt(3.65e-07*3.65e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.63/et)*(1.63/et)));
  de = 1.0/sin(2 * atan(exp(-(2.043+2.172)/2))) * det;
  }

  if(2.172<=abs(eta) && abs(eta)<2.322){
  det = et * (sqrt(2.02e-07*2.02e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.68/et)*(1.68/et)));
  de = 1.0/sin(2 * atan(exp(-(2.172+2.322)/2))) * det;
  }

  if(2.322<=abs(eta) && abs(eta)<2.500){
  det = et * (sqrt(5.27e-07*5.27e-07 + (1.12/sqrt(et))*(1.12/sqrt(et)) + (1.78/et)*(1.78/et)));
  de = 1.0/sin(2 * atan(exp(-(2.322+2.500)/2))) * det;
  }

  return(de);
}
