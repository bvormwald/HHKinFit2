#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHKinFitMasterSingleHiggsSoftLimits.h"
#include "HHKinFit2/HHKinFit2/HHKinFit.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectEConstBeta.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectMET.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectComposite.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraintLikelihood.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraintEHardM.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectEConstBeta.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2/HHFitConstraintSoftBoundary.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHLimitSettingException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterSingleHiggsSoftLimits.h"
#include "HHKinFit2/HHKinFit2/interface/HHKinFit.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstBeta.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectMET.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectComposite.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraintLikelihood.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraintEHardM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstBeta.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraintSoftBoundary.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHLimitSettingException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"
#endif

#include "TMatrixD.h"
#include "TRandom3.h"

#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <iterator>

void HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::fit()
{
  for(unsigned int i = 0; i < m_hypos.size(); ++i)
    { 

      HHFitObjectE* tau1Fit = new HHFitObjectEConstM(m_tauvis1);
      HHFitObjectE* tau2Fit = new HHFitObjectEConstM(m_tauvis2);

      //prepare MET object
      HHFitObjectMET* metFit = new HHFitObjectMET(m_MET);
      metFit->setCovMatrix(m_MET_COV);
  
      //prepare composite object: Higgs
      HHFitObject* higgs  = new HHFitObjectComposite(tau1Fit, tau2Fit, metFit);
  
      int mh = m_hypos[i];

      std::cout << mh << std::endl;
  
      try{
        tau1Fit->setFitLimitsE(tau1Fit->getInitial4Vector(),mh,tau2Fit->getInitial4Vector());
        tau2Fit->setFitLimitsE(tau2Fit->getInitial4Vector(),mh,tau1Fit->getInitial4Vector());
      }
      catch(HHLimitSettingException const& e){
        std::cout << "Exception while setting tau limits:" << std::endl;
        std::cout << e.what() << std::endl;
        std::cout << "Tau energies are not compatible with invariant mass constraint." << std::endl;

        m_map_chi2[m_hypos[i]] = -pow(10,10);
        m_map_prob[m_hypos[i]] = -pow(10,10);
        m_bestHypo =HHFitHypothesisSingleHiggs (-pow(10,10));
        m_chi2_best = -pow(10,10);
        continue;
      }
  


      //prepare constraints
      HHFitConstraint* c_invmh = new HHFitConstraintEHardM(tau1Fit, tau2Fit, mh);
      HHFitConstraint* c_balance = new HHFitConstraint4Vector(higgs, true, true, false, false);
      HHFitConstraint* c_softlimit1 = new HHFitConstraintSoftBoundary(tau1Fit,0.1);
      HHFitConstraint* c_softlimit2 = new HHFitConstraintSoftBoundary(tau2Fit,0.1);

      //fit
      HHKinFit2::HHKinFit* fitObject = new HHKinFit2::HHKinFit();

      tau1Fit->setInitStart((tau1Fit->getUpperFitLimitE()+tau1Fit->getLowerFitLimitE())/2);
      tau1Fit->setInitPrecision(0.1);
      tau1Fit->setInitStepWidth(0.1*(tau1Fit->getUpperFitLimitE() - tau1Fit->getLowerFitLimitE()));
      tau1Fit->setInitDirection(1.0);

      fitObject->addFitObjectE(tau1Fit);

      fitObject->addConstraint(c_invmh);
      fitObject->addConstraint(c_balance);
      fitObject->addConstraint(c_softlimit1);
      fitObject->addConstraint(c_softlimit2);

      tau1Fit->resetLimits();
      tau2Fit->resetLimits();
      
      fitObject->fit();

      m_map_convergence[m_hypos[i]] = fitObject->getConvergence();    
    


      double chi2 = fitObject->getChi2();
      m_map_chi2[m_hypos[i]] = chi2;
      m_map_prob[m_hypos[i]] = TMath::Prob(chi2, 1);

      if(chi2 < m_chi2_best)
        {
          m_bestHypo = m_hypos[i];
          m_chi2_best = chi2;
        }
    
      TLorentzVector fittedTau1 =  ( (TLorentzVector)tau1Fit->getFit4Vector()  );
      m_map_fittedTau1[m_hypos[i]] = fittedTau1;
      TLorentzVector fittedTau2 =  ( (TLorentzVector)tau2Fit->getFit4Vector()  );
      m_map_fittedTau2[m_hypos[i]] = fittedTau2;


      delete c_invmh;
      delete c_balance;
      delete fitObject;

      delete tau1Fit;
      delete tau2Fit;
      delete metFit;
      delete higgs;
    }
}

void HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::addHypo(HHFitHypothesisSingleHiggs hypo)
{
  m_hypos.push_back(hypo);
}


HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::HHKinFitMasterSingleHiggsSoftLimits(
                                                                TLorentzVector const& tauvis1,
                                                                TLorentzVector const& tauvis2,
                                                                TVector2 const& met, 
                                                                TMatrixD const& met_cov, 
                                                                bool istruth,
                                                                TLorentzVector const& higgsgen)
:m_MET_COV(TMatrixD(4,4))
{
  
  m_tauvis1 = HHLorentzVector(tauvis1.Px(), tauvis1.Py(), tauvis1.Pz(), tauvis1.E());  
  m_tauvis2 = HHLorentzVector(tauvis2.Px(), tauvis2.Py(), tauvis2.Pz(), tauvis2.E());
 
  m_tauvis1.SetMkeepE(1.77682);
  m_tauvis2.SetMkeepE(1.77682);
   
  m_MET = met;
  m_MET_COV = met_cov;

  m_chi2_best = pow(10,10);
  m_bestHypo = 0;

  //  if (istruth){
  //    TRandom3 r(0);
  //
  //    HHLorentzVector recoil;
  //    if(heavyhiggsgen != NULL){
  //       Double_t pxRecoil = r.Gaus(-(heavyhiggsgen->Px() ), 10.0);
  //       Double_t pyRecoil = r.Gaus(-(heavyhiggsgen->Py() ), 10.0);
  //
  //       recoil = HHLorentzVector(pxRecoil, pyRecoil, 0,
  //				sqrt(pxRecoil*pxRecoil+pyRecoil*pyRecoil));
  //    }
  //    else{
  //      recoil = HHLorentzVector(0,0,0,0);
  //      std::cout << "WARNING! Truthinput mode active but no Heavy Higgs gen-information given! Setting Recoil to Zero!" << std::endl;
  //    }
  //
  //    TMatrixD recoilCov(2,2);
  //    recoilCov(0,0)=100;  recoilCov(0,1)=0;
  //    recoilCov(1,0)=0;    recoilCov(1,1)=100;
  //
  //    HHLorentzVector recoHH = m_bjet1 + m_bjet2 + m_tauvis1 + m_tauvis2 + recoil;
  //    m_MET = TVector2(-recoHH.Px(), -recoHH.Py() );
  //
  //    m_MET_COV = TMatrixD(2,2);
  //    m_MET_COV = recoilCov + bjet1Cov + bjet2Cov;
  //
  //  }
}

HHKinFit2::HHFitHypothesisSingleHiggs HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::getBestHypothesis(){
  return(m_bestHypo);
}

double HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::getBestChi2(){
  return(m_chi2_best);
}

double HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::getChi2(HHFitHypothesisSingleHiggs hypo){
  return(m_map_chi2[hypo]);
}

double HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::getFitProb(HHFitHypothesisSingleHiggs hypo){
  return(m_map_prob[hypo]);
}

int HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::getConvergence(HHFitHypothesisSingleHiggs hypo){
  return(m_map_convergence[hypo]);
}

TLorentzVector HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::getFittedTau1(HHFitHypothesisSingleHiggs hypo){
  return(m_map_fittedTau1[hypo]);
}

TLorentzVector HHKinFit2::HHKinFitMasterSingleHiggsSoftLimits::getFittedTau2(HHFitHypothesisSingleHiggs hypo){
  return(m_map_fittedTau2[hypo]);
}
