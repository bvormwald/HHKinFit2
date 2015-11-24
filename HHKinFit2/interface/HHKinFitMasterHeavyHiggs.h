#ifndef HHKinFitMasterHeavyHiggs_H
#define HHKinFitMasterHeavyHiggs_H

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector2.h>

#ifdef HHKINFIT2
#include "HHLorentzVector.h"
#include "HHKinFit.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/HHKinFit.h"
#endif

#include <stdio.h>
#include <map>
#include <utility>
#include <vector>
#include <sstream>

namespace HHKinFit2{

typedef std::pair< int, int >                                 HHFitHypothesisHeavyHiggs;
typedef std::map< HHFitHypothesisHeavyHiggs, double >         HHFitResultD;
typedef std::map< HHFitHypothesisHeavyHiggs, int >            HHFitResultI;
typedef std::map< HHFitHypothesisHeavyHiggs, bool >           HHFitResultB;
typedef std::map< HHFitHypothesisHeavyHiggs, TLorentzVector > HHFitResultTLor;

class HHKinFitMasterHeavyHiggs{
  public:
 
  HHKinFitMasterHeavyHiggs(TLorentzVector const& bjet1,
                           TLorentzVector const& bjet2,
                           TLorentzVector const& tauvis1,
                           TLorentzVector const& tauvis2,
                           TVector2 const& met, TMatrixD const& met_cov,
                           double sigmaEbjet1 = -1.0, double sigmaEbjet2 = -1.0,
                           bool istruth=false,
                           TLorentzVector const& higgsgen=TLorentzVector(0,0,0,0));
  
  //the main action, runs over all hypotheses and performs the fit
  void fit();

  //sets the hypotheses
  void addHypo(HHFitHypothesisHeavyHiggs hypo);
  void addHypo(int mh1, int mh2);

  //Getters
  HHFitHypothesisHeavyHiggs getBestHypothesis();
  double getBestChi2();
  
  //Getters for fit results
  double getChi2(HHFitHypothesisHeavyHiggs hypo);
  double getChi2(int mh1 = 125, int mh2 = 125);

  double getChi2BJet1(HHFitHypothesisHeavyHiggs hypo);
  double getChi2BJet1(int mh1 = 125, int mh2 = 125);

  double getChi2BJet2(HHFitHypothesisHeavyHiggs hypo);
  double getChi2BJet2(int mh1 = 125, int mh2 = 125);
    
  double getChi2Balance(HHFitHypothesisHeavyHiggs hypo);
  double getChi2Balance(int mh1 = 125, int mh2 = 125);
  
  double getFitProb(HHFitHypothesisHeavyHiggs hypo);
  double getFitProb(int mh1 = 125, int mh2 = 125);
  
  double getMH(HHFitHypothesisHeavyHiggs hypo);
  double getMH(int mh1 = 125, int mh2 = 125);

  int getConvergence(HHFitHypothesisHeavyHiggs hypo);
  int getConvergence(int mh1 = 125, int mh2 = 125);

  TLorentzVector getFittedTau1(HHFitHypothesisHeavyHiggs hypo);
  TLorentzVector getFittedTau1(int mh1 = 125, int mh2 = 125);

  TLorentzVector getFittedTau2(HHFitHypothesisHeavyHiggs hypo);
  TLorentzVector getFittedTau2(int mh1 = 125, int mh2 = 125);

  TLorentzVector getFittedBJet1(HHFitHypothesisHeavyHiggs hypo);
  TLorentzVector getFittedBJet1(int mh1 = 125, int mh2 = 125);
  
  TLorentzVector getFittedBJet2(HHFitHypothesisHeavyHiggs hypo);
  TLorentzVector getFittedBJet2(int mh1 = 125, int mh2 = 125);
    
  //For Gen Studies to check smearing
  TLorentzVector getUnfittedBJet1(){return(m_bjet1);};
  TLorentzVector getUnfittedBJet2(){return(m_bjet2);};
  
  double getBJet1Resolution(){return m_sigma_bjet1;};
  double getBJet2Resolution(){return m_sigma_bjet2;};

  //DEPRECATED/////////////////////////////////////////////////////////
  //only here for backwards compatibility
  HHKinFitMasterHeavyHiggs(const TLorentzVector* bjet1,
                           const TLorentzVector* bjet2,
                           const TLorentzVector* tauvis1,
                           const TLorentzVector* tauvis2,
                           const TLorentzVector* met = 0, TMatrixD met_cov = TMatrixD(4,4), 
                           double sigmaEbjet1 = -1.0, double sigmaEbjet2 = -1.0,
                           bool istruth=false, TLorentzVector* heavyhiggsgen=0);
  void doFit();
  HHFitHypothesisHeavyHiggs getLowestChi2Hypothesis();
  void setAdvancedBalance(const TLorentzVector* met, TMatrixD met_cov);
  /////////////////////////////////////////////////////////////////////
  
  
  TLorentzVector neutrinos;
  TVector2 smearedMET;
  
  TLorentzVector initialHH;
  TLorentzVector finalHH;
  double b1Px_standardDevs;

private:
  double GetPFBJetRes(double eta, double et);

  //hypotheses
  std::vector< HHFitHypothesisHeavyHiggs > m_hypos;

  //input vectors
  HHLorentzVector m_bjet1;
  HHLorentzVector m_bjet2;
  HHLorentzVector m_tauvis1;
  HHLorentzVector m_tauvis2;

  TVector2 m_MET;
  TMatrixD m_MET_COV;

  double m_sigma_bjet1;
  double m_sigma_bjet2;

  //fit result
  HHFitResultI m_map_convergence;
  
  HHFitResultD m_map_chi2;
  HHFitResultD m_map_chi2BJet1;
  HHFitResultD m_map_chi2BJet2;
  HHFitResultD m_map_chi2Balance;
  HHFitResultD m_map_prob;
  HHFitResultD m_map_mH;

  HHFitResultTLor m_map_fittedTau1;
  HHFitResultTLor m_map_fittedTau2;
  HHFitResultTLor m_map_fittedB1;
  HHFitResultTLor m_map_fittedB2;

  HHFitHypothesisHeavyHiggs m_bestHypo;
  double m_bestChi2;

  //For ToyMC Studies
  TMatrixD m_bjet1_COV;
  TMatrixD m_bjet2_COV;
};
}
#endif
