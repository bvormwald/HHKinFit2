#ifndef HHKinFitMasterSingleHiggs_H
#define HHKinFitMasterSingleHiggs_H

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector2.h>

#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/HHKinFit.h"
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

typedef int                                 HHFitHypothesisSingleHiggs;
typedef std::map< HHFitHypothesisSingleHiggs, double >         HHFitResultD;
typedef std::map< HHFitHypothesisSingleHiggs, int >            HHFitResultI;
typedef std::map< HHFitHypothesisSingleHiggs, bool >           HHFitResultB;
typedef std::map< HHFitHypothesisSingleHiggs, TLorentzVector > HHFitResultTLor;

class HHKinFitMasterSingleHiggs{
  public:
  HHKinFitMasterSingleHiggs(TLorentzVector const& tauvis1,
                            TLorentzVector const& tauvis2,
                            TVector2 const& met, TMatrixD const& met_cov,
                            bool istruth=false,
                            TLorentzVector const& higgsgen=TLorentzVector(0,0,0,0));

  //the main action, runs over all hypotheses and performs the fit
  void fit();

  void addHypo(HHFitHypothesisSingleHiggs hypo);

  //Getters
  HHFitHypothesisSingleHiggs getBestHypothesis();
  double getBestChi2();

  //Getters for fit results
  double getChi2(HHFitHypothesisSingleHiggs hypo);
 
  double getFitProb(HHFitHypothesisSingleHiggs hypo);


  int getConvergence(HHFitHypothesisSingleHiggs hypo);

  TLorentzVector getFittedTau1(HHFitHypothesisSingleHiggs hypo);

  TLorentzVector getFittedTau2(HHFitHypothesisSingleHiggs hypo);
  
private:  
  //hypotheses
  std::vector< HHFitHypothesisSingleHiggs > m_hypos;

  //input vectors
  HHLorentzVector m_tauvis1;
  HHLorentzVector m_tauvis2;

  TVector2 m_MET;
  TMatrixD m_MET_COV;

  //fit result
  HHFitResultI m_map_convergence;
  
  HHFitResultD m_map_chi2;
  HHFitResultD m_map_prob;

  HHFitResultTLor m_map_fittedTau1;
  HHFitResultTLor m_map_fittedTau2;

  HHFitHypothesisSingleHiggs m_bestHypo;
  double m_chi2_best;
};
}
#endif
