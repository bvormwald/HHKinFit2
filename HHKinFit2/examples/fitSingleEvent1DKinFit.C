#include <iostream>
#include <iomanip>

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TMatrixD.h"
#include "HHKinFitMasterSingleHiggs.h"

#include "TGraph.h" 
#include "TFile.h" 

using HHKinFit2::HHKinFitMasterSingleHiggs;

int main(int argc, char* argv[])
{
  TLorentzVector tauvis1(-15.902388,37.502301,43.366285,59.503192);
  TLorentzVector tauvis2(16.824116,-21.066462,7.050932,27.867069);
  TVector2 met(8.283318,11.756746);
  TMatrixD met_cov(4,4);
  met_cov[0][0]= 359.15655;
  met_cov[0][1]= 0;
  met_cov[1][0]= met_cov[0][1];
  met_cov[1][1]= 359.15655;
  
  HHKinFitMasterSingleHiggs singlehiggsfit(tauvis1,tauvis2,met,met_cov);
  
  singlehiggsfit.addHypo(125);
  singlehiggsfit.fit();


  std::cout << "Single Higgs fit finished." << std::endl;

  std::cout << "chi2: " << singlehiggsfit.getChi2(125) << std::endl;
  std::cout << "Fit prob: " << singlehiggsfit.getFitProb(125) << std::endl;
  std::cout << "Convergence: " << singlehiggsfit.getConvergence(125) << std::endl;

  TLorentzVector fittedTau1 = singlehiggsfit.getFittedTau1(125);
  TLorentzVector fittedTau2 = singlehiggsfit.getFittedTau2(125);

  std::cout << "Energy of fitted tau1: " << fittedTau1.E() << std::endl;
  std::cout << "Energy of fitted tau2: " << fittedTau2.E() << std::endl;
  
  return (0);
}
