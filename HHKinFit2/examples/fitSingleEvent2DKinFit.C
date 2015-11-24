#include <iostream>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector2.h"
#include "HHKinFitMasterHeavyHiggs.h"


using HHKinFit2::HHKinFitMasterHeavyHiggs;

int main(int argc, char* argv[])
{
  TLorentzVector bjet1(-18.0706,28.8731,19.1597,39.3622);
  TLorentzVector bjet2(14.1252,71.0363,-160.683,176.315);
  TLorentzVector tauvis1(-20.9251,-24.865,-4.66771,32.8318);
  TLorentzVector tauvis2(53.6093,-24.6612,-52.1178,78.7299);
  TVector2 met(31.2118,-38.0866);
  TMatrixD met_cov(2,2);
  met_cov[0][0]= 200;
  met_cov[0][1]= 0;
  met_cov[1][0]= 0;
  met_cov[1][1]= 200;

  HHKinFitMasterHeavyHiggs heavyhiggsfit(bjet1, bjet2, tauvis1, tauvis2, met, met_cov);

  /* 
  //You can set your own bjet resolutions by using optional arguments:

  double bjet1Res = 10.0;
  double bjet2Res = 10.0;
  HHKinFitMasterHeavyHiggs heavyhiggsfit(bjet1, bjet2, tauvis1, tauvis2, met, met_cov, bjet1Res, bjet2Res);
  */

  //If no mass hypothesis is given, the MSSM case with a 125GeV light higgs is assumed.
  //heavyhiggsfit.addHypo(125,110); 

  heavyhiggsfit.fit();

  std::cout << "Heavy Higgs fit finished." << std::endl;

  std::cout << "m_H: " << heavyhiggsfit.getMH() << "GeV" << std::endl;
  std::cout << "chi2: " << heavyhiggsfit.getChi2() << std::endl;
  std::cout << "Fit prob: " << heavyhiggsfit.getFitProb() << std::endl;
  std::cout << "Convergence: " << heavyhiggsfit.getConvergence() << std::endl;
  
  TLorentzVector fittedTau1 = heavyhiggsfit.getFittedTau1();
  TLorentzVector fittedTau2 = heavyhiggsfit.getFittedTau2();
  TLorentzVector fittedBJet1 = heavyhiggsfit.getFittedBJet1();
  TLorentzVector fittedBJet2 = heavyhiggsfit.getFittedBJet2();

  std::cout << "Energy of fitted tau1: " << fittedTau1.E() << std::endl;
  std::cout << "Energy of fitted tau2: " << fittedTau2.E() << std::endl;
  std::cout << "Energy of fitted BJet1: " << fittedBJet1.E() << std::endl;
  std::cout << "Energy of fitted BJet2: " << fittedBJet2.E() << std::endl;

  return (0);
}
