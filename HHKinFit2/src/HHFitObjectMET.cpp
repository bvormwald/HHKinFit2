#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHFitObjectMET.h"
#include "HHKinFit2/HHKinFit2/HHLorentzVector.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectMET.h"
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#endif

#include <iostream>
#include <iomanip>

HHKinFit2::HHFitObjectMET::HHFitObjectMET(TVector2 const& v,TVector2 const& vfit)
  :HHFitObject(HHLorentzVector(v.X(),v.Y(),0,sqrt(v.X()*v.X()+v.Y()*v.Y()))){
  this->setFit4Vector(HHLorentzVector(vfit.X(),vfit.Y(),0,sqrt(vfit.X()*vfit.X()+vfit.Y()*vfit.Y())));
}


void
HHKinFit2::HHFitObjectMET::setCovMatrix(double xx, double yy, double xy){
  TMatrixD cov(4,4);
  cov(0,0)=xx;
  cov(1,1)=yy;
  cov(0,1)=xy;
  cov(1,0)=xy;
  m_covmatrix=cov;
}

void
HHKinFit2::HHFitObjectMET::setCovMatrix(TMatrixD const covmat){
  m_covmatrix=covmat;
}


void
HHKinFit2::HHFitObjectMET::print() const{
  std::cout << "---" << std::endl;
  std::cout << "MET fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printCovMatrix();
}


void
HHKinFit2::HHFitObjectMET::printInitial4Vector() const{
  std::cout <<  "initial vector (px,py,pt)    "
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Px()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Py()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Pt()
            << std::endl;
}

void
HHKinFit2::HHFitObjectMET::printFit4Vector() const{
  std::cout <<  "  final vector (px,py,pt)    "
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Px()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Py()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Pt()
            << std::endl;
}

void
HHKinFit2::HHFitObjectMET::printCovMatrix() const{
  std::cout <<  "covariance matrix:" << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,1) << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,1) << std::endl;
}
