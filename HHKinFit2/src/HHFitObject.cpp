#ifdef HHKINFIT2
#include "HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

#include <iostream>
#include <iomanip>

HHKinFit2::HHFitObject::HHFitObject()
  :m_fit4vector(HHLorentzVector()),
   m_initial4vector(HHLorentzVector()),
   m_covmatrix(TMatrixD(4,4)){

}

HHKinFit2::HHFitObject::HHFitObject(HHLorentzVector const initial4vector)
  :m_fit4vector(initial4vector),
   m_initial4vector(initial4vector),
   m_covmatrix(TMatrixD(4,4)){

}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObject::getInitial4Vector() const{
  return(m_initial4vector);
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObject::getFit4Vector() const{
  return(m_fit4vector);
}

TMatrixD
HHKinFit2::HHFitObject::getCovMatrix() const{
  return(m_covmatrix);
}

void
HHKinFit2::HHFitObject::setFit4Vector(HHLorentzVector const vec){
  this->m_fit4vector=vec;
}

void
HHKinFit2::HHFitObject::setCovMatrix(TMatrixD const covmatrix){
  this->m_covmatrix=covmatrix;
}

void
HHKinFit2::HHFitObject::reset(){
  this->m_fit4vector=this->m_initial4vector;
  m_covmatrix=TMatrixD(4,4);
}

void
HHKinFit2::HHFitObject::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "general fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printCovMatrix();
}

void
HHKinFit2::HHFitObject::printInitial4Vector() const{
  std::cout <<  "initial vector (px,py,pz,E,m,beta)"
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Px()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Py()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Pz()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().E()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().M()
            << std::setw(10) << std::fixed<< std::setprecision(4) << this->getInitial4Vector().Beta()
            << std::endl;
}

void
HHKinFit2::HHFitObject::printFit4Vector() const{
  std::cout <<  "  final vector (px,py,pz,E,m,beta)"
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Px()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Py()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Pz()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().E()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().M()
            << std::setw(10) << std::fixed<< std::setprecision(4) << this->getFit4Vector().Beta()
            << std::endl;
}

void
HHKinFit2::HHFitObject::printCovMatrix() const{
  std::cout <<  "covariance matrix:" << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,1)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,2)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,3) << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,1)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,2)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,3) << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(2,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(2,1)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(2,2)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(2,3) << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(3,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(3,1)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(3,2)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(3,3) << std::endl;
}
