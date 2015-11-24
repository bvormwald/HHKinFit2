#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"
#endif

#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include <iostream>
#include <sstream>
#include <cmath>

HHKinFit2::HHFitConstraint4Vector::HHFitConstraint4Vector(HHFitObject* object, bool px, bool py, bool pz, bool E)
  : HHFitConstraint(object),
    m_ncomp(0),
    m_cov(),
    m_invcov(),
    m_indices()
{
  m_components[0]=px;
  m_components[1]=py;
  m_components[2]=pz;
  m_components[3]=E;
  for (int i=0; i<4; i++) {
    if (m_components[i]) {
      m_indices.push_back(i);
      m_ncomp++;
    }
  }
  m_cov.ResizeTo(m_ncomp,m_ncomp);

  for (int i=0; i<m_ncomp; i++)
    for (int j=0; j<m_ncomp; j++)
      m_cov(i,j)=m_fitobject->getCovMatrix()(m_indices[i],m_indices[j]);


  if(m_ncomp>1){
    TMatrixDEigen eigenmatrix(m_cov);
    for (int i=0; i<m_ncomp; i++){
      if (eigenmatrix.GetEigenValues()(i,i)<0){
        std::stringstream msg;
        msg << "covariance matrix is not positive-definite, but has at least one negative eigenvalue";
        throw(HHCovarianceMatrixException(msg.str()));
      }
    }
  }

  m_invcov.ResizeTo(m_ncomp,m_ncomp);
  m_invcov = m_cov;
  m_invcov.Invert();

//  m_cov.Print();
//  m_invcov.Print();
}

double
HHKinFit2::HHFitConstraint4Vector::getChi2() const{
  HHLorentzVector res = m_fitobject->getFit4Vector()-m_fitobject->getInitial4Vector();
  double chi2sum=0;
  for(int i=0; i<m_ncomp; i++){
    for(int j=0;j<m_ncomp;j++){
      chi2sum+=res(m_indices[i])*res(m_indices[j])*m_invcov(i,j);
      //std::cout << res(i) << "*" <<res(j) << "*" << invcov(i,j) << "=" << res(i)*res(j)*invcov(i,j) << std::endl;
    }
  }
  return(chi2sum);
}

void
HHKinFit2::HHFitConstraint4Vector::printChi2() const{
  HHLorentzVector res = m_fitobject->getFit4Vector()-m_fitobject->getInitial4Vector();
  m_fitobject->getFit4Vector().Print();
  m_fitobject->getInitial4Vector().Print();
  res.Print();
  double chi2sum=0;
  for(int i=0; i<m_ncomp; i++){
    for(int j=0;j<m_ncomp;j++){
      chi2sum+=res(m_indices[i])*res(m_indices[j])*m_invcov(i,j);
      std::cout << res(m_indices[i]) << "*" <<res(m_indices[j]) << "*" << m_invcov(i,j) << "=" << res(m_indices[i])*res(m_indices[j])*m_invcov(i,j) << std::endl;
    }
  }
}

double
HHKinFit2::HHFitConstraint4Vector::getLikelihood() const{
    double likelihood=(1.0/pow(2*TMath::Pi(),m_ncomp*1.0/2.0))*(1.0/sqrt(m_cov.Determinant()))*exp(-0.5*this->getChi2());
    return(likelihood);
}
