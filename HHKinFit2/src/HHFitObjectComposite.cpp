#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectComposite.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectComposite.h"
#endif

#include <iostream>

HHKinFit2::HHFitObjectComposite::HHFitObjectComposite(std::vector<HHFitObject*> const& subobjects)
  : m_subobjects(subobjects),
    m_cov_set(false){
}

HHKinFit2::HHFitObjectComposite::HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2)
  : m_subobjects(),
    m_cov_set(false){
  this->addSubobject(subobject1);
  this->addSubobject(subobject2);
}

HHKinFit2::HHFitObjectComposite::HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2, HHFitObject* subobject3)
  : m_subobjects(),
    m_cov_set(false){
  this->addSubobject(subobject1);
  this->addSubobject(subobject2);
  this->addSubobject(subobject3);
}

HHKinFit2::HHFitObjectComposite::HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2, HHFitObject* subobject3, HHFitObject* subobject4, HHFitObject* subobject5)
  : m_subobjects(),
    m_cov_set(false){
  this->addSubobject(subobject1);
  this->addSubobject(subobject2);
  this->addSubobject(subobject3);
  this->addSubobject(subobject4);
  this->addSubobject(subobject5);
}


HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectComposite::getFit4Vector() const{
  HHLorentzVector p(0,0,0,0);
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    p += (*it)->getFit4Vector();
//  m_fit4vector = p;
  return(p);
}


HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectComposite::getInitial4Vector() const {
  HHLorentzVector p(0,0,0,0);
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it)
    p += (*it)->getInitial4Vector();
//  m_initial4vector = p;
  return(p);
}

TMatrixD
HHKinFit2::HHFitObjectComposite::getCovMatrix() const{
  if(m_cov_set){
    return(m_covmatrix);
  }
  else{
    TMatrixD cov(4,4);
    //int i=m_subobjects.size();
    for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it){
      //i--;
      //if(i==0){
      cov += (*it)->getCovMatrix();
      //      (*it)->getCovMatrix().Print();
      //}
      //else{
      //cov -= (*it)->getCovMatrix();
      //cov += (*it)->getCovMatrix();
      //(*it)->getCovMatrix().Print();
      //  m_covmatrix=cov;
      // }
    }
    //    cov.Print();
    //    std::cout << "--------" << std::endl;
    return(cov);
  }
}


void
HHKinFit2::HHFitObjectComposite::setCovMatrix(TMatrixD const cov){ 
  m_covmatrix=cov;
  m_cov_set=true;
}

void
HHKinFit2::HHFitObjectComposite::setSubobjects(std::vector<HHFitObject*> const& subobjects){
  m_subobjects=subobjects;
}


void
HHKinFit2::HHFitObjectComposite::addSubobject(HHFitObject* subobject){
  m_subobjects.push_back(subobject);
}

void
HHKinFit2::HHFitObjectComposite::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "composite fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  std::cout <<  "-----------------" << std::endl;
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it){
    (*it)->printInitial4Vector();
  }
  std::cout <<  "-----------------" << std::endl;
  for (std::vector<HHFitObject*>::const_iterator it = m_subobjects.begin() ; it != m_subobjects.end(); ++it){
    (*it)->printFit4Vector();
  }
  std::cout <<  "-----------------" << std::endl;
  this->printCovMatrix();
}
