#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyConstraintException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHLimitSettingException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyConstraintException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHLimitSettingException.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>


HHKinFit2::HHFitObjectE::HHFitObjectE(HHLorentzVector const& initial4vector)
  :HHFitObject(initial4vector),
   m_upperLimitE(pow(10,10)),
   m_lowerLimitE(0),
   m_initstart(-pow(10,10)),
   m_initprec(-pow(10,10)),
   m_initstep(-pow(10,10)),
   m_initdirection(-pow(10,10)){
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectE::scaleE(double scale) const{
  return(this->changeE(scale*this->getFit4Vector().E()));
}


void
HHKinFit2::HHFitObjectE::changeEandSave(double E, bool respectLimits){
  if (respectLimits){
    if((E<this->getLowerFitLimitE())||(E>this->getUpperFitLimitE())){
      std::stringstream msg;
      msg << "target energy is out of limits: "<<"E(set)="<<E<<" "<<"E(limits)=["<<this->getLowerFitLimitE()<<","<< this->getUpperFitLimitE() << "]";
      throw(HHEnergyRangeException(msg.str()));
    }
  }
  this->setFit4Vector(changeE(E));
}

void
HHKinFit2::HHFitObjectE::scaleEandSave(double scale, bool respectLimits){
  this->changeEandSave(scale*this->getFit4Vector().E(), respectLimits);
}

void
HHKinFit2::HHFitObjectE::constrainEtoMinvandSave(double m, HHLorentzVector const& other4vector, bool respectLimits){
  //  this->setFit4Vector(constrainEtoMinv(m, other4vector));
  double E = constrainEtoMinv(m, other4vector).E();
  this->changeEandSave(E, respectLimits);
}

double
HHKinFit2::HHFitObjectE::getE() const{
  return(m_fit4vector.E());
}


void
HHKinFit2::HHFitObjectE::setFitLimitsE(double const lowerlimit, double const upperlimit){
  this->setLowerFitLimitE(lowerlimit);
  this->setUpperFitLimitE(upperlimit);
}

void
HHKinFit2::HHFitObjectE::setFitLimitsE(HHLorentzVector const& own4vectorMin, double const minv, HHLorentzVector const& other4vectorMin){
  this->setLowerFitLimitE(own4vectorMin);
  this->setUpperFitLimitE(minv,other4vectorMin);
}


double
HHKinFit2::HHFitObjectE::getUpperFitLimitE() const{
  return(m_upperLimitE);
}

double
HHKinFit2::HHFitObjectE::getLowerFitLimitE() const{
  return(m_lowerLimitE);
}

void
HHKinFit2::HHFitObjectE::setUpperFitLimitE(double upperlimit){
  if (upperlimit<this->getLowerFitLimitE()){
    std::stringstream msg;
    msg << "Cannot set upper limit: E(upper)="<<upperlimit<<" < E(lower)="<<this->getLowerFitLimitE()<<"";
    throw(HHLimitSettingException(msg.str()));
  }

  m_upperLimitE = upperlimit;
}

void
HHKinFit2::HHFitObjectE::setUpperFitLimitE(double minv, HHLorentzVector const& other4vectorMin){
  try
  {
    this->setUpperFitLimitE(constrainEtoMinv(minv,other4vectorMin).E());
  }
  catch(HHEnergyConstraintException const& e)
  {
    std::cout << e.what() << std::endl;
    std::stringstream msg;
    msg << "Cannot set upper limit. Energy of second particle would be invalid/negative.";
    throw(HHLimitSettingException(msg.str()));
  }
}

void
HHKinFit2::HHFitObjectE::setLowerFitLimitE(double lowerlimit){
  m_lowerLimitE = lowerlimit;
}

void
HHKinFit2::HHFitObjectE::setLowerFitLimitE(double minv, HHLorentzVector const& other4vectorMin){
  try
  {
    double lowerFitLimit = calculateEConstrainedToMinv(minv,other4vectorMin);
    if(lowerFitLimit > m_lowerLimitE)
    {
      this->setLowerFitLimitE(lowerFitLimit);
    }
  }
  catch(HHEnergyConstraintException const& e)
  {
    std::cout << e.what() << std::endl;
    std::stringstream msg;
    msg << "Cannot set upper limit. Energy of second particle would be invalid/negative.";
    throw(HHLimitSettingException(msg.str()));
  }
}

void
HHKinFit2::HHFitObjectE::setLowerFitLimitE(HHLorentzVector const& lowerlimit){
  m_lowerLimitE = lowerlimit.E();
}

void
HHKinFit2::HHFitObjectE::setCovMatrix(double dE){
  TMatrixD cov(4,4);
  
  //NOTE: only dE is used for Cov calculation!
  double p=m_initial4vector.P();
  double e=m_initial4vector.E();
  double phi=m_initial4vector.Phi();
  double theta=m_initial4vector.Theta();
  double dp = e/p*dE; // error propagation p=sqrt(e^2-m^2)
  double dpt = sin(theta)*dp;

  cov(0,0) = pow(cos(phi)*dpt,2);
  cov(1,1) = pow(sin(phi)*dpt,2);
  cov(0,1) = sin(phi)*cos(phi)*dpt*dpt;
  cov(1,0) = sin(phi)*cos(phi)*dpt*dpt;

  cov(3,3)=dE*dE;
  m_covmatrix=cov;
}

void
HHKinFit2::HHFitObjectE::setCovMatrix(TMatrixD cov){
  m_covmatrix=cov;
}

void
HHKinFit2::HHFitObjectE::setInitPrecision(double prec){
  m_initprec=prec;
}

void
HHKinFit2::HHFitObjectE::setInitDirection(double daN){
  m_initdirection=daN;
}

void
HHKinFit2::HHFitObjectE::setInitStepWidth(double h){
  m_initstep=h;
}

void
HHKinFit2::HHFitObjectE::setInitStart(double start){
  m_initstart=start;
}

double
HHKinFit2::HHFitObjectE::getInitPrecision(){
  return(m_initprec);
}

double
HHKinFit2::HHFitObjectE::getInitDirection(){
  return(m_initdirection);
}

double
HHKinFit2::HHFitObjectE::getInitStepWidth(){
  return(m_initstep);
}

double
HHKinFit2::HHFitObjectE::getInitStart(){
  return(m_initstart);
}

void
HHKinFit2::HHFitObjectE::reset(){
  HHFitObject::reset();
  resetLimits();
  m_initstart = -pow(10,10);
  m_initprec = -pow(10,10);
  m_initstep = -pow(10,10);
  m_initdirection = -pow(10,10);
}

void
HHKinFit2::HHFitObjectE::resetLimits(){
  m_upperLimitE = pow(10,10);
  m_lowerLimitE = 0;
}


void
HHKinFit2::HHFitObjectE::print() const{
  std::cout << "---" << std::endl;
  std::cout << "energy component fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printLimits();
  this->printCovMatrix();
}


void
HHKinFit2::HHFitObjectE::printLimits() const{
  std::cout <<  "limits: "
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getLowerFitLimitE()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getUpperFitLimitE()
            << std::endl;
}

