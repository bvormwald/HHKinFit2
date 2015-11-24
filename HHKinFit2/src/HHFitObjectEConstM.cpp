#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyConstraintException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyConstraintException.h"
#endif

#include <cassert>
#include <iostream>
#include <sstream>
#include <math.h>

HHKinFit2::HHFitObjectEConstM::HHFitObjectEConstM(HHLorentzVector const& initial4vector)
  :HHFitObjectE(initial4vector){

}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectEConstM::constrainEtoMinv(double m, HHLorentzVector const& pset) const
{
  double E = calculateEConstrainedToMinv(m, pset);
  if (E < 0)
  {
    std::stringstream msg;
    msg << "problem in constraining energy to inv. mass in FitObjectEConstM: minv="<< m 
	<< " E(set)="<<E;
    throw(HHEnergyConstraintException(msg.str()));
  }
  return(this->changeE(E));
}

double HHKinFit2::HHFitObjectEConstM::calculateEConstrainedToMinv(double m, HHLorentzVector const& pset) const
{
  HHLorentzVector pmod = getInitial4Vector();
  HHLorentzVector combined = pset + pmod;

  double Mc = m;
  double M1c = pset.M();
  double M2c = pmod.M();
//  double M = combined.M();

  double C = 0.5*(pow(Mc,2)-pow(M1c,2)-pow(M2c,2));

//  int loopCount = 0;
//  while(fabs(M-Mc) > 0.000001){
//    loopCount++;
  double P1x = pset.Px();
  double P1y = pset.Py();
  double P1z = pset.Pz();
  double P1  = pset.P();
  double E1  = pset.E();

  double P2x = pmod.Px();
  double P2y = pmod.Py();
  double P2z = pmod.Pz();
  double P2  = pmod.P();

  double cosa = (P1x*P2x+P1y*P2y+P1z*P2z)/(P1*P2);
  double E2new = -1;

  double cp = 0.0;
  double dp = 0.0;
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;

  if(cosa==0){
    E2new=C/E1;
  }
  else{
    cp=C/(cosa*P1);
    dp=E1/(cosa*P1);
    a=pow(dp,2)-1;
    b=-2*dp*cp;
    c=pow(cp,2)+pow(M2c,2);

    if (cosa>0) E2new = (1./(2*a))*(-b+sqrt(pow(b,2)-4*a*c));
    if (cosa<0) E2new = (1./(2*a))*(-b-sqrt(pow(b,2)-4*a*c));
  }

  if (isnan(E2new)||isinf(E2new)){
    std::stringstream msg;
    msg << "problem in constraining energy to inv. mass in FitObjectEConstM: minv="<<Mc<< " E(set)="<<E2new << " cos(alpha)="<<cosa<<" P1="<<P1<<" P2="<< P2;
    std::cout << "E1:     " << E1 << std::endl;
    std::cout << "P1:     " << P1 << std::endl;
      
    std::cout << "P2:     " << P2 << std::endl;
    std::cout << "E2new:  " << E2new << std::endl;

    std::cout << "a:      " << a << std::endl;
    std::cout << "dp:     " << dp << std::endl;
    std::cout << "cosa:   " << cosa << std::endl;

     
    std::cout << "b:      " << b << std::endl;
    std::cout << "c:      " << c << std::endl;
    std::cout << "cp:     " << cp << std::endl;
    std::cout << "M2c:    " << M2c << std::endl;
      
    std::cout << "sqrt^2: " <<  pow(b,2)-4*a*c << std::endl;
    throw(HHEnergyConstraintException(msg.str()));
  }
  return E2new;
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectEConstM::changeE(double E) const{
  HHLorentzVector temp = this->getFit4Vector();
  temp.SetEkeepM(E);

  return(temp);
}

void
HHKinFit2::HHFitObjectEConstM::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "energy component fit object with constant mass:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printLimits();
//  this->printCovMatrix();
}
