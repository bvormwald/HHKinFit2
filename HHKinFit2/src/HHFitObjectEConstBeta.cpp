#ifdef HHKINFIT2
#include "HHFitObjectEConstBeta.h"
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHEnergyConstraintException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstBeta.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyConstraintException.h"
#endif

#include <cassert>
#include <iostream>
#include <sstream>
#include <math.h>

HHKinFit2::HHFitObjectEConstBeta::HHFitObjectEConstBeta(HHLorentzVector const& initial4vector)
  :HHFitObjectE(initial4vector){

}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectEConstBeta::constrainEtoMinv(double m, HHLorentzVector const& pset) const{
  
  double E = calculateEConstrainedToMinv(m, pset);
  if ( E < 0)
  {
    std::stringstream msg;
    msg << "problem in constraining energy to inv. mass in FitObjectEConstBeta: minv="
	<< m << " E(set)=" << E;
    throw(HHEnergyConstraintException(msg.str()));
  }
  return(this->changeE(E));  
}

double HHKinFit2::HHFitObjectEConstBeta::calculateEConstrainedToMinv(double m, HHLorentzVector const& pset) const
{
  HHLorentzVector pmod = getInitial4Vector();

  //energy and momenta of involved particles
  double P1x = pset.Px();
  double P1y = pset.Py();
  double P1z = pset.Pz();
  double P1  = pset.P();
  double E1  = pset.E();

  double P2x = pmod.Px();
  double P2y = pmod.Py();
  double P2z = pmod.Pz();
  double P2  = pmod.P();

  //betas and angle
  double b1c = pset.Beta();
  double b2c = pmod.Beta();
  double cosa = (P1x*P2x+P1y*P2y+P1z*P2z)/(P1*P2);

  
  //calculatae energy of this particle
  double a=1-pow(b2c,2);
  double b=2*E1*(1-b1c*b2c*cosa);
  double c=(1-pow(b1c,2))*pow(E1,2)-pow(m,2);
  double E2new = (1./(2*a))*(-b+sqrt(pow(b,2)-4*a*c));

  if (isnan(E2new)||isinf(E2new)){
    std::stringstream msg;
    msg << "problem in constraining energy to inv. mass in FitObjectEConstBeta: minv="<<m<< " E(set)="<<E2new;
    std::cout << "E1:     " << E1 << std::endl;
    std::cout << "P1:     " << P1 << std::endl;
    
    std::cout << "P2:     " << P2 << std::endl;
    std::cout << "E2new:  " << E2new << std::endl;
    
    std::cout << "a:      " << a << std::endl;
    std::cout << "b2c:    " << b2c << std::endl;
    
    std::cout << "b:      " << b << std::endl;
    std::cout << "c:      " << c << std::endl;
    std::cout << "cosa:   " << cosa << std::endl;
    std::cout << "b1c:    " << b1c << std::endl;
    
    std::cout << "sqrt^2: " <<  pow(b,2)-4*a*c << std::endl;
    throw(HHEnergyConstraintException(msg.str()));
  }
  return E2new;
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectEConstBeta::changeE(double E) const{
  HHLorentzVector temp = this->getFit4Vector();
  temp.SetEkeepBeta(E);

  return(temp);
}


void
HHKinFit2::HHFitObjectEConstBeta::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "energy component fit object with constant beta (=p/E):" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printCovMatrix();
  this->printLimits();
}
