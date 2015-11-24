#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#endif

#include <cmath>
#include <sstream>

HHKinFit2::HHLorentzVector::HHLorentzVector(double x, double y, double z, double t)
  : TLorentzVector(x,y,z,t){
}

HHKinFit2::HHLorentzVector::HHLorentzVector()
  : TLorentzVector(){
}

void
HHKinFit2::HHLorentzVector::SetEEtaPhiM(double E, double eta, double phi, double m){
  double p = sqrt(E*E-m*m);
  double sinth = sin(2*atan(exp(-eta)));
  double pt = sinth*p;
  this->SetPtEtaPhiE(pt,eta,phi,E);
}


HHKinFit2::HHLorentzVector
HHKinFit2::HHLorentzVector::operator+(HHLorentzVector const& rhs) const{
  return(HHLorentzVector(rhs.Px()+this->Px(),rhs.Py()+this->Py(),rhs.Pz()+this->Pz(),rhs.E()+this->E()));
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHLorentzVector::operator-(HHLorentzVector const& rhs) const{
  return(HHLorentzVector(this->Px()-rhs.Px(),this->Py()-rhs.Py(),this->Pz()-rhs.Pz(),this->E()-rhs.E()));
}

HHKinFit2::HHLorentzVector
HHKinFit2::HHLorentzVector::operator-() const
{
  return(HHLorentzVector(-this->Px(), -this->Py(), -this->Pz(), -this->E() ) );
}

HHKinFit2::HHLorentzVector&
HHKinFit2::HHLorentzVector::operator=(const HHLorentzVector& other)
{
  SetPxPyPzE(other.Px(), other.Py(), other.Pz(), other.E() );
  return *this;
}



void
HHKinFit2::HHLorentzVector::SetEkeepM(double E){
  if(E<M()){
    std::stringstream msg;
    msg << "target energy is smaller than the particle mass: "<<"E(set)="<<E<<" "<<"m="<<M();
    throw(HHEnergyRangeException(msg.str()));
  }
  
  double pnew = sqrt(pow(E,2)-pow(M(),2));
  if (isnan(pnew)) {
      std::cout << "WARNING: SetEkeepM(): Targeted E is smaller than m. Set P=1, E=sqrt(m**2+1**2)" << std::endl;
      std::cout << "E: " << E << std::endl;
      std::cout << "M: " << M() << std::endl;
      E=sqrt(pow(M(),2)+1.0);
      pnew = 1.0;
      std::cout << "New E: " << E << std::endl;
      std::cout << "ptnew: " << pnew * sin(2.*atan(exp(-Eta()))) << std::endl;
  }
  double ptnew = pnew * sin(2.*atan(exp(-Eta())));
  
  SetPtEtaPhiE(ptnew,Eta(),Phi(),E);

}


void
HHKinFit2::HHLorentzVector::SetMkeepE(double m){
  double energy=E();
  if(energy<m){
    std::cout << "SetMkeepE()::energy is smaller than the particle mass: "<<"m(set)="<<m<<" "<<"E="<<energy << std::endl;
    energy = 1.1*m;
  }

  double pnew = sqrt(pow(energy,2)-pow(m,2));
  double ptnew = pnew * sin(2.*atan(exp(-Eta())));
  SetPtEtaPhiE(ptnew,Eta(),Phi(),energy);
}


void
HHKinFit2::HHLorentzVector::SetEkeepBeta(double E){
  double pnew = Beta()*E;
  double ptnew = pnew * sin(2.*atan(exp(-Eta())));
  SetPtEtaPhiE(ptnew,Eta(),Phi(),E);
}
