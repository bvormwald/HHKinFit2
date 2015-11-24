#ifdef HHKINFIT2
#include "HHKinFit.h"
#include "HHFitConstraint4Vector.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectE.h"
#include "HHFitObject.h"
#include "HHFitObjectComposite.h"
#include "PSMath.h"
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHInvMConstraintException.h"
#include "exceptions/HHLimitSettingException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHKinFit.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint4Vector.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraintEHardM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraintLikelihood.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectEConstM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectComposite.h"
#include "HHKinFit2/HHKinFit2/interface/PSMath.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHInvMConstraintException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHLimitSettingException.h"
#endif

#include "TAxis.h"

#include <iomanip>
#include <iostream>
#include <sstream>

HHKinFit2::HHKinFit::HHKinFit()
: m_fitobjects(std::vector<HHFitObjectE*>()),
  m_constraints(std::vector<HHFitConstraint*>()),
  m_chi2(99999),
  m_convergence(0),
  m_printlevel(0),
  m_maxloops(500){
}

void
HHKinFit2::HHKinFit::setPrintLevel(int printlevel){
  m_printlevel=printlevel;
}

///todo: compare with old kinfit method
void 
HHKinFit2::HHKinFit::fit(){
  
  //  ----------  for PSfit -----
  const int np = m_fitobjects.size();
  double a[np];
  double astart[np];
  double alimit[np][2];
  double aprec[np];
  double daN[np];
  double h[np];
  double chi2iter[1], aMemory[np][5], g[np], H[np * np], Hinv[np * np];
  bool noNewtonShifts = false;
  
  double chi2before;
  double abefore[np];
  double daNbefore[np];
  double hbefore[np];
  double chi2iterbefore[1], aMemorybefore[np][5], gbefore[np], Hbefore[np * np], Hinvbefore[np * np];


  int iter = 0;             //  number of iterations
  int method = 1;           //  initial fit method, see PSfit()
  int mode = 1;             //  mode =1 for start of a new fit by PSfit()
  //   int icallNewton = -1;     //  init start of Newton Method
  //   int iloop = 0;            // counter for falls to fit function
  
  int iterbefore = 0;             //  number of iterations
  int methodbefore = 1;           //  initial fit method, see PSfit()
  int modebefore = 1;             //  mode =1 for start of a new fit by PSfit()
  
  for (unsigned int i=0; i<m_fitobjects.size(); i++){
    //// fill initial tau fit parameters
    //astart[i] = m_fitobjects[i]->getE();    // energy of first tau
    //aprec[i]  = 0.1;                        //0.1                 // precision for fit
    //// fill initial step width
    //h[i] = 0.2*(m_fitobjects[i]->getUpperFitLimitE()-m_fitobjects[i]->getLowerFitLimitE());   //
    //daN[i] = 1.0;   //0.0                 // initial search direction

    astart[i]    = m_fitobjects[i]->getInitStart();
    aprec[i]     = m_fitobjects[i]->getInitPrecision();
    h[i]         = m_fitobjects[i]->getInitStepWidth();
    daN[i]       = m_fitobjects[i]->getInitDirection();
    alimit[i][0] = 1.00001*m_fitobjects[i]->getLowerFitLimitE();
    alimit[i][1] = 0.99999*m_fitobjects[i]->getUpperFitLimitE();

    if( alimit[i][1] - alimit[i][0] < 0.2 ) //Allow for at least 200MeV of wiggle room
    {
      std::stringstream msg;
      msg << "Safety margin of limits for parameter " << i << " too small. Lower Limit: " << alimit[i][0] << "Upper Limit: " << alimit[i][1];
      throw(HHLimitSettingException(msg.str()));
    }

    // Check initial values against fit range
    if (astart[i] - h[i] < alimit[i][0]) {
      astart[i] = alimit[i][0] + h[i];
    }
    else if (astart[i] + h[i] > alimit[i][1]) {
      astart[i] = alimit[i][1] - h[i];
    }

    a[i] = astart[i];
 
    aMemory[i][0] = -999.0;
    aMemory[i][1] = -995.0;
    aMemory[i][2] = -990.0;
    aMemory[i][3] = -985.0;
    aMemory[i][3] = -980.0;
  }
  
  for (int iloop = 0; iloop < m_maxloops * 10 && iter < m_maxloops; iloop++) { // FIT loop
    bool respectLimits = (iter>=0) ; // do not respect limits when calculating numerical derivative
    for (unsigned int i=0; i<m_fitobjects.size();i++){
      try{
	if(isnan(a[i]))
	{
	  std::cout << "WARNING! PSMath changed E of fit object " << i 
		    << "to NAN!" << std::endl;
	  std::stringstream msg;
	  msg << "PS: PSMath changed E of fit object" << i << "to NAN!";
	  throw(HHEnergyRangeException(msg.str()));
	}
        m_fitobjects[i]->changeEandSave(a[i],respectLimits);
      }
      catch(HHKinFit2::HHEnergyRangeException const& e){
        std::cout << e.what() << std::endl;

        std::cout << "iter  ="<< iterbefore << " " << iter << std::endl;
        std::cout << "iloop =" << iloop << std::endl;
        std::cout << "fitobject ="<< i << std::endl;
        std::cout << "method="<< methodbefore << " " << method << std::endl;
        std::cout << "mode  ="<< modebefore << " " << mode << std::endl;
        std::cout << "noNewtonShifts="<< noNewtonShifts << std::endl;
        std::cout << "printlevel="<< m_printlevel << std::endl;
        std::cout << "np="<< np << std::endl;
        std::cout << "chi2="<< chi2before << " " << m_chi2 << std::endl;
        std::cout << "chi2iter[0]="<< chi2iterbefore[0] << " " << chi2iter[0] << std::endl;

        for (unsigned int i=0; i<m_fitobjects.size();i++) std::cout << "a[" << i << "]      ="<< abefore[i] << " " << a[i] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size();i++) std::cout << "astart[" << i << "] ="<< astart[i] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size();i++) std::cout << "aprec[" << i << "]  ="<< aprec[i] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size();i++) std::cout << "daN[" << i << "]    ="<< daNbefore[i] << " " << daN[i] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size();i++) std::cout << "h[" << i << "]      ="<< hbefore[i] << " " << h[i] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size();i++) std::cout << "aMemory[" << i << "]="<< aMemorybefore[i][0] << " "
            << aMemorybefore[i][1] << " "
            << aMemorybefore[i][2] << " "
            << aMemorybefore[i][3] << " "
            << aMemorybefore[i][4] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size();i++) std::cout << "after: aMemory[" << i << "]="<< aMemorybefore[i][0] << " "
            << aMemory[i][1] << " "
            << aMemory[i][2] << " "
            << aMemory[i][3] << " "
            << aMemory[i][4] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size();i++) std::cout << "g[" << i << "]      =" << gbefore[i] << " " << g[i] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size()*m_fitobjects.size();i++) std::cout << "H[" << i << "]      =" << Hbefore[i] << " " << H[i] << std::endl;
        for (unsigned int i=0; i<m_fitobjects.size()*m_fitobjects.size();i++) std::cout << "Hinv[" << i << "]   =" << Hinvbefore[i] <<  " " << Hinv[i] <<std::endl;
        throw(e);
      }
    }

    m_chi2=this->getChi2(respectLimits);

   
    chi2before = m_chi2;
    chi2iterbefore[0]=chi2iter[0];
    iterbefore = iter;
    methodbefore = method;
    modebefore = mode;
    for (unsigned int i=0; i<m_fitobjects.size();i++){
      abefore[i]=a[i];
      daNbefore[i]=daN[i];
      hbefore[i]=h[i];
      aMemorybefore[i][0]=aMemory[i][0];
      aMemorybefore[i][1]=aMemory[i][1];
      aMemorybefore[i][2]=aMemory[i][2];
      aMemorybefore[i][3]=aMemory[i][3];
      aMemorybefore[i][4]=aMemory[i][4];
      gbefore[i]=g[i];
    }
    
    for (unsigned int i=0; i<m_fitobjects.size()*m_fitobjects.size();i++){
      Hbefore[i]=H[i];
      Hinvbefore[i]=Hinv[i];
    }
    
    if (m_convergence != 0) break;
    m_convergence = PSMath::PSfit(iloop, iter, method, mode, noNewtonShifts, m_printlevel,
                                  np, a, astart, alimit, aprec,
                                  daN, h, aMemory, m_chi2, chi2iter, g, H,
                                  Hinv);
  }
  // ------ end of FIT loop
  
  int convergenceFlags = 0;

  if(m_convergence != 0) //Check for convergence at limit
  {
    for(int param = 0; param < np; ++param) //Set bitwise flags for each param
    {
      if(a[param] < (alimit[param][0] + 2*aprec[param]) )
      {
	convergenceFlags = convergenceFlags | (1 << param);
      }
      else if(a[param] > (alimit[param][1] - 2*aprec[param]) )
      {
	convergenceFlags = convergenceFlags | (1 << param);
      }
    }
    //2=at a[0] limit, 3=at a[1] limit 4=at a[0] and a[1] limit 
    m_convergence = m_convergence + convergenceFlags; 
  }
  
  //	  if (m_logLevel>1)
  //	    std::cout << "Convergence is " << m_convergence << std::endl;
  
}

double
HHKinFit2::HHKinFit::getChi2(bool respectLimits) const{
  double chi2=0;

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();
      it != m_constraints.end(); 
      ++it)
  {
    (*it)->prepare(respectLimits);
  }
  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    chi2 += (*it)->getChi2();

  return(chi2);
}

void
HHKinFit2::HHKinFit::printChi2() const{
  double chi2=0;

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    (*it)->prepare();

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it){
    std::cout << (*it)->getChi2() << std::endl;
    chi2 += (*it)->getChi2();
    (*it)->printChi2();
  }
  std::cout << chi2 << std::endl;
  std::cout << "----------------------------------------------------------------------------------------------"<<std::endl;
}

double
HHKinFit2::HHKinFit::getL(bool respectLimits) const{
  double L=1;

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    (*it)->prepare(respectLimits);

  for(std::vector<HHFitConstraint*>::const_iterator it = m_constraints.begin();it != m_constraints.end(); ++it)
    L *= (*it)->getLikelihood();

  return(L);
}


int
HHKinFit2::HHKinFit::getConvergence() const{
  return(m_convergence);
}

std::vector<HHKinFit2::HHFitObjectE*>
HHKinFit2::HHKinFit::getListOfFitObjects() const{
  return(m_fitobjects);
}

std::vector<HHKinFit2::HHFitConstraint*>
HHKinFit2::HHKinFit::getListOfConstraints() const{
  return(m_constraints);
}

void
HHKinFit2::HHKinFit::addFitObjectE(HHFitObjectE* fitobject){
  m_fitobjects.push_back(fitobject);
}

void
HHKinFit2::HHKinFit::addConstraint(HHFitConstraint* constraint){
  m_constraints.push_back(constraint);
}

TGraph*
HHKinFit2::HHKinFit::getChi2Function(int steps){
  int npoints(m_fitobjects.size()*steps);
  TGraph* gr = new TGraph(npoints);
  gr->SetName("chi2function");
  double stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/steps);
  for (int i=0; i<npoints; i++){
    double e = 1.00001*m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
    m_fitobjects[0]->changeEandSave(e);
    double chi2(this->getChi2());
    //    std::cout << i << " " << e << " " << chi2 << std::endl;
    gr->SetPoint(i,e,chi2);
  }
  gr->SetMinimum(0);
  gr->GetXaxis()->SetTitle("E_{1}[GeV]");
  gr->GetYaxis()->SetTitle("#chi^{2}");
  return(gr);
}


TGraph*
HHKinFit2::HHKinFit::getLFunction(int steps){
  int npoints(m_fitobjects.size()*steps);
  TGraph* gr = new TGraph(npoints);
  gr->SetName("Lfunction");
  double stepsize((m_fitobjects[0]->getUpperFitLimitE() - m_fitobjects[0]->getLowerFitLimitE())/steps);
  for (int i=0; i<npoints; i++){
    double e = 1.00001*m_fitobjects[0]->getLowerFitLimitE()+ i*stepsize;
    m_fitobjects[0]->changeEandSave(e);
    double L(this->getL());
    //    std::cout << i << " " << e << " " << chi2 << std::endl;
    gr->SetPoint(i,e,L);
  }
  gr->SetMinimum(0);
  gr->GetXaxis()->SetTitle("E_{1}[GeV]");
  gr->GetYaxis()->SetTitle("L");
  return(gr);
}
