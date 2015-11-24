#include "HHEventGenerator.h"
#include "HHLorentzVector.h"

#include "TF1.h"
#include <iostream>

using HHKinFit2::HHEventGenerator;

int main(int argc, char* argv[])
{
  TF1* pdf1 = new TF1("PDF1","2*x",0,2);
  TF1* pdf2 = new TF1("PDF2","2*x",0,2);
  HHEventGenerator testgenerator(300,125,125,pdf1,pdf2);

  int events = 10;

  for(unsigned int i=0; i<events; i++){
    std::cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << std::endl;
    testgenerator.generateEvent();
  }

  return(0);
}
