#ifndef HHInvMConstraintException_
#define HHInvMConstraintException_

#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace HHKinFit2{
  class HHInvMConstraintException: public HHEnergyRangeException
{
public:
 HHInvMConstraintException(const std::string message)
   : HHEnergyRangeException(message) {};
};
}
#endif
