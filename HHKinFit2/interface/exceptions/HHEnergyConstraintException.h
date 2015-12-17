#ifndef HHEnergyConstraintException_
#define HHEnergyConstraintException_

#include <iostream>
#include <stdexcept>
#include <string>

namespace HHKinFit2{
class HHEnergyConstraintException: public std::runtime_error
{
public:
  HHEnergyConstraintException(const std::string message)
   : std::runtime_error(message) {};
};
}
#endif
