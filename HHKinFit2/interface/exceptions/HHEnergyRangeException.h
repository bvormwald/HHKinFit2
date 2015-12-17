#ifndef HHEnergyRangeException_
#define HHEnergyRangeException_

#include <iostream>
#include <stdexcept>
#include <string>

namespace HHKinFit2{
class HHEnergyRangeException: public std::runtime_error
{
public:
 HHEnergyRangeException(const std::string message)
   : std::runtime_error(message) {};
};
}
#endif
