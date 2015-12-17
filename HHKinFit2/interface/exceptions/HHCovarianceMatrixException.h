#ifndef HHCovarianceMatrixException_
#define HHCovarianceMatrixException_

#include <iostream>
#include <stdexcept>
#include <string>

namespace HHKinFit2{
class HHCovarianceMatrixException: public std::runtime_error
{
public:
  HHCovarianceMatrixException(const std::string message)
   : std::runtime_error(message) {};
};
}

#endif
