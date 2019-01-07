#pragma once

#include "ABC_Model.hpp"

namespace Models
{

class SIAM_Square : public ABC_Model_2D<IO::IOSIAM>
{
  using IOModel_t = IO::IOSIAM;

public:
  SIAM_Square(const Json &jj) : ABC_Model_2D<IOModel_t>(jj){};
};
} // namespace Models