#pragma once

#include "ABC_Model.hpp"
#include "ABC_H0.hpp"

namespace Models
{

class ModelSquare8x8 : public ABC_Model_2D<IO::IOSquare8x8, ABC_H0<Nx8, Nx8>>
{

  using H0_t = ABC_H0<Nx8, Nx8>;
  using IOModel_t = IO::IOSquare8x8;

public:
  ModelSquare8x8(const Json &jj) : ABC_Model_2D<IOModel_t, H0_t>(jj){};
};
} // namespace Models