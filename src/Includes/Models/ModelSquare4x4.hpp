#pragma once

#include "ABC_Model.hpp"
#include "ABC_H0.hpp"

namespace Models
{

class ModelSquare4x4 : public ABC_Model_2D<IO::IOSquare4x4, ABC_H0<Nx4, Nx4>>
{

  using H0_t = ABC_H0<Nx4, Nx4>;
  using IOModel_t = IO::IOSquare4x4;

public:
  ModelSquare4x4(const Json &jj) : ABC_Model_2D<IOModel_t, H0_t>(jj){};
};
} // namespace Models