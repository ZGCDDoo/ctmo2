#pragma once

#include "ABC_Model.hpp"
#include "ABC_H0.hpp"

namespace Models
{

class ModelSquare4x4_DCA : public ABC_Model_2D<IO::IOSquare4x4_DCA, ABC_H0<Nx4, Nx4>>
{

  using H0_t = ABC_H0<Nx4, Nx4>;
  using IOModel_t = IO::IOSquare4x4_DCA;

public:
  static const size_t Nc = H0_t::Nc;

  ModelSquare4x4_DCA(const Json &jj) : ABC_Model_2D<IOModel_t, H0_t>(jj){};
};
} // namespace Models