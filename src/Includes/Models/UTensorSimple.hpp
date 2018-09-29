#pragma once

#include "../Utilities/Utilities.hpp"

namespace Models
{

class UTensor
{

public:
  UTensor(const Json &jjSim) : NOrb_(jjSim["NOrb"].get<size_t>()),
                               U_(jjSim["U"].get<double>()),
                               JH_(jjSim["J_H"].get<double>()),
                               UPrime_(U_ - 2 * JH_),
                               auxMu_(jjSim["mu"].get<double>() - U_ / 2.0 - (NOrb_ == 1 ? 0.0 : NOrb_ * (2.0 * UPrime_ - JH_) / 2.0))

  {
  }

  ~UTensor()
  {
  }

  double auxMu() const { return auxMu_; };
  double U() const { return U_; };
  double JH() const { return JH_; };
  double UPrime() const { return UPrime_; };

protected:
  const size_t NOrb_;

  const size_t U_;
  const size_t JH_;
  const size_t UPrime_;
  const double auxMu_;
};

} // namespace Models
