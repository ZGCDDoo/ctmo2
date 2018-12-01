#pragma once

#include "../Utilities/Utilities.hpp"

namespace Models
{

class UTensor
{

public:
  UTensor(const Json &jjSim) : NOrb_(jjSim["NOrb"].get<size_t>()),
                               U_(jjSim["U"].get<double>()),
                               JH_(jjSim["IsOrbitalDiagonal"].get<bool>() ? 0.0 : jjSim["J_H"].get<double>()),
                               UPrime_(jjSim["IsOrbitalDiagonal"].get<bool>() ? 0.0 : jjSim["UPrime"].get<double>()), //Normaly, U_ - 2 * JH_),
                               w0Phonon_(jjSim["w0Phonon"].get<double>()),
                               gPhonon_(jjSim["gPhonon"].get<double>())

  {
    if (jjSim["IsOrbitalDiagonal"])
    {
      auxMu_ = jjSim["mu"].get<double>() - U_ / 2.0;
    }
    else if (std::abs(JH_) < 1e-10)
    {
      auxMu_ = jjSim["mu"].get<double>() - U_ / 2.0 - static_cast<double>(NOrb_ - 1) * UPrime_ / 2.0;
    }
    else
    {
      auxMu_ = jjSim["mu"].get<double>() - U_ / 2.0 - static_cast<double>(NOrb_ - 1) * UPrime_ / 2.0 - static_cast<double>(NOrb_ - 1) * (UPrime_ - JH_) / 2.0;
    }
  }

  ~UTensor()
  {
  }

  double auxMu() const { return auxMu_; };

  double U() const { return U_; };

  double JH() const { return JH_; };

  double UPrime() const { return UPrime_; };
  double gPhonon() const { return gPhonon_; };
  double w0Phonon() const { return w0Phonon_; };

protected:
  const size_t NOrb_;

  const double U_;
  const double JH_;
  const double UPrime_;
  const double w0Phonon_;
  const double gPhonon_;

  double auxMu_;
};

} // namespace Models
