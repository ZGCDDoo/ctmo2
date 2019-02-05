#pragma once

#include "../Utilities/Logging.hpp"

namespace Models
{

class UTensor
{

  public:
    explicit UTensor(const Json &jjSim)
        : NOrb_(jjSim["model"]["nOrb"].get<size_t>()), U_(jjSim["model"]["U"].get<double>()),
          JH_(jjSim["solver"]["isOrbitalDiagonal"].get<bool>() ? 0.0 : jjSim["model"]["J_H"].get<double>()),
          UPrime_(jjSim["solver"]["isOrbitalDiagonal"].get<bool>() ? 0.0
                                                                   : jjSim["model"]["UPrime"].get<double>()), // Normaly, U_ - 2 * JH_),
          w0Phonon_(jjSim["model"]["w0Phonon"].get<double>()), gPhonon_(jjSim["model"]["gPhonon"].get<double>())

    {
        if (jjSim["solver"]["isOrbitalDiagonal"])
        {
            auxMu_ = jjSim["model"]["mu"].get<double>() - U_ / 2.0;
        }
        else if (std::abs(JH_) < 1e-10)
        {
            auxMu_ = jjSim["model"]["mu"].get<double>() - U_ / 2.0 - static_cast<double>(NOrb_ - 1) * UPrime_ / 2.0;
        }
        else
        {
            auxMu_ = jjSim["model"]["mu"].get<double>() - U_ / 2.0 - static_cast<double>(NOrb_ - 1) * UPrime_ / 2.0 -
                     static_cast<double>(NOrb_ - 1) * (UPrime_ - JH_) / 2.0;
        }

        Logging::Debug("UTensor Created. ");
    }

    ~UTensor() {}

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
