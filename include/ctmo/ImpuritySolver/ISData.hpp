#pragma once

#include "ctmo/Foundations/GreenTau.hpp"
#include "ctmo/Foundations/Matrix.hpp"
#include "ctmo/Model/ABC_Model.hpp"
#include "ctmo/ImpuritySolver/VerticesSimple.hpp"

namespace Markov
{

class ABC_MarkovChain;

namespace Obs
{

class FillingAndDocc;
class GreenBinning;
class Observables;

using namespace LinAlg;

class ISDataCT
{
    using GreenTau_t = GreenTau::GreenCluster0Tau;

  public:
    ISDataCT(const Json &jjSim, const std::shared_ptr<Models::ABC_Model_2D> &modelPtr)
        : modelPtr_(modelPtr),
#ifdef AFM
          green0CachedUp_(modelPtr->greenCluster0MatUp(), modelPtr_->ioModelPtr(), jjSim["solver"]["ntau"]),
          green0CachedDown_(modelPtr->greenCluster0MatDown(), modelPtr_->ioModelPtr(), jjSim["solver"]["ntau"]),
#endif
#ifndef AFM
          green0CachedUp_(modelPtr->greenCluster0MatUp(), modelPtr_->ioModelPtr(), jjSim["solver"]["ntau"]),
#endif
          MupPtr_(new Matrix_t()), MdownPtr_(new Matrix_t()), beta_(modelPtr->beta()), NOrb_(modelPtr->NOrb()), sign_(1)

    {
        Logging::Trace("ISData Created. ");
    }

    double beta() const { return beta_; };
    double NOrb() const { return NOrb_; };

  private:
    friend class Markov::Obs::Observables;
    friend class Markov::Obs::GreenBinning;
    friend class Markov::Obs::FillingAndDocc;

    friend class Markov::ABC_MarkovChain;

    std::shared_ptr<Models::ABC_Model_2D> modelPtr_;
    GreenTau_t green0CachedUp_;
#ifdef AFM
    GreenTau_t green0CachedDown_;
#endif
    std::shared_ptr<Matrix_t> MupPtr_;
    std::shared_ptr<Matrix_t> MdownPtr_;
    Diagrammatic::Vertices vertices_;
    const double beta_;
    const size_t NOrb_;
    Sign_t sign_;
};
} // namespace Obs
} // namespace Markov