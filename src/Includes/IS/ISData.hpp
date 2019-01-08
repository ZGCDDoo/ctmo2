#pragma once

#include "../Utilities/Utilities.hpp"
#include "../Utilities/GreenTau.hpp"
#include "../Utilities/Matrix.hpp"
#include "../Models/ABC_Model.hpp"
#include "../Utilities/VerticesSimple.hpp"

namespace Markov
{

#if defined SUBMATRIX
class ABC_MarkovChainSubMatrix;
#else
class ABC_MarkovChain;
#endif

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
  ISDataCT(const Json &jjSim, const Models::ABC_Model_2D &model) :
#ifdef AFM
                                                                   green0CachedUp_(model.greenCluster0MatUp(), jjSim),
                                                                   green0CachedDown_(model.greenCluster0MatDown(), jjSim),
#endif
#ifndef AFM
                                                                   green0CachedUp_(model.greenCluster0MatUp(), jjSim),
#endif
                                                                   MupPtr_(new Matrix_t()),
                                                                   MdownPtr_(new Matrix_t()),
                                                                   vertices_(),
                                                                   beta_(model.beta()),
                                                                   NOrb_(model.NOrb()),
                                                                   sign_(1)

  {
  }

  double beta() const
  {
    return beta_;
  };
  double NOrb() const
  {
    return NOrb_;
  };

private:
  friend class Markov::Obs::Observables;
  friend class Markov::Obs::GreenBinning;
  friend class Markov::Obs::FillingAndDocc;
  // friend class Markov::Obs::GreenTauMesure;

#if defined SUBMATRIX
  friend class Markov::ABC_MarkovChainSubMatrix;
#else
  friend class Markov::ABC_MarkovChain;
#endif

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