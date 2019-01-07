#pragma once

#include "../Utilities/Utilities.hpp"
#include "../Utilities/GreenTau.hpp"
#include "../Utilities/Matrix.hpp"
#include "../Models/ABC_Model.hpp"
#include "../Utilities/VerticesSimple.hpp"

namespace Markov
{

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
  ISDataCT(const Json &jj, const Models::ABC_Model_2D &model) :
#ifdef AFM
                                                                green0CachedUp_(model.greenCluster0MatUp(), jj),
                                                                green0CachedDown_(model.greenCluster0MatDown(), jj),
#endif
#ifndef AFM
                                                                green0CachedUp_(model.greenCluster0MatUp(), jj),
#endif
                                                                MupPtr_(new Matrix_t()),
                                                                MdownPtr_(new Matrix_t()),
                                                                vertices_(),
                                                                beta_(jj["beta"].get<double>()),
                                                                NOrb_(jj["NOrb"].get<size_t>()),
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

  // #if defined SUBMATRIX
  // friend class Markov::ABC_MarkovChainSubMatrix;
  // #else
  // friend class Markov::ABC_MarkovChain;
  // #endif

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