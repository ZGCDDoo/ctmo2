#pragma once

#include "ABC_MarkovChain.hpp"
#include "../Utilities/Vertices.hpp"

namespace Markov
{

template <typename TIOModel, typename TModel>
class MarkovChain : public ABC_MarkovChain<TIOModel, TModel>
{

public:
  MarkovChain(const Json &jj, const size_t &seed) : ABC_MarkovChain<TIOModel, TModel>(jj, seed), auxH_(jj["delta"].get<double>()){};

  ~MarkovChain(){};

  //Overriding
  double gammaTrad(const FermionSpin_t &spin, const AuxSpin_t &auxTo, const AuxSpin_t &auxFrom) override { return (auxH_.gamma(spin, auxTo, auxFrom)); }
  double FAux(const FermionSpin_t &spin, const AuxSpin_t &aux) override { return (auxH_.FAux(spin, aux)); }

private:
  Diagrammatic::AuxHelper auxH_;
};
} // namespace Markov
