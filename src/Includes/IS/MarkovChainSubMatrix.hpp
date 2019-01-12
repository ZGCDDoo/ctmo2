#pragma once

#include "ABC_MarkovChainSubMatrix.hpp"

namespace MarkovSub
{

class MarkovChainSub : public ABC_MarkovChainSubMatrix
{

public:
  MarkovChainSub(const Json &jj, const size_t &seed) : ABC_MarkovChainSubMatrix(jj, seed){};

  ~MarkovChainSub(){};

  // //Overriding
  // double gammaUpSubMatrix(const AuxSpin_t &auxTo, const AuxSpin_t &auxFrom) override { return (this->model_.gammaUp(auxTo, auxFrom)); }
  // double gammaDownSubMatrix(const AuxSpin_t &auxTo, const AuxSpin_t &auxFrom) override { return (this->model_.gammaDown(auxTo, auxFrom)); }
  // double FAux(const VertexPart &vp) override
  // {
  //   return 0.0; //return (auxH_.FAux(vp));
  // }
};
} // namespace Markov