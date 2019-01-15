#pragma once

#include "ABC_MarkovChainSubMatrix.hpp"

namespace MarkovSub
{

class MarkovChainSub : public ABC_MarkovChainSubMatrix
{

public:
  MarkovChainSub(const Json &jjSim, const size_t &seed) : ABC_MarkovChainSubMatrix(jjSim, seed), auxH_(jjSim["model"]["delta"].get<double>()){};

  ~MarkovChainSub(){};

  // //Overriding
  double gammaSubMatrix(const VertexPart &vpTo, const VertexPart &vpFrom) const override
  {
    return auxH_.gamma(vpTo, vpFrom);
  }

  double FAux(const VertexPart &vp) const override
  {
    return (auxH_.FAux(vp));
  }

private:
  Diagrammatic::AuxHelper auxH_;
};
} // namespace MarkovSub