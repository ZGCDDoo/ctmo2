#pragma once

#ifdef GREEN_STYLE
#include "ABC_MarkovChain_GreenStyle.hpp"
#else
#include "ABC_MarkovChain.hpp"
#endif

#include "../Utilities/VerticesSimple.hpp"

namespace Markov
{

class MarkovChain : public ABC_MarkovChain
{

public:
  MarkovChain(const Json &jj, const size_t &seed) : ABC_MarkovChain(jj, seed), auxH_(jj["model"]["delta"].get<double>()){};

  ~MarkovChain(){};

  //Overriding
  // double gammaTrad(const FermionSpin_t &spin, const AuxSpin_t &auxTo, const AuxSpin_t &auxFrom) override { return (auxH_.gamma(spin, auxTo, auxFrom)); }

#ifdef GREEN_STYLE
  double FAux(const VertexPart &vp) override
  {
    if (vp.vtype() == Diagrammatic::VertexType::Phonon)
    {
      return (auxH_.auxPh(vp.aux()));
    }
    else
    {
      return (auxH_.auxValue(vp.spin(), vp.aux()));
    }
  }

  double FAuxBar(const VertexPart &vp) override
  {
    if (vp.vtype() == Diagrammatic::VertexType::Phonon)
    {
      return (auxH_.auxPh(vp.aux()));
    }
    else
    {
      return (auxH_.auxValueBar(vp.spin(), vp.aux()));
    }
  }

#else
  double FAux(const VertexPart &vp) override
  {
    return (auxH_.FAux(vp));
  }

  double FAuxBar(const VertexPart &vp) override
  {
    return (auxH_.FAuxBar(vp));
  }
#endif

private:
  Diagrammatic::AuxHelper auxH_;
};
} // namespace Markov
