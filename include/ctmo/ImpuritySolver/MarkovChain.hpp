#pragma once

#include "ctmo/ImpuritySolver/ABC_MarkovChain.hpp"
#include "ctmo/ImpuritySolver/VerticesSimple.hpp"

namespace Markov
{

class MarkovChain : public ABC_MarkovChain
{

  public:
    MarkovChain(const Json &jjSim, const size_t &seed) : ABC_MarkovChain(jjSim, seed), auxH_(jjSim["model"]["delta"].get<double>()){};

    MarkovChain(const MarkovChain &markovChain) = default;
    MarkovChain(MarkovChain &&markovChain) = default;
    MarkovChain &operator=(const MarkovChain &markovChain) = delete;
    MarkovChain &operator=(MarkovChain &&markovChain) = delete;

    ~MarkovChain() override = default;

#ifdef GREEN_STYLE
    double FAux(const VertexPart &vp) const override
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

    double FAuxBar(const VertexPart &vp) const override
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
    double FAux(const VertexPart &vp) const override { return (auxH_.FAux(vp)); }

    double FAuxBar(const VertexPart &vp) const override { return (auxH_.FAuxBar(vp)); }

    double gamma(const VertexPart &vpI, const VertexPart &vpJ) const override { return auxH_.gamma(vpI, vpJ); }

#endif

  private:
    Diagrammatic::AuxHelper auxH_;
};
} // namespace Markov
