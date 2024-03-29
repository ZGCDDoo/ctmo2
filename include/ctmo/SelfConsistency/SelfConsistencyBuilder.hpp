#pragma once

#ifdef DCA
#include "ctmo/SelfConsistency/SelfConsistency_DCA.hpp"
#else
#include "ctmo/SelfConsistency/SelfConsistency_CDMFT.hpp"
#endif

namespace SelfCon
{

std::unique_ptr<ABC_SelfConsistency> SelfConsistencyBuilder(const Json &jjSim, const FermionSpin_t &spin)
{

    const size_t NOrb = jjSim["model"]["nOrb"].get<size_t>();
    Models::ABC_Model_2D model(jjSim);
    IO::Base_IOModel ioModel(jjSim);
    ClusterCubeCD_t greenImpurity;

    if (spin == FermionSpin_t::Up)
    {
        greenImpurity = ioModel.ReadGreenDat("greenUp.dat", NOrb);
    }
    else if (spin == FermionSpin_t::Down)
    {
        greenImpurity = ioModel.ReadGreenDat("greenDown.dat", NOrb);
    }

#ifdef DCA

    using SelfCon_t = SelfCon::SelfConsistencyDCA;
    return std::make_unique<SelfCon_t>(SelfCon_t(jjSim, model, greenImpurity, spin));

#else

    using SelfCon_t = SelfCon::SelfConsistency;
    return std::make_unique<SelfCon_t>(SelfCon_t(jjSim, model, greenImpurity, spin));

#endif
}

} // namespace SelfCon