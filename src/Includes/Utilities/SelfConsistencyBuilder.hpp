#pragma once

#include "SelfConsistency.hpp"

namespace SelfCon
{

std::unique_ptr<ABC_SelfConsistency> SelfConsistencyBuilder(const Json &jj, const FermionSpin_t &spin)
{
    const size_t NOrb = jj["NOrb"].get<size_t>();

    Models::ABC_Model_2D model(jj);
    IO::Base_IOModel ioModel(jj);

    ClusterCubeCD_t greenImpurity;
    if (spin == FermionSpin_t::Up)
    {
        greenImpurity = ioModel.ReadGreenDat("greenUp.dat", NOrb);
    }
    else if (spin == FermionSpin_t::Down)
    {
        greenImpurity = ioModel.ReadGreenDat("greenDown.dat", NOrb);
    }

    using SelfCon_t = SelfCon::SelfConsistency;
    return std::make_unique<SelfCon_t>(SelfCon_t(jj, model, greenImpurity, spin));
}

} // namespace SelfCon