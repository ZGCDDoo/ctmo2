#pragma once

#include "ctmo/MonteCarlo/MonteCarlo.hpp"
#include "ctmo/ImpuritySolver/MarkovChain.hpp"

namespace MC
{

std::unique_ptr<ABC_MonteCarlo> MonteCarloBuilder(const Json &jjSim, const size_t &seed)
{
#ifdef HAVEMPI
    mpi::environment env;
    mpi::communicator world;
#endif

    using Model_t = Models::ABC_Model_2D;
    using MarkovInt_t = Markov::MarkovChain;

    // Init a dummy model just to be sure that all files are present:
    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        const Model_t modelDummy(jjSim);
    }
#ifdef HAVEMPI
    world.barrier();
#endif

    return std::make_unique<MC::MonteCarlo<MarkovInt_t>>(std::make_shared<MarkovInt_t>(jjSim, seed), jjSim);
}

} // namespace MC
