
#include "Includes/IS/MonteCarloBuilder.hpp"
#include "Includes/Utilities/SelfConsistencyBuilder.hpp"
#include "Includes/Utilities/FS.hpp"
#include "Includes/PrintVersion.hpp"
#include "Includes/Utilities/Logging.hpp"

int main(int argc, char **argv)
{

#ifdef HAVEMPI
    mpi::environment env;
    mpi::communicator world;
#endif

    if (argc != 3)
    {
        throw std::runtime_error("Miseria: Wrong number of input parameters. Stupido !");
    }

    const std::string paramsName = argv[1];
    const int ITER = atoi(argv[2]);
    const std::string fname_params = paramsName + std::to_string(ITER) + std::string(".json");
    Json jjSim;

#ifndef HAVEMPI
    PrintVersion::PrintVersion();
    Logging::Init();
    Logging::Info("Iteration " + std::to_string(ITER));
    std::ifstream fin(fname_params);
    fin >> jjSim;
    fin.close();
    std::cout << "Iter = " << ITER << std::endl;
    const size_t seed = jjSim["monteCarlo"]["seed"].get<size_t>();

    //init a model, to make sure all the files are present and that not all proc write to the same files

    const std::unique_ptr<MC::ABC_MonteCarlo> monteCarloMachinePtr = MC::MonteCarloBuilder(jjSim, seed);

    monteCarloMachinePtr->RunMonteCarlo();

    const std::unique_ptr<SelfCon::ABC_SelfConsistency> selfconUpPtr = SelfCon::SelfConsistencyBuilder(jjSim, FermionSpin_t::Up);
    selfconUpPtr->DoSCGrid();

    IO::FS::PrepareNextIter(paramsName, ITER);

#endif

#ifdef HAVEMPI

    std::string jjSimStr;

    if (mpiUt::Rank() == mpiUt::master)
    {
        PrintVersion::PrintVersion();
        mpiUt::Print("ITER = " + std::to_string(ITER));
        std::ifstream fin(fname_params);
        fin >> jjSim;
        jjSimStr = jjSim.dump();
        fin.close();
    }

    mpi::broadcast(world, jjSimStr, mpiUt::master);
    jjSim = Json::parse(jjSimStr);
    world.barrier();
    //wait_all
    const size_t rank = world.rank();
    const size_t seed = jjSim["monteCarlo"]["seed"].get<size_t>() + 2797 * rank;

    {
        std::unique_ptr<MC::ABC_MonteCarlo> monteCarloMachinePtr = MC::MonteCarloBuilder(jjSim, seed);

        monteCarloMachinePtr->RunMonteCarlo();
    }

    world.barrier();

    const std::unique_ptr<SelfCon::ABC_SelfConsistency> selfconUpPtr = SelfCon::SelfConsistencyBuilder(jjSim, FermionSpin_t::Up);
    selfconUpPtr->DoSCGrid();

    if (mpiUt::Rank() == mpiUt::master)
    {
        IO::FS::PrepareNextIter(paramsName, ITER);
    }
#endif

    return EXIT_SUCCESS;
}
