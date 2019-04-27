#define SLMC
//#define HAVE_POSTGRES


#include "ctmo/MonteCarlo/MonteCarloBuilder.hpp"
#include "ctmo/SelfConsistency/SelfConsistencyBuilder.hpp"
#include "ctmo/Foundations/FS.hpp"
#include "ctmo/Foundations/PrintVersion.hpp"
#include "ctmo/Foundations/CMDParser.hpp"


#include <csignal>


void SIGINT_Handler(int sigNum)
{

    Logging::Info("Caught signal " + std::to_string(sigNum));
    Logging::Info("If MonteCarloBuilder returned singleton, could exit more gracefully by saving");

#ifdef HAVEMPI
    mpi::environment env;
    mpi::communicator world;
    world.barrier();
#endif

    std::cout << "Rank = " << mpiUt::Tools::Rank() << std::endl;

    exit(sigNum);
}

int main(int argc, char **argv)
{

#ifdef HAVEMPI
    mpi::environment env;
    mpi::communicator world;
#endif


    signal(SIGINT, SIGINT_Handler);

    CMDParser::CMDInfo cmdInfo;
    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        cmdInfo = CMDParser::GetProgramOptions(argc, argv);
    }

    int ITER = cmdInfo.iter();
    bool exitFromCMD = cmdInfo.exitFromCMD();
    std::string fnameParams = cmdInfo.fileName();
    Json jjSim;

#ifndef HAVEMPI

    if (exitFromCMD)
    {
        return EXIT_SUCCESS;
    }

    std::ifstream fin(fnameParams);
    fin >> jjSim;
    fin.close();

    Logging::Init(jjSim["logging"]);
    Logging::Info(PrintVersion::GetVersion());
    Logging::Info("Iteration " + std::to_string(ITER));
    const size_t seed = jjSim["monteCarlo"]["seed"].get<size_t>();

    // init a model, to make sure all the files are present and that not all proc write to the same files

    Logging::Trace("ABC_MonteCarlo Creation...");
    const std::unique_ptr<MC::ABC_MonteCarlo> monteCarloMachinePtr = MC::MonteCarloBuilder(jjSim, seed);
    Logging::Trace("ABC_MonteCarlo Created !");

    monteCarloMachinePtr->RunMonteCarlo();

#endif

#ifdef HAVEMPI

    mpi::broadcast(world, exitFromCMD, mpiUt::Tools::master);
    if (exitFromCMD)
    {
        return EXIT_SUCCESS;
    }

    std::string jjSimStr;

    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {

        std::ifstream fin(fnameParams);
        fin >> jjSim;
        jjSimStr = jjSim.dump();
        fin.close();
    }

    mpi::broadcast(world, jjSimStr, mpiUt::Tools::master);
    mpi::broadcast(world, ITER, mpiUt::Tools::master);

    jjSim = Json::parse(jjSimStr);

    Logging::Init(jjSim["logging"]);
    Logging::Info(PrintVersion::GetVersion());
    Logging::Info("Iteration " + std::to_string(ITER));
    world.barrier();
    // wait_all

    const size_t rank = world.rank();
    const size_t seed = jjSim["monteCarlo"]["seed"].get<size_t>() + 2797 * rank;

    Logging::Trace("ABC_MonteCarlo Creation...");
    std::unique_ptr<MC::ABC_MonteCarlo> monteCarloMachinePtr = MC::MonteCarloBuilder(jjSim, seed);
    Logging::Trace("ABC_MonteCarlo Created...");

    monteCarloMachinePtr->RunMonteCarlo();
    world.barrier();

#endif

    return EXIT_SUCCESS;
}
