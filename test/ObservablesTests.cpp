#define SUBMATRIX

#include <gtest/gtest.h>

#include "../src/Includes/IS/Obs/Observables.hpp"

using Model_t = Models::ABC_Model_2D;
using IOModel_t = IO::Base_IOModel;
using GreenBinning_t = Markov::Obs::GreenBinning;
using Obs_t = Markov::Obs::Observables;
using ISDataCT_t = Markov::Obs::ISDataCT;

// const double DELTA = 1e-11;
const std::string FNAME = "../test/data/cdmft_square2x2/params1.json";

Obs_t BuildObs() //for Square2x2
{

    std::ifstream fin(FNAME);
    Json jj;
    fin >> jj;
    fin.close();
    Model_t model(jj);
    std::shared_ptr<ISDataCT_t> dataCT(
        new ISDataCT_t(
            jj,
            model));

    std::shared_ptr<Model_t> modelPtr(new Model_t(jj));

    Obs_t obs(dataCT, jj);
    return obs;
}

TEST(ObsTests, Init)
{
    Obs_t Obs = BuildObs();
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
