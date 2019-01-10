#define SUBMATRIX

#include <gtest/gtest.h>

#include "../src/Includes/IS/Obs/GreenBinning.hpp"

using Model_t = Models::ABC_Model_2D;
using IOModel_t = IO::Base_IOModel;
using GreenBinning_t = Markov::Obs::GreenBinning;
using ISDataCT_t = Markov::Obs::ISDataCT;

// const double DELTA = 1e-11;
const std::string FNAME = "../test/data/cdmft_square2x2/params1.json";

GreenBinning_t BuildGreenBinning() //for Square2x2
{

    std::ifstream fin(FNAME);
    Json jj;
    fin >> jj;
    fin.close();
    std::shared_ptr<Model_t> modelPtr(new Model_t(jj));

    std::shared_ptr<ISDataCT_t> dataCT(
        new ISDataCT_t(
            jj,
            modelPtr));

    std::cout << "Here 1 " << std::endl;
    GreenBinning_t greenBinning(dataCT, jj, FermionSpin_t::Up);
    std::cout << "Here 2 " << std::endl;

    return greenBinning;
}

TEST(GreenBinningTests, Init)
{
    GreenBinning_t greenBinning = BuildGreenBinning();
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
