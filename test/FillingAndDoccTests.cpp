

#include <gtest/gtest.h>

#include "../src/Includes/IS/Obs/FillingAndDocc.hpp"

using Model_t = Models::ABC_Model_2D;
using IOModel_t = IO::Base_IOModel;
using FillingAndDocc_t = Markov::Obs::FillingAndDocc;
using ISDataCT_t = Markov::Obs::ISDataCT;

// const double DELTA = 1e-11;
const std::string FNAME = "../test/data/cdmft_square2x2/params1.json";

FillingAndDocc_t BuildFillingAndDocc() // for Square2x2
{

    std::ifstream fin(FNAME);
    Json jj;
    fin >> jj;
    fin.close();
    std::shared_ptr<Model_t> modelPtr(new Model_t(jj));

    std::shared_ptr<ISDataCT_t> dataCT(new ISDataCT_t(jj, modelPtr));

    Utilities::EngineTypeFibonacci3217_t rng(0);
    std::shared_ptr<Utilities::UniformRngFibonacci3217_t> urngPtr(
        new Utilities::UniformRngFibonacci3217_t(rng, Utilities::UniformDistribution_t(0.0, 1.0)));

    const size_t N_T_INV = 5;
    FillingAndDocc_t fillingAndDocc(dataCT, urngPtr, N_T_INV);
    return fillingAndDocc;
}

TEST(FillingAndDoccTests, Init) { FillingAndDocc_t fillingAndDocc = BuildFillingAndDocc(); }

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
