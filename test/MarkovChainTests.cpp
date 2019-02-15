#include <gtest/gtest.h>

#include "ctmo/ImpuritySolver/MarkovChain.hpp"

using namespace LinAlg;

const double DELTA = 1e-10;
const std::string FNAME = "../../test/data/DMFT/params1.json";
using Model_t = Models::ABC_Model_2D;
using IOModel_t = IO::Base_IOModel;

Markov::MarkovChain BuildMarkovChain() // for SIAM_Square
{
    std::ifstream fin(FNAME);
    Json jj;
    fin >> jj;
    fin.close();
    std::cout << "Reading in Json in BuildMarkovChain() " << std::endl;
    const size_t seed = 10224;
    Markov::MarkovChain markovchain(jj, seed);
    std::cout << "After BuildMarkovChain() " << std::endl;
    return markovchain;
}

TEST(MarkovChainTests, Init) { Markov::MarkovChain mc = BuildMarkovChain(); }

TEST(MonteCarloTest, DoStep)
{
    Markov::MarkovChain mc = BuildMarkovChain();

    for (size_t i = 0; i < 2; i++)
    {
        mc.InsertVertex();
    }

    std::cout << "After Insert " << std::endl;
    for (size_t i = 0; i < 1001; i++)
    {
        mc.InsertVertex();
        mc.RemoveVertex();
    }

    std::cout << "After Remove " << std::endl;
    size_t ii = 0;
    for (ii = 0; ii < 15000; ii++)
    {
        mc.DoStep();
    }
    std::cout << "ii = " << ii << std::endl;

    std::cout << "After DOstep " << std::endl;

    Matrix_t tmpUp;
    Matrix_t tmpDown;
    tmpUp = mc.Nup();
    tmpDown = mc.Ndown();
    std::cout << "Mup.size = " << tmpUp.n_cols() << std::endl;
    mc.CleanUpdate();

    for (size_t i = 0; i < tmpUp.n_rows(); i++)
    {
        for (size_t j = 0; j < tmpUp.n_rows(); j++)
        {
            ASSERT_NEAR(tmpUp(i, j), mc.Nup()(i, j), DELTA);
            // ASSERT_NEAR(tmpDown(i, j), mc.Ndown()(i, j), DELTA);
        }
    }

    for (size_t i = 0; i < tmpDown.n_rows(); i++)
    {
        for (size_t j = 0; j < tmpDown.n_rows(); j++)
        {
            ASSERT_NEAR(tmpDown(i, j), mc.Ndown()(i, j), DELTA);
        }
    }

    ASSERT_EQ(tmpUp.n_cols(), mc.Nup().n_rows());
    ASSERT_EQ(tmpUp.n_cols(), tmpUp.n_rows());
    std::cout << "dims = " << tmpUp.n_cols() << std::endl;
    mc.SaveTherm();
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
