

#include <gtest/gtest.h>
#include "../src/Includes/Utilities/Integrator.hpp"
#include "../src/Includes/Models/ABC_H0.hpp"
#include "TestTools.hpp"

using namespace Integrator;

const double DELTA = 1.49e-8;
const double t_diag = -1.0;
const double tp_diag = -0.4;

const double t_offdiag = -1.02;
const double tp_offdiag = 0.230;

const size_t Nc = 4;
const size_t Nx = 2;
const size_t Ny = 2;
const size_t NOrb = 2;
const double INTRA = -0.01;

struct FctTest
{
  public:
    static const size_t Nc = 4;
    static const size_t TNX = 2;
    static const size_t TNY = 2;

    FctTest(){};

    size_t n_rows() const { return Nc; };
    size_t n_cols() const { return Nc; };

    ClusterMatrixCD_t operator()(const double &x, const double &y, const double &z)
    {
        ClusterMatrixCD_t tmp(Nc, Nc);
        tmp.zeros();
        for (size_t i = 0; i < Nc; i++)
        {
            for (size_t j = 0; j < Nc; j++)
            {
                tmp(i, j) = cd_t(i * x + i * i * x * x + z, j * y);
            }
        }
        return tmp;
    }
};

TEST(IntegratorTest, TestGridKTilde)
{

    Models::ABC_H0 h0(TestTools::BuildJson());
    ClusterMatrixCD_t tLocTest = GridKTilde(h0, 100);

    ClusterMatrixCD_t diagonal_part = {{0.0, t_diag, t_diag, tp_diag},
                                       {t_diag, 0.0, 0.0, t_diag},
                                       {t_diag, 0.0, 0.0, t_diag},
                                       {tp_diag, t_diag, t_diag, 0.0}};

    ClusterMatrixCD_t off_diagonal = {{0.0, t_offdiag, t_offdiag, tp_offdiag},
                                      {t_offdiag, 0.0, tp_offdiag, t_offdiag},
                                      {t_offdiag, tp_offdiag, 0.0, t_offdiag},
                                      {tp_offdiag, t_offdiag, t_offdiag, 0.0}};

    ClusterMatrixCD_t goodResult(Nc * NOrb, Nc * NOrb);
    goodResult.zeros();

    goodResult.submat(0, 0, Nc - 1, Nc - 1) = diagonal_part;
    goodResult.submat(Nc, Nc, NOrb * Nc - 1, NOrb * Nc - 1) = diagonal_part;
    goodResult.submat(0, Nc, Nc - 1, NOrb * Nc - 1) = off_diagonal + INTRA * ClusterMatrixCD_t(Nc, Nc).eye();
    goodResult.submat(Nc, 0, NOrb * Nc - 1, Nc - 1) = off_diagonal + INTRA * ClusterMatrixCD_t(Nc, Nc).eye();

    assert(goodResult.n_rows == tLocTest.n_rows);
    for (size_t i = 0; i < goodResult.n_cols; i++)
    {
        for (size_t j = 0; j < goodResult.n_rows; j++)
        {

            ASSERT_NEAR(tLocTest(i, j).real(), goodResult(i, j).real(), DELTA);
            ASSERT_NEAR(tLocTest(i, j).imag(), goodResult(i, j).imag(), DELTA);
        }
    }
}

TEST(IntegratorTest, TestCubature)
{
    double xmin2[3] = {-1.1, 0.0, 0.0};
    double xmax2[3] = {1.1, 1.0, 1.0};
    FctTest fcttest;

    ClusterMatrixCD_t test2 = Cubature<FctTest>(fcttest, xmin2, xmax2, 0, 1.49e-8, 1.49e-8);

    for (size_t i = 0; i < test2.n_rows; i++)
    {
        for (size_t j = 0; j < test2.n_cols; j++)
        {
            cd_t tmpGood = cd_t(1.1 + 0.88733333333333 * static_cast<double>(i * i), 1.1 * static_cast<double>(j));
            ASSERT_NEAR(test2(i, j).real(), tmpGood.real(), DELTA);
            ASSERT_NEAR(test2(i, j).imag(), tmpGood.imag(), DELTA);
        }
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
