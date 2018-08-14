

#include <gtest/gtest.h>
#include "../src/Includes/Utilities/Integrator.hpp"
#include "../src/Includes/Models/ABC_H0.hpp"
#include "TestTools.hpp"

using namespace Integrator;

const double DELTA = 1.49e-8;
const double t = -1.0;
const double tp = -0.4;

struct FctTest
{
  public:
    static const size_t Nc = 4;
    static const size_t TNX = 2;
    static const size_t TNY = 2;

    FctTest(){};

    size_t n_rows() const { return Nc; };
    size_t n_cols() const { return Nc; };

    ClusterMatrixCD_t operator()(const double &x, const double &y)
    {
        ClusterMatrixCD_t tmp(Nc, Nc);
        tmp.zeros();
        for (size_t i = 0; i < Nc; i++)
        {
            for (size_t j = 0; j < Nc; j++)
            {
                tmp(i, j) = cd_t(i * x + i * i * x * x, j * y);
            }
        }
        return tmp;
    }

  private:
};

TEST(IntegratorTest, TestGridKTilde)
{

    using H0_t = Models::ABC_H0<2, 2>;
    H0_t h0(TestTools::BuildJson());
    ClusterMatrixCD_t tLocTest = GridKTilde(h0, 100);
    const ClusterMatrixCD_t goodResult = {{0.0, t, t, tp},
                                          {t, 0.0, 0.0, t},
                                          {t, 0.0, 0.0, t},
                                          {tp, t, t, 0.0}};

    for (size_t i = 0; i < goodResult.n_rows; i++)
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
    double xmin2[2] = {-1.1, 0.0};
    double xmax2[2] = {1.1, 1.0};
    FctTest fcttest;

    ClusterMatrixCD_t test2 = Cubature<FctTest>(fcttest, xmin2, xmax2, 0, 1.49e-8, 1.49e-8);

    for (size_t i = 0; i < test2.n_rows; i++)
    {
        for (size_t j = 0; j < test2.n_cols; j++)
        {
            cd_t tmpGood = cd_t(2.0 * 1.1 * 1.1 * 1.1 / 3.0 * static_cast<double>(i * i), 1.1 * static_cast<double>(j));
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
