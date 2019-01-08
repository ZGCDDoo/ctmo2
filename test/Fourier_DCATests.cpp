

#include <gtest/gtest.h>

#include "../src/Includes/Utilities/Fourier_DCA.hpp"
#include "../src/Includes/Models/ABC_H0.hpp"

using namespace FourierDCA;

Json BuildJson()
{
    Json tJson = R"(
    {
        "model": {
            "cluster": {
                "Nx": 3,
                "Ny": 3,
                "Nz": 1
            },
            "nOrb": 2,
            "nkpts": 100,
            "tParameters": {
                "00": {
                    "tIntra": 0.0,
                    "tx": -1.0,
                    "ty": -1.0,
                    "tz": 0.0,
                    "tx=y": -0.40,
                    "tx=-y": 0.0,
                    "tx=z": -0.0,
                    "tx=-z": 0.0,
                    "ty=z": 0.0,
                    "ty=-z": 0.0,
                    "t2x": 0.0,
                    "t2y": 0.0,
                    "t2z": 0.0,
                    "t3": 0.0
                },
                "01": {
                    "tIntra": -0.01,
                    "tx": -1.02,
                    "ty": -1.02,
                    "tz": 0.0,
                    "tx=y": 0.230,
                    "tx=-y": 0.230,
                    "tx=z": -0.0,
                    "tx=-z": 0.0,
                    "ty=z": 0.0,
                    "ty=-z": 0.0,
                    "t2x": -0.90,
                    "t2y": -0.90,
                    "t2z": 0.0,
                    "t3": 0.0
                },
                "11": {
                    "tIntra": 0.0,
                    "tx": -1.0,
                    "ty": -1.0,
                    "tz": 0.0,
                    "tx=y": -0.40,
                    "tx=-y": 0.0,
                    "tx=z": -0.0,
                    "tx=-z": 0.0,
                    "ty=z": 0.0,
                    "ty=-z": 0.0,
                    "t2x": 0.0,
                    "t2y": 0.0,
                    "t2z": 0.0,
                    "t3": 0.0
                }
            }
        }
    }
    )"_json;

    return tJson;
}

const double DELTA = 1e-11;
TEST(FourierDCATests, KToR)
{
    arma::arma_rng::set_seed_random();
    const size_t Nx = 3;
    const size_t NMat = 1;
    const Models::ABC_H0 h0(BuildJson());

    DataK_t greenK(Nx * Nx, Nx * Nx, NMat);
    greenK.zeros();
    for (size_t nn = 0; nn < NMat; nn++)
    {
        SiteVectorCD_t tmp(Nx * Nx);
        tmp.randu();
        // tmp(1) = tmp(2);
        greenK.slice(nn).diag() = tmp;
        greenK.slice(nn).print();
        std::cout << "\n\n";
    }

    ClusterCubeCD_t greenRTest = KtoR(greenK, h0.RSites(), h0.KWaveVectors());
    DataK_t greenKTest = RtoK(greenRTest, h0.RSites(), h0.KWaveVectors());

    assert(greenKTest.n_cols == Nx * Nx);
    assert(greenKTest.n_rows == Nx * Nx);
    assert(greenK.n_rows == Nx * Nx);
    assert(greenKTest.n_cols == Nx * Nx);

    std::cout << "\n\n"
              << " greenR " << std::endl;
    greenRTest.print();

    std::cout << "\n\n"
              << " greenKTest " << std::endl;
    greenKTest.print();

    std::cout << "\n\n"
              << " greenK " << std::endl;
    greenK.print();

    std::cout << "Here " << std::endl;
    for (size_t nn = 0; nn < NMat; nn++)
    {
        for (size_t ii = 0; ii < greenKTest.n_rows; ii++)
        {
            for (size_t jj = 0; jj < greenKTest.n_cols; jj++)
            {
                // std::cout << "ii, jj = " << ii << ", " << jj << std::endl;
                ASSERT_NEAR(greenKTest(ii, jj, nn).real(), greenK(ii, jj, nn).real(), DELTA);
                ASSERT_NEAR(greenKTest(ii, jj, nn).imag(), greenK(ii, jj, nn).imag(), DELTA);
            }
        }
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
