
#include <gtest/gtest.h>
#include "../src/Includes/Models/ABC_H0.hpp"
#include "TestTools.hpp"

const double DELTA = 1e-7;
// const double DELTA_SMALL = 1e-11;
const size_t Nc = 4;
const size_t Nx = 2;
const size_t Ny = 2;
const size_t NOrb = 2;
const double INTRA = -0.01;

Json BuildJson()
{
    Json tJson = R"(
    {   "NOrb": 2,
        "tParameters": 
            {"00": 
                {"tIntra": 0.0, "tx": -1.0, "ty": -1.0, "tx=y": -0.40, "tx=-y": 0.0, "t2x" : 0.0, "t2y": 0.0},
            "01":
                {"tIntra": -0.01, "tx": -1.02, "ty": -1.02, "tx=y": 0.230, "tx=-y": 0.230, "t2x" : -0.90, "t2y": -0.90},
            "11":
                {"tIntra": 0.0, "tx": -1.0, "ty": -1.0, "tx=y": -0.40, "tx=-y": 0.0, "t2x" : 0.0, "t2y": 0.0}

            }
    }
    )"_json;

    return tJson;
}

TEST(ABC_H0_Tests, Init)
{

    Json tJson = R"(
    {   "NOrb": 2,
        "tParameters": 
            {"00": 
                {"tIntra": 0.0, "tx": 1.0, "ty": 1.08, "tx=y": 1.07, "tx=-y": 1.06, "t2x" : 1.05, "t2y": 1.04},
            "01":
                {"tIntra": -100.0, "tx": 1.10, "ty": 1.09, "tx=y": 1.04, "tx=-y": 1.03, "t2x" : 1.02, "t2y": 1.01},
            "11":
                {"tIntra": 0.201, "tx": -1.01, "ty": 1.01, "tx=y": 1.02, "tx=-y": 1.03, "t2x" : -1.04, "t2y": 1.09}

            }
    }
    )"_json;

    std::cout << "tJson.size() = " << tJson.size() << std::endl;
    Json jj = tJson["tParameters"];
    assert(jj.size() == 3);

    Models::ABC_H0<Nx, Ny> h0(tJson);

    std::vector<double> tIntraVec = {0.0, -100.0, 0.201};
    std::vector<double> txVec = {1.0, 1.10, -1.01};
    std::vector<double> tyVec = {1.08, 1.09, 1.01};
    std::vector<double> txyVec = {1.07, 1.04, 1.02};
    std::vector<double> tx_yVec = {1.06, 1.03, 1.03};
    std::vector<double> t2xVec = {1.05, 1.02, -1.04};
    std::vector<double> t2yVec = {1.04, 1.01, 1.09};

    for (size_t ii = 0; ii < tIntraVec.size(); ii++)
    {
        ASSERT_DOUBLE_EQ(h0.tIntraOrbitalVec().at(ii), tIntraVec.at(ii));
        ASSERT_DOUBLE_EQ(h0.txVec().at(ii), txVec.at(ii));
        ASSERT_DOUBLE_EQ(h0.tyVec().at(ii), tyVec.at(ii));
        ASSERT_DOUBLE_EQ(h0.txyVec().at(ii), txyVec.at(ii));
        ASSERT_DOUBLE_EQ(h0.tx_yVec().at(ii), tx_yVec.at(ii));
        ASSERT_DOUBLE_EQ(h0.t2xVec().at(ii), t2xVec.at(ii));
        ASSERT_DOUBLE_EQ(h0.t2yVec().at(ii), t2yVec.at(ii));
    }
}

TEST(ABC_H0_Tests, Init2)
{
    Json tJson = BuildJson();

    Models::ABC_H0<Nx, Ny> h0(tJson);

    ASSERT_EQ(h0.RSites().size(), Nc);
    ASSERT_EQ(h0.KWaveVectors().size(), Nc);
    ASSERT_EQ(h0(0.0, 0.0).n_cols, Nc * NOrb);

    ASSERT_DOUBLE_EQ(h0.KWaveVectors()[0][0], 0.0);
    ASSERT_DOUBLE_EQ(h0.KWaveVectors()[0][1], 0.0);
    ASSERT_DOUBLE_EQ(h0.KWaveVectors()[1][0], M_PI);
    ASSERT_DOUBLE_EQ(h0.KWaveVectors()[1][1], 0.0);
    ASSERT_DOUBLE_EQ(h0.KWaveVectors()[2][0], 0.0);
    ASSERT_DOUBLE_EQ(h0.KWaveVectors()[2][1], M_PI);
    ASSERT_DOUBLE_EQ(h0.KWaveVectors()[3][0], M_PI);
    ASSERT_DOUBLE_EQ(h0.KWaveVectors()[3][1], M_PI);

    ASSERT_DOUBLE_EQ(h0.RSites()[0][0], 0.0);
    ASSERT_DOUBLE_EQ(h0.RSites()[0][1], 0.0);
    ASSERT_DOUBLE_EQ(h0.RSites()[1][0], 1.0);
    ASSERT_DOUBLE_EQ(h0.RSites()[1][1], 0.0);
    ASSERT_DOUBLE_EQ(h0.RSites()[2][0], 0.0);
    ASSERT_DOUBLE_EQ(h0.RSites()[2][1], 1.0);
    ASSERT_DOUBLE_EQ(h0.RSites()[3][0], 1.0);
    ASSERT_DOUBLE_EQ(h0.RSites()[3][1], 1.0);

    ASSERT_DOUBLE_EQ(-4.6847345710802575, h0.Eps0k(-0.1, 0.3, 0));
}

TEST(ABC_H0_Tests, Hopping)
{

    Json tJson = BuildJson();

    Models::ABC_H0<Nx, Nx> h0(tJson);

    ClusterMatrixCD_t GoodHoppingKTilde(Nc * NOrb, Nc * NOrb);
    GoodHoppingKTilde.zeros();
    ClusterMatrixCD_t diagonal_part = {
        {cd_t(0.0, 0.0), cd_t(-1.98006658, 0.19866933), cd_t(-1.82533561, -0.56464247), cd_t(-0.76842440, -0.15576734)},
        {cd_t(-1.98006658, -0.19866933), cd_t(0.0, 0.0), cd_t(-0.72216088, -0.30532472), cd_t(-1.82533561, -0.56464247)},
        {cd_t(-1.82533561, 0.56464247), cd_t(-0.72216088, 0.30532472), cd_t(0.0, 0.0), cd_t(-1.98006658, 0.19866933)},
        {cd_t(-0.76842440, 0.15576734), cd_t(-1.82533561, 0.56464247), cd_t(-1.98006658, -0.19866933), cd_t(0.0, 0.0)}};

    ClusterMatrixCD_t off_diagonal = {{
                                          cd_t(-3.249723947, 0),
                                          cd_t(-2.019667909, 0.2026427174),
                                          cd_t(-1.861842327, -0.5759353229),
                                          cd_t(0.857086533, 0.1737400415),
                                      },
                                      {
                                          cd_t(-2.019667909, -0.2026427174),
                                          cd_t(-3.249723947, 0),
                                          cd_t(0.8054850475, 0.3405536159),
                                          cd_t(-1.861842327, -0.5759353229),
                                      },
                                      {
                                          cd_t(-1.861842327, 0.5759353229),
                                          cd_t(0.8054850475, -0.3405536159),
                                          cd_t(-3.249723947, 0),
                                          cd_t(-2.019667909, 0.2026427174),
                                      },
                                      {
                                          cd_t(0.857086533, -0.1737400415),
                                          cd_t(-1.861842327, 0.5759353229),
                                          cd_t(-2.019667909, -0.2026427174),
                                          cd_t(-3.249723947, 0),
                                      }};

    GoodHoppingKTilde.submat(0, 0, Nc - 1, Nc - 1) = diagonal_part;
    GoodHoppingKTilde.submat(Nc, Nc, NOrb * Nc - 1, NOrb * Nc - 1) = diagonal_part;
    GoodHoppingKTilde.submat(0, Nc, Nc - 1, NOrb * Nc - 1) = off_diagonal + INTRA * ClusterMatrixCD_t(Nc, Nc).eye();
    GoodHoppingKTilde.submat(Nc, 0, NOrb * Nc - 1, Nc - 1) = off_diagonal + INTRA * ClusterMatrixCD_t(Nc, Nc).eye();

    ClusterMatrixCD_t hopping = h0(0.1, -0.3);
    for (size_t i = 0; i < Nc * NOrb; i++)
    {
        for (size_t j = 0; j < Nc * NOrb; j++)
        {
            //std::cout << "i ,j " << i << " " << j << std::endl;
            ASSERT_NEAR(hopping(i, j).real(), GoodHoppingKTilde(i, j).real(), DELTA);
            ASSERT_NEAR(hopping(i, j).imag(), GoodHoppingKTilde(i, j).imag(), DELTA);
        }
    }
}

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
