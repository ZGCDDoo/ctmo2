

#include <gtest/gtest.h>
#include "../src/Includes/Models/ABC_H0.hpp"

const size_t Nx = 2;
const size_t Nc = 4;
const size_t NOrb = 2;
const double DELTA = 1e-7;

Json BuildJson()
{
    Json tJson = R"(
    {   "NOrb": 2,
        "tParameters": 
            {"00": 
                {"tIntra": 0.0, "tx": -1.0, "ty": -1.0, "tx=y": -0.40, "tx=-y": 0.0, "t2x" : 0.0, "t2y": 0.0},
            "01":
                {"tIntra": 0.0, "tx": 0.0, "ty": 0.0, "tx=y": 0.0, "tx=-y": 0.0, "t2x" : 0.0, "t2y": 0.0},
            "11":
                {"tIntra": 0.0, "tx": -1.0, "ty": -1.0, "tx=y": -0.40, "tx=-y": 0.0, "t2x" : 0.0, "t2y": 0.0}

            }
    }
    )"_json;

    return tJson;
}

TEST(H0Triangle2DTest, Init)
{
    Json tJson = BuildJson();

    Models::ABC_H0<Nx, Nx> h0(tJson);

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

TEST(H0Triangle2DTest, Hopping)
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

    GoodHoppingKTilde.submat(0, 0, Nc - 1, Nc - 1) = diagonal_part;

    GoodHoppingKTilde.submat(Nc, Nc, NOrb * Nc - 1, NOrb * Nc - 1) = diagonal_part;

    ClusterMatrixCD_t hopping = h0(0.1, -0.3);
    for (size_t i = 0; i < Nc; i++)
    {
        for (size_t j = 0; j < Nc; j++)
        {
            //std::cout << "i ,j " << i << " " << j << std::endl;
            ASSERT_NEAR(hopping(i, j).real(), GoodHoppingKTilde(i, j).real(), DELTA);
            ASSERT_NEAR(hopping(i, j).imag(), GoodHoppingKTilde(i, j).imag(), DELTA);
        }
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
