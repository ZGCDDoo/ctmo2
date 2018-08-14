

#include <gtest/gtest.h>
#include "../src/Includes/Models/ABC_H0.hpp"
#include "TestTools.hpp"

const double DELTA = 1e-7;
const double DELTA_SMALL = 1e-11;
const size_t Nc = 4;
const size_t Nx = 2;
const size_t Ny = 2;

TEST(ModelTriangle2DTest, Init)
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

    Models::ABC_H0<Nx, Nx> h0(tJson);

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

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
