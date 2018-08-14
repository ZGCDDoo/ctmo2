

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
                {"tIntra": 0.0, "tx": 1.0, "ty": 1.0, "tx=y": 1.0, "tx=-y": 1.0, "t2x" : 1.0, "t2y": 1.0},
            "01":
                {"tIntra": 0.0, "tx": 1.0, "ty": 1.0, "tx=y": 1.0, "tx=-y": 1.0, "t2x" : 1.0, "t2y": 1.0},
            "11":
                {"tIntra": 0.0, "tx": 1.0, "ty": 1.0, "tx=y": 1.0, "tx=-y": 1.0, "t2x" : 1.0, "t2y": 1.0}

            }
    }
    )"_json;

    std::cout << "tJson.size() = " << tJson.size() << std::endl;
    Json jj = tJson["tParameters"];
    assert(jj.size() == 3);

    Models::ABC_H0<2, 2> h0_(tJson);
}

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
