
#include <gtest/gtest.h>
#include "../src/Includes/Models/ABC_Model.hpp"
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
    Json jj = R"(
{
    "monteCarlo": {
        "measurementTime": 5.0,
        "thermalizationTime": 0.5,
        "thermFromConfig": false,
        "seed": 1096954689
    },
    "model": {
        "beta": 33.3,
        "hybUpFile": "../test/data/DMFT/hybfm_SIAM_Square.dat",
        "modelFile": "../data/SIAM_Square.model",
        "J_H": 0.0,
        "nkpts": 60,
        "nOrb": 1,
        "U": 0.0,
        "UPrime": 0.0,
        "delta": 0.01,
        "gPhonon": 0.08,
        "w0Phonon": 0.08,
        "mu": 10.10,
        "cluster": {
            "Nx": 1,
            "Ny": 1,
            "Nz": 1
        },
        "tParameters": {
            "00": {
                "t2x": 0.0,
                "t2y": 0.0,
                "t2z": 0.0,
                "tIntra": 0.0,
                "tx": -0.3926,
                "tx=-y": 0.0,
                "tx=-z": 0.0,
                "tx=y": 0.0,
                "tx=z": 0.0,
                "ty": -0.3926,
                "ty=-z": 0.0,
                "ty=z": 0.0,
                "tz": -0.3926,
                "t3": 0.0
            }
        }
    },
    "solver": {
        "cleanUpdate": 100,
        "updatesMeas": 50,
        "averageOrbitals": true,
        "isOrbitalDiagonal": true,
        "eCutGreen": 100,
        "ntau": 5000,
        "n_tau_sampling": 5
    },
    "selfCon": {
        "eCutSelfCon": 200,
        "weightsR": 0.2,
        "weightsI": 0.2
    },
    "logging": {
        "level": "debug"
    }
}
    )"_json;

    return jj;
}

TEST(ABC_Model_Tests, Init)
{

    Models::ABC_Model_2D model(BuildJson());

    ASSERT_DOUBLE_EQ(model.mu(), 10.1);
    ASSERT_DOUBLE_EQ(model.beta(), 33.3);
    ASSERT_EQ(model.NOrb(), 1);
    ASSERT_EQ(model.Nc(), 1);
}

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
