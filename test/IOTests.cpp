
#include <gtest/gtest.h>
#include "../src/Includes/Utilities/IO.hpp"

const double DELTA = 1e-12;
const size_t NOrb = 5;
const std::string FNAME = "../test/data/cdmft_square4x4/hyb_5Orb.dat";

Json BuildJson()
{
    Json tJson = R"(
    {   "Nx": 4,
        "Ny": 4,
        "Nz": 1,
        "ModelFile": "../data/Square4x4.model",
        "NOrb": 1,
        "NKPTS": 100,
        "tParameters": 
            {"00": 
                {"tIntra": 0.0, "tx": -1.120, "ty": -1.120, "tz": 0.0, "tx=y": 0.1560, "tx=-y": 0.156, "tx=z": 0.0, "tx=-z": 0.0, "ty=z": 0.0, "ty=-z": 0.0, "t2x" : 0.23700, "t2y": 0.23700, "t2z": 0.0, "t3": 0.0}
            }
    }
    )"_json;

    return tJson;
}

TEST(IOModelTests, IndepToFullAndBack)
{

    IO::Base_IOModel ioModel(BuildJson());
    const ClusterCubeCD_t greenCube = ioModel.ReadGreenDat(FNAME, 5);
    const ClusterMatrixCD_t greenIndep = ioModel.FullCubeToIndep(greenCube);

    for (size_t nn = 0; nn < greenCube.n_slices; nn++)
    {
        const ClusterMatrixCD_t matrixFull = ioModel.IndepToFull(greenIndep.row(nn), NOrb);
        assert(matrixFull.n_rows == greenCube.n_rows);
        assert(matrixFull.n_cols == greenCube.n_cols);
        std::cout << "matrixFull.n_rows = " << matrixFull.n_rows << " ioModel.GetNIndepSuperSites(NOrb) = " << ioModel.GetNIndepSuperSites(NOrb) << std::endl;

        assert(matrixFull.n_rows == NOrb * ioModel.Nc);

        for (size_t ii = 0; ii < greenCube.n_rows; ii++)
        {
            for (size_t jj = 0; jj < greenCube.n_rows; jj++)
            {
                ASSERT_DOUBLE_EQ(matrixFull(ii, jj).real(), greenCube(ii, jj, nn).real());
                ASSERT_DOUBLE_EQ(matrixFull(ii, jj).imag(), greenCube(ii, jj, nn).imag());
            }
        }
    }
}

TEST(IOModelTests, ReadAndWrite)
{
    IO::Base_IOModel ioModel(BuildJson());
    const ClusterCubeCD_t greenCube = ioModel.ReadGreenDat(FNAME, 5);
    ioModel.SaveCube("tmp_test", greenCube, 10.0, NOrb, 15);
    const ClusterCubeCD_t readGreenCube = ioModel.ReadGreenDat("tmp_test.dat", 5);

    for (size_t nn = 0; nn < greenCube.n_slices; nn++)
    {

        for (size_t ii = 0; ii < greenCube.n_rows; ii++)
        {
            for (size_t jj = 0; jj < greenCube.n_rows; jj++)
            {
                ASSERT_NEAR(readGreenCube(ii, jj, nn).real(), greenCube(ii, jj, nn).real(), DELTA);
                ASSERT_NEAR(readGreenCube(ii, jj, nn).imag(), greenCube(ii, jj, nn).imag(), DELTA);
            }
        }
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
