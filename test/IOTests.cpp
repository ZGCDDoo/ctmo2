
#include <gtest/gtest.h>

#include "../src/Includes/Models/ModelSquare4x4.hpp"

using Model_t = Models::ModelSquare4x4;
using IOModel_t = IO::IOSquare4x4;

const double DELTA = 1e-12;
const size_t NOrb = 5;
const std::string FNAME = "../test/data/cdmft_square4x4/hyb_5Orb.dat";

TEST(IOModelTests, IndepToFullAndBack)
{
    IOModel_t ioModel;
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
    IOModel_t ioModel;
    const size_t NSS = ioModel.Nc * NOrb;
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
