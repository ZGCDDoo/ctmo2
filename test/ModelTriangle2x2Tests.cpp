

#include <gtest/gtest.h>
#include "../src/Includes/Models/ModelTriangle2x2.hpp"
#include "TestTools.hpp"

const double DELTA = 1e-7;
const size_t Nc = 4;
const std::string fname = "../test/data/cdmft_triangle/testtriangle.json";

using Model_t = Models::ModelTriangle2x2;

TEST(ModelTriangle2DTest, Init)
{

    std::ifstream fin(fname);
    Json jj;
    fin >> jj;
    fin.close();
    Model_t modelTriangle(jj);

}

//Test that the model gives the correct thing for a simulation done by CT-Hyb patrick for the parameters given in testtriangle.json
TEST(ModelTriangle2DTest, SimpleTriangle)
{
    std::ifstream fin(fname);
    Json jj;
    fin >> jj;
    fin.close();
    Model_t modelTriangle(jj);

    //Let us compare the values for the greenMat0_
    //std::cout << " modelTriangle.green0Mat().slice(-1) = \n"
    //        << modelTriangle.green0Mat().slice(0) << std::endl;

    ClusterMatrixCD_t greenMatTest0 = modelTriangle.greenCluster0MatUp().data().slice(0);
    ClusterMatrixCD_t greenMatTest1 = modelTriangle.greenCluster0MatUp().data().slice(1);

    ClusterMatrixCD_t goodGreenMat0 = {
        {cd_t(-0.12492550, -0.55707086), cd_t(0.19616562, 0.01962542), cd_t(0.19616562, 0.01962542), cd_t(0.13204398, 0.23970799)},
        {cd_t(0.19616562, 0.01962542), cd_t(-0.12548027, -0.56134078), cd_t(-0.00277436, 0.25388105), cd_t(0.19616562, 0.01962542)},
        {cd_t(0.19616562, 0.01962542), cd_t(-0.00277436, 0.25388105), cd_t(-0.12548027, -0.56134078), cd_t(0.19616562, 0.01962542)},
        {cd_t(0.13204398, 0.23970799), cd_t(0.19616562, 0.01962542), cd_t(0.19616562, 0.01962542), cd_t(-0.12492550, -0.55707086)}};

    cd_t goodGreenMat100(-0.03771745, -0.44068635);
    ASSERT_NEAR(goodGreenMat100.real(), greenMatTest1(0, 0).real(), DELTA);
    ASSERT_NEAR(goodGreenMat100.imag(), greenMatTest1(0, 0).imag(), DELTA);

    for (size_t i = 0; i < Nc; i++)
        for (size_t j = 0; j < Nc; j++)
        {
            std::cout << "i ,j " << i << " " << j << std::endl;
            ASSERT_NEAR(goodGreenMat0(i, j).real(), greenMatTest0(i, j).real(), DELTA);
            ASSERT_NEAR(goodGreenMat0(i, j).imag(), greenMatTest0(i, j).imag(), DELTA);
        }
}

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
