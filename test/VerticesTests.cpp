

#include <gtest/gtest.h>
#include "../src/Includes/Utilities/Vertices.hpp"
#include "TestTools.hpp"

const double DELTA = 1e-7;
const double DELTA_SMALL = 1e-11;
const double delta = 0.01;
const double U = 3.0;
const double Beta = 10.0;
const double mu = 1.8941850792671628;
const size_t Nc = 4;
const std::string fname = "../test/data/cdmft_triangle/testtriangle.json";

TEST(Vertices2DTest, Init)
{
    Vertices::AuxHelper auxHelper(delta);
    // Vertices::VertexPart vertexPart();
    // std::ifstream fin(fname);
    // Json jj;
    // fin >> jj;
    // fin.close();
    // Model_t modelTriangle(jj);

    // ASSERT_NEAR(modelTriangle.U(), U, DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.mu(), mu, DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.auxMu(), mu - U / 2.0, DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.delta(), delta, DELTA_SMALL);

    // Vertex vertex1(1.0, 1, AuxSpin_t::Up);
    // ASSERT_NEAR(modelTriangle.auxUp(vertex1.aux()), 1.0 + delta, DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.auxDown(vertex1.aux()), -delta, DELTA_SMALL);

    // Vertex vertex2(1.10, 0, AuxSpin_t::Down);
    // Vertex vertexZero(1.0, 1, AuxSpin_t::Zero);
    // ASSERT_NEAR(modelTriangle.auxUp(vertex2.aux()), -delta, DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.auxDown(vertex2.aux()), 1.0 + delta, DELTA_SMALL);

    // ASSERT_NEAR(modelTriangle.FAuxUp(vertex1.aux()), (1.0 + delta) / delta, DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.FAuxUp(vertex2.aux()), delta / (1.0 + delta), DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.FAuxUp(vertexZero.aux()), 1.0, DELTA);

    // ASSERT_NEAR(modelTriangle.FAuxDown(vertex2.aux()), (1.0 + delta) / delta, DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.FAuxDown(vertex1.aux()), delta / (1.0 + delta), DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.FAuxDown(vertexZero.aux()), 1.0, DELTA_SMALL);

    // const double KAux = -U * Beta * Nc / (((1.0 + delta) / delta - 1.0) * (delta / (1.0 + delta) - 1.0));
    // ASSERT_NEAR(modelTriangle.KAux(), KAux, DELTA_SMALL);

    // ASSERT_NEAR(modelTriangle.gammaUp(vertex1.aux(), vertex2.aux()), (modelTriangle.FAuxUp(vertex1.aux()) - modelTriangle.FAuxUp(vertex2.aux())) / modelTriangle.FAuxUp(vertex2.aux()), DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.gammaUp(vertex1.aux(), vertexZero.aux()), (modelTriangle.FAuxUp(vertex1.aux()) - modelTriangle.FAuxUp(vertexZero.aux())) / modelTriangle.FAuxUp(vertexZero.aux()), DELTA_SMALL);

    // ASSERT_NEAR(modelTriangle.gammaDown(vertex2.aux(), vertex1.aux()), (modelTriangle.FAuxDown(vertex2.aux()) - modelTriangle.FAuxDown(vertex1.aux())) / modelTriangle.FAuxDown(vertex1.aux()), DELTA_SMALL);
    // ASSERT_NEAR(modelTriangle.gammaDown(vertexZero.aux(), vertex1.aux()), (modelTriangle.FAuxDown(vertexZero.aux()) - modelTriangle.FAuxDown(vertex1.aux())) / modelTriangle.FAuxDown(vertex1.aux()), DELTA_SMALL);
}

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
