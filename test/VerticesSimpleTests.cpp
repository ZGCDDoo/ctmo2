

#include <gtest/gtest.h>
#include "../src/Includes/Utilities/VerticesSimple.hpp"
#include "TestTools.hpp"

// const double DELTA = 1e-7;
const double DELTA_SMALL = 1e-11;
const double delta = 0.01;
// const double U = 3.0;
// const double Beta = 10.0;
// const double mu = 1.8941850792671628;
const size_t Nc = 4;
const std::string FNAME = "../test/data/cdmft_square2x2/params1.json";

TEST(Vertices2DTest, AuxHelper)
{
    Diagrammatic::AuxHelper auxHelper(delta);

    ASSERT_NEAR(auxHelper.delta(), delta, DELTA_SMALL);

    ASSERT_NEAR(auxHelper.auxValue(FermionSpin_t::Up, AuxSpin_t::Up), 1.0 + delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.auxValue(FermionSpin_t::Down, AuxSpin_t::Down), 1.0 + delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.auxValue(FermionSpin_t::Up, AuxSpin_t::Down), -delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.auxValue(FermionSpin_t::Down, AuxSpin_t::Up), -delta, DELTA_SMALL);

    ASSERT_NEAR(auxHelper.FAux(FermionSpin_t::Up, AuxSpin_t::Up), (1.0 + delta) / delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.FAux(FermionSpin_t::Down, AuxSpin_t::Down), (1.0 + delta) / delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.FAux(FermionSpin_t::Up, AuxSpin_t::Down), (delta) / (1.0 + delta), DELTA_SMALL);
    ASSERT_NEAR(auxHelper.FAux(FermionSpin_t::Down, AuxSpin_t::Up), (delta) / (1.0 + delta), DELTA_SMALL);

    ASSERT_NEAR(auxHelper.gamma(FermionSpin_t::Up, AuxSpin_t::Up, AuxSpin_t::Up), 0.0, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.gamma(FermionSpin_t::Up, AuxSpin_t::Down, AuxSpin_t::Down), 0.0, DELTA_SMALL);

    ASSERT_NEAR(auxHelper.gamma(FermionSpin_t::Up, AuxSpin_t::Up, AuxSpin_t::Down), (auxHelper.FAux(FermionSpin_t::Up, AuxSpin_t::Up) - auxHelper.FAux(FermionSpin_t::Up, AuxSpin_t::Down)) / (auxHelper.FAux(FermionSpin_t::Up, AuxSpin_t::Down)), DELTA_SMALL);
    ASSERT_NEAR(auxHelper.gamma(FermionSpin_t::Up, AuxSpin_t::Down, AuxSpin_t::Up), (auxHelper.FAux(FermionSpin_t::Up, AuxSpin_t::Down) - auxHelper.FAux(FermionSpin_t::Up, AuxSpin_t::Up)) / (auxHelper.FAux(FermionSpin_t::Up, AuxSpin_t::Up)), DELTA_SMALL);
}

TEST(Vertices2DTest, InitVertices)
{
    std::ifstream fin(FNAME);
    Json jj;
    fin >> jj;
    fin.close();

    Utilities::EngineTypeMt19937_t rng_(jj["SEED"].get<size_t>());
    Utilities::UniformRngMt19937_t urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0));

    Diagrammatic::Vertices vertices;
    Diagrammatic::VertexBuilder vertexBuilder(jj, Nc);

    const auto v1 = vertexBuilder.BuildVertex(urng_);
    vertices.AppendVertex(v1);
}

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
