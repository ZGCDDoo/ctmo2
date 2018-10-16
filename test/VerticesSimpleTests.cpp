

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

    Utilities::EngineTypeMt19937_t rng_(1 + jj["SEED"].get<size_t>());
    Utilities::UniformRngMt19937_t urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0));

    Diagrammatic::Vertices vertices;
    Diagrammatic::VertexBuilder vertexBuilder(jj, Nc);

    //Try Inserting a shit load of vertices
    for (size_t ii = 0; ii < 1000000; ii++)
    {
        const auto v1 = vertexBuilder.BuildVertex(urng_);
        vertices.AppendVertex(v1);
    }

    vertices.Clear();

    for (size_t ii = 0; ii < 10; ii++)
    {
        const auto v1 = vertexBuilder.BuildVertex(urng_);
        vertices.AppendVertex(v1);
    }

    std::cout << "verticesKeys = " << std::endl;
    for (const auto &key : vertices.verticesKeysVec())
    {
        std::cout << key << std::endl;
    }

    vertices.Print();

    //This is a Same spin vertex
    const auto x5Up = vertices.atUp(4);
    const auto y5Down = vertices.atDown(6);
    const auto V5 = vertices.at(5);

    ASSERT_DOUBLE_EQ(x5Up.tau(), y5Down.tau());
    ASSERT_EQ(x5Up.orbital(), y5Down.orbital());
    ASSERT_EQ(x5Up.site(), y5Down.site());

    assert(x5Up == V5.vStart());
    assert(y5Down == V5.vEnd());

    ASSERT_EQ(x5Up.spin(), FermionSpin_t::Up);
    ASSERT_EQ(y5Down.spin(), FermionSpin_t::Down);

    //This is a Diff spin vertex
    const auto x4Down = vertices.atDown(4);
    const auto y4Down = vertices.atDown(5);
    const auto V4 = vertices.at(4);

    assert(x4Down == V4.vStart());
    assert(y4Down == V4.vEnd());
    ASSERT_EQ(x4Down.spin(), FermionSpin_t::Down);
    ASSERT_EQ(y4Down.spin(), FermionSpin_t::Down);

    // ASSERT()
    vertices.RemoveVertex(0);
    vertices.Print();

    const auto indices = vertices.GetIndicesSpins(4, FermionSpin_t::Down);
    ASSERT_EQ(indices.at(0), 4);
    ASSERT_EQ(indices.at(1), 5);

    vertices.RemoveVertex(4);
    vertices.Print();

    std::cout << "verticesKeys = " << std::endl;
    for (const auto &key : vertices.verticesKeysVec())
    {
        std::cout << key << std::endl;
    }

    const auto x7Up = vertices.atUp(6);
    const auto y7Down = vertices.atDown(4);
    const auto V7 = vertices.at(7);

    assert(x7Up == V7.vStart());
    assert(y7Down == V7.vEnd());
    ASSERT_EQ(x7Up.spin(), FermionSpin_t::Up);
    ASSERT_EQ(y7Down.spin(), FermionSpin_t::Down);
}

TEST(Vertices2DTest, InitVertices2)
{
    std::ifstream fin(FNAME);
    Json jj;
    fin >> jj;
    fin.close();

    Utilities::EngineTypeMt19937_t rng_(1 + jj["SEED"].get<size_t>());
    Utilities::UniformRngMt19937_t urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0));

    Diagrammatic::Vertices vertices;
    Diagrammatic::VertexBuilder vertexBuilder(jj, Nc);

    //Try Inserting a shit load of vertices
    for (size_t ii = 0; ii < 200; ii++)
    {
        const auto v1 = vertexBuilder.BuildVertex(urng_);
        vertices.AppendVertex(v1);
    }

    for (size_t ii = 0; ii < 9000000; ii++)
    {
        const auto v1 = vertexBuilder.BuildVertex(urng_);
        vertices.AppendVertex(v1);
        vertices.RemoveVertex(urng_() * vertices.size());
        // std::cout << "ii = " << ii << std::endl;
    }

    vertices.Clear();
}

TEST(Vertices2DTest, TestBuildVertex)
{
    std::ifstream fin(FNAME);
    Json jj;
    fin >> jj;
    fin.close();

    Utilities::EngineTypeMt19937_t rng_(1 + jj["SEED"].get<size_t>());
    Utilities::UniformRngMt19937_t urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0));

    Diagrammatic::Vertices vertices;
    Diagrammatic::VertexBuilder vertexBuilder(jj, Nc);

    //Try Inserting a shit load of vertices
    for (size_t ii = 0; ii < 200; ii++)
    {
        const auto v1 = vertexBuilder.BuildVertex(urng_);
        vertices.AppendVertex(v1);
    }

    for (size_t ii = 0; ii < 9000000; ii++)
    {
        const auto v1 = vertexBuilder.BuildVertex(urng_);
        vertices.AppendVertex(v1);
        vertices.RemoveVertex(urng_() * vertices.size());
        // std::cout << "ii = " << ii << std::endl;
    }

    vertices.Clear();
}

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
