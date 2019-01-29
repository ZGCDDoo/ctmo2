

#include <gtest/gtest.h>
#include "../src/Includes/Utilities/VerticesSimple.hpp"
#include "TestTools.hpp"

const double DELTA_SMALL = 1e-11;
const double delta = 0.01;
const std::string FNAME = "../test/data/cdmft_square2x2/params1.json";

TEST(VerticesTests, AuxHelper)
{
    Diagrammatic::AuxHelper auxHelper(delta);

    ASSERT_NEAR(auxHelper.delta(), delta, DELTA_SMALL);

    ASSERT_NEAR(auxHelper.auxValue(FermionSpin_t::Up, AuxSpin_t::Up), 1.0 + delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.auxValue(FermionSpin_t::Down, AuxSpin_t::Down), 1.0 + delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.auxValue(FermionSpin_t::Up, AuxSpin_t::Down), -delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.auxValue(FermionSpin_t::Down, AuxSpin_t::Up), -delta, DELTA_SMALL);

    using Diagrammatic::VertexPart;
    const Diagrammatic::VertexType vtype = Diagrammatic::VertexType::HubbardInter;

    const auto vp_upup = VertexPart(vtype, 1.0, 1, FermionSpin_t::Up, 1, AuxSpin_t::Up);
    const auto vp_downdown = VertexPart(vtype, 1.0, 1, FermionSpin_t::Down, 1, AuxSpin_t::Down);
    const auto vp_updown = VertexPart(vtype, 1.0, 1, FermionSpin_t::Up, 1, AuxSpin_t::Down);
    const auto vp_downup = VertexPart(vtype, 1.0, 1, FermionSpin_t::Down, 1, AuxSpin_t::Up);

    ASSERT_NEAR(auxHelper.FAux(vp_upup), (1.0 + delta) / delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.FAux(vp_downdown), (1.0 + delta) / delta, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.FAux(vp_updown), (delta) / (1.0 + delta), DELTA_SMALL);
    ASSERT_NEAR(auxHelper.FAux(vp_downup), (delta) / (1.0 + delta), DELTA_SMALL);

    ASSERT_NEAR(auxHelper.gamma(vp_upup, vp_upup), 0.0, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.gamma(vp_upup, vp_updown), (auxHelper.FAux(vp_upup) - auxHelper.FAux(vp_updown)) / (auxHelper.FAux(vp_updown)), DELTA_SMALL);
    ASSERT_NEAR(auxHelper.gamma(vp_updown, vp_upup), (auxHelper.FAux(vp_updown) - auxHelper.FAux(vp_upup)) / (auxHelper.FAux(vp_upup)), DELTA_SMALL);

    ASSERT_NEAR(auxHelper.gamma(vp_downdown, vp_downdown), 0.0, DELTA_SMALL);
    ASSERT_NEAR(auxHelper.gamma(vp_downup, vp_downdown), (auxHelper.FAux(vp_downup) - auxHelper.FAux(vp_downdown)) / (auxHelper.FAux(vp_downdown)), DELTA_SMALL);
}

TEST(VerticesTests, VertexPartInit)
{
    using namespace Diagrammatic;
    VertexPart vp;
    VertexPart vp0(VertexType::Invalid, 0.1, 1, FermionSpin_t::Down, 2, AuxSpin_t::Down);

    ASSERT_EQ(vp0.vtype(), VertexType::Invalid);
    ASSERT_DOUBLE_EQ(vp0.tau(), 0.1);
    ASSERT_EQ(vp0.site(), 1);
    ASSERT_EQ(vp0.spin(), FermionSpin_t::Down);
    ASSERT_EQ(vp0.orbital(), 2);
    ASSERT_TRUE((vp0.superSite() == std::make_pair<size_t, size_t>(1, 2)));
    ASSERT_EQ(vp0.aux(), AuxSpin_t::Down);

    //Copy constructor
    VertexPart vp1(vp0);
    ASSERT_EQ(vp1.vtype(), VertexType::Invalid);
    ASSERT_DOUBLE_EQ(vp1.tau(), 0.1);
    ASSERT_EQ(vp1.site(), 1);
    ASSERT_EQ(vp1.spin(), FermionSpin_t::Down);
    ASSERT_EQ(vp1.orbital(), 2);
    ASSERT_TRUE((vp1.superSite() == std::make_pair<size_t, size_t>(1, 2)));
    ASSERT_EQ(vp1.aux(), AuxSpin_t::Down);

    //assignement operator
    VertexPart vp2(VertexType::Phonon, 0.22, 100, FermionSpin_t::Up, 9, AuxSpin_t::Up);
    vp0 = vp2;

    ASSERT_EQ(vp0.vtype(), VertexType::Phonon);
    ASSERT_DOUBLE_EQ(vp0.tau(), 0.22);
    ASSERT_EQ(vp0.site(), 100);
    ASSERT_EQ(vp0.spin(), FermionSpin_t::Up);
    ASSERT_EQ(vp0.orbital(), 9);
    ASSERT_TRUE((vp0.superSite() == std::make_pair<size_t, size_t>(100, 9)));
    ASSERT_EQ(vp0.aux(), AuxSpin_t::Up);

    ASSERT_TRUE((vp0 == vp2));

    //FlipAux
    vp2.FlipAux();
    ASSERT_EQ(vp2.aux(), AuxSpin_t::Down);

    //SetAux
    vp2.SetAux(AuxSpin_t::Zero);
    ASSERT_EQ(vp2.aux(), AuxSpin_t::Zero);

    //SetSpin
    vp2.SetSpin(FermionSpin_t::Down);
    ASSERT_EQ(vp2.spin(), FermionSpin_t::Down);

    //now try with a vector
    std::vector<VertexPart> vecVp;
    vecVp.push_back(vp2);
    vecVp.at(0) = vp0;

    ASSERT_EQ(vp0.vtype(), VertexType::Phonon);
    ASSERT_DOUBLE_EQ(vp0.tau(), 0.22);
    ASSERT_EQ(vp0.site(), 100);
    ASSERT_EQ(vp0.spin(), FermionSpin_t::Up);
    ASSERT_EQ(vp0.orbital(), 9);
    ASSERT_TRUE((vp0.superSite() == std::make_pair<size_t, size_t>(100, 9)));
    ASSERT_EQ(vp0.aux(), AuxSpin_t::Up);
}

TEST(VerticesTests, VertexInit)
{
    using namespace Diagrammatic;
    Vertex v_null;

    VertexPart vStart(VertexType::HubbardIntra, 0.1, 1, FermionSpin_t::Up, 2, AuxSpin_t::Down);
    VertexPart vEnd(VertexType::HubbardIntra, 0.1, 1, FermionSpin_t::Down, 2, AuxSpin_t::Down);

    Vertex v0(vStart.vtype(), vStart, vEnd, 1.11);

    ASSERT_EQ(v0.vtype(), VertexType::HubbardIntra);
    ASSERT_DOUBLE_EQ(v0.vStart().tau(), 0.1);
    ASSERT_EQ(v0.vStart().site(), 1);
    ASSERT_EQ(v0.vStart().spin(), FermionSpin_t::Up);
    ASSERT_EQ(v0.vEnd().spin(), FermionSpin_t::Down);
    ASSERT_EQ(v0.aux(), AuxSpin_t::Down);

    v0.SetAux(AuxSpin_t::Up);
    ASSERT_EQ(v0.aux(), AuxSpin_t::Up);
    ASSERT_EQ(v0.vStart().aux(), AuxSpin_t::Up);
    ASSERT_EQ(v0.vEnd().aux(), AuxSpin_t::Up);

    v0.FlipAux();
    ASSERT_EQ(v0.aux(), AuxSpin_t::Down);
    ASSERT_EQ(v0.vStart().aux(), AuxSpin_t::Down);
    ASSERT_EQ(v0.vEnd().aux(), AuxSpin_t::Down);

    //test copy constructor;
    Vertex v1(v0);
    ASSERT_TRUE((v1.vStart() == v0.vStart()));
    ASSERT_TRUE((v1.vEnd() == v0.vEnd()));
    ASSERT_EQ(v1.aux(), v0.aux());
    ASSERT_FALSE((v1.vStart() == v_null.vStart()));

    //test assigment operator

    const VertexPart vStart2(VertexType::HubbardIntra, 0.33, 10, FermionSpin_t::Up, 3, AuxSpin_t::Up);
    const VertexPart vEnd2(VertexType::HubbardIntra, 0.33, 10, FermionSpin_t::Down, 3, AuxSpin_t::Up);

    Vertex v2(vStart.vtype(), vStart2, vEnd2, 2.31);
    v1 = v2;

    ASSERT_DOUBLE_EQ(vStart2.tau(), 0.33);
    ASSERT_TRUE((v1.vStart() == v2.vStart()));
    ASSERT_TRUE((v1.vEnd() == v2.vEnd()));
    ASSERT_EQ(v1.aux(), v2.aux());
    ASSERT_FALSE((v1.vStart() == v_null.vStart()));
    ASSERT_TRUE(v1.vEnd() == vEnd2);

    //test assignement with vector

    std::vector<Vertex> vecVertex;
    vecVertex.push_back(v_null);
    vecVertex.at(0) = v2;

    ASSERT_TRUE((vecVertex.at(0).vStart() == v2.vStart()));
    ASSERT_TRUE((vecVertex.at(0).vEnd() == v2.vEnd()));
    ASSERT_EQ(vecVertex.at(0).aux(), v2.aux());
    ASSERT_FALSE((vecVertex.at(0).vStart() == v_null.vStart()));

    ASSERT_DOUBLE_EQ(vStart2.tau(), vecVertex.at(0).vStart().tau());
    ASSERT_EQ(vecVertex.at(0).aux(), AuxSpin_t::Up);
    ASSERT_EQ(vecVertex.at(0).vStart().aux(), AuxSpin_t::Up);

    vecVertex.push_back(v1);
    v2.SetAux(AuxSpin_t::Zero);
    vecVertex.at(1) = v2;

    ASSERT_DOUBLE_EQ(vecVertex.at(1).vStart().tau(), 0.33);
    ASSERT_EQ(vecVertex.at(1).vStart().aux(), AuxSpin_t::Zero);
    ASSERT_TRUE((vecVertex.at(1).vStart() == v2.vStart()));
}

// TEST(Vertices2DTest, InitVertices)
// {
//     std::ifstream fin(FNAME);
//     Json jj;
//     fin >> jj;
//     fin.close();

//     Utilities::EngineTypeMt19937_t rng_(1 + jj["monteCarlo"]["seed"].get<size_t>());
//     Utilities::UniformRngMt19937_t urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0));

//     Diagrammatic::Vertices vertices;
//     Diagrammatic::VertexBuilder vertexBuilder(jj, Nc);

//     //Try Inserting a shit load of vertices
//     for (size_t ii = 0; ii < 1000000; ii++)
//     {
//         const auto v1 = vertexBuilder.BuildVertex(urng_);
//         vertices.AppendVertex(v1);
//     }

//     vertices.Clear();

//     for (size_t ii = 0; ii < 10; ii++)
//     {
//         const auto v1 = vertexBuilder.BuildVertex(urng_);
//         vertices.AppendVertex(v1);
//     }

//     std::cout << "verticesKeys = " << std::endl;
//     for (const auto &key : vertices.verticesKeysVec())
//     {
//         std::cout << key << std::endl;
//     }

//     vertices.Print();

//     //This is a Same spin vertex
//     const auto x5Up = vertices.atUp(4);
//     const auto y5Down = vertices.atDown(6);
//     const auto V5 = vertices.at(5);

//     ASSERT_DOUBLE_EQ(x5Up.tau(), y5Down.tau());
//     ASSERT_EQ(x5Up.orbital(), y5Down.orbital());
//     ASSERT_EQ(x5Up.site(), y5Down.site());

//     assert(x5Up == V5.vStart());
//     assert(y5Down == V5.vEnd());

//     ASSERT_EQ(x5Up.spin(), FermionSpin_t::Up);
//     ASSERT_EQ(y5Down.spin(), FermionSpin_t::Down);

//     //This is a Diff spin vertex
//     const auto x4Down = vertices.atDown(4);
//     const auto y4Down = vertices.atDown(5);
//     const auto V4 = vertices.at(4);

//     assert(x4Down == V4.vStart());
//     assert(y4Down == V4.vEnd());
//     ASSERT_EQ(x4Down.spin(), FermionSpin_t::Down);
//     ASSERT_EQ(y4Down.spin(), FermionSpin_t::Down);

//     // Remove vertex 0
//     const size_t vertexKey0 = vertices.GetKey(0);
//     assert(vertexKey0 == 0);
//     const auto x0 = vertices.at(0).vStart();
//     const auto y0 = vertices.at(0).vEnd();
//     const size_t pp0Up = vertices.GetKeyIndex(vertexKey0, x0.spin());
//     assert(pp0Up == 0);
//     const size_t pp0Down = vertices.GetKeyIndex(vertexKey0, y0.spin());
//     assert(pp0Down == 0);

//     const size_t kkUp = vertices.NUp();
//     const size_t kkUpm1 = kkUp - 1;
//     const size_t kkDown = vertices.NDown();
//     const size_t kkDownm1 = kkDown - 1;
//     assert(x0.spin() == FermionSpin_t::Up);
//     assert(y0.spin() == FermionSpin_t::Down);
//     vertices.SwapVertexPart(pp0Up, kkUpm1, x0.spin());
//     vertices.SwapVertexPart(pp0Down, kkDownm1, y0.spin());
//     vertices.PopBackVertexPart(x0.spin());
//     vertices.PopBackVertexPart(y0.spin());
//     vertices.RemoveVertex(0);

//     std::cout << "After remove 0 \n"
//               << std::endl;
//     vertices.Print();
//     std::cout << "verticesKeys = " << std::endl;
//     for (const auto &key : vertices.verticesKeysVec())
//     {
//         std::cout << key << std::endl;
//     }

//     const size_t vertexKey4New = vertices.GetKey(4);
//     assert(vertexKey4New == 12);

//     const auto x4New = vertices.at(4).vStart();
//     const auto y4New = vertices.at(4).vEnd();
//     assert(x4New.spin() == y4New.spin());
//     const size_t i2 = vertices.GetKeyIndex(vertexKey4New + 1, x4New.spin());
//     assert(i2 == 5);
//     vertices.SwapVertexPart(i2, vertices.NDown() - 1, x4New.spin());

//     const size_t i1 = vertices.GetKeyIndex(vertexKey4New, y4New.spin());
//     assert(i1 == 4);
//     vertices.SwapVertexPart(i1, vertices.NDown() - 2, y4New.spin());

//     vertices.PopBackVertexPart(x4New.spin());
//     vertices.PopBackVertexPart(y4New.spin());
//     vertices.RemoveVertex(4);

//     std::cout << "After remove 4 \n"
//               << std::endl;

//     vertices.Print();

//     std::cout << "verticesKeys = " << std::endl;
//     for (const auto &key : vertices.verticesKeysVec())
//     {
//         std::cout << key << std::endl;
//     }

//     const auto x7Up = vertices.atUp(6);
//     const auto y7Down = vertices.atDown(4);
//     const auto V7 = vertices.at(7);

//     assert(x7Up == V7.vStart());
//     assert(y7Down == V7.vEnd());
//     ASSERT_EQ(x7Up.spin(), FermionSpin_t::Up);
//     ASSERT_EQ(y7Down.spin(), FermionSpin_t::Down);
// }

// TEST(Vertices2DTest, TestBuildVertex)
// {
//     std::ifstream fin(FNAME);
//     Json jj;
//     fin >> jj;
//     fin.close();

//     Utilities::EngineTypeMt19937_t rng_(1 + jj["monteCarlo"]["seed"].get<size_t>());
//     Utilities::UniformRngMt19937_t urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0));

//     Diagrammatic::Vertices vertices;
//     Diagrammatic::VertexBuilder vertexBuilder(jj, Nc);

//     //Try Inserting a shit load of vertices
//     for (size_t ii = 0; ii < 200; ii++)
//     {
//         const auto v1 = vertexBuilder.BuildVertex(urng_);
//         vertices.AppendVertex(v1);
//     }

//     for (size_t ii = 0; ii < 9000000; ii++)
//     {
//         const auto v1 = vertexBuilder.BuildVertex(urng_);
//         vertices.AppendVertex(v1);

//         const size_t pp = urng_() * vertices.size();
//         const size_t vertexKey = vertices.GetKey(pp);
//         const auto x = vertices.at(pp).vStart();
//         const auto y = vertices.at(pp).vEnd();

//         if (x.spin() != y.spin())
//         {
//             assert(std::abs(x.tau() - y.tau()) < 1e-10);
//             const size_t ppUp = vertices.GetKeyIndex(vertexKey, x.spin());
//             const size_t ppDown = vertices.GetKeyIndex(vertexKey, y.spin());
//             const size_t kkUp = vertices.NUp();
//             const size_t kkUpm1 = kkUp - 1;
//             const size_t kkDown = vertices.NDown();
//             const size_t kkDownm1 = kkDown - 1;
//             assert(x.spin() == FermionSpin_t::Up);
//             vertices.SwapVertexPart(ppUp, kkUpm1, x.spin());
//             vertices.SwapVertexPart(ppDown, kkDownm1, y.spin());
//             vertices.RemoveVertex(pp);
//             vertices.PopBackVertexPart(x.spin());
//             vertices.PopBackVertexPart(y.spin());
//         }
//         else
//         {
//             assert(x.spin() == y.spin());

//             assert(std::abs(x.tau() - y.tau()) < 1e-10);
//             assert(x.orbital() != y.orbital());
//             assert(x.site() == y.site());

//             const size_t pp2Spin = vertices.GetKeyIndex(vertexKey + 1, y.spin());
//             const size_t kkSpin = (x.spin() == FermionSpin_t::Up) ? vertices.NUp() : vertices.NDown();
//             const size_t kkSpinm1 = kkSpin - 1;
//             const size_t kkSpinm2 = kkSpin - 2;

//             vertices.SwapVertexPart(pp2Spin, kkSpinm1, x.spin());

//             const size_t pp1SpinNew = vertices.GetKeyIndex(vertexKey, y.spin());
//             vertices.SwapVertexPart(pp1SpinNew, kkSpinm2, y.spin());

//             vertices.RemoveVertex(pp);
//             vertices.PopBackVertexPart(x.spin());
//             vertices.PopBackVertexPart(y.spin());
//         }
//     }

//     vertices.Clear();
// }

int main(int argc, char **argv)
{
    TestTools::RemoveFilesForTests();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
